// STL
#include <cmath>
#include <array>
#include <exception>
#include <algorithm>
#include <utility>
#include <fstream>
#include <regex>
#include <map>
// Boost
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/multi_array.hpp>
// Eigen
#include <Eigen/Eigen>
// Teem
#include <teem/nrrd.h>
// spurt
#include <misc/option_parse.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/progress.hpp>
#include <spline/spline.h>

#include <crease/crease_mc_cases.hpp>
#include <data/raster.hpp>

#include <vtkDoubleArray.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkPointData.h>

// idea:
// 1. Use Peikert's determinant idea to find ridge pouints along edges
// 2. Do no rely on linear variation of the determinant, instead
//    2a. Sample determinant along edge 
//    2b. Fit a cubic spline to those samples
//    2c. Find possible local extrema in between samples that could correspond to zero crossings
//    2d. Compute zero crossings
// 3. Check for MC configuration corresponding to found edges
// 4. If no MC case or multiple ridge pouints per edge, mark for subdivision
// 5. For each subdivision voxel, go back to 1.

typedef unsigned int uint;

uint verbose = 1;

constexpr int criterion_underflow = -10;
constexpr int determinant_underflow = -20;
constexpr int weak_ridge = -30;
constexpr int low_ridge = -40;
constexpr int several_ridge_points = -50;
constexpr int no_ridge_point = -1;
constexpr int unknown_error = -100;
constexpr int one_ridge_point = 0;
constexpr int odd_face_crossings = -60;

constexpr int one_triangle = 1;
constexpr int valid_mc_case = 2;
constexpr int invalid_mc_case = 3;
constexpr int exotic_case = 4;
constexpr int not_enough_edges = 5;

typedef Eigen::Matrix<double, 2, 1> vec2;
typedef Eigen::Matrix<double, 3, 1> vec3;
typedef Eigen::Matrix<double, 4, 1> vec4;
typedef Eigen::Matrix<double, 3, 3> mat3;
typedef Eigen::Matrix<double, 7, 1> vec7;
typedef Eigen::Matrix<double, 7, 2> mat72;
typedef Eigen::Matrix<uint, Eigen::Dynamic, Eigen::Dynamic> dyn_mat_type;
typedef Eigen::SelfAdjointEigenSolver<mat3> eigensolver_type;

typedef spurt::rgrid3d grid_type;
typedef grid_type::point_type point_type;
typedef grid_type::coord_type coord_type;
typedef grid_type::vec_type vec_type;
typedef grid_type::bounds_type bounds_type;
typedef grid_type::size_type size_type;
typedef grid_type::scalar_type scalar_type;
typedef std::pair<coord_type, point_type> invoxel_point_type;

typedef spurt::image3d<double> scalar_image;
typedef spurt::image3d<vec3> vector_image;
typedef spurt::image3d<mat3> matrix_image;

typedef spurt::nrrd_utils::nrrd_data_wrapper<double> wrapper_type;

struct Less
{
    bool operator()(const coord_type& a, const coord_type& b) const
    {
        if (a[0] < b[0]) return true;
        else if (b[0] < a[0]) return false;
        else if (a[1] < b[1]) return true;
        else if (b[1] < a[1]) return false;
        else return a[2] < b[2];
    }
};

template<typename V=vec3>
struct PosLexico
{
    typedef V vec_type;
    static constexpr double eps = 1.0e-9;
    bool operator()(const vec_type& p0, const vec_type& p1) const
    {
        if (p0[0] < p1[0]-eps)
            return true;
        else if (p1[0] < p0[0]-eps)
            return false;
        else if (p0[1] < p1[1]-eps)
            return true;
        else if (p1[1] < p0[1]-eps)
            return false;
        else
            return p0[2] < p1[2]-eps;
    }
};

template<typename V>
point_type to_point(const V& v)
{
    return point_type(static_cast<scalar_type>(v[0]), 
                      static_cast<scalar_type>(v[1]), 
                      static_cast<scalar_type>(v[2]));
}

constexpr double invalid = std::numeric_limits<double>::max();

struct EdgeParam {
    uint resolution;
    double min_distance;
    double min_value;
    double min_strength;
    double epsilon;
    double min_criterion;
    double min_determinant;
};

template<typename V>
struct ValueFiller {};

template<>
struct ValueFiller<double> {
    typedef spurt::nrrd_utils::nrrd_data_wrapper<double> wrapper_type;
    double operator()(const wrapper_type& w, size_t i, size_t stride) {
        return w[i];
    }
};

template<>
struct ValueFiller<vec3> {
    typedef spurt::nrrd_utils::nrrd_data_wrapper<double> wrapper_type;
    vec3 operator()(const wrapper_type& w, size_t i, size_t stride) {
        return vec3(w[i], w[i+stride], w[i+2*stride]);
    } 
};

template <>
struct ValueFiller<mat3>
{
    typedef spurt::nrrd_utils::nrrd_data_wrapper<double> wrapper_type;
    mat3 operator()(const wrapper_type &w, size_t i, size_t stride)
    {
        return mat3({ {w[i],          w[i+stride],   w[i+2*stride]},
                      {w[i+3*stride], w[i+4*stride], w[i+5*stride]},
                      {w[i+6*stride], w[i+7*stride], w[i+8*stride]} });
    }
};

template<typename Value>
void import_nrrd(spurt::image3d<Value>& out, const std::string& filename)
{
    typedef Value value_type;
    typedef ValueFiller<value_type> filler_type;
    Nrrd* nin = nrrdNew();
    typedef spurt::image3d<value_type> dataset_type;

    if (nrrdLoad(nin, filename.c_str(), NULL))
    {
        char* err = biffGetDone(NRRD);
        std::cerr << "Peikert ridge method: " << err << std::endl;
        exit(-1);
    }

    wrapper_type wrap(nin);
    bool is_scalar = (nin->dim == 3);

    // figure out how the data is stored
    int outerstride = 1;
    int innerstride = 1;
    int spacemin=0;
    if (!is_scalar) {
        if (nin->dim == 4) {
            if (nin->axis[0].size <= 9) {
                // individual scalar entries are in innermost dim
                spacemin = 1;
                outerstride = nin->axis[0].size;
                innerstride = 1;
            }
            else if (nin->axis[3].size <= 9) {
                // individual scalar entries are in outermost dim
                spacemin = 0;
                innerstride = nin->axis[0].size*nin->axis[1].size*nin->axis[2].size;
                outerstride = 1;
            }
            else {
                throw std::runtime_error("Unrecognized Nrrd shape");
            }
        }
        else if (nin->dim == 5) {
            if (nin->axis[0].size*nin->axis[1].size <= 9) {
                // individual scalar entries are in innermost dim
                spacemin = 2;
                outerstride = nin->axis[0].size*nin->axis[1].size;
                innerstride = 1;
            }
            else if (nin->axis[3].size*nin->axis[4].size <= 9) {
                // individual scalar entries are in outermost dim
                spacemin = 0;
                innerstride = nin->axis[0].size*nin->axis[1].size*nin->axis[2].size;
                outerstride = 1;
            }
            else {
                throw std::runtime_error("Unrecognized Nrrd shape");
            }
        }
    }
    coord_type res(nin->axis[spacemin].size, 
                   nin->axis[spacemin+1].size, 
                   nin->axis[spacemin+2].size);
    point_type origin(0,0,0);
    vec_type spacings(1,1,1);
    for (int d=0; d<3; ++d)
    {
        if (!std::isnan(nin->axis[spacemin+d].min))
        {
            origin[d] = nin->axis[spacemin+d].min;
        }
        if (!std::isnan(nin->axis[spacemin+d].spacing))
        {
            spacings[d] = nin->axis[spacemin+d].spacing;
        }
    }

    grid_type agrid(res, origin, spacings);
    std::vector<value_type> data;
    size_t nvalues = res[0] * res[1] * res[2];
    data.resize(nvalues);
    filler_type filler;
    for (size_t i=0; i<nvalues; ++i) {
        size_t idx = outerstride*i;
        data[i] = filler(wrap, idx, innerstride);
    }
    out = dataset_type(agrid, data);
}

template<typename T>
T sign(const T& value) {
    if (value >= 0) return T(1);
    else return T(-1);
}

struct EdgeID: public std::pair< coord_type, coord_type >
{
    typedef std::pair< coord_type, coord_type > base_type;

    EdgeID(const coord_type& c0, const coord_type& c1)
        : base_type(c0, c1) {
        Less less;
        if (less(c1, c0)) {
            coord_type tmp = base_type::first;
            base_type::first = base_type::second;
            base_type::second = tmp;
        }
    }

    bool operator<(const EdgeID& other) const {
        Less less;
        if (less(base_type::first, other.first)) return true;
        else if (less(other.first, base_type::first)) return false;
        else return less(base_type::second, other.second);
    }
};

template<typename T>
T get_value(const spurt::image3d<T>& volume, const coord_type& coord, int depth)
{
    if (depth == 0) return volume(coord[0], coord[1], coord[2]);
    int factor = 1 << depth;
    point_type p = point_type(coord) / static_cast<double>(factor);
    point_type q(std::trunc(p[0]), std::trunc(p[1]), std::trunc(p[2]));
    p -= q;
    if (nvis::norm(p)==0) return volume(q[0], q[1], q[2]);
    else return volume.value_in_voxel(coord_type(q), p);
}

point_type get_point(const grid_type& grid, const coord_type& coord, int depth)
{
    if (depth == 0) return grid(coord[0], coord[1], coord[2]);
    int factor = 1 << depth;
    point_type p = point_type(coord)/static_cast<double>(factor);
    point_type q(std::trunc(p[0]), std::trunc(p[1]), std::trunc(p[2]));
    p -= q;
    return grid(q[0], q[1], q[2]) + p*grid.spacing();
}

void voxel_to_edges(std::vector<EdgeID>& edges, const coord_type& voxelid)
{
    edges.clear();
    for (int i=0; i<12; ++i)
    {
        auto c = spurt::marching_cubes::canonical_edge_coordinates[i];
        coord_type shift1(c[0][0], c[0][1], c[0][2]);
        coord_type shift2(c[1][0], c[1][1], c[1][2]);
        edges.push_back(EdgeID(voxelid + shift1, voxelid+shift2));
    }
}

std::ostream& operator<<(std::ostream& os, const coord_type& id) {
    os << "[" << id[0] << "," << id[1] << "," << id[2] << "]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const vec3& id) {
    os << "[" << id[0] << "," << id[1] << "," << id[2] << "]";
    return os;
}

std::ostream& operator<<(std::ostream& os, const mat3& m) {
    os << "[" << "[" << m(0,0) << "," << m(0,1) << "," << m(0,2) << "],"
       <<        "[" << m(1,0) << "," << m(1,1) << "," << m(1,2) << "]," 
       <<        "[" << m(2,0) << "," << m(2,1) << "," << m(2,2) << "] ]";
    return os;
}

template<typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& a)
{
    os << "[";
    std::for_each(a.begin(), a.end(), [&](T v) {
        os << v << ", ";
    });
    os << "]";
    return os;
}

std::ostream &
operator<<(std::ostream &os, const EdgeID &e)
{
    os << "[" << e.first << " - " << e.second << "]";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    os << "[";
    if (v.size() > 0) 
    {
        for (uint i=0; i<v.size()-1; ++i) os << v[i] << ", ";
        os << v.back();
    }
    os << "]";
    return os;
}

template<typename T>
T invlinear(const T& f0, const T& f1, double umin=0, double umax=1) {
    // f = f0 + (u - umin) / (umax - umin) * (f1 - f0)
    // f = 0 <=> -f0 / (f1 - f0) = (u - umin) / (umax - umin)
    // f = 0 <=> u = -f0 / (f1 - f0) * (umax - umin) + umin    
    return umin - f0/(f1-f0)*(umax-umin);
}

template<typename T>
T linear(double u, const T& v0, const T& v1) {
    return (1.-u)*v0 + u*v1;
}

std::pair<double, vec3> evmin(const mat3& H) {
    eigensolver_type solver;
    // std::cout << "About to compute eigenvalues for " << H << '\n';
    mat3 M = H;
    solver.compute(H);
    // std::cout << "done.\n";
    if (H.maxCoeff() > 1000.) {
        std::cout << "matrix H is " << H << '\n';
        std::cout << "copy of H is " << M << '\n';
        throw std::runtime_error("Invalid hessian value: " + std::to_string(H.maxCoeff()));
    }
    return std::make_pair(solver.eigenvalues()[0], solver.eigenvectors().col(0));
}

std::pair<vec3, mat3> evall(const mat3& H) {
    eigensolver_type solver;
    solver.compute(H);
    return std::make_pair(solver.eigenvalues(), solver.eigenvectors());
}

double determinant(const vec3& g, const mat3& H, bool normalize=false) {
    mat3 A;
    A.col(0) = g;
    A.col(1) = H * g;
    A.col(2) = H * A.col(1);

    if (normalize) {
        A.col(0) /= g.norm();
        A.col(1) /= A.col(1).norm();
        A.col(2) /= A.col(2).norm();
    }
    return A.determinant();
};

template<typename T>
class LinearFunction {
public:
    typedef T value_type;

    LinearFunction() : m_v0(), m_v1() {}

    LinearFunction(const value_type& val0, const value_type& val1) :
        m_v0(val0), m_v1(val1) {}

    value_type operator()(double u) const {
        return linear(u, m_v0, m_v1);
    }

private:
    value_type m_v0, m_v1;
};

typedef LinearFunction<double> scalar_function;
typedef LinearFunction<vec3> vector_function;
typedef LinearFunction<mat3> matrix_function;

class C0Criterion {
public:
    C0Criterion(const vector_function& gfunc, const matrix_function& hfunc) {
        m_gfunc = gfunc;
        m_hfunc = hfunc;
    }

    double operator()(double u) const {
        auto g = m_gfunc(u);
        auto h = m_hfunc(u);
        return determinant(g, h);
    }
private:
    vector_function m_gfunc;
    matrix_function m_hfunc;
};

template<typename Func>
double find_zero_crossing(const Func& func, double umin = 0, double umax = 1, double eps=1.0e-12, uint maxn=100) {
    // pruint('find zero crossing')
    double fmin = func(umin);
    double fmax = func(umax);
    if (fmin * fmax >= 0) return -1;
    // det(0) * det(1) < 0

    double f = 1;  // default value
    double u;
    for (uint n=0; ((umax - umin) > 1.0e-9 or std::fabs(f) > eps) and n < maxn; ++n) {
        u = 0.5 * (umin + umax);
        f = func(u);
        if (fmin*f < 0) {
            umax = u;
            fmax = f;
        } 
        else {
            umin = u; 
            fmin = f;
        } 
    }
    return u;
}

std::pair<vec3, double> project(const vec3& g, const mat3& h) {
    auto r = evall(h);
    vec3 evals = r.first;
    mat3 evecs = r.second;
    vec3 coords = (evecs.transpose()*g).array().abs().matrix();
    return std::make_pair(coords, evals[0]);
}

bool check_solution(double u, const vector_function& gfunc, 
                    const matrix_function& hfunc,
                    bool verbose=false) {
    auto r = evall(hfunc(u));
    vec3 evals = r.first;
    mat3 evecs = r.second;
    vec3 coords = (evecs.transpose()*gfunc(u)).array().abs().matrix();
    if (verbose) {
        std::cout << "eigenvalues: " << evals << '\n'
        << "eigenvectors: " << evecs << '\n'
        << "coordinates: " << coords << '\n';
    }

    Eigen::Index which;
    double minval = coords.minCoeff(&which);
    return (which == 0) && minval < 1.0e-6*coords.mean();
}

typedef std::array<vec3, 3> triangle_type;
int triangulate(std::vector<triangle_type>& out, 
                std::map<int, std::vector<vec3>> &edges, 
                std::ostream& os)
{
    os << "triangles: edges contains " << edges.size() << " edges with ridge pouints\n";
    out.clear();

    std::map<int, int> edge_to_nbsols;
    int nbsols = 0;
    for (auto iter=edges.begin(); iter!=edges.end(); ++iter) 
    {
        edge_to_nbsols[iter->first] = iter->second.size();
        nbsols += iter->second.size();
    }

    os << "There are " << nbsols << " ridge points over " << edges.size() << " in total\n";

    // if we only have 3 or less pouints, things are easy
    if (edges.size() == 3 && nbsols == 3) 
    {
        triangle_type T;
        int i=0;
        for (auto iter=edges.begin(); iter!=edges.end(); ++iter, ++i) 
        {
            T[i] = iter->second[0];
        }
        out.push_back(T);
        os << "3 points on 3 edges: success!\n";
        return one_triangle;
    }
    else if (edges.size() < 3)
    {
        os << "We have only 2 (or less) edges in input of triangulation.\n"
           << "Giving up (Case N<3)\n";
        return not_enough_edges;
    }

    if (nbsols == edges.size()) 
    {
        // Calculate edge case number
        int edge_case = 0;
        for (auto iter=edges.begin(); iter!=edges.end(); ++iter) 
        {
            edge_case += 1 << iter->first;
        }
        os << "edge_case = " << edge_case << '\n';
        uint triangle_case = spurt::marching_cubes::edge_code_to_case_id[edge_case];
        os << "triangle_case = " << triangle_case << '\n';

        if (triangle_case == -1) // invalid
        {
            os << "the edges do not match a valid MC case...\n Giving up. (Case NoMC)\n";
            return invalid_mc_case;
        }

        auto indices = spurt::marching_cubes::triTable[triangle_case];
        for (int i=0; i<15 && indices[i]!=-1; i+=3) 
        {
            out.push_back(triangle_type({{edges[indices[i]][0], edges[indices[i+1]][0], edges[indices[i+2]][0]}}));
            const triangle_type& t = out.back();
            os << "added triangle: " << t << '\n';
        }
        os << "A valid MC case was found and " << out.size() << " triangles "
           << "were created.\n";
        return valid_mc_case;
    }
    else
    {
        // 3 or more edges, some of which have several pouints
        // for now, do nothing
        return exotic_case;
    }
}

void fit_spline(tk::spline& spline, const std::vector<double>& u, C0Criterion& f) {
    std::vector<double> values;
    std::for_each(u.begin(), u.end(), [&](double x){ values.push_back(f(x)); });
    spline.set_points(u, values);
}

std::vector<double> bisection_spline_local_extremum(tk::spline &, double, double, bool);

std::vector<double> newton_spline_local_extremum(tk::spline &spline, double u0, double u1, bool verbose)
{
    // Newton:
    // f(x+dx) = 0 = f(x) + f'(x)dx <=> dx = -f(x)/f'(x)
    double f0 = spline(u0);
    double f1 = spline(u1);
    double u = -f0/(f1-f0); // initial guess (linear solution)
    double f = spline.deriv(1, u);
    uint iter = 0;
    std::vector<double> all_us;
    if (verbose)
    {
        std::cout << "Looking for local extremum\n";
    }
    while (std::abs(f) > 1.0e-6 && iter < 10)
    {
        double d2f = spline.deriv(2, u);
        all_us.push_back(u);
        if (d2f == 0)
        {
            if (verbose) std::cout << "vanishing 2nd derivative: resorting to bisection\n";
            return bisection_spline_local_extremum(spline, u0, u1, verbose);
        }
        double dx = -f / d2f;
        double uu = u + dx;
        if (uu > u1) u = u1;
        else if (uu < u0) u = u0;
        else u = uu;
        if (verbose)
        {
            std::cout << "u=" << u << ", df/dt=" << f << ", d2f/dt2=" << d2f << ", dx=" << dx << '\n';
        }
        f = spline.deriv(1, u);
        if (verbose)
        {
            std::cout << "df/dt=" << f << '\n';
        }
        iter += 1;
    }
    return all_us;
}

std::vector<double> bisection_spline_local_extremum(tk::spline& spline, double u0, double u1, bool verbose) 
{
    double f0 = spline.deriv(1, u0);
    double f1 = spline.deriv(1, u1);
    double u, t;
    std::vector<double> all_us;
    for (int i=0; i<10; ++i) 
    {
        t = -f0/(f1-f0);
        u = (1-t)*u0 + t*u1;
        all_us.push_back(u);
        double f = spline.deriv(1, u);
        if (verbose)
        {
            std::cout << "bisection search: iteration " << i << ", df/dt=" << f << ", u0=" << u0 << ", u1=" << u1 << '\n';
        }
        if (f*f0 > 0)
        {
            u0 = u;
            f0 = f;
        }
        else
        {
            u1 = u;
            f1 = f;
        }
    }
    t = -f0/(f1-f0);
    u = (1 - t) * u0 + t * u1;
    all_us.push_back(u);
    double f = spline.deriv(1, u);
    if (verbose)
    {
        std::cout << "final derivative: " << f << '\n';
    }
    return all_us;
}

std::pair<int, std::vector<double>>
analyze_edge(const EdgeID& anedge, const vector_image& gradient,
             const matrix_image& hessian, const EdgeParam &params,
             int depth = 0,
             bool verbose = false)
{   
    vector_function gfunc(get_value(gradient, anedge.first, depth), 
                          get_value(gradient, anedge.second, depth));
    matrix_function hfunc(get_value(hessian, anedge.first, depth), 
                          get_value(hessian, anedge.second, depth));
    C0Criterion cfunc(gfunc, hfunc);

    // 0 --- 1/3 --- 2/3 --- 3/3: N=4, du=1/3 
    
    double du = 1./static_cast<double>(params.resolution-1);
    std::map<double, double> values;
    std::vector<double> xs;
    
    for (int i=0; i<params.resolution; ++i) {
        xs.push_back(std::min(i*du, 1.));
    }
    tk::spline init_spline;
    fit_spline(init_spline, xs, cfunc);
    std::vector<double> derivs;
    for (int i=0; i<xs.size(); ++i) {
        derivs.push_back(init_spline.deriv(1, xs[i]));
        values[xs[i]] = init_spline(xs[i]);
        if (std::isnan(derivs.back()) || std::isnan(values[xs[i]]))
        {
            std::cout << "ERROR in analyze_edge\n";
            std::cout << "edgeid was " << anedge << '\n';
        }
    }

    // check if derivative has constant sign in which case we are good
    double asign = sign(derivs[0]);
    std::vector<std::pair<double, double>> where;
    if (verbose)
    {
        std::cout << "sampling spline derivative along edge\n";
        std::cout << "sample #0: "  << derivs[0] << '\n';
    }
    // find root of derivative iff:
    // 1. downward then upward + both values positive: local min may be negative
    // 2. upward then downward + both values negative: local max may be positive
    for (int i=1; i<derivs.size(); ++i) {
        if (verbose)
        {
            std::cout << "sample #" << i << ": " << derivs[i] << '\n';
            std::cout << "value #" << i << ": " << values[xs[i]] << '\n';
        }
        double newsign = sign(derivs[i]);
        if (asign * newsign < 0) {
            if (verbose)
            {
                std::cout << "sign change\n";
            }
            if ((asign > 0 && values[i-1] < 0 && values[i] < 0) || // local maximum may be positive
                (asign < 0 && values[i-1] > 0 && values[i] > 0)) // local minimum may be negative
                where.push_back(std::make_pair(xs[i-1], xs[i]));
            asign = newsign;
        }
    }
    if (verbose)
        std::cout << where.size() << " potential hidden zero crossings detected\n";

    for (int i=0; i<where.size(); ++i)
    {
        double u0 = where[i].first;
        double u1 = where[i].second;
        // double u = bisection_spline_local_extremum(init_spline, u0, u1, verbose);
        std::vector<double> new_us = newton_spline_local_extremum(init_spline, u0, u1, verbose);
        // found zero crossing of derivative
        std::for_each(new_us.begin(), new_us.end(), [&](double u)
        {
            values[u] = cfunc(u);
        });
    }

    if (verbose) {
        std::cout << "After spline-based refinement, samples are:\n";
        for (auto iter=values.begin(); iter!=values.end(); ++iter) {
            std::cout << std::setprecision(32) << iter->first << ": " << iter->second << '\n';
        }
    }
    std::vector<double> ridge_points;

    for (auto iter=values.begin(); iter!=values.end(); ++iter)
    {
        auto next = iter;
        ++next;
        if (next == values.end()) break;
        double u0 = iter->first;
        double f0 = iter->second;
        double u1 = next->first;
        double f1 = next->second;
        if (f0*f1 < 0)
        {
            if (verbose) {
                std::cout << "zero crossing found between " << u0 << " and " << u1 << '\n';
            }
            double u = find_zero_crossing(cfunc, u0, u1);
            if (check_solution(u, gfunc, hfunc, verbose))
            {
                if (verbose)
                    std::cout << "ridge point found\n";
                ridge_points.push_back(u);
            }
            else if (verbose) {
                std::cout << "no ridge point found\n";
            }
        }
    }

    if (ridge_points.size() == 1)
        return std::make_pair(one_ridge_point, ridge_points);
    else if (ridge_points.size() > 1)
        return std::make_pair(several_ridge_points, ridge_points);
    else
        return std::make_pair(no_ridge_point, std::vector<double>());
}

std::pair<double, double> evaluate(double u, const EdgeID, const scalar_image&, 
                                   const matrix_image&, int);

int check_for_singularity(const std::map<int, std::vector<vec3>>& found)
{
    // loop over faces to identify broken ones
    for (int faceid=0; faceid<6; ++faceid) 
    {
        auto face_edges = spurt::marching_cubes::canonical_face_edge_indices[faceid];
        std::vector<std::pair<int, vec3>> edge_crossings;
        for (int i=0; i<4; ++i)
        {
            int e = face_edges[i];
            auto iter = found.find(e);
            if (iter != found.end()) {
                if (iter->second.size() == 1) 
                    edge_crossings.push_back(std::make_pair(i, iter->second[0]));
                else if (iter->second.size() > 1)
                    return several_ridge_points;
            }
        }
        if (edge_crossings.size() % 2)
        {
            // odd number of face edge crossings: Hessian singularity
            return odd_face_crossings;
        }
    }
    return 0;
}

int process_voxel(std::vector<triangle_type> &all_triangles,
                  const coord_type& avoxel,
                  const vector_image& gradient, 
                  const matrix_image& hessian,
                  std::map<EdgeID, std::vector<vec3>> &known_edges,
                  const EdgeParam &params,
                  int depth=0,
                  bool verbose=false)
{
    std::map<int, std::vector<vec3>> found;
    std::vector<EdgeID> edges;
    voxel_to_edges(edges, avoxel);
    std::ostringstream os;
    int ncrossings = 0;
    bool multiple_xings = false;
    bool failed = false;
    for (int i = 0; i < 12; ++i)
    {
        EdgeID edgeid=edges[i];
        auto iter = known_edges.find(edgeid);
        if (iter != known_edges.end())
        {
            if (verbose) std::cout << "skipping known edge\n";
            if (!iter->second.empty())
                found[i] = iter->second;
        }
        else
        {
            auto result = analyze_edge(edgeid, gradient, hessian, params, depth, verbose);
            auto coordinates = result.second;

            known_edges[edgeid] = std::vector<vec3>();
            if (!coordinates.empty())
            {
                found[i] = std::vector<vec3>();
                for (int k = 0; k < coordinates.size(); ++k)
                {
                    double u = coordinates[k];
                    auto p = (1.-u)*get_point(gradient.grid(), edgeid.first, depth) + 
                              u*get_point(gradient.grid(), edgeid.second, depth);
                    found[i].push_back(vec3(p[0], p[1], p[2]));
                    // os << "point with local coordinate " << u << " on edge " << i << " (" << edgeid.first << ", " << edgeid.second << ") = " << get_point(gradient.grid(), edgeid.first, depth) << " - " << get_point(gradient.grid(), edgeid.second, depth) << " of voxel " << avoxel << " has coordinates " << p << '\n';
                }
                known_edges[edgeid] = found[i];
                if (coordinates.size() > 1)
                {
                    multiple_xings = true;
                }
                ncrossings += coordinates.size();
            }
        }
    }

    int status = check_for_singularity(found);
    if (status != 0)
    {
        return status;
    }

    if (!multiple_xings && found.size() >= 3)
    {
        std::vector<triangle_type> tris;
        std::ostringstream log;
        int tri_case;
        if (verbose)
            tri_case = triangulate(tris, found, std::cout);
        else
            tri_case = triangulate(tris, found, log);

        if (tris.size() > 0)
        {
            std::for_each(tris.begin(), tris.end(), [&](auto T)
                          { all_triangles.push_back(T); });
            // std::cout << os.str();
            return tris.size();
        }
        else
        {
            if (verbose)
                std::cout << "Triangulation failed\n";
            return -ncrossings;
        }
    }
    else if (failed)
        return -ncrossings;
    else
        return 0;
}

std::pair<double, double> evaluate(double u, const EdgeID edgeid,
                                   const scalar_image& values, 
                                   const matrix_image& hessian,
                                   int depth=0) 
{
    const coord_type& c0 = edgeid.first;
    const coord_type& c1 = edgeid.second;
    mat3 H0 = get_value<mat3>(hessian, c0, depth);
    double v0 = get_value<double>(values, c0, depth);
    mat3 H1 = get_value<mat3>(hessian, c1, depth);
    double v1 = get_value<double>(values, c1, depth);
    scalar_function v(v0, v1);
    matrix_function h(H0, H1);
    return std::make_pair(v(u), evmin(h(u)).first);
}

std::pair<int, std::vector<double>> 
process_edge(const EdgeID& edgeid, const scalar_image& value,
             const vector_image& gradient, const matrix_image& hessian,
             const EdgeParam &params, int depth=0, bool verbose = false)
{
    if (verbose) {
        std::cout << "\n\nprocessing edge " << edgeid << '\n';
    }
    std::vector<double> solutions; // initially empty

    const coord_type& c0 = edgeid.first;
    const coord_type& c1 = edgeid.second;
    
    scalar_function sfunc(get_value<double>(value, c0, depth), 
                          get_value<double>(value, c1, depth));
    vector_function gfunc(get_value<vec3>(gradient, c0, depth), 
                          get_value<vec3>(gradient, c1, depth));
    matrix_function hfunc(get_value<mat3>(hessian, c0, depth), 
                          get_value<mat3>(hessian, c1, depth));
    double lmin0 = evmin(hfunc(0)).first;
    double lmin1 = evmin(hfunc(1)).first;
    C0Criterion cfunc(gfunc, hfunc);

    // check if we are dealing with a degenerate case or one that
    // does not meet prescribed filtering criteria
    if (std::max(std::fabs(cfunc(0)), std::fabs(cfunc(1))) < params.min_criterion) 
    {
        return std::make_pair(criterion_underflow, solutions);
    }
    else if (std::max(lmin0, lmin1) > -std::fabs(params.min_strength)) 
    {
        return std::make_pair(weak_ridge, solutions);
    }
    else if (std::max(sfunc(0), sfunc(1)) < params.min_value) 
    {
        return std::make_pair(low_ridge, solutions);
    }
    
    // initial uniform sampling
    std::vector<double> dets;
    double step = 1./double(params.resolution-1);
    for (double u=0; u<1+step; u+=step) 
    {
        dets.push_back(cfunc(u));
        if (verbose) 
        {
            std::cout << "det[" << u << "] = " << dets.back() << '\n';
        }
    }
    
    // adaptive zero-crossing search
    for (int i=0; i<dets.size()-1; ++i) 
    {
        double u0 = i*step;
        double u1 = u0 + step;
        const double& det0 = dets[i];
        const double& det1 = dets[i+1];
        if (det0 * det1 < 0) 
        {
            if (verbose) 
            {
                std::cout << "zero crossing found between " << u0 << " and " << u1 << '\n';
            }
            double u = find_zero_crossing(cfunc, u0, u1);
            if (verbose) 
            {
                std::cout << "zero crossing found at u=" << u << '\n';
            }
            if (u>=0 && check_solution(u, gfunc, hfunc, verbose)) 
            {
                solutions.push_back(u);
                if (verbose)
                    std::cout << "zero crossing found to be a ridge pouint\n";
            }
            else if (verbose) 
            {
                std::cout << "zero crossing is not a valid ridge pouint\n";
            }
        }
    }

    if (solutions.size() == 1) 
    {
        return std::make_pair(one_ridge_point, solutions);
    }
    else if (solutions.size() == 0) 
    {
        return std::make_pair(no_ridge_point, solutions); // code for nothing found
    }
    else 
    {
        return std::make_pair(several_ridge_points, solutions);
    }
}

template<typename T>
struct parser_traits;

template<>
struct parser_traits<int> 
{
    static std::regex get()
    {
        return std::regex("([-0-9]+)");
    }
    static int cast(const std::string& s) 
    {
        return std::stoi(s);
    }
};

template<>
struct parser_traits<double> 
{
    static std::regex get()
    {
        return std::regex("([+-]? *[0-9]+(\\.[0-9]+)?)");
    }
    static uint cast(const std::string &s)
    {
        return std::stod(s);
    }
};

template<typename T=uint>
void parse_values(std::vector<T>& out, const std::string& str, size_t n) 
{
    std::regex myregex = parser_traits<T>::get();
    auto begin = std::sregex_iterator(str.begin(), str.end(), myregex);
    auto end = std::sregex_iterator();
    // std::cout << "Found " << std::distance(begin, end) << " values\n";
    if (std::distance(begin, end) == n)
    {
        uint i = 0;
        for (std::sregex_iterator iter = begin; iter != end; ++iter, ++i)
        {
            std::smatch match = *iter;
            // std::cout << "value=" << match.str() << std::endl;
            out.push_back(parser_traits<T>::cast(match.str()));
        }
    }
    else 
    {
        throw std::runtime_error("invalid input");
    }
}

void save_sampled_voxel(const std::string& filename, coord_type& voxel_id, 
                        const vector_image& gradient, 
                        const matrix_image& hessian,
                        const size_t resolution=101)
{
    std::vector<double> det_values(2 * resolution * resolution * resolution);
    double dh = 1./static_cast<double>(resolution-1);
    for (int k = 0; k < resolution; ++k)
    {
        double z = k * dh;
        for (int j = 0; j < resolution; ++j)
        {
            double y = j * dh;
            for (int i = 0; i < resolution; ++i)
            {
                double x = i * dh;
                vec3 g = gradient.value_in_voxel(voxel_id, point_type(x,y,z));
                mat3 H = hessian.value_in_voxel(voxel_id, point_type(x,y,z));
                det_values[2 * (i + resolution * (j + resolution * k))] = determinant(g, H);
                det_values[2 * (i + resolution * (j + resolution * k)) + 1] = (project(g, H).first)[0];
            }
        }
    }
    size_t sizes[4] = {2, resolution, resolution, resolution};
    std::ostringstream os;
    os << filename << "_voxel_" << voxel_id[0] << "_" << voxel_id[1] << "_" << voxel_id[2] << ".nrrd";
    spurt::nrrd_utils::writeNrrd((void *)&det_values[0], os.str(), nrrdTypeDouble, 4, sizes);
}

int check_voxel(const coord_type &voxelid, const scalar_image& value,
                const matrix_image& hessian, double minval, double minstr,
                int depth=0)
{
    // check if current voxel satisfies threshold requirements
    std::vector<double> strs(8);
    std::vector<double> vals(8);
    for (int i=0; i<8; ++i)
    {
        auto shift = spurt::marching_cubes::vertices[i];
        coord_type pid = voxelid + coord_type(shift[0], shift[1], shift[2]);
        mat3 h = get_value<mat3>(hessian, pid, depth);
        strs[i] = evmin(h).first;
        vals[i] = get_value(value, pid, depth);
    }
    if (*std::max_element(vals.begin(), vals.end()) < minval) return -1;
    if (*std::min_element(strs.begin(), strs.end()) > minstr) return -2;
    else return 1;
}

void export_triangles(const std::vector<triangle_type> &triangles,  
                      const std::vector<std::array<double, 2>>& attributes,
                      const std::string &filename)
{
    std::map<vec3, size_t, PosLexico<vec3>> pos_to_id;
    std::map<size_t, size_t> new_to_old_id;
    std::vector<std::array<size_t, 3>> tris;
    size_t id = 0;
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        std::array<size_t, 3> atri;
        for (int j = 0; j < 3; ++j)
        {
            const vec3 &p = triangles[i][j];
            auto iter = pos_to_id.find(p);
            if (iter != pos_to_id.end())
            {
                atri[j] = iter->second;
            }
            else
            {
                atri[j] = id;
                new_to_old_id[id] = 3*i+j;
                pos_to_id[p] = id++;
            }
        }
        tris.push_back(atri);
        // std::cout << "atri=" << atri << '\n';
        // std::cout << "points were " << triangles[i][0] << ", " << triangles[i][1] << ", " << triangles[i][2] << '\n';
    }
    vtkDoubleArray *coords = vtkDoubleArray::New();
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(pos_to_id.size());
    vtkDoubleArray *vals = vtkDoubleArray::New();
    vals->SetName("Values");
    vals->SetNumberOfComponents(1);
    vals->SetNumberOfTuples(pos_to_id.size());
    vtkDoubleArray *strs = vtkDoubleArray::New();
    strs->SetNumberOfComponents(1);
    strs->SetNumberOfTuples(pos_to_id.size());
    strs->SetName("Strength");
    for (auto iter = pos_to_id.begin(); iter != pos_to_id.end(); ++iter)
    {
        coords->SetTuple3(iter->second, iter->first[0], iter->first[1], iter->first[2]);
        vals->SetTuple1(iter->second, attributes[new_to_old_id[iter->second]][0]);
        strs->SetTuple1(iter->second, attributes[new_to_old_id[iter->second]][1]);
    }
    vtkPoints *points = vtkPoints::New();
    points->SetData(coords);
    vtkCellArray *cells = vtkCellArray::New();
    for (size_t i = 0; i < triangles.size(); ++i)
    {
        cells->InsertNextCell(3);
        cells->InsertCellPoint(tris[i][0]);
        cells->InsertCellPoint(tris[i][1]);
        cells->InsertCellPoint(tris[i][2]);
    }
    vtkPolyData *poly = vtkPolyData::New();
    poly->SetPoints(points);
    poly->SetPolys(cells);
    poly->GetPointData()->AddArray(vals);
    poly->GetPointData()->AddArray(strs);
    vtkXMLPolyDataWriter *writer = vtkXMLPolyDataWriter::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(poly);
    writer->Write();
}

int main(int argc, const char *argv[])
{
    std::string data_name, gradient_name, hessian_name, output_name;
    double minval, minstr, eps, mind;
    uint res, niter;
    bool verbose;
    coord_type voxel_id;
    std::pair<coord_type, coord_type> bounds;
    std::string voxel_str = "(-1, -1, -1)";
    std::string bounds_str = "(0, -1; 0, -1; 0, -1)";
    int maxdepth = 3;
    

    namespace cl = spurt::command_line;
    cl::option_traits
        required_group(true, false, "Required Options"),
        positional_group(true, true, "Positional Group"),
        optional_group(false, false, "Optional Group");

    cl::option_parser parser(argv[0],
                             "Extract ridge surfaces from scalar volume using\nPeikert and Sadlo's level set method");
    try 
    {
        parser.use_short_symbols(false);
        parser.use_brackets(true);

        parser.add_value("value", data_name, "Scalar raster filename", required_group);
        parser.add_value("gradient", gradient_name, "Gradient raster filename", required_group);
        parser.add_value("hessian", hessian_name, "Hessian raster filename", required_group);
        parser.add_value("output", output_name, "Output filename", required_group);
        parser.add_value("minval", minval, 0., "Min scalar value", optional_group);
        parser.add_value("minstr", minstr, 0., "Min ridge strength", optional_group);
        parser.add_value("eps", eps, 1.0e-9, "Numerical precision", optional_group);
        parser.add_value("mindist", mind, 0.02, "Min sampling distance", optional_group);
        parser.add_value("verbose", verbose, false, "Verbose output", optional_group);
        parser.add_value("res", res, 10, "Initial (coarse) sampling resolution to find ridge pouints", optional_group);
        parser.add_value("niter", niter, 100, "Max number of iterations of the zero crossing algorithm", optional_group);
        parser.add_value("voxel", voxel_str, voxel_str, "Coordinates of single voxel to process");
        parser.add_value("bounds", bounds_str, bounds_str, "Bounds of domain to consider");
        parser.add_value("maxdepth", maxdepth, maxdepth, "Maximum refinement depth");
        parser.parse(argc, argv);
    }
    catch (std::runtime_error &e)
    {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    catch (std::exception &e)
    {
        std::cerr << "ERROR(2): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options enteredso far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }

    bool valid_voxel_id = false;
    bool valid_bounds = false;
    {
        std::vector<int> iv;
        parse_values<int>(iv, voxel_str, 3);
        std::cout << iv << '\n';
        if (iv[0] != -1)
        {
            voxel_id = { iv[0], iv[1], iv[2] };
            valid_voxel_id = true;
        }

        std::vector<int> dv;
        parse_values<int>(dv, bounds_str, 6);
        std::cout << dv << '\n';
        if (dv[1] != -1)
        {
            bounds.first = { dv[0], dv[2], dv[4] };
            bounds.second = { dv[1], dv[3], dv[5] };
            valid_bounds = true;
        }
    }

    scalar_image values;
    vector_image gradients;
    matrix_image hessians;
    import_nrrd(values, data_name);
    import_nrrd(gradients, gradient_name);
    import_nrrd(hessians, hessian_name);

    EdgeParam params;
    params.resolution = res;
    params.epsilon = eps;
    params.min_distance = mind;
    params.min_value = minval;
    params.min_strength = minstr;
    params.min_determinant = 1.0e-9;
    params.min_criterion = 1.0e-9;

    auto shape = values.grid().resolution();
    if (verbose) std::cout << "grid shape: " << shape << '\n';
    std::vector< coord_type > voxel_ids;

    if (valid_voxel_id)
    {
        std::cout << "valid voxel id\n";
        valid_voxel_id = true;
        voxel_ids.clear();
        voxel_ids.push_back(voxel_id);
        verbose = true;
    }
    else if (valid_bounds)
    {
        std::cout << "valid bounds\n";
        voxel_ids.clear();
        const coord_type& l = bounds.first;
        const coord_type& h = bounds.second;
        for (int k=std::max<int>(0, l[2]); k<std::min<int>(shape[2]-1, h[2]-1); k++)
        {
            for (int j=std::max<int>(0, l[1]); j<std::min<int>(shape[1]-1, h[1]-1); ++j)
            {
                for (int i=std::max<int>(0, l[0]); i<std::min<int>(shape[0]-1, h[0]-1); ++i)
                {
                    voxel_ids.push_back(coord_type(i,j,k));
                }
            }
        }
    }
    else 
    {
        std::cout << "including every voxel\n";
        int _n = 0;
        voxel_ids.resize((shape[0]-1)*(shape[1]-1)*(shape[2]-1));
        for (int k = 0; k < shape[2] - 1; ++k)
        {
            for (int j = 0; j < shape[1] - 1; ++j)
            {
                for (int i = 0; i < shape[0] - 1; ++i, ++_n)
                {
                    voxel_ids[_n] = coord_type(i, j, k);
                }
            }
        }
    }
    std::cout << "There are " << voxel_ids.size() << " voxels in input for a total of " << 12*voxel_ids.size() << " (redundant) edges\n";

    int nskipped = 0;
    int nfailed_to_triangulate = 0;
    int nnot_enough_points = 0;
    int nsucceeded = 0;
    int nprocessed = 0;
    int nexcluded = 0;
    int nnone = 0;
    int neven = 0;
    int nweak = 0;
    int nlow = 0;
    int nunderflow = 0;
    int nfiltered = 0;
    int nfailed_even = 0;
    int nfailed_odd = 0;
    int nfailed_multi = 0;
    std::vector<triangle_type> all_triangles;
    std::vector<vec4> rejected;

    if (valid_voxel_id)
        save_sampled_voxel(output_name, voxel_id, gradients, hessians);
    
    spurt::ProgressDisplay progress;
    std::vector<coord_type> failed_voxels;
    std::vector<uint> voxel_to_edge;
    std::map<coord_type, std::vector<uint>, Less> voxel_counter;
    int skipped_voxels=0;
    int nedges = 0;
    int nfound = 0;

    srand48(130819751900);
    for (int depth=0; depth<maxdepth; ++depth)
    {
        std::cout << "\nProcessing depth=" << depth << ", " << voxel_ids.size() << " voxels\n";
        // we reset the edge information at each new depth since we have not encountered
        // those half edges yet.
        std::map<EdgeID, std::vector<vec3>> all_processed_edges;
        progress.begin(voxel_ids.size(), "Extract ridges", 1000000, "tris: 0, failed: 0, skipped: 0, edges: 0, found: 0");
        for (int n=0; n<voxel_ids.size(); ++n) 
        {
            std::string update_str = "tris: " + std::to_string(all_triangles.size()) + 
            ", failed: " + std::to_string(failed_voxels.size()) + ", skipped: " + std::to_string(skipped_voxels)
            + ", edges: " + std::to_string(all_processed_edges.size()) + ", found: " + std::to_string(nfound);
            progress.update(n, update_str);
            coord_type id = voxel_ids[n];
            bool voxel_failed = false;

            // check if current voxel satisfies threshold requirements
            int ret = check_voxel(id, values, hessians, minval, minstr, depth);
            if (ret < 0)
            {
                skipped_voxels++;
                // if (verbose) {
                //     if (ret == -1)
                //         std::cout << "voxel was discarded by value\n";
                //     else if (ret == -2)
                //         std::cout << "voxel was discarded by ridge strength\n";
                // }
                continue;
            }
            int ntris = process_voxel(all_triangles, id, gradients, hessians, all_processed_edges, params, depth, verbose);
            if (ntris < 0) 
            {
                if (verbose) std::cout << "failed voxel\n";
                failed_voxels.push_back(id);
                if (ntris == several_ridge_points)
                {
                    nfailed_multi++;
                }
                else if (ntris == odd_face_crossings)
                {   
                    nfailed_odd++;
                }
            }
            else if (ntris > 0) 
            {   
                nfound++;
                if (verbose) std::cout << ntris << " triangles found\n";

            }
        }
        std::cout << "\nThere were " << failed_voxels.size() << " failed voxels\n";
        std::cout << "of those:\n" 
                  << "* " << nfailed_multi << " voxels had edges with multiple crossings (" << 100.*float(nfailed_multi)/float(failed_voxels.size()) << "%)\n"
                  << "* " << nfailed_odd << " voxels had odd number of face crossings ("  <<  100.*float(nfailed_odd)/float(failed_voxels.size()) << "%)\n";

        std::vector<vec3> all_edge_points;
        std::string suffix = "_d=" + std::to_string(depth);
        for (auto iter=all_processed_edges.begin(); iter!=all_processed_edges.end(); ++iter)
        {
            auto pts = iter->second;
            if (pts.empty())
                continue;
            all_edge_points.insert(all_edge_points.end(), pts.begin(), pts.end());
        }
        if (all_edge_points.size() > 0)
        {
            size_t sizes[3] = { 3, all_edge_points.size() };
            spurt::nrrd_utils::writeNrrd((void*)&all_edge_points[0], output_name + "_all_points" + suffix + ".nrrd", 
                                        nrrdTypeDouble, 2, sizes);
        }

        if (!all_triangles.empty())
        {
            size_t sizes[3] = {3, 3, all_triangles.size()};
            spurt::nrrd_utils::writeNrrd((void *)&all_triangles[0], output_name + "_mesh" + suffix + ".nrrd",
                                         nrrdTypeDouble, 3, sizes);
            std::vector<std::array<double, 2>> attributes(3 * all_triangles.size());
            std::cout << "interpolating value and ridge strength...\n";
            std::cout << "there are " << all_triangles.size() << " triangles\n";
            for (size_t n = 0; n < all_triangles.size(); ++n)
            {
                auto tri = all_triangles[n];
                for (uint k = 0; k < 3; ++k)
                {
                    auto p = tri[k];
                    // compute value and ridge strength at that position
                    auto v = values.value(to_point(p));
                    auto h = hessians.value(to_point(p));
                    auto l3_ev3 = evmin(h);
                    attributes[3 * n + k] = {v, l3_ev3.first};
                }
            }
            std::cout << "\ndone\n";
            sizes[0] = 2;
            spurt::nrrd_utils::writeNrrd((void *)&attributes[0], output_name + "_attributes" + suffix + ".nrrd",
                                         nrrdTypeDouble, 3, sizes);
            export_triangles(all_triangles, attributes, output_name + "_polydata" + suffix + ".vtp");
        }
        else
        {
            std::cout << "No triangles found\n";
        }

        // subdivide failed voxels
        std::vector<coord_type> new_voxels;
        for (int n=0; n<failed_voxels.size(); ++n)
        {
            coord_type id = failed_voxels[n];
            coord_type base = 2*id;
            new_voxels.push_back(base);
            new_voxels.push_back(base + coord_type(1, 0, 0));
            new_voxels.push_back(base + coord_type(1, 1, 0));
            new_voxels.push_back(base + coord_type(0, 1, 0));
            new_voxels.push_back(base + coord_type(0, 0, 1));
            new_voxels.push_back(base + coord_type(1, 0, 1));
            new_voxels.push_back(base + coord_type(1, 1, 1));
            new_voxels.push_back(base + coord_type(0, 1, 1));
        }
        voxel_ids.swap(new_voxels);
        all_triangles.clear();
        all_edge_points.clear();
        all_processed_edges.clear();
        failed_voxels.clear();
        skipped_voxels = 0;
        nedges = 0;
        nfound = 0;
    }
    if (!rejected.empty()) 
    {
        size_t sizes[2] = { 4, rejected.size() };
        spurt::nrrd_utils::writeNrrd((void*)&rejected[0], "rejected.nrrd", nrrdTypeDouble, 2, sizes);
    }
    else 
    {
        std::cout << "No rejected points\n";
    }

    return 0;
}
