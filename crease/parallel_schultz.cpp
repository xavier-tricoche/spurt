// STL
#include <cmath>
#include <array>
#include <exception>
#include <algorithm>
#include <utility>
#include <fstream>
#include <regex>
// Teem
#include <teem/nrrd.h>
// spurt
#include <math/types.hpp>
#include <data/image.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <spline/spline.h>

#include <crease/crease_mc_cases.hpp>

#include <vtk/vtk_utils.hpp>

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
std::atomic<size_t> tbb_progress_counter;
std::atomic<size_t> tbb_failed_to_triangulate;
std::atomic<size_t> tbb_none, tbb_skipped, tbb_succeed, tbb_processed, tbb_not_enough_points;
std::mutex add_triangles_mutex, add_point_mutex, update_progress_mutex;

using namespace spurt;

int verbose = 1;

constexpr int criterion_underflow = -10;
constexpr int weak_ridge = -30;
constexpr int low_ridge = -40;
constexpr int several_ridge_points = -50;
constexpr int no_ridge_points = -1;
constexpr int unknown_error = -100;
constexpr int one_ridge_points = 0;

constexpr int one_triangle = 1;
constexpr int valid_mc_case = 2;
constexpr int invalid_mc_case = 3;
constexpr int exotic_case = 4;
constexpr int not_enough_edges = 5;


typedef lvec3 coord_type;
typedef vec3 pos_type;
typedef bounding_box<pos_type> bounds_type;
typedef double scalar_type;
typedef long size_type;
typedef vec3 vector_type;
typedef mat3 matrix_type;
typedef std::array<pos_type, 3> triangle_type;

const scalar_type invalid_scalar = std::numeric_limits<scalar_type>::max();
const size_type invalid_size = std::numeric_limits<size_type>::max();

const pos_type invalid_pos = pos_type(invalid_scalar);
const coord_type invalid_coord = coord_type(-1);

// scalar, vector, and tensor field with C2 B-spline interpolation
typedef image<size_type, scalar_type, 3, scalar_type, 
              kernels::MitchellNetravaliBC, coord_type, pos_type> 
              scalar_image_type;
typedef image<size_type, scalar_type, 3, vector_type, 
              kernels::MitchellNetravaliBC, coord_type, pos_type> 
              vector_image_type;
typedef image<size_type, scalar_type, 3, matrix_type, 
              kernels::MitchellNetravaliBC, coord_type, pos_type> 
              matrix_image_type;
typedef scalar_image_type::grid_type grid_type;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> dyn_mat_type;

template<typename Value_>
using _image = image<long, double, 3, Value_, kernels::MitchellNetravaliBC, lvec3, vec3>;

template<typename Value_>
void import_nrrd(_image<Value_>& out, const std::string& filename)
{
    Nrrd* nin = nrrdNew();
    std::cout << "importing " << filename << '\n';
    if (nrrdLoad(nin, filename.c_str(), NULL))
    {
        char* err = biffGetDone(NRRD);
        std::cerr << "Thomas Schultz's ridge method: " << err << std::endl;
        exit(-1);
    }
    std::cout << "NRRD imported\n";
    typedef typename _image<Value_>::base_type raster_type;
    raster_type _raster = spurt::nrrd_utils::to_raster<size_type, scalar_type, 3, Value_>(nin, std::is_scalar<Value_>::value);
    std::cout << "Raster created\n";
    out = _image<Value_>(_raster);
    std::cout << "image created\n";
}

template<typename T>
T sign(const T& value) {
    if (value >= 0) return T(1);
    else return T(-1);
}

int distance(const coord_type& c0, const coord_type& c1) {
    return linf_norm(c1-c0);
}

struct EdgeID : public std::pair<coord_type, coord_type>
{
    typedef coord_type value_type;
    typedef std::pair<coord_type, coord_type> base_type;
    typedef spurt::lexicographical_order less_type;

    EdgeID(const coord_type& c0, const coord_type& c1)
        : base_type(c0, c1) {
        less_type less;
        if (less(base_type::second, base_type::first)) {
            std::swap(base_type::first, base_type::second);
        }
    }

    bool operator<(const EdgeID& other) const {
        less_type less;
        if (less(base_type::first, other.first)) return true;
        else if (less(other.first, base_type::first)) return false;
        else return less(base_type::second, other.second);
    }

    coord_type& operator[](int i) {
        if (i==0) {
            return base_type::first;
        }
        else if (i==1) {
            return base_type::second;
        }
        else {
            throw std::runtime_error("Invalid index for EdgeID");
        }
    }

    const coord_type& operator[](int i) const {
        if (i == 0)
        {
            return base_type::first;
        }
        else if (i == 1)
        {
            return base_type::second;
        }
        else
        {
            throw std::runtime_error("Invalid index for EdgeID");
        }
    }
};

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
        for (int i=0; i<v.size()-1; ++i) os << v[i] << ", ";
        os << v.back();
    }
    os << "]";
    return os;
}

template<typename T>
T invlinear(const T& f0, const T& f1, scalar_type umin=0, scalar_type umax=1) {
    // f = f0 + (u - umin) / (umax - umin) * (f1 - f0)
    // f = 0 <=> -f0 / (f1 - f0) = (u - umin) / (umax - umin)
    // f = 0 <=> u = -f0 / (f1 - f0) * (umax - umin) + umin    
    return umin - f0/(f1-f0)*(umax-umin);
}

template<typename T>
T linear(scalar_type u, const T& v0, const T& v1) {
    return (1.-u)*v0 + u*v1;
}

std::pair<scalar_type, vector_type> evmin(const matrix_type& H) {
    vector_type evals;
    matrix_type evecs;
    sym_eigensystem(evals, evecs, H);
    return std::make_pair(evals[2], vector_type(evecs.column(2)));
}

std::pair<vector_type, matrix_type> evall(const matrix_type& H) {
    vector_type evals;
    matrix_type evecs;
    sym_eigensystem(evals, evecs, H);
    return std::make_pair(evals, evecs);
}

scalar_type ridge_strength(const matrix_type& H) 
{
    vector_type evals;
    matrix_type evecs;
    sym_eigensystem(evals, evecs, H);
    return evals[2];
}

scalar_type delta23(const vector_type& evals) 
{
    return evals[1]-evals[2];
}

scalar_type l3(const vector_type& evals, double theta)
{
    scalar_type delta = delta23(evals);
    if (delta > theta) return 0;
    else return (1. - delta/theta)*(1. - delta/theta);
}

scalar_type dl3(const vector_type& evals, const vector_type& devals, 
                double theta)
{
    scalar_type _l3 = l3(evals, theta);
    if (_l3 == 0) return 0;
    else
    {
        return -2/theta * (devals[1] - devals[2])*(1 - delta23(evals)/theta);
    }
}

matrix_type tensor_T(const matrix_type& H, scalar_type theta=0.01)
{
    matrix_type Q;
    vector_type l;
    sym_eigensystem(l, Q, H);
    matrix_type lambda = matrix_type::identity();
    lambda(2,2) = l3(l, theta);
    return Q * lambda * transpose(Q);
}

// compute partial derivative of tensor T from partial derivative
// of Hessian
matrix_type tensor_dT(const matrix_type& dH, const matrix_type& H, 
                      scalar_type theta=0.01)
{
    /*
        cf "Matrix Differential Calculus with Applications in Statistics and Econometrics" by J.R. Magnus and H. Neudecker, Wiley, 2019
        dH_e = d\Lambda = diag(d\lambda_i's) (1)
        (1) follows from: 
            a) H = Q \Lambda Q^t
            b) dH = dQ \Lambda Q^t + Q d\Lambda Q^t + Q \Lambda dQ^t (from a)
            c) dH_e = Q^t dQ \Lambda Q^t Q + Q^t Q d\Lambda Q^t Q + Q^t Q \Lambda dQ^t Q (from b)
            d) dH_e = Q^t dQ \Lambda + d\Lambda + \Lambda dQ^t Q (from c)
            e) dH_e = (Q^t dQ + dQ^t Q) \Lambda + d\Lambda (from d)
            f) d(Q^tQ) = d(||Q||^2) = dQ^tQ + QdQ^t = 0
            g) dH_e = d\Lambda (from e and f)

        Similarly and by definition dT_e = diag(0, 0, dl_3)
    */
    matrix_type Q;
    vector_type l;
    sym_eigensystem(l, Q, H);
    scalar_type _l3 = l3(l, theta);
    if (_l3 == 0) return matrix_type(0); // dT_e = 0 => dT = 0
    else {
        matrix_type dH_e = transpose(Q) * dH * Q; // (c)
        vector_type dl(dH_e(0,0), dH_e(1,1), dH_e(2,2)); // (g)
        scalar_type _dl3 = dl3(l, dl, theta);
        matrix_type dT_e = 0;
        dT_e(2,2) = _dl3;
        return Q * dT_e * transpose(Q); // back to dataset coordinates
    }
}

vector_type vector_h(const vector_type& g, const matrix_type& H, 
                     scalar_type theta=0.01)
{
    matrix_type T = tensor_T(H, theta);
    return T*g - g;
}

vector_type vector_dh(const vector_type& g, const matrix_type& H, 
                      const vector_type& dg, const matrix_type& dH, 
                      scalar_type theta)
{
    matrix_type T = tensor_T(H, theta);
    matrix_type dT = tensor_dT(dH, H, theta);
    return dT*g + T*dg - dg;
}

vector_type ridge_normal(const pos_type& point,
                         const vector_type& g,  const matrix_type& H, 
                         const vector_image_type& gradient,
                         const matrix_image_type& hessian,
                         scalar_type theta=0.01)
{
    auto H_prime = hessian.derivative(point); // 3rd order tensor
    matrix_type h_prime;
    h_prime.column(0) = vector_dh(g, H, H.column(0), H_prime.layer(0), theta);
    h_prime.column(1) = vector_dh(g, H, H.column(1), H_prime.layer(1), theta);
    h_prime.column(2) = vector_dh(g, H, H.column(2), H_prime.layer(2), theta);

    // h' = outer(n,n);
    vector_type hevals;
    matrix_type hevecs;
    sym_eigensystem(hevals, hevecs, h_prime);
    vector_type hevals_plus = spurt::abs(hevals);
    return hevecs.column(spurt::argmax(hevals_plus));
}

std::pair<pos_type, scalar_type>
project(const vector_type& g, const matrix_type& h) {
    vector_type evals;
    matrix_type evecs;
    sym_eigensystem(evals, evecs, h);
    pos_type coords = spurt::abs(transpose(evecs)*g);
    return std::make_pair(coords, evals[2]);
}

template<typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T,N>& a) {
    os << "[";
    for (int i=0; i<N-1; ++i) {
        os << a[i] << ", ";
    }
    os << a[N-1] << "]";
    return os;
}

int case_number(const std::array<bool, 12>& edges) {
    std::cout << "edges=" << edges << '\n';
    std::array<int, 8> states({{0, 0, 0, 0, 0, 0, 0, 0}});
    // 0: not set
    // -1 or 1: positive or negative value
    states[0] = -1; // arbitrary initialization
    // The following order ensures that vertices are visited in order
    int edges_order[] = {0, 1, 2, 8, 9, 10, 11, 3, 4, 5, 6, 7};
    for (int e = 0; e < 12; ++e)
    {
        int edgeid = edges_order[e];
        auto v0 = spurt::marching_cubes::canonical_edge_indices[edgeid][0];
        auto v1 = spurt::marching_cubes::canonical_edge_indices[edgeid][1];
        if (v0 > v1)
        {
            std::swap(v0, v1);
        }
        if (e < 7)
        {
            if (edges[edgeid])
                states[v1] = -states[v0];
            else
                states[v1] = states[v0];
        }
        else
        {
            if ((edges[edgeid] && edges[v1] * edges[v0] > 0) ||
                (!edges[edgeid] && edges[v1] * edges[v0] < 0))
                return -1; // failed
        }
        std::cout << "e=" << e << ", edgeid=" << edgeid << ", states=" << states << '\n';
    }

    int cn = 0;
    for (int i = 0; i < 8; ++i)
    {
        if (states[i] > 0)
        {
            cn += 1 << i;
        }
    }
    return cn;

}

std::pair<int, scalar_type>
process_edge(const EdgeID &e,
             const scalar_image_type &data,
             const vector_image_type &gradient,
             const matrix_image_type &hessian,
             bool verbose = false,
             scalar_type theta = 0.01)
{
    vector_type h0 = vector_h(gradient(e[0]), hessian(e[0]));
    vector_type h1 = vector_h(gradient(e[1]), hessian(e[1]));
    if (inner(h0, h1) < 0) 
    {
        scalar_type n0 = h0.norm();
        scalar_type n1 = h1.norm();
        // model zero crossing as: (1-u) n0 - u n1 = 0 <=> u (n0 + n1) = n0
        scalar_type u = n0/(n1+n0);
        return std::pair(1, u);
    }
    else
    {
        return std::pair(0, 0.);
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
struct parser_traits<scalar_type> 
{
    static std::regex get()
    {
        return std::regex("([+-]? *[0-9]+(\\.[0-9]+)?)");
    }
    static int cast(const std::string &s)
    {
        return std::stod(s);
    }
};

template<typename T=int>
void parse_values(std::vector<T>& out, const std::string& str, size_t n) 
{
    std::regex myregex = parser_traits<T>::get();
    auto begin = std::sregex_iterator(str.begin(), str.end(), myregex);
    auto end = std::sregex_iterator();
    std::cout << "Found " << std::distance(begin, end) << " values\n";
    if (std::distance(begin, end) == n)
    {
        int i = 0;
        for (std::sregex_iterator iter = begin; iter != end; ++iter, ++i)
        {
            std::smatch match = *iter;
            std::cout << "value=" << match.str() << std::endl;
            out.push_back(parser_traits<T>::cast(match.str()));
        }
    }
    else 
    {
        throw std::runtime_error("invalid input");
    }
}

void edge_neighbors(std::vector<coord_type>& neighbors, const EdgeID& eid) 
{
    coord_type i0 = eid[0];
    coord_type i1 = eid[1];
    int dim=0;
    for (int i=0; i<3; ++i) 
    {
        if (i0[i] != i1[i]) 
        {
            dim = i;
            break;
        }
    }
    neighbors.resize(4);
    int low = std::min(i0[dim], i1[dim]);
    coord_type ref = eid[0];
    ref[dim] = low;
    std::fill(neighbors.begin(), neighbors.end(), ref);
    neighbors[1][(dim+1)%3]--;
    neighbors[2][(dim+2)%3]--;
    neighbors[3][(dim+1)%3]--;
    neighbors[3][(dim+2)%3]--;
}

void grow_cluster(std::set<int>& cluster, int start, const dyn_mat_type& dist) 
{
    if (cluster.find(start) == cluster.end())
    {
        cluster.insert(start);
        for (int i=start+1; i<dist.cols(); ++i) {
            if (dist(start, i) == 1)
            {
                if (cluster.find(i) == cluster.end())
                {
                    grow_cluster(cluster, i, dist);
                }
            }
        }
    }
}

void find_neighbors(std::vector<std::vector<coord_type>>& neighbors, 
                    const std::map<coord_type, std::vector<int>, spurt::lexicographical_order> voxels) 
{
    std::vector<coord_type> all_voxels;
    std::cout << "creating an array of voxels\n";
    for (auto iter=voxels.begin(); iter!=voxels.end(); ++iter)
    {
        all_voxels.push_back(iter->first);
    }
    int nvoxels = all_voxels.size();
    std::cout<< "done. creating a distance matrix\n";
    dyn_mat_type dist = dyn_mat_type::Zero(nvoxels, nvoxels);
    for (int i=0; i<nvoxels-1; ++i)
    {
        for (int j=i+1; j<nvoxels; ++j)
        {
            dist(i,j) = dist(j,i) = distance(all_voxels[i], all_voxels[j]);
        }
    }
    std::cout << "done.\n";
    std::vector<int> cluster_id(nvoxels, -1);
    for (int i=0; i<nvoxels; ++i)
    {
        if (cluster_id[i] >= 0) continue;
        std::set<int> acluster;
        acluster.insert(i);
        for (int j=i+1; j<nvoxels; ++j)
        {
            if (dist(i,j) == 1)
            {
                grow_cluster(acluster, j, dist);
            }
        }
        std::vector<coord_type> ids;
        std::for_each(acluster.begin(), acluster.end(), [&](int n) {
            ids.push_back(all_voxels[n]);
        });
        neighbors.push_back(ids);
        std::for_each(acluster.begin(), acluster.end(), [&](int n) {
            cluster_id[n] = neighbors.size()-1;
        });
    }
}

std::pair<scalar_type, scalar_type> 
evaluate(const pos_type& point, const scalar_image_type& values, 
         const matrix_image_type& hessian) 
{
    scalar_type f = values.value(point);
    matrix_type H = hessian.value(point);
    return std::make_pair(f, ridge_strength(H));
}

int triangulate(std::vector<triangle_type>& out, 
                std::map<int, pos_type> &edges, 
                std::ostream& os)
{
    os << "triangles: edges contains " << edges.size() << " edges with ridge points\n";
    out.clear();

    // if we only have 3 or less points, things are easy
    if (edges.size() == 3) 
    {
        triangle_type T;
        int i=0;
        for (auto iter=edges.begin(); iter!=edges.end(); ++iter, ++i) 
        {
            T[i] = iter->second;
        }
        out.push_back(T);
        os << "3 points on 3 edges: " << T << ": success!\n";
        return one_triangle;
    }
    else if (edges.size() < 3)
    {
        os << "We have only 2 (or less) edges in input of triangulation.\n"
           << "Giving up (Case N<3)\n";
        return not_enough_edges;
    }
    else 
    {
        // Calculate edge case number
        int edge_case = 0;
        for (auto iter=edges.begin(); iter!=edges.end(); ++iter) 
        {
            edge_case += 1 << iter->first;
        }
        os << "edge_case = " << edge_case << '\n';
        int triangle_case = spurt::marching_cubes::edge_code_to_case_id[edge_case];
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
}

void select_voxels(std::vector<coord_type>& voxels, const bounds_type& bounds, 
                   const coord_type& shape)
{
    if (spurt::all(bounds.min() < bounds.max()))
    {
        voxels.clear();
        const pos_type& l = bounds.min();
        const pos_type& h = bounds.max();
        coord_type low = floor(l);
        coord_type high = ceil(h);
        for (int k=low[2]; k<=high[2]; ++k) 
        {
            for (int j=low[1]; j<=high[1]; ++j)
            {
                for (int i=low[0]; i<=high[0]; ++i) {
                    voxels.push_back(coord_type(i,j,k));
                }
            }
        }
    }
    else 
    {
        voxels.resize((shape[0]-1)*(shape[1]-1)*(shape[2]-1));
        size_type n = 0;
        for (int k = 0; k < shape[2] - 1; ++k)
        {
            for (int j = 0; j < shape[1] - 1; ++j)
            {
                for (int i = 0; i < shape[0] - 1; ++i, ++n)
                {
                    voxels[n] = coord_type(i, j, k);
                }
            }
        }
    }
}

void visualize(const std::vector<pos_type>& points, 
               const std::vector<triangle_type>& triangles)
{
    vtkPolyData* pd = vtk_utils::make_points(points);
    vtkPolyData* glyphs = vtk_utils::make_spheres(pd, 0.1);
    VTK_CREATE(vtkPolyDataMapper, pts_mapper);
    pts_mapper->SetInputData(glyphs);
    VTK_CREATE(vtkActor, pts_actor);
    pts_actor->SetMapper(pts_mapper);

    std::vector<vec3> tri_pts;
    VTK_CREATE(vtkCellArray, cells);
    for (auto it=triangles.begin(); it!=triangles.end() ; ++it) {
        tri_pts.push_back((*it)[0]);
        tri_pts.push_back((*it)[1]);
        tri_pts.push_back((*it)[2]);
        cells->InsertNextCell(3);
        cells->InsertCellPoint(tri_pts.size()-3);
        cells->InsertCellPoint(tri_pts.size()-2);
        cells->InsertCellPoint(tri_pts.size()-1);
    }
    vtkPolyData* pdtris = vtk_utils::make_points(tri_pts);
    pdtris->SetPolys(cells);
    VTK_CREATE(vtkPolyDataMapper, tri_mapper);
    tri_mapper->SetInputData(pdtris);
    VTK_CREATE(vtkActor, tri_actor);
    tri_actor->SetMapper(tri_mapper);

    VTK_CREATE(vtkRenderer, renderer);
    renderer->AddActor(pts_actor);
    renderer->AddActor(tri_actor);
    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
    window->SetSize(1024, 768);
    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    interactor->Initialize();
    window->Render();
    interactor->Start();
}

int main(int argc, const char* argv[]) 
{
    std::string scalar_name, gradient_name, hessian_name, output_name;
    scalar_type minval, minstr, eps, mind;
    int res, niter;
    int verbose;
    coord_type voxel_id = invalid_coord;
    bounds_type bounds(pos_type(0.), pos_type(-1.));
    spurt::vec4 dv(0.1, 0.2, 0.3, 0.4);
    
    namespace cl = spurt::command_line;
    cl::option_traits
        required_group(true, false, "Required Options"),
        positional_group(true, true, "Positional Group"),
        optional_group(false, false, "Optional Group");

    cl::option_parser parser(argv[0],
                             "Extract ridge surfaces from scalar volume using\n Thomas Schultz's 2009 method");
    try 
    {
        parser.use_short_symbols(false);
        parser.use_brackets(true);

        parser.add_value("value", scalar_name, "Scalar volume", required_group);
        parser.add_value("gradient", gradient_name, "Gradient volume", required_group);
        parser.add_value("hessian", hessian_name, "Hessian volume", required_group);
        parser.add_value("output", output_name, "Output filename", required_group);
        parser.add_value("minval", minval, 0., "Min scalar value", optional_group);
        parser.add_value("minstr", minstr, 0., "Min ridge strength (<0)", optional_group);
        // parser.add_value("eps", eps, 1.0e-9, "Numerical precision", optional_group);
        parser.add_value("verbose", verbose, 0, "Verbose output", optional_group);
        parser.add_tuple<3>("voxel", voxel_id, voxel_id, "Coordinates of single voxel to process");
        parser.add_tuple<3>("blower", bounds.min(), bounds.min(), "Lower bounds of domain to consider");
        parser.add_tuple<3>("bupper", bounds.max(), bounds.max(), "Upper bounds of domain to consider");
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

    scalar_image_type values;
    vector_image_type gradient;
    matrix_image_type hessian;
    import_nrrd(values, scalar_name);
    import_nrrd(gradient, gradient_name);
    import_nrrd(hessian, hessian_name);

    auto shape = values.grid().resolution();
    std::vector< coord_type > voxels;
    if (spurt::any(voxel_id != invalid_coord))
    {
        voxels.clear();
        voxels.push_back(voxel_id);
        verbose = 2;
    }
    else
    {
        std::cout << "selected voxels..." << std::flush;
        select_voxels(voxels, bounds, shape);
        std::cout << " done\n";
    }

    if (verbose) std::cout << "There are " << voxels.size() << " voxels in input for a total of " << 12*voxels.size() << " (redundant) edges\n";

    tbb_skipped = 0;
    int nnot_enough_points = 0;
    tbb_succeed = 0;
    int nprocessed = 0;
    int nexcluded = 0;
    int neven = 0;
    int nweak = 0;
    int nlow = 0;
    int nunderflow = 0;
    int nfiltered = 0;
    // keep track of what triangles stem from what voxels to remove them 
    // if needed in case their containing voxels are subdivided. 
    std::map<coord_type, std::vector<int>, spurt::lexicographical_order> voxel_to_triangles;
    std::vector<triangle_type> all_triangles;
    std::vector<vec4> rejected;
    std::vector<vec3> all_points;

    spurt::ProgressDisplay progress;

    struct broken_voxel {
        coord_type id;
        std::array<scalar_type, 8> values;
        std::array<scalar_type, 8> strengths;
        std::array<vector_type, 8> gradients;
        std::array<matrix_type, 8> hessians;
        std::map<int, pos_type>  edge_points;
    };

    std::vector<broken_voxel> broken_voxels;
    std::vector<EdgeID> double_edges;
    std::vector<coord_type> to_subdivide;
    std::vector<int> voxel_to_edge;
    std::map<coord_type, std::vector<int>, spurt::lexicographical_order> voxel_counter;

    srand48(130819751900);
    progress.begin(voxels.size(), "Extract ridges", 500, "tris: 0, done: 0, ok: 0, skip: 0, underflow: 0, weak: 0, low: 0, none: 0, even: 0, failed: 0");
    tbb_progress_counter = 0;
    tbb_failed_to_triangulate = 0;
    tbb::parallel_for(tbb::blocked_range<int>(0,voxels.size()),
                       [&](tbb::blocked_range<int> r) {
    for (int n=r.begin(); n!=r.end(); ++n)
    {
        coord_type voxel_id = voxels[n];

        // check if current voxel satisfies threshold requirements
        std::array<scalar_type, 8> v_values;
        std::array<scalar_type, 8> s_values;
        for (int vertex_id=0; vertex_id<8; ++vertex_id)
        {
            auto shift = spurt::marching_cubes::vertices[vertex_id];
            coord_type v = voxel_id + coord_type(shift[0], shift[1], shift[2]);
            v_values[vertex_id] = values(v[0], v[1], v[2]);
            matrix_type h = hessian(v[0], v[1], v[2]);
            s_values[vertex_id] = evmin(h).first;
        }
        if ((*std::min_element(s_values.begin(), s_values.end()) > minstr) ||
            (*std::max_element(v_values.begin(), v_values.end()) < minval))
        {
            tbb_none++;
            continue;
        }

        std::map<int, pos_type> found;
        std::string update_str = "tris: " + std::to_string(all_triangles.size()) + 
            ", ok: " + std::to_string(tbb_succeed) + 
            ", skip: " + std::to_string(tbb_skipped) + 
            ", none: " + std::to_string(tbb_none) +
            ", failed: " + std::to_string(tbb_failed_to_triangulate) +
            ", insuf:" + std::to_string(tbb_not_enough_points);
        std::ostringstream log;
        for (int edge_id = 0; edge_id < 12; ++edge_id) 
        {
            coord_type v0 = voxel_id + coord_type(ivec3(spurt::marching_cubes::canonical_edge_coordinates[edge_id][0]));
            coord_type v1 = voxel_id + coord_type(ivec3(spurt::marching_cubes::canonical_edge_coordinates[edge_id][1]));
            EdgeID global_edge_id(v0, v1);
            {
                auto result = process_edge(global_edge_id, values, gradient, hessian, verbose, 0.01);
                if (result.first == 1)
                {       
                    scalar_type u = result.second;
                    found[edge_id] = (1. - u) * global_edge_id[0] + u * global_edge_id[1];
                    std::lock_guard<std::mutex> guard(add_point_mutex);
                    all_points.push_back(found[edge_id]);
                }
                tbb_processed++;
            }
        }
        {
            std::lock_guard<std::mutex> guard(update_progress_mutex);
            progress.update(tbb_progress_counter, update_str);
        }
        ++tbb_progress_counter;

        if (found.size() >= 3) 
        {
            std::vector<triangle_type> tris;
            int tri_case = triangulate(tris, found, log);

            if (tri_case == one_triangle || tri_case == valid_mc_case) 
            {
                std::for_each(tris.begin(), tris.end(), [&](auto T) 
                { 
                    std::lock_guard<std::mutex> guard(add_triangles_mutex);
                    all_triangles.push_back(T);
                    // std::cout << "\nTriangle is " << T << '\n';
                });
                tbb_succeed++;
            }
            else 
            {
                tbb_failed_to_triangulate++;
            }
        }
    }
    });

    progress.end();

    if (!all_triangles.empty()) 
    {
        size_t sizes[3] = { 3, 3, all_triangles.size() };
        spurt::nrrd_utils::writeNrrd((void*)&all_triangles[0], output_name + "_mesh.nrrd", 
                                     nrrdTypeDouble, 3, sizes);
        std::vector<std::array<scalar_type, 2>> attributes(3*all_triangles.size());
        std::cout << "interpolating value and ridge strength...\n";
        std::cout << "there are " << all_triangles.size() << " triangles\n";
        for (size_t n=0; n<all_triangles.size(); ++n)
        {
            auto tri = all_triangles[n];
            for (int k=0; k<3; ++k)
            {
                auto p = tri[k];
                // compute value and ridge strength at that position
                auto v = values.value(p);
                auto h = hessian.value(p);
                auto l3_ev3 = evmin(h);
                attributes[3*n+k] = {v, l3_ev3.first};
            }
        }
        std::cout << "\ndone\n";
        sizes[0] = 2;
        spurt::nrrd_utils::writeNrrd((void*)&attributes[0], output_name + "_attributes.nrrd",
                                     nrrdTypeDouble, 3, sizes);
        visualize(all_points, all_triangles);
    }
    else 
    {
        std::cout << "No triangles found\n";
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
