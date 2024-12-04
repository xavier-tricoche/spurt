// STL
#include <cmath>
#include <array>
#include <exception>
#include <algorithm>
#include <utility>
#include <fstream>
#include <regex>
// Boost
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>
#include <boost/multi_array.hpp>
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

using namespace spurt;

// idea: store all ridge points, weak or even numbered with a
// corresponding code to indicate why they were not included for
// triangulation and upon running into a pathological case
// figure out if some of these points should have been used anyway,
// and, if so, how to make that determination ahead of time, if at
// all possible.

// Edges containing multiple ridge points need to be refined, yielding 2 edges that are shared by 4 
// regular voxels or 8 refined voxels. 
// Algorithm:
// - If an edge is found to have 2 or more rige points:
//    all its surrounding voxels are marked for subdivision.
// - Afterwards: 
//    all unique subdivision voxels are processed and subdivided as needed

int verbose = 1;

constexpr int criterion_underflow = -10;
constexpr int determinant_underflow = -20;
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


// scalar, vector, and tensor field with C2 B-spline interpolation
typedef image<long, double, 3, double, kernels::MitchellNetravaliBC, lvec3, vec3> scalar_image_type;
typedef image<long, double, 3,   vec3, kernels::MitchellNetravaliBC, lvec3, vec3> vector_image_type;
typedef image<long, double, 3,   mat3, kernels::MitchellNetravaliBC, lvec3, vec3> matrix_image_type;
typedef scalar_image_type::grid_type grid_type;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> dyn_mat_type;

template<typename Value_>
using _image = image<long, double, 3, Value_, kernels::MitchellNetravaliBC, lvec3, vec3>;

typedef lvec3 CoordID;

template<typename Value_>
void import_nrrd(_image<Value_>& out, const std::string& filename)
{
    Nrrd* nin = nrrdNew();

    if (nrrdLoad(nin, filename.c_str(), NULL))
    {
        char* err = biffGetDone(NRRD);
        std::cerr << "Thomas Schultz's ridge method: " << err << std::endl;
        exit(-1);
    }
    
    typedef typename _image<Value_>::base_type raster_type;
    raster_type _raster;
    spurt::nrrd_utils::to_raster(_raster, nin, std::is_scalar<Value_>::value);
    out = _image<Value_>(_raster);
}

template<typename T>
T sign(const T& value) {
    if (value >= 0) return T(1);
    else return T(-1);
}

int distance(const CoordID& c0, const CoordID& c1) {
    return linf_norm(c1-c0);
}

struct EdgeID : public std::pair<CoordId, CoordId>
{
    typedef CoordId value_type;
    typedef std::pair<CoordId, CoordId> base_type;
    typedef spurt::lexicographical_order less_type;

    EdgeID(const CoordID& c0, const CoordID& c1)
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

    CoordID& operator[](int i) {
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

    const CoordID& operator[](int i) const {
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
    vec3 evals;
    mat3 evecs;
    sym_eigensystem(evals, evecs, H);
    return std::make_pair(evals[2], vec3(evecs.column(2)));
}

std::pair<vec3, mat3> evall(const mat3& H) {
    vec3 evals;
    mat3 evecs;
    sym_eigensystem(evals, evecs, H);
    return std::make_pair(evals, evecs);
}

mat3 schultz_tensorT(const mat3& H, double theta=0.01)
{
    /*
    
    A = Q^t L Q

    A Q = L Q

    aij qjk = lk qik 
    
    ELL_3M_TRANSPOSE(evecTrans,evecs);

    (\nabla H)_E = Q (\nabla H) Q^t
    diff = lmed - lmin
    if (diff < theta):
        l3 = (1-diff/theta)^2
    else:
        l3 = 0

    (\nabla H)_E (xz) = (\nabla H)_E (zx) = (\nabla H)_E (xz) * (1-l3)/(lmax-lmin)
    (\nabla H)_E (yz) = (\nabla H)_E (zy) = (\nabla H)_E (yz) * (1-l3)/(lmed-lmin)

    if (diff < theta):
        \nabla l3 = 2/theta * (1-diff/theta) * ((\nabla H)_E (zz) - (\nabla H)_E (yy))
    else:
        \nabla l3 = 0

    (\nabla H)_E (zz) = \nabla l3
    (\nabla H)_E (xx) = (\nabla H)_E (xy) = (\nabla H)_E (yx) = (\nabla H)_E (yy) = 0

  if (ridge) {
    double diff=evals[1]-evals[2];
    double l3;
    double l3der;
    if (diff<evalDiffThresh)
      l3=(1.0-diff/evalDiffThresh)*(1.0-diff/evalDiffThresh);
    else l3=0.0;
    hessderE[2]*=(1.0-l3)/(evals[0]-evals[2]);
    hessderE[6]=hessderE[2];
    hessderE[5]*=(1.0-l3)/(evals[1]-evals[2]);
    hessderE[7]=hessderE[5];

    if (diff<evalDiffThresh)
      l3der = 2/evalDiffThresh*(1-diff/evalDiffThresh)*
        (hessderE[8]-hessderE[4]);
    else
      l3der=0;
    hessderE[8]=l3der;
    hessderE[0]=hessderE[1]=hessderE[3]=hessderE[4]=0.0;



    
    */

    mat3 Q;
    vec3 l;
    sym_eigensystem(l, Q, H);
    mat3 lambda = mat3::identity();
    lambda(2,2) = 0;
    double dlambda = l[1] - l[2];
    if (dlambda < theta)
    {
        lambda(2,2) = (1. - dlambda/theta)*(1. - dlambda/theta);
    }
    return Q * lambda * transpose(Q);
}

vec3 schultz_vectorh(const mat3& H, const vec3& g, double theta=0.01)
{
    mat3 T = schultz_tensorT(H, theta);
    return T*g - g;
}

vec3 ridge_normal(const vec3& point, const mat3& H, const vec3& g, 
                  const matrix_image& hessian)
{

}


std::pair<vec3, double> project(const vec3& g, const mat3& h) {
    vec3 evals;
    mat3 evecs;
    sym_eigensystem(evals, evecs, h);
    vec3 coords = spurt::abs(transpose(evecs)*g);
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

std::pair<int, double>
analyze_edge_schultz(const EdgeID &e,
                     const vector_image_type &gradient,
                     const matrix_image_type &hessian,
                     const scalar_image_type &data,
                     bool verbose = false,
                     double theta = 0.01)
{
    vec3 h0 = schultz_vectorh(hessian(e[0]), gradient(e[0]));
    vec3 h1 = schultz_vectorh(hessian(e[1]), gradient(e[1]));
    if (inner(h0, h1) < 0) 
    {
        double n0 = h0.norm();
        double n1 = h1.norm();
        // model zero crossing as: (1-u) n0 - u n1 = 0 <=> u (n0 + n1) = n0
        double u = n0/(n1+n0);
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
struct parser_traits<double> 
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

void edge_neighbors(std::vector<CoordID>& neighbors, const EdgeID& eid) 
{
    CoordID i0 = eid[0];
    CoordID i1 = eid[1];
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
    CoordID ref = eid[0];
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

void find_neighbors(std::vector<std::vector<CoordID>>& neighbors, 
                    const std::map<CoordID, std::vector<int>, spurt::lexicographical_order> voxels) 
{
    std::vector<CoordID> all_voxels;
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
        std::vector<CoordID> ids;
        std::for_each(acluster.begin(), acluster.end(), [&](int n) {
            ids.push_back(all_voxels[n]);
        });
        neighbors.push_back(ids);
        std::for_each(acluster.begin(), acluster.end(), [&](int n) {
            cluster_id[n] = neighbors.size()-1;
        });
    }
}

int main(int argc, const char* argv[]) 
{
    std::string data_name, gradient_name, hessian_name, output_name;
    double minval, minstr, eps, mind;
    int res, niter;
    bool verbose;
    CoordID voxel_id;
    std::pair<vec3, vec3> bounds;
    std::string voxel_str = "(-1, -1, -1)";
    std::string bounds_str = "(0, -1; 0, -1; 0, -1)";
    

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

        parser.add_value("value", data_name, "Scalar volume", required_group);
        parser.add_value("gradient", gradient_name, "Gradient volume", required_group);
        parser.add_value("hessian", hessian_name, "Hessian volume", required_group);
        parser.add_value("output", output_name, "Output filename", required_group);
        parser.add_value("minval", minval, 0., "Min scalar value", optional_group);
        parser.add_value("minstr", minstr, 0., "Min ridge strength", optional_group);
        parser.add_value("eps", eps, 1.0e-9, "Numerical precision", optional_group);
        parser.add_value("verbose", verbose, false, "Verbose output", optional_group);
        parser.add_value("voxel", voxel_str, voxel_str, "Coordinates of single voxel to process");
        parser.add_value("bounds", bounds_str, bounds_str, "Bounds of domain to consider");
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

    {
        std::vector<int> iv;
        parse_values<int>(iv, voxel_str, 3);
        std::cout << iv << '\n';
        voxel_id = { iv[0], iv[1], iv[2] };

        std::vector<double> dv;
        parse_values<double>(dv, bounds_str, 6);
        std::cout << dv << '\n';
        bounds.first = { dv[0], dv[2], dv[4] };
        bounds.second = { dv[1], dv[3], dv[5] };
    }

    scalar_image_type values;
    vector_image_type gradient;
    matrix_image_type hessian;
    import_nrrd(values, data_name);
    import_nrrd(gradient, gradient_name);
    import_nrrd(hessian, hessian_name);

    auto shape = values.grid().resolution();
    std::vector< CoordID > voxels;
    if (voxel_id[0] != -1)
    {
        voxels.clear();
        voxels.push_back(voxel_id);
        verbose = true;
    }
    else if (bounds.first[0] < bounds.second[0])
    {
        voxels.clear();
        const vec3& l = bounds.first;
        const vec3& h = bounds.second;
        CoordID low = floor(l);
        CoordID high = ceil(h);
        for (int k=low[2]; k<=high[2]; ++k) 
        {
            for (int j=low[1]; j<=high[1]; ++j)
            {
                for (int i=low[0]; i<=high[0]; ++i) {
                    voxels.push_back(CoordID(i,j,k));
                }
            }
        }
    }
    else 
    {
        int _n = 0;
        voxels.resize((shape[0]-1)*(shape[1]-1)*(shape[2]-1));
        for (int k = 0; k < shape[2] - 1; ++k)
        {
            for (int j = 0; j < shape[1] - 1; ++j)
            {
                for (int i = 0; i < shape[0] - 1; ++i, ++_n)
                {
                    voxels[_n] = CoordID({{i, j, k}});
                }
            }
        }
    }
    std::cout << "There are " << voxels.size() << " voxels in input for a total of " << 12*voxels.size() << " (redundant) edges\n";

    std::map<EdgeID, std::vector<vec3> > all_processed_edges;
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
    // keep track of what triangles stem from what voxels to remove them 
    // if neededin case their containing voxels are subdivided. 
    std::map<CoordID, std::vector<int>, spurt::lexicographical_order> voxel_to_triangles;
    std::vector<triangle_type> all_triangles;
    std::vector<vec4> rejected;

    if (voxel_id[0] != -1)
    {
        std::vector<mat3> H(8);
        std::vector<vec3> g(8);
        for (int i=0; i<8; ++i)
        {
            auto shift = spurt::marching_cubes::vertices[i];
            CoordID pid = voxel_id + CoordID(shift[0], shift[1], shift[2]);
            H[i] = hessian(pid[0], pid[1], pid[2]);
            g[i] = gradient(pid[0], pid[1], pid[2]);
        }
        std::vector<double> det_values(2*101*101*101);
        for (int k=0; k<=100; ++k)
        {
            double z = k*0.01;
            double Z = 1.-z;
            for (int j=0; j<=100; ++j)
            {
                double y = j*0.01;
                double Y = 1.-y;
                for (int i=0; i<=100; ++i)
                {
                    double x = i*0.01;
                    double X = 1.-x;
                    mat3 theH = X * Y * Z * H[0] + x * Y * Z * H[1] + x * y * Z * H[2] + X * y * Z * H[3] +
                                X * Y * z * H[4] + x * Y * z * H[5] + x * y * z * H[6] + X * y * z * H[7];
                    vec3 theg = X * Y * Z * g[0] + x * Y * Z * g[1] + x * y * Z * g[2] + X * y * Z * g[3] +
                                X * Y * z * g[4] + x * Y * z * g[5] + x * y * z * g[6] + X * y * z * g[7];
                    det_values[2 * (i + 101 * (j + 101 * k))    ] = determinant(theg, theH);
                    det_values[2 * (i + 101 * (j + 101 * k)) + 1] = (project(theg, theH).first)[0];
                }
            }
        }
        size_t sizes[4] = {2, 101, 101, 101};
        std::ostringstream os;
        os << output_name << "_voxel_" << voxel_id[0] << "_" << voxel_id[1] << "_" << voxel_id[2] << ".nrrd";
        spurt::nrrd_utils::writeNrrd((void *)&det_values[0], os.str(), nrrdTypeDouble, 4, sizes);
    }

    spurt::ProgressDisplay progress;

    struct broken_voxel {
        CoordID id;
        std::vector<double> values;
        std::vector<double> strengths;
        std::vector<vec3> gradients;
        std::vector<mat3> hessians;
        std::map<int, std::vector<vec3> > edge_points;
    };
    std::vector<broken_voxel> broken_voxels;
    std::vector<EdgeID> double_edges;
    std::vector<CoordID> to_subdivide;
    std::vector<int> voxel_to_edge;
    std::map<CoordID, std::vector<int>, spurt::lexicographical_order> voxel_counter;

    srand48(130819751900);
    progress.begin(voxels.size(), "Extract ridges", 10000, "tris: 0, done: 0, ok: 0, skip: 0, underflow: 0, weak: 0, low: 0, none: 0, even: 0, failed: 0");
    for (int n=0; n<voxels.size(); ++n) 
    {
        CoordID id = voxels[n];
        // verbose = (id == CoordID(6, 14, 13));
        // verbose = (id[0] >= 150 && id[0] <= 160 && id[1] >= 60 && id[1] <= 75 && id[2] >= 60 && id[2] <=65);
        if (verbose)
        {
            std::cout << "processing voxel: " << id << '\n';
        }

        // check if current voxel satisfies threshold requirements
        std::vector<double> v_values(8);
        std::vector<double> s_values(8);
        for (int i=0; i<8; ++i)
        {
            auto shift = spurt::marching_cubes::vertices[i];
            CoordID v = id + CoordID(shift[0], shift[1], shift[2]);
            v_values[i] = values(id[0], id[1], id[2]);
            auto h = hessian(id[0], id[1], id[2]);
            s_values[i] = evmin(h).first;
        }
        if ((*std::min_element(&s_values[0], &s_values[8]) > minstr) ||
            (*std::max_element(&v_values[0], &v_values[8]) < minval))
        {
            nnone++;
            continue;
        }

        // verbose = (id[0] >= 8 && id[0]<=10 && id[1]>=14 && id[1]<=19 && id[2]>=14 && id[2]<=17);
        std::map<int, std::vector<vec3>> found;
        std::string update_str = "tris: " + std::to_string(all_triangles.size()) + 
            ", ok: " + std::to_string(nsucceeded) + 
            ", skip: " + std::to_string(nskipped) + 
            // ", weak: " + std::to_string(nweak) +
            // ", low: " + std::to_string(nlow) +
            ", none: " + std::to_string(nnone) +
            ", failed: " + std::to_string(nfailed_to_triangulate) +
            ", insuf:" + std::to_string(nnot_enough_points);
        progress.update(n, update_str);
        std::ostringstream log;
        for (int i = 0; i < 12; ++i) 
        {
            CoordID v0 = id + CoordID(ivec3(spurt::marching_cubes::canonical_edge_coordinates[i][0]));
            CoordID v1 = id + CoordID(ivec3(spurt::marching_cubes::canonical_edge_coordinates[i][1]));
            EdgeID edgeid(v0, v1);
            auto iter = all_processed_edges.find(edgeid);
            if (iter != all_processed_edges.end())
            {   
                if (verbose) {
                    std::cout << "entry found for edge #" << i << " between " << v0 << " and " << v1 << '\n';
                    std::cout << "edge #" << i << " of voxel " << id << ": " << edgeid << " has already been processed\n";
                    std::cout << "map contains " << all_processed_edges.size() << " elements\n";
                }
                nskipped += 1;
                auto solutions = iter->second;
                if (!solutions.empty()) 
                {
                    found[i] = solutions;
                    if (verbose) {
                        std::cout << "This entry contains one or several ridge points\n";
                        std::cout << "Found now contains " << found.size() << " entries\n";
                    }
                }
                else {
                    if (verbose) {
                        std::cout << "This entry did not contain a ridge point\n";
                        std::cout << "Found now contains " << found.size() << " entries\n";
                    } 
                    continue;
                }
            }
            else 
            {
                if (verbose) 
                    std::cout << "\n\nedge #" << i << " of voxel " << id << ": " << edgeid << " is being processed" << '\n';
                auto result = analyze_edge_schultz(edgeid, gradient, hessian, values, verbose, 0.01);
                if (result.first == 1)
                {                
                    found[i] = std::vector<vec3>();
                    if (verbose) std::cout << "process_edge returned 1 solution\n";
                    double u = result.second;
                    found[i].push_back((1. - u) * edgeid[0] + u * edgeid[1]);
                    if (verbose) std::cout << "Found now contains " << found.size() << " entries\n";
                    all_processed_edges[edgeid] = found[i];
                }
                else
                {
                    if (verbose)
                        std::cout << "no ridge point\n";
                    nnone++;
                    all_processed_edges[edgeid] = std::vector<vec3>();
                }
                nprocessed++;
            }
        }

        if (found.size() >= 3) 
        {
            if (verbose) 
            {
                std::cout << "After looking at all edges found contains " << found.size() << " entries\n";
                for (auto iter = found.begin(); iter!=found.end(); ++iter)
                {
                    std::cout << "Edge #" << iter->first << " contains " << iter->second << '\n';
                    for (auto jter=iter->second.begin(); jter!=iter->second.end(); ++jter) {
                        auto r = evaluate(*jter, values, hessian);
                        std::cout << "ridge point at " << *jter << " has value " << r.first << " and ridge strength " << r.second << '\n';
                    }
                }
            }
            std::vector<triangle_type> tris;
            int tri_case;
            if (verbose)
                tri_case = triangulate(tris, found, std::cout);
            else
                tri_case = triangulate(tris, found, log);

            if (tris.size() > 0) 
            {
                std::for_each(tris.begin(), tris.end(), [&](auto T) 
                {
                    all_triangles.push_back(T);
                    // std::cout << "\nTriangle is " << T << '\n';
                });
                nsucceeded++;
            }
            else 
            {
                if (verbose)
                    std::cout << "Triangulation failed\n";
                nfailed_to_triangulate++;
                // store rejected points with triangulation diagnostics
                for (auto iter=found.begin(); iter!=found.end(); ++iter) 
                {
                    for (int l=0; l<iter->second.size(); ++l)
                    {
                        const vec3& p = iter->second[l];
                        vec4 q({p[0], p[1], p[2], tri_case});
                        rejected.push_back(q);
                    }
                }
                broken_voxel voxel;
                voxel.id = id;
                voxel.edge_points = found;
                for (int vid=0; vid<8; ++vid) {
                    CoordID pid = id + CoordID(spurt::marching_cubes::vertices[vid]);
                    voxel.values.push_back(values(pid[0], pid[1], pid[2]));
                    voxel.gradients.push_back(gradient(pid[0], pid[1], pid[2]));
                    voxel.hessians.push_back(hessian(pid[0], pid[1], pid[2]));
                    voxel.strengths.push_back(evmin(voxel.hessians.back()).first);
                }
                broken_voxels.push_back(voxel);
            }
        }
        else if (found.size() > 0) 
        {
            nnot_enough_points++;
            for (auto iter = found.begin(); iter != found.end(); ++iter)
            {
                for (int l=0; l<iter->second.size(); ++l)
                {
                    const vec3& p = iter->second[l];
                    vec4 q({p[0], p[1], p[2], not_enough_edges});
                    rejected.push_back(q);
                }
            }

            broken_voxel voxel;
            voxel.id = id;
            voxel.edge_points = found;
            for (int vid = 0; vid < 8; ++vid)
            {
                CoordID pid = id + CoordID(ivec3(spurt::marching_cubes::vertices[vid]));
                voxel.values.push_back(values(pid[0], pid[1], pid[2]));
                voxel.gradients.push_back(gradient(pid[0], pid[1], pid[2]));
                voxel.hessians.push_back(hessian(pid[0], pid[1], pid[2]));
                voxel.strengths.push_back(evmin(voxel.hessians.back()).first);
            }
            broken_voxels.push_back(voxel);
        }
    }
    std::cout << "There were " << double_edges.size() << " edges containing more than one ridge point\n";
    std::cout << "Of the " << to_subdivide.size() << " associated voxels, " << voxel_counter.size() << " voxels are unique\n";
    // std::vector<std::vector<CoordID>> all_clusters;
    // find_neighbors(all_clusters, voxel_counter);
    // std::cout << "These " << voxel_counter.size() << " voxels form " << all_clusters.size() << " clusters\n";
    // std::ofstream o("clusters.txt");
    // for (int i=0; i<all_clusters.size(); ++i)
    // {
    //     std::cout << "Cluster #" << i << " contains " << all_clusters[i] << "voxels\n";
    //     auto cluster = all_clusters[i];
    //     std::set<EdgeID> edge_ids;
    //     for (int i=0; i<cluster.size(); ++i)
    //     {
    //         auto edges = voxel_counter[cluster[i]];
    //         std::for_each(edges.begin(), edges.end(), [&](int k)
    //         {
    //             edge_ids.insert(double_edges[k]);
    //         });
    //     }
    //     std::cout << "Cluster #" << i << " contains " << edge_ids.size() << " problematic edges\n";
    //     o << cluster.size() << " " << edge_ids.size() << std::endl;
    //     for (int j=0; j<cluster.size(); ++j) {
    //         o << cluster[j] << '\n';
    //     }
    //     for (auto it=edge_ids.begin(); it!=edge_ids.end(); ++it)
    //     {
    //         o << *it << '\n';
    //     }
    //     o << std::endl;
    // }
    // o.close();

    std::vector<vec3> all_edge_points;
    for (auto iter=all_processed_edges.begin(); iter!=all_processed_edges.end(); ++iter)
    {
        auto pts = iter->second;
        if (pts.empty())
            continue;
        all_edge_points.insert(all_edge_points.end(), pts.begin(), pts.end());
    }
    {
        size_t sizes[3] = { 3, all_edge_points.size() };
        spurt::nrrd_utils::writeNrrd((void*)&all_edge_points[0], output_name + "_all_points.nrrd", 
                                     nrrdTypeDouble, 2, sizes);
    }

    progress.end();
    int nfailed = 0;
    int nfixed = 0;
    progress.begin(broken_voxels.size(), "Repair ridges", 10000, "tris: 0, done: 0, ok: 0, skip: 0, underflow: 0, weak: 0, low: 0, none: 0, even: 0, failed: 0");
    for (int n = 0; n < broken_voxels.size(); ++n)
    {
        CoordID id = broken_voxels[n].id;
        std::map<int, std::vector<vec3>> found;
        std::string update_str = "fixed: " + std::to_string(nfixed) +
                                 ", failed: " + std::to_string(nfailed);
        progress.update(n, update_str);
        std::ostringstream log;
        for (int i = 0; i < 12; ++i)
        {
            CoordID v0 = id + CoordID(ivec3(spurt::marching_cubes::canonical_edge_coordinates[i][0]));
            CoordID v1 = id + CoordID(ivec3(spurt::marching_cubes::canonical_edge_coordinates[i][1]));
            EdgeID edgeid(v0, v1);
            auto iter = all_processed_edges.find(edgeid);
            assert(iter != all_processed_edges.end());
            auto solutions = iter->second;
            if (!verbose && !solutions.empty()) {
                found[i] = solutions;
            }
            else
            {
                if (verbose)
                    std::cout << "\n\nedge #" << i << " of voxel " << id << ": " << edgeid << " is being reprocessed" << '\n';
                auto result = process_edge(edgeid, gradient, hessian, values, verbose, 10*res, minval, minstr, eps, eps);
                auto coordinates = result.second;
                if (!coordinates.empty())
                {
                    found[i] = std::vector<vec3>();
                    if (verbose)
                        std::cout << "process_edge returned " << coordinates.size() << " solutions\n";
                    for (int k = 0; k < coordinates.size(); ++k)
                    {
                        double u = coordinates[k];
                        found[i].push_back((1. - u) * edgeid[0] + u * edgeid[1]);
                    }
                    if (verbose)
                        std::cout << "Found now contains " << found.size() << " entries\n";
                    all_processed_edges[edgeid] = found[i];
                }
                else if (result.first == criterion_underflow)
                {
                    if (verbose)
                        std::cout << "criterion underflow\n";
                    all_processed_edges[edgeid] = std::vector<vec3>();
                }
                else if (result.first == weak_ridge)
                {
                    if (verbose)
                        std::cout << "weak ridge\n";
                    all_processed_edges[edgeid] = std::vector<vec3>();
                }
                else if (result.first == low_ridge)
                {
                    if (verbose)
                        std::cout << "low ridge\n";
                    all_processed_edges[edgeid] = std::vector<vec3>();
                }
                else if (result.first == no_ridge_points)
                {
                    if (verbose)
                        std::cout << "no ridge point\n";
                    all_processed_edges[edgeid] = std::vector<vec3>();
                }
                else
                {
                    if (verbose)
                        std::cout << "unrecognized case\n";
                    all_processed_edges[edgeid] = std::vector<vec3>();
                }
                nprocessed++;
            }
        }

        if (found.size() >= 3)
        {
            if (verbose)
            {
                std::cout << "After looking at all edges found contains " << found.size() << " entries\n";
                for (auto iter = found.begin(); iter != found.end(); ++iter)
                {
                    std::cout << "Edge #" << iter->first << " contains " << iter->second << '\n';
                }
            }
            std::vector<triangle_type> tris;
            int tri_case;
            if (verbose)
                tri_case = triangulate(tris, found, std::cout);
            else
                tri_case = triangulate(tris, found, log);

            if (tris.size() > 0)
            {
                std::for_each(tris.begin(), tris.end(), [&](auto T)
                              {
                                  all_triangles.push_back(T);
                                  // std::cout << "\nTriangle is " << T << '\n';
                              });
                nfixed++;
            }
            else
            {
                if (verbose)
                    std::cout << "Triangulation failed\n";
                nfailed++;
            }
        }
        else if (found.size() > 0)
        {
            nfailed++;
        }
        else
        {
            nfailed++;
        }
    }
    progress.end();
    std::cout << nfixed << " / " << broken_voxels.size() << " broken voxels were fixed\n";

    if (!all_triangles.empty()) 
    {
        size_t sizes[3] = { 3, 3, all_triangles.size() };
        spurt::nrrd_utils::writeNrrd((void*)&all_triangles[0], output_name + "_mesh.nrrd", 
                                     nrrdTypeDouble, 3, sizes);
        std::vector<std::array<double, 2>> attributes(3*all_triangles.size());
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

    if (false && !broken_voxels.empty()) {
        std::cout << "There were " << broken_voxels.size() << " broken voxels (" << nfailed_to_triangulate << ")\n";
        std::ofstream output("broken_voxels.txt");
        for (int n=0; n<broken_voxels.size(); ++n) 
        {
            const broken_voxel& voxel = broken_voxels[n];
            output << "Broken voxel #" << voxel.id << '\n';
            output << "    * values:\n";
            for (int k=0; k<8; ++k) {
                output << "        [" << k << "]: " << std::scientific << voxel.values[k] << '\n';
            }
            output << "    * gradients:\n";
            for (int k=0; k<8; ++k) {
                const vec3& g = voxel.gradients[k];
                output << "        [" << k << "]: [" << std::scientific << g[0] << ", " << g[1] << ", " << g[2] << "]\n";
            }
            output << "    * hessians:\n";
            for (int k = 0; k < 8; ++k)
            {
                const mat3 &h = voxel.hessians[k];
                output << "        [" << k << "]: [" << std::scientific
                       << "[" << h(0, 0) << ", " << h(0, 1) << ", " << h(0, 2) << "], "
                       << "[" << h(1, 0) << ", " << h(1, 1) << ", " << h(1, 2) << "], "
                       << "[" << h(2, 0) << ", " << h(2, 1) << ", " << h(2, 2) << "]]\n";

            }
            output << "    * strengths:\n";
            for (int k=0; k<8; ++k)
            {
                const double& s = voxel.strengths[k];
                output << "        [" << k << "]: " << std::scientific << s << '\n'; 
            }
            output << "    * edge points:\n";
            for (auto iter=voxel.edge_points.begin(); iter!=voxel.edge_points.end(); ++iter) {
                output << "        [" << iter->first << "]: [";
                for (int k=0; k<iter->second.size(); ++k) {
                    const vec3& p = iter->second[k];
                    output << std::scientific
                        << "[" << p[0] << ", " << p[1] << ", " << p[2] << "], ";
                } 
                output << "]\n";
            }
            output << "\n";
        }
        output.close();
    }

    return 0;
}
