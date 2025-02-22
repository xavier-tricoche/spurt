/*
    This program implements a modified version of the ridge surface extraction 
    algorithm proposed by T. Schultz et al. in their 2009 IEEE VIS paper.

    T. Schultz, H. Theisel and H.-P. Seidel, 
        "Crease Surfaces: From Theory to Extraction and Application to 
         Diffusion Tensor MRI," 
    in IEEE Transactions on Visualization and Computer Graphics, vol. 16, no. 1, pp. 109-119, Jan.-Feb. 2010, doi: 10.1109/TVCG.2009.44.

    The modifications include:
        - Crease point location on an edge is refined by Newton search,
          enabled by explicit knowledge of vector h's derivatives
        - Triangulation per voxel identifies one or more topological disks
          that are subsequently triangulated
        - The implementation is embarrasingly parallel and the extraction
          is very efficient compared to alternative approaches
        - Scale space extraction is enabled by supplying scalar, gradient, 
          hessian, and ridge strength volumes that correspond to ridge strength 
          maximizing submanifold in 4D scale space (which is somewhat of an 
          open problem though).
*/

// STL
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <exception>
#include <filesystem>
#include <locale>
#include <mutex>
#include <regex>
#include <tuple>
#include <utility>
// Teem: to import value, gradient, hessian, and ridge strength as NRRD volumes
#include <teem/nrrd.h>
// spurt
#include <data/image.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <spline/spline.h>
#include <format/filename.hpp>
// spurt/crease
#include <crease/crease_mc_cases.hpp>
#include <crease/zheng.hpp>
#include <vtk/vtk_utils.hpp>
#include <format/filename.hpp>

// VTK: to export crease mesh as a VTK polydata object
#include <vtkCellArrayIterator.h>
#include <vtkVoxel.h>

// TBB
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>

// Boost C++
#include <boost/json/src.hpp>

using namespace spurt;

int verbose = 0;

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

typedef double scalar_type;
typedef long size_type;
typedef spurt::small_vector<scalar_type, 3> vector_type;
typedef spurt::small_matrix<scalar_type, 3, 3> matrix_type;
typedef spurt::small_vector<scalar_type, 2> face_pos_type;
typedef spurt::small_vector<scalar_type, 2> face_step_type;
typedef spurt::small_vector<long, 3> coord_type;
typedef spurt::small_vector<scalar_type, 3> pos_type;
typedef spurt::bounding_box<pos_type> bounds_type;
typedef spurt::small_vector<scalar_type, 6> sym_type;
typedef spurt::small_tensor<scalar_type, 3, 3, 3> tensor_type;
typedef spurt::small_vector<scalar_type, 7> CF_type;
typedef spurt::small_matrix<scalar_type, 7, 6> dCFdT_type;
typedef spurt::small_matrix<scalar_type, 7, 2> dCFdX_type;
typedef spurt::small_matrix<scalar_type, 6, 2> dTdX_type;
typedef std::array<pos_type, 3> triangle_type;

const scalar_type invalid_scalar = std::numeric_limits<scalar_type>::max();
const size_type invalid_size = -1;

const pos_type invalid_pos = pos_type(invalid_scalar);
const coord_type invalid_coord = coord_type(invalid_size);

// scalar, vector, and tensor field with C2 B-spline interpolation
template <typename Value>
using image3 = image<size_type, scalar_type, 3, Value,
                     kernels::MitchellNetravaliBC,
                     coord_type, pos_type>;

typedef image3<scalar_type> scalar_image_type;
typedef image3<vector_type> vector_image_type;
typedef image3<matrix_type> matrix_image_type;
typedef scalar_image_type::grid_type grid_type;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> dyn_mat_type;

template <typename Value>
image3<Value> import_nrrd(const std::string &filename)
{
    Nrrd *nin = nrrdNew();
    if (nrrdLoad(nin, filename.c_str(), NULL))
    {
        char *err = biffGetDone(NRRD);
        std::cerr << "Thomas Schultz's ridge method: " << err << std::endl;
        exit(-1);
    }
    typedef typename image3<Value>::base_type raster_type;
    raster_type raster = spurt::nrrd_utils::to_raster<size_type, scalar_type, 3, Value>(nin, std::is_scalar<Value>::value);
    return image3<Value>(raster);
}

// Class to uniquely represent a small set of objects of the same type
template <typename T, size_t N, typename Order = std::less<T>>
struct Connectivity : public std::array<T, N>
{
    typedef T value_type;
    typedef std::array<T, N> base_type;
    typedef Order less_type;
    typedef Connectivity<T, N, Order> self_type;
    //
    Connectivity(const T &c0 = T(invalid_scalar), const T &c1 = T(invalid_scalar), const T &c2 = T(invalid_scalar), const T c3 = T(invalid_scalar))
        : base_type()
    {
        T cs[] = {c0, c1, c2, c3};
        for (size_t i = 0; i < N; ++i)
            (*this)[i] = cs[i];
        std::sort(base_type::begin(), base_type::end(), less_type());
    }

    bool operator<(const self_type &other) const
    {
        less_type less;
        for (size_t i = 0; i < N; ++i)
        {
            if (less((*this)[i], other[i]))
                return true;
            else if (less(other[i], (*this)[i]))
                return false;
        }
        return false;
    }
};

template <typename T, size_t N, typename Order>
std::ostream &operator<<(std::ostream &os, const Connectivity<T, N, Order> &c)
{
    os << "[ ";
    std::copy(c.begin(), c.end(), std::ostream_iterator<T>(os, " "));
    os << " ]";
    return os;
}

// super ugly: made to facilitate json export...
namespace std
{
    template <typename T, size_t N, typename Order>
    std::string to_string(const Connectivity<T, N, Order> &c)
    {
        std::ostringstream os;
        os << c;
        return os.str();
    }
}

typedef Connectivity<coord_type, 2, spurt::lexicographical_order> edge_index_t;
typedef Connectivity<coord_type, 4, spurt::lexicographical_order> face_index_t;
typedef Connectivity<size_type,  2> segment_index_t;

const edge_index_t invalid_edge_index =
    edge_index_t({invalid_coord, invalid_coord});
const face_index_t invalid_face_index =
    face_index_t({invalid_coord, invalid_coord, invalid_coord, invalid_coord});
const segment_index_t invalid_segment_index =
    segment_index_t({invalid_size, invalid_size});

std::array<edge_index_t, 4> faceid_to_edgeids(const face_index_t &f)
{
    std::array<edge_index_t, 4> r;
    // f { 0, 3, 1, 2 }
    // actual order: f[0] f[2] f[3] f[1]
    r[0] = edge_index_t(f[0], f[2]);
    r[1] = edge_index_t(f[2], f[3]);
    r[2] = edge_index_t(f[3], f[1]);
    r[3] = edge_index_t(f[1], f[0]);
    return r;
}

std::vector<coord_type> faceid_to_voxelids(const face_index_t &f,
                                           const coord_type &res)
{
    coord_type diag = f[3] - f[0];
    int external_dim = -1;
    for (int d = 0; d < 3; ++d)
    {
        if (diag[d] == 0)
        {
            external_dim = d;
            break;
        }
    }
    if (external_dim == -1)
        throw std::runtime_error("invalid face computation");

    std::vector<coord_type> voxels;
    vector_type e_ext(0);
    e_ext[external_dim] = 1;
    if (f[0][external_dim] > 0)
        voxels.push_back(f[0] - e_ext);
    if (f[0][external_dim] < res[external_dim] - 1)
        voxels.push_back(f[0]);
    return voxels;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v)
{
    os << "[";
    if (v.size() > 0)
    {
        for (int i = 0; i < v.size() - 1; ++i)
            os << v[i] << ", ";
        os << v.back();
    }
    os << "]";
    return os;
}

std::pair<scalar_type, vector_type> evmin(const matrix_type &H)
{
    vector_type evals;
    matrix_type evecs;
    sym_eigensystem(evals, evecs, H);
    return std::make_pair(evals[2], vector_type(evecs.column(2)));
}

std::pair<vector_type, matrix_type> evall(const matrix_type &H)
{
    vector_type evals;
    matrix_type evecs;
    sym_eigensystem(evals, evecs, H);
    return std::make_pair(evals, evecs);
}

scalar_type ridge_strength(const matrix_type &H)
{
    vector_type evals;
    matrix_type evecs;
    sym_eigensystem(evals, evecs, H);
    return evals[2];
}

/*

The next few functions are the core ingredients of Schultz's method:
- tensor T
- vector h
- their directional derivatives

*/

scalar_type delta23(const vector_type &evals)
{
    return evals[1] - evals[2];
}

scalar_type l3(const vector_type &evals, scalar_type theta)
{
    scalar_type delta = delta23(evals);
    if (delta > theta)
        return 0;
    else
        return (1. - delta / theta) * (1. - delta / theta);
}

scalar_type dl3(const vector_type &evals, const vector_type &devals,
                double theta)
{
    // l3' = (1 - (lambda2-lambda3)/theta)^2
    // dl3' = 2*d(1 - (lambda2-lambda3)/theta)*(1 - (lambda2-lambda3)/theta)
    // dl3' = 2*(-(dlambda2-dlambda3)/theta)*(1 - (lambda2-lambda3)/theta)
    scalar_type _l3 = l3(evals, theta);
    if (_l3 == 0)
        return 0;
    else
    {
        return -2 / theta * (devals[1] - devals[2]) * (1 - delta23(evals) / theta);
    }
}

matrix_type tensor_T(const matrix_type &H, scalar_type theta)
{
    matrix_type Q;
    vector_type l;
    sym_eigensystem(l, Q, H);
    matrix_type lambda = matrix_type::identity();
    lambda(2, 2) = l3(l, theta);
    return Q * lambda * transpose(Q);
}

// compute partial derivative of tensor T from partial derivative
// of Hessian
matrix_type tensor_dT(const matrix_type &dH, const matrix_type &H,
                      scalar_type theta)
{
    /*
        cf "Matrix Differential Calculus with Applications in Statistics and Econometrics" by J.R. Magnus and H. Neudecker, Wiley, 2019
        dH_e = d\Lambda = diag(d\lambda_i's) (1)
        Notation: "_e": in the basis of eigenvectors
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
    if (_l3 == 0)
        return matrix_type(0); // dT_e = 0 => dT = 0
    else
    {
        matrix_type dH_e = transpose(Q) * dH * Q;           // (c)
        vector_type dl(dH_e(0, 0), dH_e(1, 1), dH_e(2, 2)); // (g)
        scalar_type _dl3 = dl3(l, dl, theta);
        matrix_type dT_e = 0;
        dT_e(2, 2) = _dl3;
        return Q * dT_e * transpose(Q); // back to dataset coordinates
    }
}

vector_type vector_h(const vector_type &g, const matrix_type &H,
                     scalar_type theta)
{
    // h = Tg - g
    matrix_type T = tensor_T(H, theta);
    return T * g - g;
}

vector_type vector_dh(const vector_type &g, const matrix_type &H,
                      const vector_type &dg, const matrix_type &dH,
                      scalar_type theta)
{
    // ∇h = ∇Tg + T∇g − ∇g
    matrix_type T = tensor_T(H, theta);
    matrix_type dT = tensor_dT(dH, H, theta);
    return dT * g + T * dg - dg;
}

vector_type ridge_normal(const pos_type &point,
                         const vector_type &g, const matrix_type &H,
                         const vector_image_type &gradient,
                         const matrix_image_type &hessian,
                         scalar_type theta,
                         bool verbose = false)
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

scalar_type edge_linearity_value(const matrix_type &T1, const matrix_type &T2)
{
    return trace(transpose(T1) * T2) / sqrt(trace(transpose(T1) * T1)) / sqrt(trace(transpose(T2) * T2));
}

scalar_type theta_threshold(const matrix_image_type &hessian, double eps = 0.005)
{
    vector_type evals;
    matrix_type evecs;
    std::vector<scalar_type> _l23;
    for (const matrix_type &H : hessian)
    {
        sym_eigensystem(evals, evecs, H);
        _l23.push_back(evals[1] - evals[2]);
    }
    auto iters = std::minmax_element(_l23.begin(), _l23.end());
    return eps * (*iters.second - *iters.first);
}

// end of Schultz core section

template <typename T, size_t N>
std::ostream &operator<<(std::ostream &os, const std::array<T, N> &a)
{
    os << "[";
    for (int i = 0; i < N - 1; ++i)
    {
        os << a[i] << ", ";
    }
    os << a[N - 1] << "]";
    return os;
}

struct h_root_finder
{
    int edge_dim(const edge_index_t &e)
    {
        coord_type de = spurt::abs(e[1] - e[0]);
        return std::distance(de.begin(), std::max_element(de.begin(), de.end()));
    }

    vector_type compute_h(const pos_type &p)
    {
        vector_type gv = m_gradient.value(p);
        matrix_type hess = m_hessian.value(p);
        return vector_h(gv, hess, m_theta);
    }

    std::pair<vector_type, vector_type> compute_h_dh(const pos_type &p, int dim)
    {
        vector_type dp(0);
        dp[dim] = 1.0e-6;
        vector_type gv = m_gradient.value(p);
        matrix_type hess = m_hessian.value(p);
        tensor_type hess_prime = m_hessian.derivative(p);
        vector_type dg = m_gradient.derivative(p).column(dim);
        vector_type h = vector_h(gv, hess, m_theta);
        vector_type dh = vector_dh(gv, hess, dg, hess_prime.layer(dim), m_theta);
        return std::make_pair(h, dh);
    }

    h_root_finder(const vector_image_type &grad, const matrix_image_type &hess,
                  scalar_type theta, int maxiter = 10)
        : m_gradient(grad), m_hessian(hess), m_theta(theta), m_maxiter(maxiter) {};

    bool lnsearch(scalar_type &u, vector_type &h, const scalar_type du,
                  const pos_type &p0, const pos_type &p1,
                  bool verbose = false)
    {
        scalar_type lambda = 1.0;
        const scalar_type alpha = 1e-4;

        scalar_type u0 = u;
        vector_type h0 = h;
        scalar_type d;
        if (u0 + du < 0)
            d = -u0;
        else if (u0 + du > 1)
            d = (1 - u0);
        else
            d = du;
        if (verbose)
        {
            std::cout << std::fixed << std::setprecision(12) << "linesearch: u=" << u << ", h=" << h << ", du=" << du << ", d=" << d << '\n';
        }
        scalar_type v0 = spurt::norm(h0);
        pos_type p;
        for (unsigned int i = 0; i < 20; ++i)
        {
            u = u0 + lambda * d;
            p = (1 - u) * p0 + u * p1;
            h = compute_h(p);
            if (spurt::norm(h) < (1 - alpha * lambda) * v0)
            {
                return true;
            }
            // if (verbose) {
            //     std::cout << std::fixed << std::setprecision(12) << "u=" << u << ", h=" << h << '\n';
            // }

            lambda *= 0.5;
        }

        return false;
    }

    bool search(scalar_type &u, vector_type &h, const edge_index_t &e,
                int maxiter = -1, scalar_type eps = 1.0e-6,
                bool verbose = false)
    {
        scalar_type u0 = u;
        const pos_type p0 = m_gradient.grid()(e[0]);
        const pos_type p1 = m_gradient.grid()(e[1]);
        int dim = edge_dim(e);
        if (verbose)
        {
            std::cout << "e=" << e << ", dim=" << dim << '\n';
        }
        vector_type dh;
        vector_type dp(0);
        dp[dim] = eps;
        if (maxiter < 0)
            maxiter = m_maxiter;
        for (int i = 0; i < maxiter; ++i)
        {
            pos_type p = (1 - u) * p0 + u * p1;
            std::tie(h, dh) = compute_h_dh(p, dim);
            auto h0 = h;
            auto dv = -(h / dh);
            auto du = spurt::mean(dv);
            auto u1 = u;
            if (!lnsearch(u, h, du, p0, p1, verbose))
            {
                break;
            }
            if (spurt::norm(h) < eps)
            {
                return true;
            }
        }
        return false;
    }

    const vector_image_type &m_gradient;
    const matrix_image_type &m_hessian;
    scalar_type m_theta;
    int m_maxiter;
};

std::pair<int, scalar_type>
process_edge(const edge_index_t &e,
             const scalar_image_type &data,
             const vector_image_type &gradient,
             const matrix_image_type &hessian,
             scalar_type theta,
             int verbose = 0,
             bool improve = true)
{
    typedef Eigen::Array<scalar_type, 3, 1> array_type;
    vector_type h0 = vector_h(gradient(e[0]), hessian(e[0]), theta);
    vector_type h1 = vector_h(gradient(e[1]), hessian(e[1]), theta);
    if (inner(h0, h1) < 0)
    {
        scalar_type n0 = h0.norm();
        scalar_type n1 = h1.norm();

        vector_type dv = -h0 / (h1 - h0);
        scalar_type ubest = 0;
        scalar_type bestnorm = std::numeric_limits<scalar_type>::max();
        for (int i = 0; i < 3; ++i)
        {
            scalar_type u = dv[i];
            if (u < 0 || u > 1)
                continue;
            pos_type p = (1 - u) * gradient.grid()(e[0]) + u * gradient.grid()(e[1]);
            scalar_type n = spurt::norm(vector_h(gradient.value(p), hessian.value(p), theta));
            if (n < bestnorm)
            {
                bestnorm = n;
                ubest = u;
            }
        }
        scalar_type u = ubest;

        if (improve)
        {
            pos_type p = (1 - u) * e[0] + u * e[1];
            vector_type h = vector_h(gradient.value(p), hessian.value(p), theta);
            scalar_type hinit = spurt::norm(h);
            h_root_finder root(gradient, hessian, theta, 100);
            auto ignored = root.search(u, h, e);
        }
        return std::pair(1, u);
    }
    return std::pair(0, 0);
}

std::pair<scalar_type, scalar_type>
evaluate(const pos_type &point, const scalar_image_type &values,
         const matrix_image_type &hessian)
{
    scalar_type f = values.value(point);
    matrix_type H = hessian.value(point);
    return std::make_pair(f, ridge_strength(H));
}

size_t nfour = 0;
size_t nthree = 0;
size_t nmore_than_four = 0;
size_t ntwo = 0;
size_t n_one = 0;

// voxel ridge information
struct CreasePoint
{
    pos_type position;
    edge_index_t edge_index; // valid if on an edge
    face_index_t face_index; // valid if L point on face
    scalar_type value;
    scalar_type strength;
    bool singular;

    CreasePoint(const pos_type &p = invalid_pos,
                const edge_index_t &e = invalid_edge_index,
                const face_index_t &f = invalid_face_index,
                const scalar_type &v = invalid_scalar,
                const scalar_type &str = invalid_scalar,
                bool s = false)
        : position(p), edge_index(e), face_index(f), value(v),
          strength(str), singular(s) {}
};
typedef CreasePoint crease_point_t;

std::ostream &operator<<(std::ostream &os, const crease_point_t &cp)
{
    os << "Crease Point: {\n"
       << " * position: " << cp.position << '\n'
       << " * edge_index: " << cp.edge_index << '\n'
       << " * face_index: " << cp.face_index << '\n'
       << " * value: " << cp.value << '\n'
       << " * strength: " << cp.strength << '\n'
       << " * singular: " << (cp.singular ? "true" : "false") << '\n'
       << "}";
    return os;
}

struct FaceInfo
{
    face_index_t face_id;
    size_type Lpoint_index;
    std::vector<segment_index_t> segments;

    FaceInfo(const face_index_t &f = invalid_face_index, size_type lptid = invalid_size)
        : face_id(f), Lpoint_index(lptid), segments() {}
};
typedef FaceInfo face_info_t;

typedef std::array<size_type, 3> triangle_index_t;
struct VoxelInfo
{
    coord_type voxel_id;
    std::vector<face_info_t> ridge_faces;
    std::vector<triangle_index_t> ridge_triangles;

    VoxelInfo(const coord_type &v = invalid_coord)
        : voxel_id(v), ridge_faces(), ridge_triangles() {}
};
typedef VoxelInfo voxel_info_t;

// all the unique edges of a given set of voxel coordinates
void unique_edge_indices(std::vector<edge_index_t> &edge_indices,
                         const std::vector<coord_type> &voxel_indices)
{
    const coord_type ex = coord_type(1, 0, 0);
    const coord_type ey = coord_type(0, 1, 0);
    const coord_type ez = coord_type(0, 0, 1);
    coord_type basis[] = {ex, ey, ez};

    std::set<coord_type, spurt::lexicographical_order> unique_edge_start_indices;
    for (int dim = 0; dim < 3; ++dim)
    {
        const coord_type &b0 = basis[dim];
        const coord_type &b1 = basis[(dim + 1) % 3];
        const coord_type &b2 = basis[(dim + 2) % 3];
        unique_edge_start_indices.clear();
        for (const coord_type &vc : voxel_indices)
        {
            unique_edge_start_indices.insert(vc);
            unique_edge_start_indices.insert(vc + b1);
            unique_edge_start_indices.insert(vc + b2);
            unique_edge_start_indices.insert(vc + b1 + b2);
        }
        for (const coord_type &pc : unique_edge_start_indices)
        {
            edge_indices.push_back(edge_index_t(pc, pc + b0));
        }
    }
}

// edge to varying dimension
// prerequisite: e is valid edge index:
// first and second index are distinct and differ only by a single coordinate
int edge_to_dim(const edge_index_t &e)
{
    const coord_type ex = coord_type(1, 0, 0);
    const coord_type ey = coord_type(0, 1, 0);
    const coord_type ez = coord_type(0, 0, 1);

    coord_type step = spurt::abs(e[1] - e[0]);
    if (spurt::all(step == ex))
        return 0;
    else if (spurt::all(step == ey))
        return 1;
    else if (spurt::all(step == ez))
        return 2;
    else
        throw std::runtime_error("invalid edge index in edge_to_dim");
}

// all the unique faces adjacents to the given set of edges
void unique_faces(std::vector<face_info_t> &faces,
                  const std::map<edge_index_t, size_type> &edge_index_to_crease_point_index,
                  const std::vector<coord_type> &valid_voxel_indices)
{
    const coord_type ex = coord_type(1, 0, 0);
    const coord_type ey = coord_type(0, 1, 0);
    const coord_type ez = coord_type(0, 0, 1);
    coord_type basis[] = {ex, ey, ez};

    // std::cout << "edges in input are: \n";
    // for (auto e : edge_index_to_crease_point_index) {
    //     std::cout << e.first << '\n';
    // }

    coord_type dummy_res(10000000, 10000000, 10000000);

    std::set<coord_type, spurt::lexicographical_order> valid_voxels;
    valid_voxels.insert(valid_voxel_indices.begin(), valid_voxel_indices.end());

    std::set<face_index_t> _unique_face_indices;
    for (auto v : edge_index_to_crease_point_index)
    {
        const edge_index_t &e = v.first;
        // std::cout << "currently processing edge " << e << '\n';
        int dim = edge_to_dim(e);
        const coord_type &a = e[0];
        const coord_type &b = e[1];
        for (int d = 0; d < 3; ++d)
        {
            if (d == dim)
                continue;
            face_index_t face_id0(a, b, b + basis[d], a + basis[d]);
            face_index_t face_id1(a, b, b - basis[d], a - basis[d]);
            // std::cout << "Two candidate face indices are:\n"
            //     << face_id0 << '\n'
            //     << face_id1 << '\n';
            std::vector<coord_type> voxels0 = faceid_to_voxelids(face_id0, dummy_res);
            std::vector<coord_type> voxels1 = faceid_to_voxelids(face_id1, dummy_res);
            for (auto v0 : voxels0)
            {
                // std::cout << "Voxel for face " << face_id0 << " is " << v0 << '\n';
                if (valid_voxels.find(v0) != valid_voxels.end())
                {
                    _unique_face_indices.insert(face_id0);
                }
            }
            for (auto v1 : voxels1)
            {
                // std::cout << "Voxel for face " << face_id1 << " is " << v1 << '\n';
                if (valid_voxels.find(v1) != valid_voxels.end())
                {
                    _unique_face_indices.insert(face_id1);
                }
            }
        }
    }

    faces.clear();
    faces.reserve(_unique_face_indices.size());
    for (auto fid : _unique_face_indices)
    {
        faces.push_back(face_info_t(fid));
    }
}

// all the unique voxels adjacent to the given set of faces
void unique_voxels(std::vector<voxel_info_t> &voxels,
                   const std::vector<face_info_t> &faces,
                   const std::vector<coord_type> &valid_voxels)
{
    std::set<coord_type, spurt::lexicographical_order> _valid_voxels;
    _valid_voxels.insert(valid_voxels.begin(), valid_voxels.end());
    std::map<coord_type, voxel_info_t, spurt::lexicographical_order> _unique_voxels;

    coord_type dummy_res(10000000, 10000000, 10000000);
    for (auto f : faces)
    {
        std::vector<coord_type> voxels = faceid_to_voxelids(f.face_id, dummy_res);
        for (auto v : voxels)
        {
            if (_valid_voxels.find(v) == _valid_voxels.end())
                continue;
            auto iter = _unique_voxels.find(v);
            if (iter != _unique_voxels.end())
            {
                iter->second.ridge_faces.push_back(f);
            }
            else
            {
                _unique_voxels[v] = voxel_info_t(v);
                _unique_voxels[v].ridge_faces.push_back(f);
            }
        }
    }
    voxels.clear();
    voxels.reserve(_unique_voxels.size());
    for (auto v2i : _unique_voxels)
    {
        voxels.push_back(v2i.second);
    }
}

void select_voxel_indices(std::vector<coord_type> &voxel_indices,
                          const bounds_type &bounds,
                          const coord_type &shape)
{
    if (spurt::all(bounds.min() < bounds.max()))
    {
        bounds_type b = bounds;
        for (int i = 0; i < 3; ++i)
        {
            if (b.min()[i] < 0)
                b.min()[i] = 0;
            if (b.max()[i] > shape[i])
                b.max()[i] = shape[i];
        }
        voxel_indices.clear();
        const pos_type &l = b.min();
        const pos_type &h = b.max();
        coord_type low = floor(l);
        coord_type high = ceil(h);
        for (int k = low[2]; k < high[2]; ++k)
        {
            for (int j = low[1]; j < high[1]; ++j)
            {
                for (int i = low[0]; i < high[0]; ++i)
                {
                    voxel_indices.push_back(coord_type(i, j, k));
                }
            }
        }
    }
    else
    {
        voxel_indices.resize((shape[0] - 1) * (shape[1] - 1) * (shape[2] - 1));
        size_type n = 0;
        for (int k = 0; k < shape[2] - 1; ++k)
        {
            for (int j = 0; j < shape[1] - 1; ++j)
            {
                for (int i = 0; i < shape[0] - 1; ++i, ++n)
                {
                    voxel_indices[n] = coord_type(i, j, k);
                }
            }
        }
    }
}

struct value_set
{
    scalar_type value;
    scalar_type strength;
    vector_type gradient;
    matrix_type hessian;
    vector_type normal;
};

template <typename T>
struct container_wrapper : public std::vector<T>
{
    typedef std::vector<T> base_type;
    typedef typename base_type::size_type idx_type;
    container_wrapper() : base_type() {};

    idx_type add(const T &val)
    {
        base_type::push_back(val);
        return base_type::size() - 1;
    }
};

struct GraphMap
{
    typedef std::array<size_type, 2> neighbor_t;
    typedef std::map<size_type, neighbor_t> map_t;
    map_t m_edges;

    GraphMap(const std::vector<segment_index_t> &segments)
        : m_edges()
    {
        for (auto s : segments)
        {
            for (int i = 0; i < 2; ++i)
            {
                size_type v = s[i];
                size_type w = s[(i + 1) % 2];
                auto iter = m_edges.find(v);
                if (iter == m_edges.end())
                {
                    m_edges[v] = neighbor_t({w, invalid_size});
                }
                else
                {
                    iter->second[1] = w;
                }
            }
        }
    }

    std::vector<size_type> unique_ids() const
    {
        std::vector<size_type> r;
        for (auto kv : m_edges)
        {
            r.push_back(kv.first);
        }
        return r;
    }

    size_type next_id(size_type at, size_type from = invalid_size) const
    {
        const neighbor_t &ns = const_cast<map_t &>(m_edges)[at];
        if (from == invalid_size)
            return ns[0];
        else if (ns[0] == from)
            return ns[1];
        else
            return ns[0];
    }
};

void compute_cycles(std::vector<std::vector<size_type>> &cycles,
                    const std::vector<segment_index_t> &segments,
                    bool verbose = false)
{

    typedef std::vector<size_type> cycle_t;
    GraphMap gmap(segments);
    std::set<size_type> checked;

    std::vector<size_type> unique_ids = gmap.unique_ids();
    for (size_type id : unique_ids)
    {
        if (checked.find(id) != checked.end())
        {
            continue;
        }
        checked.insert(id);
        cycles.push_back(cycle_t());
        cycles.back().push_back(id);
        size_type from = id;
        size_type at = gmap.next_id(id);
        while (at != invalid_size)
        {
            cycles.back().push_back(at);
            if (at == id)
            {
                break;
            }
            checked.insert(at);
            size_type next = gmap.next_id(at, from);
            from = at;
            at = next;
        }
    }
}

bool find_Lpoint_on_face(crease_point_t &found,
                         const face_index_t &face_index,
                         const scalar_image_type &values,
                         const matrix_image_type &hessian)
{
    pos_type Lpt;
    pos_type orig = values.grid()(face_index[0]);
    vector_type e0 = values.grid()(face_index[2]) - orig;
    vector_type e1 = values.grid()(face_index[1]) - orig;
    spurt::zheng::Face<scalar_type> f(orig, e0, e1);
    if (spurt::zheng::findLPoint(Lpt, f, hessian))
    {
        found = crease_point_t(Lpt, invalid_edge_index, face_index, 0, 0, true);
        return true;
    }

    return false;
}

bool connect_crease_points(face_info_t &faceinfo,
                           const std::vector<size_type> &crease_point_ids,
                           const container_wrapper<crease_point_t> &crease_points,
                           const vector_image_type &gradient,
                           const matrix_image_type &hessian,
                           const scalar_type &theta0)
{
    if (crease_point_ids.size() == 1)
    {
        return false;
    }
    else if (crease_point_ids.size() == 2)
    {
        faceinfo.segments.push_back(segment_index_t(crease_point_ids[0], crease_point_ids[1]));
        return true;
    }
    else
    {
        std::vector<int> valid_indices;
        size_type L_index = invalid_size;
        std::vector<vector_type> normals;
        for (int i = 0; i < crease_point_ids.size(); ++i)
        {
            auto id = crease_point_ids[i];
            if (!crease_points[id].singular)
            {
                valid_indices.push_back(i);
                auto p = crease_points[id].position;
                auto g = gradient.value(p);
                auto h = hessian.value(p);
                normals.push_back(ridge_normal(p, g, h, gradient, hessian, theta0));
            }
            else
            {
                L_index = id;
            }
        }

        std::vector<scalar_type> scores;
        std::vector<std::array<int, 2>> pairings;
        for (int i = 0; i < valid_indices.size() - 1; ++i)
        {
            for (int j = i + 1; j < valid_indices.size(); ++j)
            {
                pairings.push_back({valid_indices[i], valid_indices[j]});
            }
        }
        for (auto apair : pairings)
        {
            auto id0 = crease_point_ids[apair[0]];
            auto id1 = crease_point_ids[apair[1]];
            const pos_type &p0 = crease_points[id0].position;
            const pos_type &p1 = crease_points[id1].position;
            vector_type seg = p1 - p0;
            seg /= spurt::norm(seg);
            scores.push_back(std::abs(spurt::inner(seg, normals[apair[0]])) +
                             std::abs(spurt::inner(seg, normals[apair[1]])));
        }
        int best = std::distance(
            scores.begin(),
            std::min_element(scores.begin(), scores.end()));

        const size_type &best_id0 = crease_point_ids[pairings[best][0]];
        const size_type &best_id1 = crease_point_ids[pairings[best][1]];
        faceinfo.segments.push_back(segment_index_t(best_id0, best_id1));

        // Identify remaining crease points
        std::vector<size_type> others;
        for (int i = 0; i < crease_point_ids.size() && others.size() < 2; ++i)
        {
            auto id = crease_point_ids[i];
            if (id != best_id0 && id != best_id1)
            {
                others.push_back(id);
            }
        }
        if (others.size() == 2)
            faceinfo.segments.push_back(segment_index_t(others[0], others[1]));

        return true;
    }
}

void visualize(const container_wrapper<crease_point_t> &cpoints,
               const std::vector<triangle_index_t> &all_triangles,
               const std::vector<std::vector<size_type>> &all_open_cycles,
               const std::string &fname,
               bool vis)
{
    std::cout << "There are " << cpoints.size() << " points in input\n";
    std::cout << "There are " << all_triangles.size() << " triangles in input\n";
    std::vector<pos_type> points(cpoints.size(), pos_type(0));
    std::vector<scalar_type> values(cpoints.size(), 0);
    std::vector<scalar_type> strengths(cpoints.size(), 0);

    for (size_type i = 0; i < cpoints.size(); ++i)
    {
        points[i] = cpoints[i].position;
        values[i] = cpoints[i].value;
        strengths[i] = cpoints[i].strength;
    }
    vtkSmartPointer<vtkPolyData> pd = vtk_utils::make_points(points);

    VTK_CREATE(vtkCellArray, triangles);
    for (size_type i = 0; i < all_triangles.size(); ++i)
    {
        const triangle_index_t &T = all_triangles[i];
        triangles->InsertNextCell(3);
        triangles->InsertCellPoint(T[0]);
        triangles->InsertCellPoint(T[1]);
        triangles->InsertCellPoint(T[2]);
    }
    pd->SetPolys(triangles);

    VTK_CREATE(vtkCellArray, lines);
    for (size_type i = 0; i < all_open_cycles.size(); ++i)
    {
        const std::vector<size_type> &aline = all_open_cycles[i];
        lines->InsertNextCell(aline.size());
        for (size_type j : aline)
            lines->InsertCellPoint(j);
    }
    pd->SetLines(lines);

    pd = vtk_utils::add_scalars(pd, values, true, "values", false);
    pd = vtk_utils::add_scalars(pd, strengths, true, "ridge_strength", true);

    VTK_CREATE(vtkXMLPolyDataWriter, writer);
    writer->SetInputData(pd);
    auto name = filename::replace_extension(fname, "vtp");
    writer->SetFileName(name.c_str());
    writer->Write();

    if (!vis)
        return;

    VTK_CREATE(vtkPolyDataMapper, mesh_mapper);
    mesh_mapper->SetInputData(pd);
    mesh_mapper->ScalarVisibilityOn();
    VTK_CREATE(vtkActor, mesh_actor);
    mesh_actor->SetMapper(mesh_mapper);
    mesh_actor->GetProperty()->SetLineWidth(2);
    mesh_actor->GetProperty()->VertexVisibilityOn();

    VTK_CREATE(vtkRenderer, renderer);
    renderer->AddActor(mesh_actor);
    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
    window->SetSize(1920, 1080);
    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    interactor->Initialize();
    window->Render();
    interactor->Start();
}

struct json_traits
{
    typedef boost::json::value jvalue_t;
    typedef boost::json::array jarray_t;

    template <typename T, typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
    jvalue_t map(const T &t)
    {
        return jvalue_t(t);
    }

    template <typename T, size_t N>
    jvalue_t map(const std::array<T, N> &a)
    {
        jarray_t jv;
        for (const T &t : a)
        {
            jv.push_back(map(t));
        }
        return jv;
    }

    template <typename T, size_t N>
    jvalue_t map(const spurt::small_vector<T, N> &v)
    {
        jarray_t jv;
        for (const T &t : v)
        {
            jv.push_back(map(t));
        }
        return jv;
    }

    template <typename T>
    jvalue_t map(const std::vector<T> &v)
    {
        jarray_t jv;
        for (const T &t : v)
        {
            jv.push_back(map(t));
        }
        return jv;
    }

    template <typename T, size_t N, typename Order>
    jvalue_t map(const Connectivity<T, N, Order> &c)
    {
        jarray_t jv;
        for (const T &t : c)
        {
            jv.push_back(map(t));
        }
        return jv;
    }

    jvalue_t map(const face_info_t &f)
    {
        jvalue_t jv = {
            {"face_id", map(f.face_id)},
            {"Lpoint_index", map(f.Lpoint_index)},
            {"segments", map(f.segments)}};
        return jv;
    }

    template <typename T>
    jvalue_t operator()(const T &t)
    {
        return map(t);
    }
};

void export_crease_points(const container_wrapper<crease_point_t> &unique_crease_points, const std::string &filename)
{
    using namespace boost::json;
    std::ofstream out(filename, std::ios::out);
    out << "[";
    for (size_type i = 0; i < unique_crease_points.size(); ++i)
    {
        const crease_point_t &cp = unique_crease_points[i];
        const pos_type &p = cp.position;
        const edge_index_t &e = cp.edge_index;
        const face_index_t &f = cp.face_index;
        json_traits jt;
        value jv = {
            {"index", i},
            {"position", jt(p)},
            {"edge_index", jt(e)},
            {"face_index", jt(f)},
            {"value", cp.value},
            {"strength", cp.strength},
            {"singular", cp.singular}};
        out << jv;
        if (i < unique_crease_points.size() - 1)
            out << ",";
    }
    out << "]";
    out.close();
}

void export_faces(const std::vector<face_info_t> &all_faces,
                  const std::string &filename)
{
    std::ofstream out(filename, std::ios::out);
    out << "[";
    using namespace boost::json;
    json_traits jt;
    for (size_type i = 0; i < all_faces.size(); ++i)
    {
        const face_info_t &f = all_faces[i];
        out << jt(f);
        if (i < all_faces.size() - 1)
            out << ",";
    }
    out << "]";
    out.close();
}

void export_voxels(const std::vector<voxel_info_t> &all_voxels,
                   const std::string &filename)
{
    std::ofstream out(filename, std::ios::out);
    out << "[";
    using namespace boost::json;
    json_traits jt;
    for (size_type i = 0; i < all_voxels.size(); ++i)
    {
        const voxel_info_t &v = all_voxels[i];
        value jv = {
            {"voxel_id", jt(v.voxel_id)},
            {"ridge_faces", jt(v.ridge_faces)},
            {"ridge_triangles", jt(v.ridge_triangles)}};
        out << jv;
        if (i < all_voxels.size() - 1)
            out << ",";
    }
    out << "]";
    out.close();
}

int main(int argc, const char *argv[])
{
    std::string scalar_name, gradient_name, hessian_name, output_name, strength_name;
    scalar_type minval, minstr, eps, mind;
    std::string minval_str, minstr_str;
    size_type minsize;
    int res, niter;
    int verbose;
    bool vis, export_data;
    coord_type voxel_id = invalid_coord;
    bounds_type bounds(pos_type(0.), pos_type(-1.));
    spurt::vec4 dv(0.1, 0.2, 0.3, 0.4);

    namespace cl = spurt::command_line;
    cl::option_traits
        required_group(true, false, "Required Options"),
        positional_group(true, true, "Positional Group"),
        optional_group(false, false, "Optional Group");

    cl::option_parser parser(argv[0],
                             "Extract ridge surfaces from scalar volume using Thomas Schultz's 2009 method");
    try
    {
        parser.use_short_symbols(false);
        parser.use_brackets(true);

        parser.add_value("value", scalar_name, "Scalar volume", required_group);
        parser.add_value("gradient", gradient_name, "Gradient volume", required_group);
        parser.add_value("hessian", hessian_name, "Hessian volume", required_group);
        parser.add_value("output", output_name, "Output filename", required_group);
        parser.add_value("strength", strength_name, "Ridge strength volume", required_group);
        parser.add_value("minval", minval_str, "0", "Min scalar value", optional_group);
        parser.add_value("minstr", minstr_str, "0", "Min ridge strength (<=0)", optional_group);
        parser.add_flag("vis", vis, "Visualize resulting surfaces", optional_group);
        parser.add_tuple<3>("voxel", voxel_id, voxel_id, "Coordinates of single voxel to process");
        parser.add_tuple<3>("blower", bounds.min(), bounds.min(), "Lower bounds of domain to consider");
        parser.add_tuple<3>("bupper", bounds.max(), bounds.max(), "Upper bounds of domain to consider");
        parser.add_flag("export", export_data, "Export all intermeidate data", optional_group);
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

    std::locale::global(std::locale(""));
    std::cout.imbue(std::locale());

    spurt::timer load_time;
    load_time.start();
    scalar_image_type values = import_nrrd<scalar_type>(scalar_name);
    vector_image_type gradient = import_nrrd<vector_type>(gradient_name);
    matrix_image_type hessian = import_nrrd<matrix_type>(hessian_name);
    scalar_image_type strength = import_nrrd<scalar_type>(strength_name);
    load_time.stop();
    std::cout << "loading data took " << load_time.wall_time() << '\n';
    double theta0 = theta_threshold(hessian);
    std::cout << "theta threshold is " << theta0 << '\n';

    if (minval_str.back() == '%')
    {
        int pct = std::stoi(minval_str.substr(0, minval_str.size() - 1));
        std::vector<double> vals(values.begin(), values.end());
        std::sort(vals.begin(), vals.end());
        minval = vals[std::floor(pct * vals.size() / 100.)];
    }
    else
        minval = std::stof(minval_str);

    if (minstr_str.back() == '%')
    {
        int pct = 100 - std::stoi(minstr_str.substr(0, minstr_str.size() - 1));
        std::vector<double> vals(strength.begin(), strength.end());
        std::sort(vals.begin(), vals.end());
        minstr = vals[std::floor(pct * vals.size() / 100.)];
    }
    else
        minstr = std::stof(minstr_str);

    std::cout << "Filtering thresholds set to: value: " << minval << ", ridge strength: " << minstr << '\n';

    auto shape = values.grid().resolution();
    std::vector<coord_type> all_voxel_indices;
    if (spurt::any(voxel_id != invalid_coord))
    {
        all_voxel_indices.clear();
        all_voxel_indices.push_back(voxel_id);
        verbose = 2;
    }
    else
    {
        select_voxel_indices(all_voxel_indices, bounds, shape);
    }
    std::vector<edge_index_t> all_edge_indices;
    unique_edge_indices(all_edge_indices, all_voxel_indices);
    /*

        Algorithm:

        1: process all edges of all voxels (skipping redundant ones),
           looking for crease points
    */

    // edge -> globalcrease point id
    std::map<edge_index_t, size_type> edge_to_crease_point_index;

    // information about crease points
    container_wrapper<crease_point_t> unique_crease_points;

    std::atomic<int> nskipped = 0;
    std::atomic<int> nfailed_to_triangulate = 0;
    std::atomic<int> nnot_enough_points = 0;
    std::atomic<int> nsucceeded = 0;
    std::atomic<int> nprocessed = 0;
    std::atomic<int> nexcluded = 0;
    std::atomic<int> nnone = 0;
    std::atomic<int> neven = 0;
    std::atomic<int> nweak = 0;
    std::atomic<int> nlow = 0;
    std::atomic<int> nunderflow = 0;
    std::atomic<int> nfiltered = 0;
    std::atomic<int> noddfaces = 0;

    spurt::ProgressDisplay progress;

    VTK_CREATE(vtkCellArray, vtk_cells);
    VTK_CREATE(vtkCellArray, vtk_lines);
    VTK_CREATE(vtkCellArray, vtk_verts);
    std::vector<edge_index_t> double_edges;
    std::vector<coord_type> to_subdivide;
    std::vector<int> voxel_to_edge;

    srand48(130819751900);
    std::atomic<int> edge_counter = 0;
    progress.begin(all_edge_indices.size(), "Extract ridge points", 10000, "done: 0, found: 0, weak: 0");
    std::mutex update_mtx;
    std::mutex edge_add_mtx;
    tbb::parallel_for(tbb::blocked_range<int>(0, all_edge_indices.size()),
                      [&](tbb::blocked_range<int> r)
                      {
                          for (int n = r.begin(); n < r.end(); ++n)
                          {
                              const edge_index_t &the_edge_index = all_edge_indices[n];
                              ++edge_counter;
                              std::unique_lock<std::mutex> lock(update_mtx, std::defer_lock);
                              if (lock.try_lock())
                              {
                                  std::ostringstream update_oss;
                                  update_oss.imbue(std::locale());
                                  update_oss << "done: " << edge_counter
                                             << ", found: " << unique_crease_points.size()
                                             << ", weak: " << nweak;
                                  progress.update(edge_counter, update_oss.str());
                              }

                              // check if current edge satisfies threshold requirements
                              const coord_type &the_coord_0 = the_edge_index[0];
                              const coord_type &the_coord_1 = the_edge_index[1];
                              if (values(the_coord_0) < minval || values(the_coord_1) < minval)
                              {
                                  nnone++;
                                  nweak++;
                                  continue;
                              }
                              if (evmin(hessian(the_coord_0)).first > minstr ||
                                  evmin(hessian(the_coord_1)).first > minstr)
                              {
                                  nweak++;
                                  nnone++;
                                  continue;
                              }
                              const pos_type &the_point_0 = values.grid()(the_coord_0);
                              const pos_type &the_point_1 = values.grid()(the_coord_1);

                              auto result = process_edge(the_edge_index, values, gradient, hessian, theta0, verbose, true);
                              if (result.first == 1)
                              {
                                  scalar_type u = result.second;
                                  pos_type p = (1. - u) * the_point_0 + u * the_point_1;
                                  auto v = evaluate(p, values, hessian);

                                  {
                                      std::lock_guard<std::mutex> lock(edge_add_mtx);
                                      edge_to_crease_point_index[the_edge_index] = unique_crease_points.add(crease_point_t(p, the_edge_index, invalid_face_index, v.first, v.second, false));
                                  }
                              }
                              else
                              {
                                  nnone++;
                              }
                              nprocessed++;
                          }
                      });
    // all valid edges have now been processed
    std::cout << "\n\nEdge processing complete:\n";
    std::cout << progress << '\n';
    std::cout << "  * processed: " << all_edge_indices.size() << '\n';
    std::cout << "  * found: " << unique_crease_points.size() << '\n';
    std::cout << "  * weak (skipped): " << nweak << '\n';

    /*


        2: process all faces adjacent to edges with crease points looking for
           crease surface segments, handling ambiguous/pathological cases

    */

    // export unique_crease_points

    // now processing active faces
    std::vector<face_info_t> all_faces;
    unique_faces(all_faces, edge_to_crease_point_index, all_voxel_indices);

    std::atomic<size_type> n_pathological = 0;
    std::atomic<size_type> n_singular = 0;
    std::atomic<size_type> n_failed = 0;
    std::atomic<int> face_counter = 0;
    std::atomic<int> voxel_counter = 0;
    std::mutex debug_mtx;
    std::mutex lpoint_add_mtx;
    std::cout << "\nProcessing all " << all_faces.size() << " unique faces\n";
    progress.begin(all_faces.size(), "Extract ridge points", 10000, "done: 0, found: 0, pathological: 0, singular: 0");
    std::atomic<size_type> n_faces_found = 0;

    tbb::parallel_for(tbb::blocked_range<int>(0, all_faces.size()),
                      [&](tbb::blocked_range<int> r)
                      {
                          for (int n = r.begin(); n < r.end(); ++n)
                          {
                              face_info_t &the_face = all_faces[n];
                              const face_index_t &the_face_index = the_face.face_id;
                              ++face_counter;
                              std::unique_lock<std::mutex> lock(update_mtx, std::defer_lock);
                              if (lock.try_lock())
                              {
                                  std::ostringstream update_oss;
                                  update_oss << "done: " << face_counter
                                             << ", found: " << n_faces_found
                                             << ", pathological: " << n_pathological
                                             << ", singular: " << n_singular;
                                  progress.update(face_counter, update_oss.str());
                              }

                              std::map<int, pos_type> Lpoints;
                              std::vector<pos_type> epoints;
                              std::map<int, std::pair<int, int>> segments; // face to edge pair
                              std::map<int, int> boundary_segments;        // face to edge (edge - L point pair)
                              bool pathological = false;

                              std::array<edge_index_t, 4> the_edge_indices = faceid_to_edgeids(the_face_index);
                              std::vector<edge_index_t> active_edges;
                              std::vector<size_type> crease_point_indices;
                              for (int e = 0; e < 4; ++e)
                              {
                                  auto iter = edge_to_crease_point_index.find(the_edge_indices[e]);
                                  if (iter != edge_to_crease_point_index.end())
                                  {
                                      active_edges.push_back(the_edge_indices[e]);
                                      crease_point_indices.push_back(iter->second);
                                  }
                              }

                              if (active_edges.size() % 2)
                              {
                                  crease_point_t Lpoint;
                                  if (find_Lpoint_on_face(Lpoint, the_face_index, values, hessian))
                                  {
                                      std::lock_guard<std::mutex> lock(lpoint_add_mtx);
                                      the_face.Lpoint_index = unique_crease_points.add(Lpoint);
                                      crease_point_indices.push_back(the_face.Lpoint_index);
                                      ++n_singular;
                                  }
                                  else
                                      ++n_pathological;
                              }

                              connect_crease_points(the_face, crease_point_indices, unique_crease_points, gradient, hessian, theta0);
                              ++n_faces_found;
                          }
                      });
    // all active faces have now been processed
    std::cout << "\n\nFace processing complete.\n";
    std::cout << progress << '\n';
    std::cout << "  * found: " << n_faces_found << '\n';
    std::cout << "  * singular: " << n_singular << '\n';
    std::cout << "  * pathological: " << n_pathological << '\n';

    // export all_faces

    /*
        Algorithm

        3: process all voxels with active faces and mesh their segments

    */

    std::vector<voxel_info_t> active_voxels;
    unique_voxels(active_voxels, all_faces, all_voxel_indices);
    // information about voxels
    std::vector<triangle_index_t> all_triangles;
    std::vector<std::vector<size_type>> all_open_cycles;
    std::mutex add_triangles_mtx;
    std::mutex add_cycles_mtx;
    std::atomic<size_type> nmorepts = 0;
    n_failed = 0;
    size_type n_triangulated = 0;
    std::cout << "\nProcessing all " << active_voxels.size() << " unique active voxels\n";
    progress.begin(active_voxels.size(), "Mesh active voxels", 10000, "done: 0, solved: 0, triangles: 0");
    tbb::parallel_for(tbb::blocked_range<int>(0, active_voxels.size()),
                      [&](tbb::blocked_range<int> r)
                      {
                          for (int n = r.begin(); n < r.end(); ++n)
                          {
                              voxel_info_t &the_voxel = active_voxels[n];
                              const coord_type &the_voxel_index = the_voxel.voxel_id;
                              ++voxel_counter;

                              std::unique_lock<std::mutex> lock(update_mtx, std::defer_lock);
                              if (lock.try_lock())
                              {
                                  std::ostringstream update_oss;
                                  update_oss.imbue(std::locale());
                                  update_oss << "done: " << voxel_counter
                                             << ", solved: " << n_triangulated
                                             << ", triangles: " << all_triangles.size();
                                  progress.update(voxel_counter, update_oss.str());
                              }

                              const std::vector<face_info_t> &faces = the_voxel.ridge_faces;
                              std::vector<segment_index_t> voxel_face_segments;
                              std::vector<size_type> voxel_Lpoints;
                              for (auto f : faces)
                              {
                                  voxel_face_segments.insert(voxel_face_segments.end(), f.segments.begin(), f.segments.end());
                                  if (f.Lpoint_index != invalid_size)
                                  {
                                      voxel_Lpoints.push_back(f.Lpoint_index);
                                  }
                              }

                              // if 2 boundary points were found, create a boundary segment
                              // connecting them to account for singular line forming crease
                              // boundary
                              if (voxel_Lpoints.size() == 2)
                              {
                                  voxel_face_segments.push_back(segment_index_t(voxel_Lpoints[0], voxel_Lpoints[1]));
                              }

                              std::vector<std::vector<size_type>> cycles;
                              if (!verbose)
                              {
                                  compute_cycles(cycles, voxel_face_segments, false);
                              }
                              else
                              {
                                  std::scoped_lock<std::mutex> lock(debug_mtx);
                                  compute_cycles(cycles, voxel_face_segments, true);
                              }
                              for (auto c : cycles)
                              {
                                  if (c.size() < 4 || c.front() != c.back())
                                  {
                                      std::scoped_lock<std::mutex> alock(add_cycles_mtx);
                                      all_open_cycles.push_back(c);
                                  }
                                  else if (c.size() == 4)
                                  {
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]})); // c[3] == c[0]
                                  }
                                  else if (c.size() == 5)
                                  {
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[2], c[3]})); // c[4] == c[0]
                                  }
                                  else if (c.size() == 6)
                                  {
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[2], c[4]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[3], c[4]}));
                                  }
                                  else if (c.size() == 7)
                                  {
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[2], c[4]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[3], c[4]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[4], c[5]}));
                                  }
                                  else if (c.size() == 8)
                                  {
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[2], c[6]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[3], c[4]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[4], c[6]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[4], c[5], c[6]}));
                                  }
                                  else if (c.size() == 9)
                                  {
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[2], c[7]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[3], c[5]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[3], c[4], c[5]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[5], c[6], c[7]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[5], c[7]}));
                                  }
                                  else if (c.size() == 10)
                                  {
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[3], c[4]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[4], c[5], c[6]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[6], c[7], c[8]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[4], c[6]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[2], c[6]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[6], c[8]}));
                                  }
                                  else if (c.size() == 11)
                                  {
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[3], c[4]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[4], c[5], c[6]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[6], c[7], c[8]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[2], c[9]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[4], c[9]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[4], c[6], c[9]}));
                                      the_voxel.ridge_triangles.push_back(triangle_index_t({c[6], c[8], c[9]}));
                                  }
                                  else if (c.size() > 11)
                                  {
                                      ++nmorepts;
                                      n_failed++;
                                  }
                              }
                              if (the_voxel.ridge_triangles.size() == 0)
                                  ++n_failed;
                              else
                              {
                                  std::scoped_lock<std::mutex> alock(add_triangles_mtx);
                                  all_triangles.insert(all_triangles.end(), the_voxel.ridge_triangles.begin(), the_voxel.ridge_triangles.end());
                                  ++n_triangulated;
                              }
                          }
                      }
    );
    progress.end();
    std::cout << "All voxels processed.\n";
    std::cout << progress << std::endl;
    std::cout << "There were " << nmorepts << " cycles with more than 12 points\n";

    std::cout << "done: " << active_voxels.size() << ", triangulated: " << n_triangulated << ", triangles: " << all_triangles.size() << ", failed: " << n_failed << '\n';

    if (export_data)
    {
        typedef std::filesystem::path path_type;
        std::string stem = path_type(output_name).replace_extension().string();
        export_crease_points(unique_crease_points, stem + "_crease_points.json");
        export_faces(all_faces, stem + "_faces.json");
        export_voxels(active_voxels, stem + "_voxels.json");
    }

    visualize(unique_crease_points, all_triangles, all_open_cycles, output_name, vis);

    return 0;
}
