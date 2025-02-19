// STL
#include <cmath>
#include <array>
#include <exception>
#include <algorithm>
#include <utility>
#include <fstream>
#include <regex>
#include <locale>
#include <list>
// Teem
#include <teem/nrrd.h>
// spurt
#include <math/types.hpp>
#include <data/image.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <spline/spline.h>
#include <format/filename.hpp>
// #include <graph/union_find.hpp>

#include <crease/crease_mc_cases.hpp>
#include <crease/zheng.hpp>

#include <vtk/vtk_utils.hpp>
#include <format/filename.hpp>
// #include <vtkConnectivityFilter.h>
#include <vtkCellArrayIterator.h>
#include <vtkVoxel.h>

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

typedef double scalar_type;
typedef long size_type;
typedef spurt::small_vector<scalar_type, 3> vector_type;
typedef spurt::small_matrix<scalar_type, 3, 3> matrix_type;
typedef spurt::small_vector<scalar_type, 2> face_pos_type;
typedef spurt::small_vector<scalar_type, 2> face_step_type;
typedef spurt::small_vector<long, 3> coord_type;
typedef spurt::small_vector<scalar_type, 3> pos_type;
typedef spurt::bounding_box<pos_type> bounds_type;
typedef spurt::small_vector<scalar_type, 6>  sym_type;
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
template<typename Value>
using image3 = image<size_type, scalar_type, 3, Value, 
                    //  kernels::MitchellNetravaliBC, 
                    kernels::Linear,
                    coord_type, pos_type>;

typedef image3<scalar_type> scalar_image_type;
typedef image3<vector_type> vector_image_type;
typedef image3<matrix_type> matrix_image_type;
typedef scalar_image_type::grid_type grid_type;
typedef Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> dyn_mat_type;

/*
    Ridge surface extraction after T. Schultz's 2009 TVCG paper 

    Algorithm:
    1. Identify all unique edges
    2. Check each edge in parallel for a ridge point


Observations/Comments:
- The current implementation does not handle non MC cases (except for those involving only 3 vertices)
- The number of broken voxels (according to the categorization above) is overwhelming and on par with the number of found ridge triangles.
- The spurious hyperbolae in radar images appear to correspond to ridge lines rather than ridge surfaces.
- The extraction could be modified to exclude ridge line-like feature points, e.g. with a threshold on the cross product between gradient and major eigenvector
- The current implementation does not handle exotic cases (except for those involving only 2 vertices per face
*/


template<typename Iterator>
std::vector<int> arg_sort(Iterator begin, Iterator end, bool fwd=true) {
    std::vector<int> ranks(std::distance(begin, end));
    std::iota(ranks.begin(), ranks.end(), 0);
    std::sort(ranks.begin(), ranks.end(), [&](int a, int b){
        if (fwd) return *(begin+a) < *(begin+b);
        else return *(begin+a) > *(begin+b);
    });
    return ranks;
}

void add_voxel_faces(VTK_SMART(vtkFloatArray) coords, 
                     VTK_SMART(vtkCellArray) cells,
                     const pos_type& min, const pos_type& max)
{
    using namespace spurt::marching_cubes;
    pos_type bounds[] = { min, max };
    long id0 = coords->GetNumberOfTuples();
    for (int v=0; v<8; ++v) {
        auto vertex = vertices[v];
        coords->InsertNextTuple3(bounds[vertex[0]][0], 
                                 bounds[vertex[1]][1], 
                                 bounds[vertex[2]][2]);

        std::cout << "Point " << v << " has coordinates [" << bounds[vertex[0]][0] << ", " << bounds[vertex[1]][1] << ", " << bounds[vertex[2]][2] << "]\n"; 
    }
    // triangulate faces
    for (int f=0; f<6; ++f) {
        auto face = canonical_face_indices[f];
        cells->InsertNextCell(3);
        cells->InsertCellPoint(id0 + face[0]);
        cells->InsertCellPoint(id0 + face[1]);
        cells->InsertCellPoint(id0 + face[2]);

        std::cout << "Adding triangle with indices [" << id0 + face[0] << ", " << id0 + face[1] << ", " << id0 + face[2] << "]\n";

        cells->InsertNextCell(3);
        cells->InsertCellPoint(id0 + face[0]);
        cells->InsertCellPoint(id0 + face[2]);
        cells->InsertCellPoint(id0 + face[3]);    

        std::cout << "Adding triangle with indices [" << id0 + face[0] << ", " << id0 + face[2] << ", " << id0 + face[3] << "]\n";
    }
}

template<typename Value>
image3<Value> import_nrrd(const std::string& filename)
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
    typedef typename image3<Value>::base_type raster_type;
    raster_type raster = spurt::nrrd_utils::to_raster<size_type, scalar_type, 3, Value>(nin, std::is_scalar<Value>::value);
    std::cout << "Raster created\n";
    return image3<Value>(raster);
}

template<typename T>
T sign(const T& value) {
    if (value >= 0) return T(1);
    else return T(-1);
}

int distance(const coord_type& c0, const coord_type& c1) {
    return linf_norm(c1-c0);
}

template<typename T, size_t N, typename Order = std::less<T> >
struct Connectivity : public std::array<T, N>
{
    typedef T value_type;
    typedef std::array<T, N> base_type;
    typedef Order less_type;
    typedef Connectivity<T, N, Order> self_type;
    Connectivity(const T& c0=T(invalid_scalar), const T& c1=T(invalid_scalar), const T& c2=T(invalid_scalar), const T c3=T(invalid_scalar)) 
        : base_type() {
        T cs[] = { c0, c1, c2, c3 };
        for (size_t i=0; i<N; ++i) (*this)[i] = cs[i];
        std::sort(base_type::begin(), base_type::end(), less_type());
    }

    bool operator<(const self_type& other) const {
        less_type less;
        for (size_t i=0; i<N; ++i) {
            if (less((*this)[i], other[i])) return true;
            else if (less(other[i], (*this)[i])) return false;
        }
        return false;
    }
};

template<typename T, size_t N, typename Order>
std::ostream& operator<<(std::ostream& os, const Connectivity<T, N, Order>& c) {
    os << "[ ";
    std::copy(c.begin(), c.end(), std::ostream_iterator<T>(std::cout, " "));
    os << " ]";
    return os;
}

typedef Connectivity<coord_type, 2, spurt::lexicographical_order> edge_index_t;
typedef Connectivity<coord_type, 4, spurt::lexicographical_order> face_index_t;
typedef Connectivity<size_type, 2> segment_index_t;

const edge_index_t invalid_edge_index = 
    edge_index_t({invalid_coord, invalid_coord});
const face_index_t invalid_face_index = 
    face_index_t({invalid_coord, invalid_coord, invalid_coord, invalid_coord});
const segment_index_t invalid_segment_index = 
    segment_index_t({invalid_size, invalid_size});

std::array<edge_index_t, 4> faceid_to_edgeids(const face_index_t& f) {
    std::array<edge_index_t, 4> r;
    // f { 0, 3, 1, 2 }
    // actual order: f[0] f[2] f[3] f[1]
    r[0] = edge_index_t(f[0], f[2]);
    r[1] = edge_index_t(f[2], f[3]);
    r[2] = edge_index_t(f[3], f[1]);
    r[3] = edge_index_t(f[1], f[0]);
    return r;
}

std::vector<coord_type> faceid_to_voxelids(const face_index_t& f, 
                                           const coord_type& res) {
    coord_type diag = f[3] - f[0];
    int external_dim = -1;
    for (int d=0; d<3; ++d) {
        if (diag[d] == 0) {
            external_dim = d;
            break;
        }
    }
    if (external_dim == -1) 
        throw std::runtime_error("invalid face computation");

    std::cout << "external dim for " << f << " is " << external_dim << '\n';

    std::vector<coord_type> voxels;
    vector_type e_ext(0);
    e_ext[external_dim] = 1;
    if (f[0][external_dim] > 0) voxels.push_back(f[0]-e_ext);
    if (f[0][external_dim] < res[external_dim]-1) voxels.push_back(f[0]);
    return voxels;
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

face_index_t face_local_to_global(const coord_type& voxel, int localfaceid) {
    std::array<coord_type, 4> cd;
    for (int i=0; i<4; ++i) {
        int vid = spurt::marching_cubes::canonical_face_indices[localfaceid][i];
        int* c = spurt::marching_cubes::vertices[vid];
        coord_type delta(c[0], c[1], c[2]);
        cd[i] = voxel + delta;
    }
    return face_index_t(cd[0], cd[1], cd[2], cd[3]);
}

face_index_t face_edge_neighbor_to_global(const coord_type& voxel, 
                                    int face, int edge_id) {
    int f = spurt::marching_cubes::face_edge_neighbor(face, edge_id);
    return face_local_to_global(voxel, f);
}

namespace peikert {
    // Peikert-Sadlo Eurovis 2008 method
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

    std::pair<pos_type, scalar_type>
    project(const vector_type& g, const matrix_type& h) {
    vector_type evals;
    matrix_type evecs;
    sym_eigensystem(evals, evecs, h);
    pos_type coords = spurt::abs(transpose(evecs)*g);
    return std::make_pair(coords, evals[2]);
    }

    double determinant(const vector_type& g, const matrix_type& H, bool     normalize=false) {
    matrix_type A;
    A.column(0) = g;
    A.column(1) = H * g;
    A.column(2) = H * A.column(1);

    if (normalize) {
        A.column(0) /= g.norm();
        A.column(1) /= spurt::norm(A.column(1));
        A.column(2) /= spurt::norm(A.column(2));
    }
    return spurt::determinant(A);
    };
} // Peikert-Sadlo method stuff ends

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

scalar_type l3(const vector_type& evals, scalar_type theta)
{
    scalar_type delta = delta23(evals);
    if (delta > theta) return 0;
    else {
        // std::cout << "delta underflow\n";
        return (1. - delta/theta)*(1. - delta/theta);
    }
}

scalar_type dl3(const vector_type& evals, const vector_type& devals, 
                double theta)
{
    // l3' = (1 - (lambda2-lambda3)/theta)^2
    // dl3' = 2*d(1 - (lambda2-lambda3)/theta)*(1 - (lambda2-lambda3)/theta)
    // dl3' = 2*(-(dlambda2-dlambda3)/theta)*(1 - (lambda2-lambda3)/theta)
    scalar_type _l3 = l3(evals, theta);
    if (_l3 == 0) return 0;
    else
    {
        return -2/theta * (devals[1] - devals[2])*(1 - delta23(evals)/theta);
    }
}

matrix_type tensor_T(const matrix_type& H, scalar_type theta)
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
                     scalar_type theta)
{
    // h = Tg - g
    matrix_type T = tensor_T(H, theta);
    return T*g - g;
}

vector_type vector_dh(const vector_type& g, const matrix_type& H, 
                      const vector_type& dg, const matrix_type& dH, 
                      scalar_type theta)
{
    // ∇h = ∇Tg + T∇g − ∇g
    matrix_type T = tensor_T(H, theta);
    matrix_type dT = tensor_dT(dH, H, theta);
    return dT*g + T*dg - dg;
}

vector_type ridge_normal(const pos_type& point,
                         const vector_type& g, const matrix_type& H, 
                         const vector_image_type& gradient,
                         const matrix_image_type& hessian,
                         scalar_type theta,
                         bool verbose = false)
{
    if (verbose) std::cout << "value of h at ridge point: " << vector_h(g, H, theta) << '\n';
    auto H_prime = hessian.derivative(point); // 3rd order tensor
    matrix_type h_prime;
    h_prime.column(0) = vector_dh(g, H, H.column(0), H_prime.layer(0), theta);
    h_prime.column(1) = vector_dh(g, H, H.column(1), H_prime.layer(1), theta);
    h_prime.column(2) = vector_dh(g, H, H.column(2), H_prime.layer(2), theta);

    // h' = outer(n,n);
    vector_type hevals;
    matrix_type hevecs;
    sym_eigensystem(hevals, hevecs, h_prime);
    if (verbose) {
        std::cout << "eigenvalues of h_prime: " << hevals << '\n';
        std::cout << "eigenvectors of h_prime: " << hevecs << '\n';
    }
    vector_type hevals_plus = spurt::abs(hevals);
    if (verbose) std::cout << "Return eigenvector #" << spurt::argmax(hevals_plus) << '\n';
    return hevecs.column(spurt::argmax(hevals_plus));
}

scalar_type edge_linearity_value(const matrix_type& T1, const matrix_type& T2) 
{
    return trace(transpose(T1)*T2) 
        /sqrt(trace(transpose(T1)*T1))
        /sqrt(trace(transpose(T2)*T2));
}

scalar_type theta_threshold(const matrix_image_type& hessian, double eps=0.005) {
    vector_type evals;
    matrix_type evecs;
    std::vector<scalar_type> _l23;
    for (const matrix_type& H : hessian) {
        sym_eigensystem(evals, evecs, H);
        _l23.push_back(evals[1] - evals[2]);
    }
    auto iters = std::minmax_element(_l23.begin(), _l23.end());
    return eps*(*iters.second - *iters.first);
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

// void save_sampled_voxel(const std::string& filename, coord_type& voxel_id, 
//                         const scalar_image_type& values,
//                         const vector_image_type& gradient, 
//                         const matrix_image_type& hessian,
//                         const scalar_type theta,
//                         const size_t resolution=101)
// {
//     std::vector<double> det_values(5 * resolution * resolution * resolution);
//     double dh = 1./static_cast<double>(resolution-1);
//     for (int k = 0; k < resolution; ++k)
//     {
//         double z = k * dh;
//         for (int j = 0; j < resolution; ++j)
//         {
//             double y = j * dh;
//             for (int i = 0; i < resolution; ++i)
//             {
//                 double x = i * dh;
//                 vec3 g = gradient.value_in_voxel(voxel_id, pos_type(x,y,z));
//                 mat3 H = hessian.value_in_voxel(voxel_id, pos_type(x,y,z));
//                 det_values[3 * (i + resolution * (j + resolution * k))] = determinant(g, H);
//                 det_values[3 * (i + resolution * (j + resolution * k)) + 1] = (project(g, H).first)[2];
//                 scalar_type v = values.value_in_voxel(voxel_id, pos_type(x,y,z));
//                 det_values[3 * (i + resolution * (j + resolution * k)) + 2] = v;
//             }
//         }
//     }
//     size_t sizes[4] = {3, resolution, resolution, resolution};
//     std::ostringstream os;
//     os << filename << "_voxel_" << voxel_id[0] << "_" << voxel_id[1] << "_" << voxel_id[2] << ".nrrd";
//     spurt::nrrd_utils::writeNrrd((void *)&det_values[0], os.str(), nrrdTypeDouble, 4, sizes);
// }

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
process_edge(const edge_index_t &e,
             const scalar_image_type &data,
             const vector_image_type &gradient,
             const matrix_image_type &hessian,
             scalar_type theta,
             int verbose = 0,
             bool improve=false)
{
    vector_type h0 = vector_h(gradient(e[0]), hessian(e[0]), theta);
    vector_type h1 = vector_h(gradient(e[1]), hessian(e[1]), theta);
    if (inner(h0, h1) < 0) 
    {
        scalar_type n0 = h0.norm();
        scalar_type n1 = h1.norm();
        // model zero crossing as: (1-u) n0 - u n1 = 0 <=> u (n0 + n1) = n0
        scalar_type u = n0/(n1+n0);
        if (verbose >= 2) {
            std::cout << "zero crossing detected: h0=" << h0 << ", h1=" << h1 << ", u=" << u << ", n0=" << n0 << ", n1=" << n1 << '\n';
            auto T1 = tensor_T(hessian(e[0]), theta);
            auto T2 = tensor_T(hessian(e[1]), theta);
            // std::cout << "T1=" << T1 << '\n';
            // std::cout << "T2=" << T2 << '\n';
            std::cout << "edge linearity test: " << edge_linearity_value(T1, T2) << '\n';
        }
        if (improve) {
            pos_type p;
            pos_type p0 = gradient.grid()(e[0]);
            pos_type p1 = gradient.grid()(e[1]);
            scalar_type u0 = 0;
            scalar_type u1 = 1;
            int maxiter = 20;
            bool failed = true;
            for (int i=0; i<maxiter; ++i) {
                p = (1-u)*p0 + u*p1;
                vector_type h2 = vector_h(gradient.value(p), hessian.value(p), theta);
                scalar_type n2 = h2.norm();
                // std::cout << "current h value at u=" << u << " is " << h2 << ", norm =" << n2 << '\n';
                if (n2 < 1.0e-6) {
                    // std::cout << "solution found.\n";
                    failed = false;
                    break;
                }
                std::array<scalar_type, 3> ns({n0, n1, n2});
                auto ranks = arg_sort(ns.begin(), ns.end());
                if (ranks[2] == 2) {
                    // std::cout << "linear estimate increased h norm\n";
                    break;
                }
                else if (ranks[1] == 2) {
                    if (ranks[0] == 0) {
                        u1 = u;
                        n1 = n2;
                        u = 0.5*(u0+u);
                    }
                    else {
                        u0 = u;
                        n0 = n2;
                        u = 0.5*(u+u1);
                    }
                }
                else {
                    if (ranks[0] == 1) {
                        u1 = u;
                        n1 = n2;
                        u = 0.5*(u0+u);
                    }
                    else {
                        u0 = u;
                        n0 = n2;
                        u = 0.5*(u+u1);
                    }
                }
            }
            if (failed) {
                int nsteps = 100;
                std::vector<scalar_type> hnorms(nsteps+1);
                hnorms[0] = h0.norm();
                hnorms[nsteps] = h1.norm();
                for (int i=1; i<nsteps; ++i) {
                    u = i/scalar_type(nsteps);
                    pos_type p = (1-u)*p0 + u*p1;
                    vector_type h = vector_h(gradient.value(p), hessian.value(p), theta);
                    hnorms[i] = h.norm();
                }
                auto ranks = arg_sort(hnorms.begin(), hnorms.end());
                // std::cout << "min norm found (" << hnorms[ranks[0]] << ") found at u = " << ranks[0]/scalar_type(nsteps) << '\n';
                // std::cout << "hnorms were ";
                // for (auto anorm : hnorms) std::cout << anorm << ", ";
                // std::cout << '\n';
                u = ranks[0]/scalar_type(nsteps);
            }
        }

        return std::pair(1, u);
    }
    return std::pair(0, 0);
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

size_t nfour=0;
size_t nthree=0;
size_t nmore_than_four=0;
size_t ntwo=0;
size_t n_one=0;
int crude_meshing(const std::map<int, size_type>& found,
                   vtkCellArray* tris, 
                   vtkCellArray* edges, 
                   vtkCellArray* verts,
                   std::ostream& os = std::cout)
{
    std::vector<size_type> ids;
    for (auto it=found.begin(); it!=found.end(); ++it) {
        ids.push_back(it->second);
    }
    if (ids.size() == 3) {
        ++nthree;
        tris->InsertNextCell(3);
        tris->InsertCellPoint(ids[0]);
        tris->InsertCellPoint(ids[1]);
        tris->InsertCellPoint(ids[2]);
        return one_triangle;
    }
    else if (ids.size() > 3) {
        if (ids.size() == 4) ++nfour;
        else ++nmore_than_four;
        // Calculate edge case number
        int edge_case = 0;
        for (auto iter=found.begin(); iter!=found.end(); ++iter) 
        {
            edge_case += 1 << iter->first;
        }
        os << "edge_case = " << edge_case << '\n';
        int triangle_case = spurt::marching_cubes::edge_code_to_case_id[edge_case];
        os << "triangle_case = " << triangle_case << '\n';

        if (triangle_case == -1) // invalid
        {
            for (int i=0; i<ids.size()-1; ++i) {
                for (int j=i+1; j<ids.size(); ++j) {
                    edges->InsertNextCell(2);
                    edges->InsertCellPoint(ids[i]);
                    edges->InsertCellPoint(ids[j]);
                }
            } 
            os << "the edges do not match a valid MC case...\n Giving up. (Case NoMC)\n";
            return invalid_mc_case;
        }
        else {
            auto indices = spurt::marching_cubes::triTable[triangle_case];
            for (int i=0; i<15 && indices[i]!=-1; i+=3) 
            {
                tris->InsertNextCell(3);
                tris->InsertCellPoint(found.at(indices[i]));               tris->InsertCellPoint(found.at(indices[i+1]));
                tris->InsertCellPoint(found.at(indices[i+2]));
            }
            os << "A valid MC case was found\n";
            return valid_mc_case;
        }
    }
    else if (ids.size() == 2) {
        ++ntwo;
        edges->InsertNextCell(2);
        edges->InsertCellPoint(ids[0]);
        edges->InsertCellPoint(ids[1]);
        return not_enough_edges;
    }
    else if (ids.size() == 1) {
        ++n_one;
        verts->InsertNextCell(1);
        verts->InsertCellPoint(ids[0]);
        return not_enough_edges;
    }
    else {
        return not_enough_edges;
    }
}

// voxel ridge information
struct CreasePoint {
    pos_type position;
    edge_index_t edge_index; // valid if on an edge
    face_index_t face_index; // valid if L point on face
    scalar_type value;
    scalar_type strength;
    bool singular;

    CreasePoint(const pos_type& p=invalid_pos, 
                const edge_index_t& e=invalid_edge_index,
                const face_index_t& f=invalid_face_index,
                const scalar_type& v=invalid_scalar,
                const scalar_type& str=invalid_scalar,
                bool s=false)
        : position(p), edge_index(e), face_index(f), value(v), 
          strength(str), singular(s) {} 
};
typedef CreasePoint crease_point_t;

std::ostream& operator<<(std::ostream& os, const crease_point_t& cp) {
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

struct FaceInfo {
    face_index_t face_id;
    size_type Lpoint_index;
    std::vector<segment_index_t> segments;

    FaceInfo(const face_index_t& f=invalid_face_index, size_type lptid=invalid_size)
        : face_id(f), Lpoint_index(lptid), segments() {}
};
typedef FaceInfo face_info_t;

typedef std::array<size_type, 3> triangle_index_t;
struct VoxelInfo {
    coord_type voxel_id;
    std::vector<face_info_t> ridge_faces;
    std::vector<triangle_index_t> ridge_triangles;

    VoxelInfo(const coord_type& v = invalid_coord) 
        : voxel_id(v), ridge_faces(), ridge_triangles() {}
};
typedef VoxelInfo voxel_info_t;

// all the unique edges of a given set of voxel coordinates
void unique_edge_indices(std::vector<edge_index_t>& edge_indices, 
                         const std::vector<coord_type>& voxel_indices) {
    const coord_type ex = coord_type(1,0,0);
    const coord_type ey = coord_type(0,1,0);
    const coord_type ez = coord_type(0,0,1);
    coord_type basis[] = { ex, ey, ez };

    std::set<coord_type, spurt::lexicographical_order> unique_edge_start_indices;
    for (int dim=0; dim<3; ++dim) {
        const coord_type& b0 = basis[dim];
        const coord_type& b1 = basis[(dim+1)%3];
        const coord_type& b2 = basis[(dim+2)%3];
        unique_edge_start_indices.clear();
        for (const coord_type& vc : voxel_indices) {
            unique_edge_start_indices.insert(vc);
            unique_edge_start_indices.insert(vc+b1);
            unique_edge_start_indices.insert(vc+b2);
            unique_edge_start_indices.insert(vc+b1+b2);
        }
        for (const coord_type& pc : unique_edge_start_indices) {
            edge_indices.push_back(edge_index_t(pc,pc+b0));
        }
    }
}

// edge to varying dimension
// prerequisite: e is valid edge index:
// first and second index are distinct and differ only by a single coordinate
int edge_to_dim(const edge_index_t& e) {
    const coord_type ex = coord_type(1,0,0);
    const coord_type ey = coord_type(0,1,0);
    const coord_type ez = coord_type(0,0,1);

    coord_type step = spurt::abs(e[1]-e[0]);
    if (spurt::all(step == ex)) return 0;
    else if (spurt::all(step == ey)) return 1;
    else if (spurt::all(step == ez)) return 2;
    else throw std::runtime_error("invalid edge index in edge_to_dim");
}

// all the unique faces adjacents to the given set of edges
void unique_faces(std::vector<face_info_t>& faces, 
                  const std::map<edge_index_t, size_type>& edge_index_to_crease_point_index,
                  const std::vector<coord_type>& valid_voxel_indices)
{
    const coord_type ex = coord_type(1,0,0);
    const coord_type ey = coord_type(0,1,0);
    const coord_type ez = coord_type(0,0,1);
    coord_type basis[] = { ex, ey, ez };

    std::cout << "edges in input are: \n";
    for (auto e : edge_index_to_crease_point_index) {
        std::cout << e.first << '\n';
    }

    coord_type dummy_res(10000000, 10000000, 10000000);

    std::set<coord_type, spurt::lexicographical_order> valid_voxels;
    valid_voxels.insert(valid_voxel_indices.begin(), valid_voxel_indices.end());

    std::set<face_index_t> _unique_face_indices;
    for (auto v : edge_index_to_crease_point_index) {
        const edge_index_t& e = v.first;
        std::cout << "currently processing edge " << e << '\n';
        int dim = edge_to_dim(e);
        const coord_type& a = e[0];
        const coord_type& b = e[1];
        for (int d=0; d<3; ++d) {
            if (d==dim) continue;
            face_index_t face_id0(a, b, b+basis[d], a+basis[d]);
            face_index_t face_id1(a, b, b-basis[d], a-basis[d]);
            std::cout << "Two candidate face indices are:\n"
                << face_id0 << '\n'
                << face_id1 << '\n'; 
            std::vector<coord_type> voxels0 = faceid_to_voxelids(face_id0, dummy_res);
            std::vector<coord_type> voxels1 = faceid_to_voxelids(face_id1, dummy_res);
            for (auto v0 : voxels0) {
                std::cout << "Voxel for face " << face_id0 << " is " << v0 << '\n';
                if (valid_voxels.find(v0) != valid_voxels.end()) {
                    _unique_face_indices.insert(face_id0);
                }
            }
            for (auto v1 : voxels1) {
                std::cout << "Voxel for face " << face_id1 << " is " << v1 << '\n';
                if (valid_voxels.find(v1) != valid_voxels.end()) {
                    _unique_face_indices.insert(face_id1);
                }
            }
        }
    }

    faces.clear();
    faces.reserve(_unique_face_indices.size());
    for (auto fid : _unique_face_indices) {
        faces.push_back(face_info_t(fid));
    }
} 


// all the unique voxels adjacent to the given set of faces
void unique_voxels(std::vector<voxel_info_t>& voxels, 
                   const std::vector<face_info_t>& faces,
                   const std::vector<coord_type>& valid_voxels)
{
    std::set<coord_type, spurt::lexicographical_order> _valid_voxels;
    _valid_voxels.insert(valid_voxels.begin(), valid_voxels.end());
    std::map<coord_type, voxel_info_t, spurt::lexicographical_order> _unique_voxels;

    coord_type dummy_res(10000000, 10000000, 10000000);
    for (auto f : faces) {
        std::vector<coord_type> voxels = faceid_to_voxelids(f.face_id, dummy_res);
        for (auto v : voxels) {
            if (_valid_voxels.find(v) == _valid_voxels.end()) 
                continue;
            auto iter = _unique_voxels.find(v);
            if (iter != _unique_voxels.end()) {
                iter->second.ridge_faces.push_back(f);
            }
            else {
                _unique_voxels[v] = voxel_info_t(v);
                _unique_voxels[v].ridge_faces.push_back(f);
            }
        }
    }
    voxels.clear();
    voxels.reserve(_unique_voxels.size());
    for (auto v2i : _unique_voxels) {
        voxels.push_back(v2i.second);
    }
}

void select_voxel_indices(std::vector<coord_type>& voxel_indices, 
                          const bounds_type& bounds, 
                          const coord_type& shape)
{
    if (spurt::all(bounds.min() < bounds.max()))
    {
        voxel_indices.clear();
        const pos_type& l = bounds.min();
        const pos_type& h = bounds.max();
        coord_type low = floor(l);
        coord_type high = ceil(h);
        for (int k=low[2]; k<=high[2]; ++k) 
        {
            for (int j=low[1]; j<=high[1]; ++j)
            {
                for (int i=low[0]; i<=high[0]; ++i) {
                    voxel_indices.push_back(coord_type(i,j,k));
                }
            }
        }
    }
    else 
    {
        voxel_indices.resize((shape[0]-1)*(shape[1]-1)*(shape[2]-1));
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

struct value_set {
    scalar_type value;
    scalar_type strength;
    vector_type gradient;
    matrix_type hessian;
    vector_type normal;
};

struct broken_voxel {
    coord_type id;
    pos_type min, max;
    std::array<value_set, 8> vox_values;
    std::vector<value_set> edg_values;
    std::vector<pos_type> edg_points;
    std::vector<pos_type> L_points;
    std::vector<std::pair<int, int>> face_edges;
};

void export_broken_voxels(const std::vector<broken_voxel>& voxels, 
                          const std::string filename) {

    std::cout << "filename = " << filename << '\n';
    VTK_CREATE(vtkFloatArray, vox_v);
    VTK_CREATE(vtkFloatArray, vox_s);
    VTK_CREATE(vtkFloatArray, vox_g);
    VTK_CREATE(vtkFloatArray, vox_h);
    VTK_CREATE(vtkFloatArray, vox_c);
    vox_c->SetNumberOfComponents(3);
    vox_v->SetNumberOfComponents(1);
    vox_s->SetNumberOfComponents(1);
    vox_g->SetNumberOfComponents(3);
    vox_h->SetNumberOfComponents(9);
    VTK_CREATE(vtkCellArray, faces);

    VTK_CREATE(vtkFloatArray, edg_v);
    VTK_CREATE(vtkFloatArray, edg_s);
    VTK_CREATE(vtkFloatArray, edg_g);
    VTK_CREATE(vtkFloatArray, edg_h);
    VTK_CREATE(vtkFloatArray, edg_c);
    VTK_CREATE(vtkFloatArray, edg_n);
    edg_c->SetNumberOfComponents(3);
    edg_v->SetNumberOfComponents(1);
    edg_s->SetNumberOfComponents(1);
    edg_g->SetNumberOfComponents(3);
    edg_h->SetNumberOfComponents(9);
    edg_n->SetNumberOfComponents(3);
    VTK_CREATE(vtkCellArray, verts);
    VTK_CREATE(vtkCellArray, segments);

    for ( const broken_voxel& v : voxels ) {
        // Insert voxel geometry
        add_voxel_faces(vox_c, faces, v.min, v.max);
        for (int i=0; i<8; ++i) {
            auto data = v.vox_values[i];
            vox_v->InsertNextTuple1(data.value);
            vox_s->InsertNextTuple1(data.strength);
            vox_g->InsertNextTuple3(
                data.gradient[0], data.gradient[1], data.gradient[2]
            );
            vox_h->InsertNextTuple9(
                data.hessian(0,0), data.hessian(0,1), data.hessian(0,2), 
                data.hessian(1,0), data.hessian(1,1), data.hessian(1,2),
                data.hessian(2,0), data.hessian(2,1), data.hessian(2,2)
            );
        }

        std::map<int, int> edge_point_indexing;
        for (int i=0; i<v.edg_points.size(); ++i) {
            auto p = v.edg_points[i];
            auto d = v.edg_values[i];
            int id = edg_c->GetNumberOfTuples();
            edge_point_indexing[i] = id;
            edg_c->InsertNextTuple3(p[0], p[1], p[2]);
            edg_v->InsertNextTuple1(d.value);
            edg_s->InsertNextTuple1(d.strength);
            edg_g->InsertNextTuple3(
                d.gradient[0], d.gradient[1], d.gradient[2]
            );
            edg_n->InsertNextTuple3(d.normal[0], d.normal[1], d.normal[2]);
            edg_h->InsertNextTuple9(
                d.hessian(0,0), d.hessian(0,1), d.hessian(0,2), 
                d.hessian(1,0), d.hessian(1,1), d.hessian(1,2),
                d.hessian(2,0), d.hessian(2,1), d.hessian(2,2)
            );
            verts->InsertNextCell(1);
            verts->InsertCellPoint(id);
        }
        for (int i=0; i<v.L_points.size(); ++i) {
            auto p = v.L_points[i];
            int id = edg_c->GetNumberOfTuples();
            edg_c->InsertNextTuple3(p[0], p[1], p[2]);
            verts->InsertNextCell(1);
            verts->InsertCellPoint(id);
        }
        for (int i=0; i<v.face_edges.size(); ++i) {
            auto idpair = v.face_edges[i];
            int id0 = edge_point_indexing[idpair.first];
            int id1 = edge_point_indexing[idpair.second];
            segments->InsertNextCell(2);
            segments->InsertCellPoint(id0);
            segments->InsertCellPoint(id1);
        }
    }

    VTK_CREATE(vtkPolyData, vpoly);
    VTK_CREATE(vtkPoints, vpoints);
    vpoints->SetData(vox_c);
    vpoly->SetPoints(vpoints);
    vpoly->SetPolys(faces);
    vox_v->SetName("value");
    vox_s->SetName("ridge_strength");
    vox_g->SetName("gradient");
    vox_h->SetName("hessian");
    vpoly->GetPointData()->AddArray(vox_v);
    vpoly->GetPointData()->AddArray(vox_s);
    vpoly->GetPointData()->AddArray(vox_g);
    vpoly->GetPointData()->AddArray(vox_h);

    VTK_CREATE(vtkXMLPolyDataWriter, writer);
    writer->SetFileName((filename + "_voxel_points.vtp").c_str());
    writer->SetInputData(vpoly);
    writer->Write();

    VTK_CREATE(vtkPolyData, cpoly);
    std::cout << "cpoints contains " << edg_c->GetNumberOfTuples() << " points\n";
    VTK_CREATE(vtkPoints, cpoints);
    cpoints->SetData(edg_c);
    cpoly->SetPoints(cpoints);
    cpoly->SetVerts(verts);
    cpoly->SetLines(segments);
    edg_v->SetName("value");
    std::cout << "values contain " << edg_v->GetNumberOfTuples() << " values\n";    
    edg_s->SetName("ridge_strength");
    edg_g->SetName("gradient");
    edg_h->SetName("hessian");
    edg_n->SetName("normal");
    cpoly->GetPointData()->AddArray(edg_v);
    cpoly->GetPointData()->AddArray(edg_s);
    cpoly->GetPointData()->AddArray(edg_g);
    cpoly->GetPointData()->AddArray(edg_h);
    cpoly->GetPointData()->AddArray(edg_n);

    writer->SetFileName((filename + "_crease_points.vtp").c_str());
    writer->SetInputData(cpoly);
    writer->Write();
}

template<typename T>
struct container_wrapper : public std::vector<T> {
    typedef std::vector<T> base_type;
    typedef typename base_type::size_type idx_type;
    container_wrapper(): base_type() {};

    idx_type add(const T& val) {
        base_type::push_back(val);
        return base_type::size()-1;
    }
};

struct GraphMap {
    typedef std::array<size_type, 2> neighbor_t;
    typedef std::map<size_type, neighbor_t> map_t;
    map_t m_edges;
    
    GraphMap(const std::vector<segment_index_t>& segments) 
        : m_edges()
    {
        for ( auto s : segments) {
            for (int i=0; i<2; ++i) {
                size_type v = s[i];
                size_type w = s[(i+1)%2];
                auto iter = m_edges.find(v);
                if (iter == m_edges.end()) {
                    m_edges[v] = neighbor_t({w, invalid_size});
                }
                else {
                    iter->second[1] = w;
                }
            }
        }
    }

    std::vector<size_type> unique_ids() const {
        std::vector<size_type> r;
        for (auto kv : m_edges) {
            r.push_back(kv.first);
        }
        return r;
    }

    size_type next_id(size_type at, size_type from=invalid_size) const {
        const neighbor_t& ns = const_cast<map_t&>(m_edges)[at];
        if (from == invalid_size) return ns[0];
        else if (ns[0] == from) return ns[1];
        else return ns[0];
    }
};

void compute_cycles(std::vector< std::vector< size_type > >& cycles, 
                    const std::vector<segment_index_t>& segments,
                    bool verbose=false) {

    typedef std::vector<size_type> cycle_t;
    GraphMap gmap(segments);
    std::set<size_type> checked;

    std::vector<size_type> unique_ids = gmap.unique_ids();
    if (verbose) {
        std::cout << "input segments are:\n";
        for (auto s : segments) {
            std::cout << s[0] << " <-> " << s[1] << '\n';
        }

        std::cout << "unique ids:\n";
        std::copy(unique_ids.begin(), unique_ids.end(), std::ostream_iterator<size_type>(std::cout, " "));
        std::cout << '\n';
    }
    for (size_type id : unique_ids) {
        if (verbose) std::cout << "entering new loop at " << id << '\n';
        if (checked.find(id) != checked.end()) {
            if (verbose) std::cout << "skipping already visited id\n";
            continue;
        }
        cycles.push_back(cycle_t());
        cycles.back().push_back(id);
        if (verbose) std::cout << "new cycle started at " << id << '\n';

        size_type from = id;
        size_type at = gmap.next_id(id);
        while (at != invalid_size) {
            if (verbose) std::cout << "In loop: from=" << from << ", at=" << at << '\n';
            cycles.back().push_back(at);
            if (at == id) {
                if (verbose) std::cout << "cycle closed. done\n";
                break;
            }
            size_type next = gmap.next_id(at, from);
            from = at;
            at = next;
        }
    }
}

void visualize(const container_wrapper<crease_point_t>& cpoints, 
               const std::vector<triangle_index_t>& all_triangles, 
               const std::vector<std::vector<size_type>>& all_open_cycles,
               const std::string& fname,
               bool novis)
{
    std::cout << "There are " << cpoints.size() << " points in input\n";
    std::cout << "There are " << all_triangles.size() << " triangles in input\n";
    std::cout << "crease points\n";
    for (size_type i=0; i<cpoints.size(); ++i) {
        std::cout << i << ": " << cpoints[i] << '\n';
    }
    for (size_type i=0; i<all_triangles.size(); ++i) {
        std::cout << i << ": " << all_triangles[i] << '\n';
    }
    std::vector<pos_type> points(cpoints.size(), pos_type(0));
    std::vector<scalar_type> values(cpoints.size(), 0);
    std::vector<scalar_type> strengths(cpoints.size(), 0);

    for (size_type i=0; i<cpoints.size(); ++i) {
        points[i] = cpoints[i].position;
        values[i] = cpoints[i].value;
        strengths[i] = cpoints[i].strength;
    }
    vtkSmartPointer<vtkPolyData> pd = vtk_utils::make_points(points);

    VTK_CREATE(vtkCellArray, triangles);
    for (size_type i=0; i<all_triangles.size(); ++i) {
        const triangle_index_t& T = all_triangles[i];
        triangles->InsertNextCell(3);
        triangles->InsertCellPoint(T[0]);
        triangles->InsertCellPoint(T[1]);
        triangles->InsertCellPoint(T[2]);
    }
    pd->SetPolys(triangles);

    VTK_CREATE(vtkCellArray, lines);
    for (size_type i=0; i<all_open_cycles.size(); ++i) {
        const std::vector<size_type>& aline = all_open_cycles[i];
        lines->InsertNextCell(aline.size());
        for (size_type j : aline) lines->InsertCellPoint(j);
    }
    pd->SetLines(lines);

    pd = vtk_utils::add_scalars(pd, values, true, "values", false);
    pd = vtk_utils::add_scalars(pd, strengths, true, "ridge_strength", true);

    pd->PrintSelf(std::cout, vtkIndent(1));

    VTK_CREATE(vtkXMLPolyDataWriter, writer);
    writer->SetInputData(pd);
    auto name = filename::replace_extension(fname, "vtp");
    writer->SetFileName(name.c_str());
    writer->Write();

    if (novis) return;

    VTK_CREATE(vtkPolyDataMapper, mesh_mapper);
    mesh_mapper->SetInputData(pd);
    mesh_mapper->ScalarVisibilityOn();
    VTK_CREATE(vtkActor, mesh_actor);
    mesh_actor->SetMapper(mesh_mapper);
    mesh_actor->GetProperty()->SetLineWidth(2);
    mesh_actor->GetProperty()->VertexVisibilityOn();
    // mesh_actor->GetProperty()->EdgeVisibilityOn();

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

int main(int argc, const char* argv[]) 
{
    std::string scalar_name, gradient_name, hessian_name, output_name, strength_name;
    scalar_type minval, minstr, eps, mind;
    std::string minval_str, minstr_str;
    size_type minsize;
    int res, niter;
    int verbose;
    bool novis;
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
        parser.add_value("strength", strength_name, "Ridge strength volume", required_group);
        parser.add_value("minval", minval_str, "0", "Min scalar value", optional_group);
        parser.add_value("minstr", minstr_str, "0", "Min ridge strength (<=0)", optional_group);
        // parser.add_value("minsize", minsize, 10, "Min size of connected component on output");
        // parser.add_value("eps", eps, 1.0e-9, "Numerical precision", optional_group);
        parser.add_value("verbose", verbose, 0, "Verbose output", optional_group);
        parser.add_flag("novis", novis, "Skip visualization and exit", optional_group);
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

    if (minval_str.back() == '%') {
        int pct = std::stoi(minval_str.substr(0, minval_str.size()-1));
        std::vector<double> vals(values.begin(), values.end());
        std::sort(vals.begin(), vals.end());
        minval = vals[std::floor(pct*vals.size()/100.)];
    }
    else minval = std::stof(minval_str);

    if (minstr_str.back() == '%') {
        int pct = 100-std::stoi(minstr_str.substr(0, minstr_str.size()-1));
        std::vector<double> vals(strength.begin(), strength.end());
        std::sort(vals.begin(), vals.end());
        minstr = vals[std::floor(pct*vals.size()/100.)];
    }
    else minstr = std::stof(minstr_str);

    std::cout << "Filtering thresholds set to: value: " << minval << ", ridge strength: " << minstr << '\n';

    auto shape = values.grid().resolution();
    std::vector< coord_type > all_voxel_indices;
    if (spurt::any(voxel_id != invalid_coord))
    {
        all_voxel_indices.clear();
        all_voxel_indices.push_back(voxel_id);
        verbose = 2;
        // save_sampled_voxel(output_name, voxel_id, values, gradient, hessian, theta0, 101);
    }
    else
    {
        std::cout << "selected voxels..." << std::flush;
        select_voxel_indices(all_voxel_indices, bounds, shape);
        std::cout << " done\n";
    }
    std::vector<edge_index_t> all_edge_indices;
    unique_edge_indices(all_edge_indices, all_voxel_indices);

    if (verbose) std::cout << "There are " << all_voxel_indices.size() << " voxels in input for a total of " << all_edge_indices.size() << " unique edges\n";

    /*
    
        Algorithm:
        
        1: process all edges of all voxels (skipping redundant ones),
           looking for crease points
    */

    // edge -> globalcrease point id
    std::map< edge_index_t, size_type > edge_to_crease_point_index;

    // information about crease points
    container_wrapper<crease_point_t> unique_crease_points;

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
    int noddfaces = 0;

    spurt::ProgressDisplay progress;

    VTK_CREATE(vtkCellArray, vtk_cells);
    VTK_CREATE(vtkCellArray, vtk_lines);
    VTK_CREATE(vtkCellArray, vtk_verts);

    std::vector<broken_voxel> broken_voxels;
    std::vector<edge_index_t> double_edges;
    std::vector<coord_type> to_subdivide;
    std::vector<int> voxel_to_edge;
    std::map<coord_type, std::vector<int>, spurt::lexicographical_order> voxel_counter;

    srand48(130819751900);
    std::cout << "Processing all " << all_edge_indices.size() << " unique edges\n";
    if (verbose > 1) {
        std::cout << "These edges are:\n";
        for (int i=0; i<all_edge_indices.size(); ++i) {
            std::cout << all_edge_indices[i] << '\n';
        }
    }
    progress.begin(all_edge_indices.size(), "Extract ridge points", 10000, "done: 0, found: 0, weak: 0");
    for (int n=0; n<all_edge_indices.size(); ++n) 
    {
        const edge_index_t& the_edge_index = all_edge_indices[n];
        std::ostringstream update_oss;
        update_oss.imbue(std::locale());
        update_oss << "done: " << n
            << ", found: "  << unique_crease_points.size() 
            << ", weak: "  << nweak;
        progress.update(n, update_oss.str());

        // check if current edge satisfies threshold requirements
        const coord_type& the_coord_0 = the_edge_index[0];
        const coord_type& the_coord_1 = the_edge_index[1];
        if (values(the_coord_0) < minval || values(the_coord_1) < minval) {
            nnone++;
            nweak++;
            continue;
        }
        if (evmin(hessian(the_coord_0)).first > minstr || 
            evmin(hessian(the_coord_1)).first > minstr) {
            nweak++;
            nnone++;
            continue;
        }
        const pos_type& the_point_0 = values.grid()(the_coord_0);
        const pos_type& the_point_1 = values.grid()(the_coord_1);

        std::ostringstream log;
        if (verbose >= 2) 
            std::cout << "\n\nedge #" << the_edge_index << " is being processed" << '\n';
        auto result = process_edge(the_edge_index, values, gradient, hessian, theta0, verbose, true);
        if (result.first == 1)
        {        
            if (verbose) 
                std::cout << "process_edge returned 1 solution\n";
            scalar_type u = result.second;
            pos_type p = (1. - u) * the_point_0 + u * the_point_1;
            auto v = evaluate(p, values, hessian);
            edge_to_crease_point_index[the_edge_index] = unique_crease_points.add(crease_point_t(p, the_edge_index, invalid_face_index, v.first, v.second, false));
        }
        else
        {
            if (verbose >= 2)
                std::cout << "no ridge point\n";
            nnone++;
        }
        nprocessed++;
    }
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

    // now processing active faces
    std::vector<face_info_t> all_faces;
    unique_faces(all_faces, edge_to_crease_point_index, all_voxel_indices);
    // information about faces
    // std::map< face_index_t, size_type > face_to_Lpoint;
    // std::map< face_index_t, segment_index_t > face_to_segment;
    // std::map< face_index_t, segment_index_t > face_to_boundary_segment;

    if (verbose > 1) {
        std::cout << "The unique faces are:\n";
        for (int i=0; i<all_faces.size(); ++i) {
            std::cout << all_faces[i].face_id << '\n';
        }
    }

    size_type n_pathological = 0;
    size_type n_singular = 0;
    size_type n_failed = 0;
    std::cout << "\nProcessing all " << all_faces.size() << " unique faces\n";
    progress.begin(all_faces.size(), "Extract ridge points", 10000, "done: 0, found: 0, pathological: 0, singular: 0, failed: 0");
    size_type n_faces_found = 0;
    for (int n=0; n<all_faces.size(); ++n) 
    {
        face_info_t& the_face = all_faces[n];
        const face_index_t& the_face_index = the_face.face_id;

        std::ostringstream update_oss;
        update_oss.imbue(std::locale());
        update_oss << "done: " << n
            << ", found: "  << n_faces_found
            << ", pathological: "  << n_pathological
            << ", singular: " << n_singular
            << ", failed: " << n_failed;
        progress.update(n, update_oss.str());

        std::map<int, pos_type> Lpoints;
        std::vector<pos_type> epoints;
        std::map<int, std::pair<int, int> > segments; // face to edge pair
        std::map<int, int> boundary_segments; // face to edge (edge - L point pair)
        bool pathological = false;

        std::array<edge_index_t, 4> the_edge_indices = faceid_to_edgeids(the_face_index);
        std::vector<edge_index_t> active_edges;
        std::vector<size_type> crease_point_indices;
        for (int e=0; e<4; ++e) {
            auto iter = edge_to_crease_point_index.find(the_edge_indices[e]);
            if (iter != edge_to_crease_point_index.end()) {
                active_edges.push_back(the_edge_indices[e]);
                crease_point_indices.push_back(iter->second);
            }
        }

        if (active_edges.size() == 1) {
            n_failed++;
            continue;
        }
        else if (active_edges.size() == 2) {
            the_face.segments.push_back(segment_index_t(crease_point_indices[0], crease_point_indices[1]));
            ++n_faces_found;
        }
        else if (active_edges.size() == 3)
        {
            pathological = true;
            n_pathological++;
            if (verbose)
                std::cout << "Face " << the_face_index << " contains 3 crease points\n";
            bool Lfound = false;
            pos_type Lpt;
            pos_type orig = values.grid()(the_face_index[0]);
            vector_type e0 = values.grid()(the_face_index[2]) - orig;
            vector_type e1 = values.grid()(the_face_index[1]) - orig;
            spurt::zheng::Face<scalar_type> f(orig, e0, e1);
            if (spurt::zheng::findLPoint(Lpt, f, hessian)) {
                if (verbose) std::cout << "L point found at " << Lpt << '\n';
                Lfound = true;
                the_face.Lpoint_index = unique_crease_points.add(crease_point_t(Lpt, invalid_edge_index, the_face_index, 0, 0, true));
                ++n_singular;
            }
            else {
                if (verbose) std::cout << "L point not found\n";
            }

            // compute ridge normals at ridge points
            pos_type p[3];
            vector_type g[3];
            matrix_type h[3];
            vector_type n[3];
            for (int l=0; l<3; ++l) {
                p[l] = unique_crease_points[crease_point_indices[l]].position;
                g[l] = gradient.value(p[l]);
                h[l] = hessian.value(p[l]);
                n[l] = ridge_normal(p[l], g[l], h[l], gradient, hessian, theta0);
            }
            std::array<scalar_type, 3> scores;
            for (int l=0; l<3; ++l) {
                int ll = (l+1)%3;
                vector_type seg = p[(l+1)%3]-p[l];
                seg /= spurt::norm(seg);
                scores[l] = std::abs(spurt::inner(seg, n[l])) + 
                            std::abs(spurt::inner(seg, n[ll]));
            }
            int l = std::distance(
                    scores.begin(), 
                    std::min_element(scores.begin(), scores.end())
            );
            int ll = (l+1)%3;
            std::cout << "The best connection found is edge " << active_edges[l] << " - " << active_edges[ll] << '\n';
            the_face.segments.push_back(segment_index_t(crease_point_indices[l], crease_point_indices[ll]));
            if (Lfound) {
                int lll = (l+2)%3;
                the_face.segments.push_back(segment_index_t(crease_point_indices[lll], the_face.Lpoint_index));
            }
        }
        else if (active_edges.size() == 4) {
            std::cout << "Face " << the_face_index << " contains 4 crease points\n";
            n_failed++;
        }
        else {
            std::cout << "Face " << the_face_index << " contains " << active_edges.size() << " crease points\n";
            n_failed++;
        }
    }
    // all active faces have now been processed
    std::cout << "\n\nFace processing complete.\n";
    std::cout << progress << '\n';
    std::cout << "  * found: " << n_faces_found << '\n';
    std::cout << "  * pathological: " << n_pathological << '\n';
    std::cout << "  * failed: " << n_failed << '\n';
    
    /*
        Algorithm

        3: process all voxels with active faces and mesh their segments
    
    */

    std::vector<voxel_info_t> active_voxels;
    unique_voxels(active_voxels, all_faces, all_voxel_indices);
    // information about voxels
    std::vector< triangle_index_t > all_triangles;
    std::vector< std::vector< size_type > > all_open_cycles;
    n_failed = 0;
    size_type n_triangulated = 0;
    std::cout << "\nProcessing all " << active_voxels.size() << " unique active voxels\n";
    progress.begin(active_voxels.size(), "Mesh active voxels", 10000, "done: 0, solved: 0, boundary: 0, failed: 0, triangles: 0");
    for (int n=0; n<active_voxels.size(); ++n) 
    {
        voxel_info_t& the_voxel = active_voxels[n];
        const coord_type& the_voxel_index = the_voxel.voxel_id;

        std::ostringstream update_oss;
        update_oss.imbue(std::locale());
        update_oss << "done: " << n
            << ", solved: "  << n_triangulated
            << ", failed: " << n_failed 
            << ", triangles: " << all_triangles.size();
        progress.update(n, update_oss.str());

        const std::vector<face_info_t>& faces = the_voxel.ridge_faces;
        std::vector<segment_index_t> all_face_segments;
        std::vector<size_type> all_Lpoints;
        for (auto f : faces) {
            all_face_segments.insert(all_face_segments.end(), f.segments.begin(), f.segments.end());
            if (f.Lpoint_index != invalid_size) {
                all_Lpoints.push_back(f.Lpoint_index);
            }
        }

        // if 2 boundary points were found, create a boundary segment 
        // connecting them
        if (all_Lpoints.size() == 2) {
            all_face_segments.push_back(segment_index_t(all_Lpoints[0], all_Lpoints[1]));
        }

        bool show_all = false;
        if (all_triangles.size()>50 && all_triangles.size() < 60) 
            show_all = true;

        std::vector< std::vector< size_type > > cycles;
        compute_cycles(cycles, all_face_segments, show_all);

        if (show_all) {
            std::cout << "There are " << all_triangles.size() << " triangles so far\n";
            std::cout << "face segments are:\n";
            for (int k=0; k<all_face_segments.size(); ++k) {
                std::cout << all_face_segments[k] << '\n';
            }
            std::cout << "cycles are:\n";
            for (auto c : cycles) {
                std::copy(c.begin(), c.end(), std::ostream_iterator<size_type>(std::cout, ", "));
                std::cout << '\n';
                if (c[0] != c.back()) {
                    std::cout << "this is not a cycle\n";
                }
            }
        }

        for (auto c : cycles) {
            if (c.size() < 4 || c.front() != c.back()) {
                std::cout << "ERROR: one returned cycle is ";
                std::copy(c.begin(), c.end(), std::ostream_iterator<size_type>(std::cout, ", "));
                std::cout << '\n';
                all_open_cycles.push_back(c);
            }
            else if (c.size() == 4) {
                if (show_all) std::cout << "there are 3 unique points in cycle\n";
                the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]})); // c[3] == c[0]
            }  
            else if (c.size() == 5) {
                if (show_all) std::cout << "there are 4 unique points in cycle\n";
                the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]})); 
                the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[2], c[3]})); // c[4] == c[0]
            }
            else if (c.size() == 6) {
                if (show_all) std::cout << "there are 5 unique points in cycle\n";
                the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]})); 
                the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[2], c[4]})); 
                the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[3], c[4]})); 
            }
            else if (c.size() == 7) {
                if (show_all) std::cout << "there are 6 unique points in cycle\n";
                the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[1], c[2]})); 
                the_voxel.ridge_triangles.push_back(triangle_index_t({c[0], c[2], c[4]})); 
                the_voxel.ridge_triangles.push_back(triangle_index_t({c[2], c[3], c[4]})); 
                the_voxel.ridge_triangles.push_back(triangle_index_t({c[3], c[4], c[5]})); 
            }
            else if (c.size() > 7) {
                std::cout << "ERROR: the current cycle contains more than 6 points\n";
                std::copy(c.begin(), c.end(), std::ostream_iterator<size_type>(std::cout, ", "));
                std::cout << '\n';
                n_failed++;
            }
        }
        if (the_voxel.ridge_triangles.size() == 0) ++n_failed;
        else {
            std::cout << "There were " << the_voxel.ridge_triangles.size() << " triangles for this voxel\n";

            all_triangles.insert(all_triangles.end(), the_voxel.ridge_triangles.begin(), the_voxel.ridge_triangles.end());
            ++n_triangulated;
        }
    }
    progress.end();
    std::cout << "All voxels processed.\n";
    std::cout << progress << std::endl;
    
    std::cout << "done: " << active_voxels.size() << ", triangulated: " << n_triangulated << ", triangles: " << all_triangles.size() << ", failed: " << n_failed << '\n';

    visualize(unique_crease_points, all_triangles, all_open_cycles, output_name, novis);

    return 0;
}
