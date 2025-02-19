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

#include <crease/crease_mc_cases.hpp>

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

inline scalar_type sq(const scalar_type& a) {
    return a*a;
}

inline scalar_type sqdiff(const scalar_type& a, const scalar_type& b) {
    return sq(a) - sq(b);
}

inline vector_type dsqdiff(const scalar_type& a, const scalar_type& b,
                           const vector_type& da, const vector_type& db)
{
    /* d(a*a - b*b) = 2ada - 2bdb */
    return 2*a*da - 2*b*db;
}

inline scalar_type dsqdiffd0(const scalar_type& a, const scalar_type& b) {
    /* d(a^2 - b^2)/da = 2*a */
    return 2*a;
}

inline scalar_type sqsum(const scalar_type& a, const scalar_type& b) {
    return sq(a) + sq(b);
}

inline vector_type dsqsum(const scalar_type& a, const scalar_type& b,
                           const vector_type& da, const vector_type& db)
{
    /* d(a*a + b*b) = 2ada + 2bdb */
    return 2*a*da + 2*b*db;
}

CF_type zheng_CF(const matrix_type& H) {
    const scalar_type& x = H(0,0);
    const scalar_type& y = H(1,1);
    const scalar_type& z = H(2,2);
    const scalar_type& u = H(0,1);
    const scalar_type& v = H(0,2);
    const scalar_type& w = H(1,2);

    CF_type CF(0);

    CF[0] = 
        x * ((y*y - z*z) + (u*u - v*v)) + 
        y * ((z*z - x*x) + (w*w - u*u)) + 
        z * ((x*x - y*y) + (v*v - w*w));

    CF[1] = 
        w *(2*(w*w-x*x) - (v*v + u*u) + 2*(y*x + z*x - y*z)) + 
        u*v*(2*x - z - y);

    CF[2] = 
        v * (2*(v*v - y*y) - (u*u + w*w) + 2*(z*y + x*y - z*x)) + 
        w*u*(2*y - x - z);

    CF[3] = 
        u*(2*(u*u - z*z) - (w*w + v*v) + 2*(x*z + y*z - x*y)) + 
        v*w*(2*z - y - x);

    CF[4] = w*(v*v-u*u) + u*v*(y-z);
    CF[5] = v*(u*u-w*w) + w*u*(z-x);
    CF[6] = u*(w*w-v*v) + v*w*(x-y);

    return CF;
}

sym_type dC1dT(const matrix_type& H, int Cid) {
    const scalar_type& x = H(0,0);
    const scalar_type& y = H(1,1);
    const scalar_type& z = H(2,2);
    const scalar_type& u = H(0,1);
    const scalar_type& v = H(0,2);
    const scalar_type& w = H(1,2);

    sym_type r(0);
    switch (Cid) {
        case 0: {
            r(0) = u*u - v*v + 2*x*(z - y) + y*y - z*z;
            r(1) = w*w - u*u + 2*y*(x - z) - x*x + z*z;
            r(2) = v*v - w*w + 2*z*(y - x) - y*y + x*x;
            r(3) = 2 * u * (x - y);
            r(4) = 2 * v * (z - x);;
            r(5) = 2 * w * (y - z);
        }
        break;
        case 1: {
            r(0) = 2*u*v + 2*w*(-2*x + y + z);
            r(1) = w*(2*x - 2*z) - u*v;
            r(2) = w*(2*x - 2*y) - u*v;
            r(3) = -2*u*w + v*(2*x - y - z);
            r(4) = -2*v*w + u*(2*x - y - z);
            r(5) = -u*u - v*v + 6*w*w - 2*x*x + 2*x*y + 2*x*z - 2*y*z;
        }
        break;
        case 2: {
            r(0) = 2*v*(y-z) - u*w;
            r(1) = 2*u*w + 2*v*(x - 2*y + z);
            r(2) = -(u*w) + 2*v*(-x + y);
            r(3) = -2*u*v - w*(x - 2*y + z);
            r(4) = -u*u + 6*v*v - w*w + 2*(x - y)*(y - z);
            r(5) = -2*v*w - u*(x - 2*y + z);
        }
        break;
        case 3: { 
            r(0) = -(v*w) + 2*u*(-y + z);
            r(1) = -(v*w) + 2*u*(-x + z);
            r(2) = 2*v*w + 2*u*(x + y - 2*z);
            r(3) = 6*u*u - v*v - w*w - 2*x*y + 2*x*z + 2*y*z - 2*z*z;
            r(4) = -2*u*v - w*(x + y - 2*z);
            r(5) = -2*u*w - v*(x + y - 2*z);
        }
        break;
        case 4: {
            r(0) = 0;
            r(1) = u*v;
            r(2) = -u*v;
            r(3) = -2*u*w + v*(y-z);
            r(4) = 2*v*w + u*(y - z);
            r(5)  = v*v - u*u;
        }
        break;
        case 5: {
            r(0) = -u*w;
            r(1) = 0;
            r(2) = u*w;
            r(3) = 2*u*v + w*(-x + z);
            r(4) = u*u - w*w;
            r(5) = -2*v*w + u*(-x + z);
        }
        break;
        case 6: {
            r(0) = v*w;
            r(1) = -v*w;
            r(2) = 0;
            r(3) = w*w - v*v;
            r(4) = -2*u*v + w*(x - y);
            r(5) = 2*u*w + v*(x - y);
        }  
        break;
        default:
            throw std::runtime_error("Invalid Cid");
    }
    return r;
}



dCFdX_type zheng_dCFdX(const matrix_type& H, 
                       const vec6& dHdx, const vec6& dHdy) {

    const scalar_type& x = H(0,0);
    const scalar_type& y = H(1,1);
    const scalar_type& z = H(2,2);
    const scalar_type& u = H(0,1);
    const scalar_type& v = H(0,2);
    const scalar_type& w = H(1,2);

    CF_type CF(0);
    dCFdT_type dCFdT(0);

    // first compute derivative of constraint function w.r.t 
    // tensor coefficients using notations and order above 
    // dCF/dT: 7 x 6 matrix

    // dC0/d(x,y,z,u,v,w)
    dCFdT(0,0) = u*u - v*v + 2*x*(z - y) + y*y - z*z;
    dCFdT(0,1) = w*w - u*u + 2*y*(x - z) - x*x + z*z;
    dCFdT(0,2) = v*v - w*w + 2*z*(y - x) - y*y + x*x;
    dCFdT(0,3) = 2 * u * (x - y);
    dCFdT(0,4) = 2 * v * (z - x);
    dCFdT(0,5) = 2 * w * (y - z);

    // dC1/d(x,y,z,u,v,w)
    dCFdT(1,0) = 2*u*v + 2*w*(-2*x + y + z);
    dCFdT(1,1) = w*(2*x - 2*z) - u*v;
    dCFdT(1,2) = w*(2*x - 2*y) - u*v;
    dCFdT(1,3) = -2*u*w + v*(2*x - y - z);
    dCFdT(1,4) = -2*v*w + u*(2*x - y - z);
    dCFdT(1,5) = -u*u - v*v + 6*w*w - 2*x*x + 2*x*y + 2*x*z - 2*y*z;

    // dC2/d(x,y,z,u,v,w)
    dCFdT(2,0) = 2*v*(y-z) - u*w;
    dCFdT(2,1) = 2*u*w + 2*v*(x - 2*y + z);
    dCFdT(2,2) = -(u*w) + 2*v*(-x + y);
    dCFdT(2,3) = -2*u*v - w*(x - 2*y + z);
    dCFdT(2,4) = -u*u + 6*v*v - w*w + 2*(x - y)*(y - z);
    dCFdT(2,5) = -2*v*w - u*(x - 2*y + z);
    
    // dC3/d(x,y,z,u,v,w)
    dCFdT(3,0) = -(v*w) + 2*u*(-y + z);
    dCFdT(3,1) = -(v*w) + 2*u*(-x + z);
    dCFdT(3,2) = 2*v*w + 2*u*(x + y - 2*z);
    dCFdT(3,3) = 6*u*u - v*v - w*w - 2*x*y + 2*x*z + 2*y*z - 2*z*z;
    dCFdT(3,4) = -2*u*v - w*(x + y - 2*z);
    dCFdT(3,5) = -2*u*w - v*(x + y - 2*z);

    // dC4/d(x,y,z,u,v,w)
    dCFdT(4,0) = 0;
    dCFdT(4,1) = u*v;
    dCFdT(4,2) = -u*v;
    dCFdT(4,3) = -2*u*w + v*(y-z);
    dCFdT(4,4) = 2*v*w + u*(y - z);
    dCFdT(4,5)  = v*v - u*u;

    // dC5/d(x,y,z,u,v,w)
    dCFdT(5,0) = -u*w;
    dCFdT(5,1) = 0;
    dCFdT(5,2) = u*w;
    dCFdT(5,3) = 2*u*v + w*(-x + z);
    dCFdT(5,4) = u*u - w*w;
    dCFdT(5,5) = -2*v*w + u*(-x + z);

    // dC6/d(x,y,z,u,v,w)
    dCFdT(6,0) = v*w;
    dCFdT(6,1) = -v*w;
    dCFdT(6,2) = 0;
    dCFdT(6,3) = w*w - v*v;
    dCFdT(6,4) = -2*u*v + w*(x - y);
    dCFdT(6,5) = 2*u*w + v*(x - y);

    // std::cout << "dCFdT: " << dCFdT << '\n';

    dTdX_type dTdX;
    dTdX(0,0) = dHdx(0); // dT00/dx
    dTdX(0,1) = dHdy(0); // dT00/dy
    dTdX(1,0) = dHdx(1); // dT11/dx
    dTdX(1,1) = dHdy(1); // dT11/dy
    dTdX(2,0) = dHdx(2); // dT22/dx
    dTdX(2,1) = dHdy(2); // dT22/dy
    dTdX(3,0) = dHdx(3); // dT01/dx
    dTdX(3,1) = dHdy(3); // dT01/dy
    dTdX(4,0) = dHdx(4); // dT02/dx
    dTdX(4,1) = dHdy(4); // dT02/dy
    dTdX(5,0) = dHdx(5); // dT12/dx
    dTdX(5,1) = dHdy(5); // dT12/dy

    // dCFdT: 7 x 6
    // dTdX: 6 x 2
    return dCFdT*dTdX;
}

face_step_type zheng_tangent(const dCFdX_type& dCFdX, const CF_type& CF) {
    // dCFdX: 7 x 2
    // CF: 7 x 1
    // dCFdX^T: 2 x 7
    // dCFdX^T * CF: 2 x 1
    return spurt::moore_penrose_pseudoinverse(dCFdX) * CF;
}

struct face_pos {
    scalar_type m_cst;
    int m_dim0, m_dim1;

    static pos_type lift(scalar_type x, scalar_type y, 
                         int dim0, int dim1, scalar_type cst) {
        pos_type q(cst);
        q[dim0] = x;
        q[dim1] = y;
        return q;
    }

    static void initialize(int faceid, int& dim0, int& dim1, scalar_type& cst) {
        switch (faceid) {
            case 0: dim0 = 0; dim1 = 1; cst = 0; break;
            case 1: dim0 = 0; dim1 = 1; cst = 1; break;
            case 2: dim0 = 0; dim1 = 2; cst = 0; break;
            case 3: dim0 = 1; dim1 = 2; cst = 1; break;
            case 4: dim0 = 0; dim1 = 2; cst = 1; break;
            case 5: dim0 = 1; dim1 = 2; cst = 0; break;
            default: throw std::runtime_error("Invalid face id");
        }
    }

    face_pos(int faceid) {
        initialize(faceid, m_dim0, m_dim1, m_cst);
    }

    pos_type operator()(scalar_type x, scalar_type y) const {
        return lift(x, y, m_dim0, m_dim1, m_cst);
    }

    pos_type operator()(const face_pos_type& p) const {
        return this->operator()(p[0], p[1]);
    }

    int dim0() const { return m_dim0; }

    int dim1() const { return m_dim1; }
};

static bool is_inside_face(const face_pos_type& p) {
    if (p[0] < 0 || p[0] > 1) return false;
    if (p[1] < 0 || p[1] > 1) return false;
    return true;
}

static scalar_type max_length(const face_pos_type& x, 
                              const face_step_type& dir,
                              scalar_type amax=1)
{
    face_step_type d = dir / spurt::norm(dir);
    scalar_type maxl = amax;
    for (int n=0; n<2; ++n) {
        scalar_type u = d[n];
        scalar_type l;
        if (u > 0) l = (1-x[n])/u;
        else if (u < 0) l = x[n]/(-u);
        else l=1;
        maxl = std::min(maxl, l);
    }
    return maxl;
}

/*
    numerical search for L point on unit square 
    function f: (u,v) -> f(u,v) = 7x1 vector
    function f': df/duv = 7 x 2 matrix

    Goal: find x1 such that f(x1) = 0
    Current state:
    f(x) = f0
    f'(x) = f0'
    f(x+dx) = f0 + f0'*dx = 0
    f0'dx = -f0 <=> dx = -f0'^+f0, 
        with f0'^+ : Moore-Penrose pseudo inverse
        f0'^+ = (f'0^T f'0)^{-1} f'0^T

    x_1 = x_0 - (f'^T f')^{-1} f^T CF

    (2x7 . 7x2)^{-1} . 1x7 

*/

struct CF_on_face {
    CF_on_face(const matrix_image_type& hessian, const coord_type& voxel_id,
               const face_pos& face)
        : m_hessian(hessian), m_face(face), m_voxel(voxel_id) {}

    CF_type operator()(const face_pos_type& p) const {
        pos_type x = m_face(p);
        matrix_type H = m_hessian.value_in_voxel(m_voxel, x);
        return zheng_CF(H);
    }

    face_step_type tangent(const face_pos_type& p) const {
        // const scalar_type eps = 1.0e-6;
        pos_type x = m_face(p);
        // std::cout << "m_face(" << p << ")=" << x << '\n';
        matrix_type H = m_hessian.value_in_voxel(m_voxel, x);
        // std::cout << "H(p)=H(" << p << ")=H(" << x << ")=" << H << '\n';
        CF_type CF = zheng_CF(H);
        // debugging
        /*
        pos_type x_ = m_face(p[0]+eps, p[1]);
        pos_type _x = m_face(p[0]-eps, p[1]);
        pos_type y_ = m_face(p[0], p[1]+eps);
        pos_type _y = m_face(p[0], p[1]-eps);
        matrix_type Hdx = m_hessian.value_in_voxel(m_voxel, x_);
        std::cout << "H(p+dx) = H(" << p+face_pos_type(eps,0) << ")=H(" << x_ << ")=" << Hdx << '\n';
        matrix_type Hdx_ = m_hessian.value_in_voxel(m_voxel, _x);
        std::cout << "H(p-dx) = H(" << p-face_pos_type(eps,0) << ")=H("<< _x << ")=" << Hdx_ << '\n';
        matrix_type Hdy = m_hessian.value_in_voxel(m_voxel, y_);
        std::cout << "H(p+dy) = H(" << p+face_pos_type(0, eps) << ")=H(" << y_ << ")=" << Hdy << '\n';
        matrix_type Hdy_ = m_hessian.value_in_voxel(m_voxel, _y);
        std::cout << "H(p-dy) = H(" << p-face_pos_type(0, eps) << ")=H(" << _y << ")=" << Hdy_ << '\n';
        matrix_type dHdx_alt = (Hdx - Hdx_)/(2*eps);
        matrix_type dHdy_alt = (Hdy - Hdy_)/(2*eps);
        vec6 dHdx2, dHdy2;
        dHdx2[0] = dHdx_alt(0,0);
        dHdx2[1] = dHdx_alt(1,1);
        dHdx2[2] = dHdx_alt(2,2);
        dHdx2[3] = dHdx_alt(0,1);
        dHdx2[4] = dHdx_alt(0,2);
        dHdx2[5] = dHdx_alt(1,2);
        dHdy2[0] = dHdy_alt(0,0);
        dHdy2[1] = dHdy_alt(1,1);
        dHdy2[2] = dHdy_alt(2,2);
        dHdy2[3] = dHdy_alt(0,1);
        dHdy2[4] = dHdy_alt(0,2);
        dHdy2[5] = dHdy_alt(1,2);
        */
        
        tensor_type dH = m_hessian.derivative_in_voxel(m_voxel, x);
        vec6 dHdx, dHdy;
        dHdx[0] = dH(0,0,m_face.dim0());
        dHdx[1] = dH(1,1,m_face.dim0());
        dHdx[2] = dH(2,2,m_face.dim0());
        dHdx[3] = dH(0,1,m_face.dim0());
        dHdx[4] = dH(0,2,m_face.dim0());
        dHdx[5] = dH(1,2,m_face.dim0());
        dHdy[0] = dH(0,0,m_face.dim1());
        dHdy[1] = dH(1,1,m_face.dim1());
        dHdy[2] = dH(2,2,m_face.dim1());
        dHdy[3] = dH(0,1,m_face.dim1());
        dHdy[4] = dH(0,2,m_face.dim1());
        dHdy[5] = dH(1,2,m_face.dim1());

        /*
        std::cout << "dHdx=" << dHdx << '\n';
        std::cout << "dHdx_alt=" << dHdx2 << '\n';
        std::cout << "dHdy=" << dHdy << '\n';
        std::cout << "dHdy_alt=" << dHdy2 << '\n';

        vec6 delta_x = dHdx2 - dHdx;
        vec6 delta_y = dHdy2 - dHdy; 
        std::cout << "error in dHdx: " << spurt::norm(delta_x) << '\n';
        std::cout << "error in dHdy: " << spurt::norm(delta_y) << '\n';
        */

        dCFdX_type dCFdX = zheng_dCFdX(H, dHdx, dHdy);
        /*
        std::cout << "CF(p) = " << CF << '\n';

        sym_type H0;
        H0(0) = H(0,0);
        H0(1) = H(0,1);
        H0(2) = H(0,2);
        H0(3) = H(1,1);
        H0(4) = H(1,2);
        H0(5) = H(2,2);
        for (int i=0; i<6; ++i) {
            sym_type h_ = H0;
            sym_type _h = H0;
            h_(i) += eps;
            _h(i) -= eps;
            matrix_type hh_, _hh;
            hh_(0,0) = h_[0];
            hh_(0,1) = hh_(1,0) = h_[1];
            hh_(0,2) = hh_(2,0) = h_[2];
            hh_(1,1) = h_[3];   
            hh_(1,2) = hh_(2,1) = h_[4];
            hh_(2,2) = h_[5];
            _hh(0,0) = _h[0];
            _hh(0,1) = _hh(1,0) = _h[1];
            _hh(0,2) = _hh(2,0) = _h[2];
            _hh(1,1) = _h[3];
            _hh(1,2) = _hh(2,1) = _h[4];
            _hh(2,2) = _h[5];
            CF_type CF_ = zheng_CF(hh_);
            CF_type _CF = zheng_CF(_hh);
            CF_type dCFdX_ = (CF_ - _CF) / (2*eps);
            std::cout << "dCFdT[" << i << "] = " << dCFdX_ << '\n';           
        }
        std::cout << '\n';
        for (int c=0; c<7; ++c) {
            std::cout << "dC[" << c << "]/dT = " << dC1dT(H, c) << '\n';
        }
    

        CF_type CFdx = zheng_CF(Hdx);
        std::cout << "CF(p+dx) = " << CFdx << '\n';
        CF_type CFdx_ = zheng_CF(Hdx_);
        std::cout << "CF(p-dx) = " << CFdx_ << '\n';
        CF_type CFdy = zheng_CF(Hdy);
        std::cout << "CF(p+dy) = " << CFdy << '\n';
        CF_type CFdy_ = zheng_CF(Hdy_);
        std::cout << "CF(p-dy) = " << CFdy_ << '\n';
        CF_type dCFdx = (CFdx - CFdx_) / (2*eps);   
        CF_type dCFdy = (CFdy - CFdy_) / (2*eps);

        std::cout << "dCFdX = " << dCFdX << '\n';
        std::cout << "dCFdx alt = " << dCFdx << '\n';
        std::cout << "dCFdy alt = " << dCFdy << '\n';
        */
        return zheng_tangent(dCFdX, CF);
    }

    const matrix_image_type& m_hessian;
    face_pos m_face;
    coord_type m_voxel;
};


bool lnsearch(const CF_on_face& CF, // value function
              face_pos_type& x, // current state
              CF_type& f, // value at current state
              const face_step_type& dd, // tangent at current state
              scalar_type ml=0.1) 
{
    scalar_type lambda = 1.0;
    const scalar_type alpha = 1e-4;
    const scalar_type maxl = max_length(x, dd, ml);

    // std::cout << "maxl = " << maxl << '\n';
    // std::cout << "tangent = " << dd << '\n';
    
    face_pos_type x0 = x;
    CF_type f0 = f;
    face_step_type d = (norm(dd) > maxl) ? dd * maxl / norm(dd) : dd;
    scalar_type v0 = spurt::norm(f0);
    // std::cout << "v0 = " << v0 << '\n';
    // std::cout << "d=" << d << '\n';
    for (unsigned int i = 0; i < 20; ++i) {
        x = x0 + lambda * d;
        f = CF(x);
        // std::cout << "f=" << f << ", x=" << x << ", lambda=" << lambda << ", norm=" << spurt::norm(f) << '\n';
        if (spurt::norm(f) < (1 - alpha*lambda)*v0) {
            return true;
        }
        
        lambda *= 0.5;
    }
    
    return false;
}

bool findLPoint(face_pos_type& Lpt, 
                const coord_type& voxel_id, int faceid, 
                const matrix_image_type& hessian,
                int maxiter=100)
{   
    // we will work in local voxel coordinates
    face_pos to3d(faceid); 
    CF_on_face CF(hessian, voxel_id, to3d);

    // create a dense sampling of norm(f) for debugging

    VTK_CREATE(vtkImageData, image);
    image->SetDimensions(100, 100, 1);
    VTK_CREATE(vtkFloatArray, data);
    data->SetNumberOfComponents(1);
    VTK_CREATE(vtkFloatArray, verts);
    verts->SetNumberOfComponents(3);
    for (int i=0; i<100; ++i) {
        for (int j=0; j<100; ++j) {
            face_pos_type X(i/100., j/100.);
            CF_type f = CF(X);
            data->InsertNextTuple1(spurt::norm(f));
        }
    }
    std::ostringstream os;
    os << "normf_" << voxel_id << "_" << faceid << ".vti";
    image->GetPointData()->SetScalars(data);
    VTK_CREATE(vtkXMLImageDataWriter, writer);
    writer->SetInputData(image);
    writer->SetFileName(os.str().c_str());
    writer->Write();

    face_pos_type X(0.5, 0.5);
    verts->InsertNextTuple3(X[0]*100., X[1]*100., 0);    
    CF_type f = CF(X);
    VTK_CREATE(vtkCellArray, cells);
    // std::cout << "initial norm is " << spurt::norm(f) << '\n';
    for (int i=0 ; i<maxiter; ++i) {
        // std::cout << "iteration " << i << '\n';
        face_step_type d = -CF.tangent(X);
        if (!lnsearch(CF, X, f, d)) {
            break;
        }
        verts->InsertNextTuple3(X[0]*100., X[1]*100., 0);
        cells->InsertNextCell(2);
        cells->InsertCellPoint(i);
        cells->InsertCellPoint(i+1);
        std::cout << "current norm is " << spurt::norm(f) << '\n';
        if (spurt::norm(f) < 1.0e-6) {
            std::cout << "Solution found\n";
            Lpt = X;
            return true;
        }
    }

    VTK_CREATE(vtkPoints, points);
    points->SetData(verts);
    VTK_CREATE(vtkPolyData, poly);
    poly->SetPoints(points);
    poly->SetLines(cells);
    VTK_CREATE(vtkXMLPolyDataWriter, writer2);
    writer2->SetInputData(poly);
    os.clear();
    os << "path_" << voxel_id << "_" << faceid << ".vtp";
    writer2->SetFileName(os.str().c_str());
    writer2->Write();
    return false;
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
    else {
        std::cout << "delta underflow\n";
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

matrix_type approx_nabla_h(const pos_type& p, 
                           const vector_image_type& gradient,
                           const matrix_image_type& hessian,
                           scalar_type theta)
{
    const scalar_type eps = 0.0001;
    matrix_type nabla_h(0);
    for (int dim=0; dim<3; ++dim) {
        pos_type _p = p;
        pos_type p_ = p;
        _p[dim] -= eps;
        if (!gradient.grid().bounds().inside(_p)) _p[dim] = p[dim];
        p_[dim] += eps;
        if (!gradient.grid().bounds().inside(p_)) p_[dim] = p[dim];
        vector_type _g, g_;
        _g = gradient.value(_p);
        g_ = gradient.value(p_);
        matrix_type _hess, hess_;
        _hess = hessian.value(_p);
        hess_ = hessian.value(p_);
        vector_type _h, h_;
        _h = vector_h(_g, _hess, theta);
        h_ = vector_h(g_, hess_, theta);
        nabla_h.column(dim) = 1/(p_[dim]-_p[dim])*(h_-_h);
    }
    return nabla_h;
}

vector_type ridge_normal(const pos_type& point,
                         const vector_type& g,  const matrix_type& H, 
                         const vector_image_type& gradient,
                         const matrix_image_type& hessian,
                         scalar_type theta)
{
    std::cout << "value of h at ridge point: " << vector_h(g, H, theta) << '\n';
    auto H_prime = hessian.derivative(point); // 3rd order tensor
    matrix_type h_prime;
    h_prime.column(0) = vector_dh(g, H, H.column(0), H_prime.layer(0), theta);
    h_prime.column(1) = vector_dh(g, H, H.column(1), H_prime.layer(1), theta);
    h_prime.column(2) = vector_dh(g, H, H.column(2), H_prime.layer(2), theta);

    matrix_type h_prime_alt = approx_nabla_h(point, gradient, hessian, theta);
    std::cout << "Derived h_prime:\n" << h_prime << '\n';
    std::cout << "Approx h_prime:\n" << h_prime_alt << '\n';
    std::cout << "Difference: " << (h_prime-h_prime_alt).norm() << '\n';

    {
        vector_type hevals;
        matrix_type hevecs;
        sym_eigensystem(hevals, hevecs, H);
        std::cout << "eigenvalues of H: " << hevals << '\n';
        std::cout << "eigenvectors of H: " << hevecs << '\n';
        sym_eigensystem(hevals, hevecs, h_prime_alt);
        std::cout << "eigenvalues of approx h_prime: " << hevals << '\n';
        std::cout << "eigenvectors of approx h_prime: " << hevecs << '\n';
    }

    // h' = outer(n,n);
    vector_type hevals;
    matrix_type hevecs;
    sym_eigensystem(hevals, hevecs, h_prime);
    std::cout << "eigenvalues of h_prime: " << hevals << '\n';
    std::cout << "eigenvectors of h_prime: " << hevecs << '\n';
    vector_type hevals_plus = spurt::abs(hevals);
    std::cout << "Return eigenvector #" << spurt::argmax(hevals_plus) << '\n';
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

std::pair<pos_type, scalar_type>
project(const vector_type& g, const matrix_type& h) {
    vector_type evals;
    matrix_type evecs;
    sym_eigensystem(evals, evecs, h);
    pos_type coords = spurt::abs(transpose(evecs)*g);
    return std::make_pair(coords, evals[2]);
}

double determinant(const vector_type& g, const matrix_type& H, bool normalize=false) {
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


template<typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T,N>& a) {
    os << "[";
    for (int i=0; i<N-1; ++i) {
        os << a[i] << ", ";
    }
    os << a[N-1] << "]";
    return os;
}

void save_sampled_voxel(const std::string& filename, coord_type& voxel_id, 
                        const scalar_image_type& values,
                        const vector_image_type& gradient, 
                        const matrix_image_type& hessian,
                        const scalar_type theta,
                        const size_t resolution=101)
{
    std::vector<double> det_values(5 * resolution * resolution * resolution);
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
                vec3 g = gradient.value_in_voxel(voxel_id, pos_type(x,y,z));
                mat3 H = hessian.value_in_voxel(voxel_id, pos_type(x,y,z));
                det_values[3 * (i + resolution * (j + resolution * k))] = determinant(g, H);
                det_values[3 * (i + resolution * (j + resolution * k)) + 1] = (project(g, H).first)[2];
                scalar_type v = values.value_in_voxel(voxel_id, pos_type(x,y,z));
                det_values[3 * (i + resolution * (j + resolution * k)) + 2] = v;
            }
        }
    }
    size_t sizes[4] = {3, resolution, resolution, resolution};
    std::ostringstream os;
    os << filename << "_voxel_" << voxel_id[0] << "_" << voxel_id[1] << "_" << voxel_id[2] << ".nrrd";
    spurt::nrrd_utils::writeNrrd((void *)&det_values[0], os.str(), nrrdTypeDouble, 4, sizes);
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
                std::cout << "current h value at u=" << u << " is " << h2 << ", norm =" << n2 << '\n';
                if (n2 < 1.0e-6) {
                    std::cout << "solution found.\n";
                    failed = false;
                    break;
                }
                std::array<scalar_type, 3> ns({n0, n1, n2});
                auto ranks = arg_sort(ns.begin(), ns.end());
                if (ranks[2] == 2) {
                    std::cout << "linear estimate increased h norm\n";
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
                std::cout << "min norm found (" << hnorms[ranks[0]] << ") found at u = " << ranks[0]/scalar_type(nsteps) << '\n';
                std::cout << "hnorms were ";
                for (auto anorm : hnorms) std::cout << anorm << ", ";
                std::cout << '\n';
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

void visualize(const std::vector<pos_type>& points, 
               const std::vector<scalar_type>& values,
               const std::vector<scalar_type>& ridge_strength,
               vtkCellArray* triangles, vtkCellArray* edges, 
               vtkCellArray* verts, const std::string& fname,
               bool novis,
               size_type minsize=-1)
{
    vtkSmartPointer<vtkPolyData> pd = vtk_utils::make_points(points);
    pd->SetPolys(triangles);
    pd->SetLines(edges);
    pd->SetVerts(verts);
    pd = vtk_utils::add_scalars(pd, values, true, "values", false);
    pd = vtk_utils::add_scalars(pd, ridge_strength, true, "ridge_strength", true);

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
    std::vector< coord_type > voxels;
    if (spurt::any(voxel_id != invalid_coord))
    {
        voxels.clear();
        voxels.push_back(voxel_id);
        verbose = 2;
        save_sampled_voxel(output_name, voxel_id, values, gradient, hessian, theta0, 101);
    }
    else
    {
        std::cout << "selected voxels..." << std::flush;
        select_voxels(voxels, bounds, shape);
        std::cout << " done\n";
    }

    if (verbose) std::cout << "There are " << voxels.size() << " voxels in input for a total of " << 12*voxels.size() << " (redundant) edges\n";

    typedef size_type glob_point_id_t;
    typedef int loc_face_id_t; // [0, 6)
    typedef int face_edge_id_t; // [0, 4) 
    typedef int voxel_edge_id_t; // [0, 12)
    typedef std::pair<voxel_edge_id_t, voxel_edge_id_t> loc_edge_pair_t;
    std::map< EdgeID, glob_point_id_t > edge_to_point_index;
    container_wrapper<pos_type> unique_points;
    container_wrapper<scalar_type> unique_values;
    container_wrapper<scalar_type> unique_strengths;

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
    // keep track of what triangles stem from what voxels to remove them 
    // if needed in case their containing voxels are subdivided. 
    std::map<coord_type, std::vector<int>, spurt::lexicographical_order> voxel_to_triangles;
    std::vector<triangle_type> all_triangles;
    std::vector<vec4> rejected;

    spurt::ProgressDisplay progress;

    VTK_CREATE(vtkCellArray, vtk_cells);
    VTK_CREATE(vtkCellArray, vtk_lines);
    VTK_CREATE(vtkCellArray, vtk_verts);

    std::vector<broken_voxel> broken_voxels;
    std::vector<EdgeID> double_edges;
    std::vector<coord_type> to_subdivide;
    std::vector<int> voxel_to_edge;
    std::map<coord_type, std::vector<int>, spurt::lexicographical_order> voxel_counter;

    srand48(130819751900);
    progress.begin(voxels.size(), "Extract ridges", 10000, "tris: 0, done: 0, ok: 0, skip: 0, underflow: 0, weak: 0, low: 0, none: 0, even: 0, failed: 0");
    for (int n=0; n<voxels.size(); ++n) 
    {
        coord_type voxel_id = voxels[n];
        if (verbose >= 2)
        {
            std::cout << "processing voxel: " << voxel_id << '\n';
        }

        std::map<int, size_type> found;
        std::ostringstream update_oss;
        update_oss.imbue(std::locale());
        update_oss << "tris: " << vtk_cells->GetNumberOfCells()  
            << ", ok: "  << nsucceeded 
            << ", skip: "  << nskipped 
            << ", none: " << nnone
            << ", failed: "  << nfailed_to_triangulate
            << ", insuf: " << nnot_enough_points;
        progress.update(n, update_oss.str());

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
            nnone++;
            continue;
        }
        std::ostringstream log;
        for (voxel_edge_id_t edge_id = 0; edge_id < 12; ++edge_id) 
        {
            coord_type v0 = voxel_id + coord_type(ivec3(spurt::marching_cubes::canonical_edge_coordinates[edge_id][0]));
            coord_type v1 = voxel_id + coord_type(ivec3(spurt::marching_cubes::canonical_edge_coordinates[edge_id][1]));
            EdgeID global_edge_id(v0, v1);
            auto iter = edge_to_point_index.find(global_edge_id);
            if (iter != edge_to_point_index.end())
            {   
                if (verbose) {
                    std::cout << "entry found for edge #" << edge_id << " between " << v0 << " and " << v1 << '\n';
                    std::cout << "edge #" << edge_id << " of voxel " << voxel_id << ": " << global_edge_id << " has already been processed\n";
                    std::cout << "map contains " << edge_to_point_index.size() << " elements\n";
                }
                nskipped += 1;
                size_type id = iter->second;
                if (id >= 0)
                {
                    found[edge_id] = id;
                    if (verbose) {
                        std::cout << "This entry contains a ridge point\n";
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
                if (verbose >= 2) 
                    std::cout << "\n\nedge #" << edge_id << " of voxel " << voxel_id << ": " << global_edge_id << " is being processed" << '\n';
                auto result = process_edge(global_edge_id, values, gradient, hessian, theta0, verbose, true);
                if (result.first == 1)
                {        
                    if (verbose) 
                        std::cout << "process_edge returned 1 solution\n";
                    scalar_type u = result.second;
                    pos_type p = (1. - u) * global_edge_id[0] + u * global_edge_id[1];
                    found[edge_id] = unique_points.add(p);
                    edge_to_point_index[global_edge_id] = found[edge_id];
                    auto v = evaluate(p, values, hessian);
                    unique_values.add(v.first);
                    unique_strengths.add(v.second);
                    if (verbose >= 2) 
                        std::cout << "Found now contains " << found.size() << " entries\n";
                }
                else
                {
                    if (verbose >= 2)
                        std::cout << "no ridge point\n";
                    nnone++;
                    edge_to_point_index[global_edge_id] = -1;
                }
                nprocessed++;
            }
        }

        std::map<loc_face_id_t, glob_point_id_t> loc_face_to_Lpoints;
        std::map<loc_face_id_t, loc_edge_pair_t > loc_face_to_segments;
        std::map<loc_face_id_t, voxel_edge_id_t> face_to_edge_segments;
        bool pathological = false;
        if (true) 
        {
            // check if we have a pathological case on our hands...
            for (loc_face_id_t f=0; f<6; ++f) {
                auto face = spurt::marching_cubes::canonical_face_edge_indices[f];
                face_pos to3d(f);
                std::cout << "processing face " << f << '\n';
                std::vector<voxel_edge_id_t> edge_ids;
                for (face_edge_id_t e=0; e<4; ++e) {
                    if (found.find(face[e]) != found.end()) {
                        std::cout << "point found in edge " << face[e] << '\n';
                        edge_ids.push_back(face[e]); 
                    }
                }
                if (edge_ids.size() == 2) {
                    loc_face_to_segments[f] = std::make_pair(edge_ids[0], edge_ids[1]);
                }
                else if (edge_ids.size() % 2) {
                    pathological = true;
                    std::cout << "Face " << f << " contains an odd number of crease points\n";
                    // do something clever about it...
                    bool Lfound = false;
                    if (edge_ids.size() == 3) {
                        face_pos_type Lpt;
                        if (findLPoint(Lpt, voxel_id, f, hessian)) {
                            std::cout << "L point found at " << Lpt << '\n';
                            auto p = to3d(Lpt);
                            Lfound = true;
                            loc_face_to_Lpoints[f] = unique_points.add(p);
                        }
                        else {
                            std::cout << "L point not found\n";
                        }

                        // compute ridge normals at ridge points
                        pos_type p[3];
                        vector_type g[3];
                        matrix_type h[3];
                        vector_type n[3];
                        for (int k=0; k<3; ++k) {
                            p[k] = unique_points[found[edge_ids[k]]];
                            g[k] = gradient.value(p[k]);
                            h[k] = hessian.value(p[k]);
                            n[k] = ridge_normal(p[k], g[k], h[k], gradient, hessian, theta0);
                        }
                        std::array<scalar_type, 3> scores;
                        for (int l=0; l<3; ++l) {
                            vector_type lk = p[(l+1)%3]-p[l];
                            lk /= spurt::norm(lk);
                            scores[l] = std::abs(spurt::inner(lk, n[l])) + std::abs(spurt::inner(lk, n[(l+1)%3]));
                        }
                        int maxid = std::distance(scores.begin(), std::min_element(scores.begin(), scores.end()));
                        std::cout << "The best connection found is edge " << edge_ids[maxid] << " - " << edge_ids[(maxid+1)%3] << '\n';
                        loc_face_to_segments[f] = 
                            std::make_pair(edge_ids[maxid], 
                                           edge_ids[(maxid+1)%3]);
                        // if (Lfound) {
                        //     boundary_segments[f] = edge_ids[(maxid+2)%3]
                        // }
                    }
                }
                else if (edge_ids.size() == 4) {
                    std::cout << "Face " << f << " contains 4 crease points\n";
                }
                else {
                    std::cout << "Face " << f << " contains " << edge_ids.size() << " crease points\n";
                }
            }
        }


        /*
        if (found.size() >= 3) 
        {
            if (verbose >= 2) 
            {
                std::cout << "After looking at all edges found contains " << found.size() << " entries\n";
                for (auto iter = found.begin(); iter!=found.end(); ++iter)
                {
                    std::cout << "Edge #" << iter->first << " contains " << iter->second << '\n';
                    auto r = evaluate(iter->second, values, hessian);
                    std::cout << "ridge point at " << iter->second << " has value " << r.first << " and ridge strength " << r.second << '\n';
                }
            }
            std::vector<triangle_type> tris;
            int tri_case;
            if (verbose)
                tri_case = triangulate(tris, found, std::cout);
            else
                tri_case = triangulate(tris, found, log);

            if (tri_case == one_triangle || tri_case == valid_mc_case) 
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
                        const pos_type& p = iter->second[l];
                        vec4 q({p[0], p[1], p[2], tri_case});
                        rejected.push_back(q);
                    }
                }
                broken_voxel broken;
                broken.id = voxel_id;
                broken.edge_points = found;
                for (int vid=0; vid<8; ++vid) {
                    coord_type pid = voxel_id + coord_type(spurt::marching_cubes::vertices[vid]);
                    broken.values[vid] = values(pid[0], pid[1], pid[2]);
                    broken.gradients[vid] = gradient(pid[0], pid[1], pid[2]);
                    broken.hessians[vid] = hessian(pid[0], pid[1], pid[2]);
                    broken.strengths[vid] = evmin(broken.hessians[vid]).first;
                }
                broken_voxels.push_back(broken);
            }
        }
        else if (found.size() > 0) 
        {
            nnot_enough_points++;
            for (auto iter = found.begin(); iter != found.end(); ++iter)
            {
                const pos_type& p = iter->second;
                vec4 q({p[0], p[1], p[2], not_enough_edges});
                rejected.push_back(q);
            }

            broken_voxel broken;
            broken.id = voxel_id;
            broken.edge_points = found;
            for (int vid = 0; vid < 8; ++vid)
            {
                coord_type pid = voxel_id + coord_type(ivec3(spurt::marching_cubes::vertices[vid]));
                broken.values[vid] = values(pid[0], pid[1], pid[2]);
                broken.gradients[vid] = gradient(pid[0], pid[1], pid[2]);
                broken.hessians[vid] = hessian(pid[0], pid[1], pid[2]);
                broken.strengths[vid] = evmin(broken.hessians[vid]).first;
            }
            broken_voxels.push_back(broken);
        }*/

        int mesh_case;
        if (verbose >= 2) 
            mesh_case = crude_meshing(found, vtk_cells, vtk_lines, vtk_verts, std::cout);
        else {
            mesh_case = crude_meshing(found, vtk_cells, vtk_lines, vtk_verts, log);
        }
        if (mesh_case == one_triangle || mesh_case == valid_mc_case) {
            ++nsucceeded;
        }
        else {
            if (mesh_case == invalid_mc_case) 
                ++nfailed_to_triangulate;
            else
                ++nnot_enough_points;

            broken_voxel broken;
            broken.id = voxel_id;
            std::map<int, int> edge_to_rank;
            for (auto iter = found.begin(); iter!=found.end(); ++iter) {
                value_set vs;
                edge_to_rank[iter->first] =  broken.edg_points.size();
                const pos_type& p = unique_points[iter->second];
                vs.value = values.value(p);
                vs.gradient = gradient.value(p);
                vs.hessian = hessian.value(p);
                vs.strength = strength.value(p);
                vs.normal = ridge_normal(p, vs.gradient, vs.hessian, gradient, hessian, theta0);
                std::cout << "normal is " << vs.normal << '\n';
                broken.edg_points.push_back(p);
                broken.edg_values.push_back(vs);
            }
            for (auto p : loc_face_to_Lpoints) {
                broken.L_points.push_back(unique_points[p.second]);
            }
            for (auto p : loc_face_to_segments) {
                broken.face_edges.push_back(std::make_pair(edge_to_rank[p.second.first], edge_to_rank[p.second.second]));
            }
            broken.min = values.grid()(voxel_id);
            broken.max = values.grid()(voxel_id + coord_type(1,1,1));
            for (int vid = 0; vid < 8; ++vid)
            {
                coord_type pid = voxel_id + coord_type(ivec3(spurt::marching_cubes::vertices[vid]));

                value_set& data = broken.vox_values[vid];
                
                data.value = values(pid[0], pid[1], pid[2]);
                data.gradient = gradient(pid[0], pid[1], pid[2]);
                data.hessian = hessian(pid[0], pid[1], pid[2]);
                data.strength = strength(pid[0], pid[1], pid[2]);
            }
            broken_voxels.push_back(broken);
        }
    }

    std::cout << "\n\nstats: n3 = " << nthree << ", n4 = " << nfour << ", n>4 = " << nmore_than_four << ", n2 = " << ntwo << ", none = " << n_one << '\n';

    std::string output_base = filename::remove_extension(output_name);
    export_broken_voxels(broken_voxels, output_base + "_broken.vtp");

    /*
    std::vector<pos_type> all_edge_points;
    for (auto iter=all_processed_edges.begin(); iter!=all_processed_edges.end(); ++iter)
    {
        auto pt = iter->second;
        if (pt != invalid_pos)
            continue;
        all_edge_points.push_back(pt);
    }
    {
        size_t sizes[3] = { 3, all_edge_points.size() };
        spurt::nrrd_utils::writeNrrd((void*)&all_edge_points[0], output_name + "_all_points.nrrd", 
                                     nrrdTypeDouble, 2, sizes);
    }

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
        visualize(all_processed_edges, all_triangles);
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
    */

    visualize(unique_points, unique_values, unique_strengths, vtk_cells, vtk_lines, vtk_verts, output_name, novis, minsize);

    return 0;
}
