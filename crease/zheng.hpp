#pragma once

#include <exception>
#include <iostream>
#include <math/small_vector.hpp>
#include <math/small_matrix.hpp>
#include <math/small_tensor.hpp>

namespace spurt {

namespace zheng {

template<typename T>
using vector2 = spurt::small_vector<T, 2>;

template<typename T>
using vector3 = spurt::small_vector<T, 3>;

template<typename T>
using matrix3 = spurt::small_matrix<T, 3, 3>;

template<typename T>
using vector6 = spurt::small_vector<T, 6>;

template<typename T>
using vector7 = spurt::small_vector<T, 7>;

template<typename T>
using matrix62 = spurt::small_matrix<T, 6, 2>;

template<typename T>
using matrix72 = spurt::small_matrix<T, 7, 2>;

template<typename T>
using matrix76 = spurt::small_matrix<T, 7, 6>;

template<typename T>
using tensor3 = spurt::small_tensor<T, 3, 3, 3>;

// 7 constraint functions
template<typename T>
vector7<T> compute_CF7(const matrix3<T>& H) {
    const T& x = H(0,0);
    const T& y = H(1,1);
    const T& z = H(2,2);
    const T& u = H(0,1);
    const T& v = H(0,2);
    const T& w = H(1,2);

    vector7<T> _CF(0);

    _CF[0] = 
        x * ((y*y - z*z) + (u*u - v*v)) + 
        y * ((z*z - x*x) + (w*w - u*u)) + 
        z * ((x*x - y*y) + (v*v - w*w));

    _CF[1] = 
        w *(2*(w*w-x*x) - (v*v + u*u) + 2*(y*x + z*x - y*z)) + 
        u*v*(2*x - z - y);

    _CF[2] = 
        v * (2*(v*v - y*y) - (u*u + w*w) + 2*(z*y + x*y - z*x)) + 
        w*u*(2*y - x - z);

    _CF[3] = 
        u*(2*(u*u - z*z) - (w*w + v*v) + 2*(x*z + y*z - x*y)) + 
        v*w*(2*z - y - x);

    _CF[4] = w*(v*v-u*u) + u*v*(y-z);
    _CF[5] = v*(u*u-w*w) + w*u*(z-x);
    _CF[6] = u*(w*w-v*v) + v*w*(x-y);

    return _CF;
}

// derivative of one of the 7 constraint functions w.r.t. tensor coefficients
template<typename T>
vector6<T> compute_dCF1dT(const matrix3<T>& H, int Cid) {
    const T& x = H(0,0);
    const T& y = H(1,1);
    const T& z = H(2,2);
    const T& u = H(0,1);
    const T& v = H(0,2);
    const T& w = H(1,2);

    vector6<T> r(0);
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

// derivative of 7 constraint functions w.r.t. 2 spatial dimensions
template<typename T>
matrix72<T> compute_dCF7dX(const matrix3<T>& H, 
                           const vector6<T>& dHdx, const vector6<T>& dHdy) {

    const T& x = H(0,0);
    const T& y = H(1,1);
    const T& z = H(2,2);
    const T& u = H(0,1);
    const T& v = H(0,2);
    const T& w = H(1,2);

    matrix76<T> dCF7dT(0);

    // first compute derivative of constraint function w.r.t 
    // tensor coefficients using notations and order above 
    // dCF/dT: 7 x 6 matrix

    // dC0/d(x,y,z,u,v,w)
    dCF7dT(0,0) = u*u - v*v + 2*x*(z - y) + y*y - z*z;
    dCF7dT(0,1) = w*w - u*u + 2*y*(x - z) - x*x + z*z;
    dCF7dT(0,2) = v*v - w*w + 2*z*(y - x) - y*y + x*x;
    dCF7dT(0,3) = 2 * u * (x - y);
    dCF7dT(0,4) = 2 * v * (z - x);
    dCF7dT(0,5) = 2 * w * (y - z);

    // dC1/d(x,y,z,u,v,w)
    dCF7dT(1,0) = 2*u*v + 2*w*(-2*x + y + z);
    dCF7dT(1,1) = w*(2*x - 2*z) - u*v;
    dCF7dT(1,2) = w*(2*x - 2*y) - u*v;
    dCF7dT(1,3) = -2*u*w + v*(2*x - y - z);
    dCF7dT(1,4) = -2*v*w + u*(2*x - y - z);
    dCF7dT(1,5) = -u*u - v*v + 6*w*w - 2*x*x + 2*x*y + 2*x*z - 2*y*z;

    // dC2/d(x,y,z,u,v,w)
    dCF7dT(2,0) = 2*v*(y-z) - u*w;
    dCF7dT(2,1) = 2*u*w + 2*v*(x - 2*y + z);
    dCF7dT(2,2) = -(u*w) + 2*v*(-x + y);
    dCF7dT(2,3) = -2*u*v - w*(x - 2*y + z);
    dCF7dT(2,4) = -u*u + 6*v*v - w*w + 2*(x - y)*(y - z);
    dCF7dT(2,5) = -2*v*w - u*(x - 2*y + z);
    
    // dC3/d(x,y,z,u,v,w)
    dCF7dT(3,0) = -(v*w) + 2*u*(-y + z);
    dCF7dT(3,1) = -(v*w) + 2*u*(-x + z);
    dCF7dT(3,2) = 2*v*w + 2*u*(x + y - 2*z);
    dCF7dT(3,3) = 6*u*u - v*v - w*w - 2*x*y + 2*x*z + 2*y*z - 2*z*z;
    dCF7dT(3,4) = -2*u*v - w*(x + y - 2*z);
    dCF7dT(3,5) = -2*u*w - v*(x + y - 2*z);

    // dC4/d(x,y,z,u,v,w)
    dCF7dT(4,0) = 0;
    dCF7dT(4,1) = u*v;
    dCF7dT(4,2) = -u*v;
    dCF7dT(4,3) = -2*u*w + v*(y-z);
    dCF7dT(4,4) = 2*v*w + u*(y - z);
    dCF7dT(4,5)  = v*v - u*u;

    // dC5/d(x,y,z,u,v,w)
    dCF7dT(5,0) = -u*w;
    dCF7dT(5,1) = 0;
    dCF7dT(5,2) = u*w;
    dCF7dT(5,3) = 2*u*v + w*(-x + z);
    dCF7dT(5,4) = u*u - w*w;
    dCF7dT(5,5) = -2*v*w + u*(-x + z);

    // dC6/d(x,y,z,u,v,w)
    dCF7dT(6,0) = v*w;
    dCF7dT(6,1) = -v*w;
    dCF7dT(6,2) = 0;
    dCF7dT(6,3) = w*w - v*v;
    dCF7dT(6,4) = -2*u*v + w*(x - y);
    dCF7dT(6,5) = 2*u*w + v*(x - y);

    matrix62<T> dTdX;
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

    // dCF7dT: 7 x 6
    // dTdX: 6 x 2
    return dCF7dT*dTdX;
}

template<typename T>
vector2<T> tangent(const matrix72<T>& dCF7dX, const vector7<T>& CF7) {
    // dCF7dX: 7 x 2
    // CF7: 7 x 1
    // dCF7dX^T: 2 x 7
    // dCF7dX^T * CF7: 2 x 1
    return spurt::moore_penrose_pseudoinverse(dCF7dX) * CF7;
}

// A data structure for a unit square 2D region embedded in a 3D mesh
template<typename T>
struct Face {
    typedef T scalar_t;
    typedef vector3<T> vec3_t;
    typedef vector2<T> vec2_t;

    std::array<int, 3> m_dims;
    vec3_t m_orig, m_e0, m_e1;

    Face(const vec3_t& origin, const vec3_t& e0, const vec3_t& e1) 
        : m_orig(origin), m_e0(e0), m_e1(e1) {
        for (int d=0; d<3; ++d) {
            if (m_e0[d] != 0) m_dims[0] = d;
            else if (m_e1[d] != 0) m_dims[1] = d;
            else m_dims[2] = d;
        }
    }

    vector3<scalar_t> operator()(scalar_t x, scalar_t y) const {
        return m_orig + x*m_e0 + y*m_e1;
    }

    vector3<scalar_t> operator()(const vec2_t& p) const {
        return this->operator()(p[0], p[1]);
    }

    const std::array<int, 3>& dims() const { return m_dims; }
};

template<typename T>
static bool is_inside_face(const vector2<T>& p) {
    if (p[0] < 0 || p[0] > 1) return false;
    if (p[1] < 0 || p[1] > 1) return false;
    return true;
}

template<typename T>
static T max_length(const vector2<T>& x, 
                    const vector2<T>& dir,
                    T amax=1)
{
    vector2<T> d = dir / spurt::norm(dir);
    T maxl = amax;
    for (int n=0; n<2; ++n) {
        T u = d[n];
        T l;
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
    f'(x) = G(x) = G0
    f(x+dx) = f0 + G0*dx = 0
    G0dx = -f0 <=> dx = -G0^+f0, 
        with G0^+ : Moore-Penrose pseudo inverse
        G0^+ = (G0^T G0)^{-1} G0^T

    x_1 = x_0 - (G^T G)^{-1} f^T CF

    (2x7 . 7x2)^{-1} . 1x7 

*/
template<typename T, typename MatrixFunction>
struct CF_on_face {
    typedef Face<T> face_t;
    typedef typename face_t::vec2_t vec2_t;
    typedef typename face_t::vec3_t vec3_t;
    typedef matrix3<T> mat3_t;
    typedef vector6<T> partial_t;
    typedef vector7<T> cfuncs_t;
    typedef tensor3<T> tensor_t;
    typedef matrix72<T> dcfuncs_t;

    CF_on_face(const MatrixFunction& hessian, const face_t& face)
        : m_hessian(hessian), m_face(face) {}

    cfuncs_t operator()(const vec2_t& p) const {
        return compute_CF7(m_hessian.value(m_face(p)));
    }

    vec2_t tangent(const vec2_t& p) const {
        mat3_t H = m_hessian.value(m_face(p));
        cfuncs_t cf_ = compute_CF7(H);
        
        tensor_t dH = m_hessian.derivative(m_face(p));
        partial_t dHdx, dHdy;
        dHdx[0] = dH(0,0,m_face.dims()[0]);
        dHdx[1] = dH(1,1,m_face.dims()[0]);
        dHdx[2] = dH(2,2,m_face.dims()[0]);
        dHdx[3] = dH(0,1,m_face.dims()[0]);
        dHdx[4] = dH(0,2,m_face.dims()[0]);
        dHdx[5] = dH(1,2,m_face.dims()[0]);
        dHdy[0] = dH(0,0,m_face.dims()[1]);
        dHdy[1] = dH(1,1,m_face.dims()[1]);
        dHdy[2] = dH(2,2,m_face.dims()[1]);
        dHdy[3] = dH(0,1,m_face.dims()[1]);
        dHdy[4] = dH(0,2,m_face.dims()[1]);
        dHdy[5] = dH(1,2,m_face.dims()[1]);

        dcfuncs_t dcfdx_ = compute_dCF7dX(H, dHdx, dHdy);
        return zheng::tangent(dcfdx_, cf_);
    }

    const MatrixFunction& m_hessian;
    face_t m_face;
};


template<typename T, typename MatrixFunction>
bool lnsearch(const CF_on_face<T, MatrixFunction>& CF, // value function
              vector2<T>& x, // current state
              vector7<T>& f, // value at current state
              const vector2<T>& dd, // tangent at current state
              T ml=0.1) 
{
    T lambda = 1.0;
    const T alpha = 1e-4;
    const T maxl = max_length(x, dd, ml);
    
    vector2<T> x0 = x;
    vector7<T> f0 = f;
    vector2<T> d = (norm(dd) > maxl) ? dd * maxl / norm(dd) : dd;
    T v0 = spurt::norm(f0);
    for (unsigned int i = 0; i < 20; ++i) {
        x = x0 + lambda * d;
        f = CF(x);
        if (spurt::norm(f) < (1 - alpha*lambda)*v0) {
            return true;
        }
        
        lambda *= 0.5;
    }
    
    return false;
}

template<typename T, typename MatrixFunction>
bool findLPoint(vector3<T>& Lpt, const Face<T>& face,
                const MatrixFunction& hessian, int maxiter=100)
{   
    // we will work in local voxel coordinates
    CF_on_face<T, MatrixFunction> CFf(hessian, face);

    vector2<T> X(0.5, 0.5); 
    vector7<T> f = CFf(X);
    for (int i=0 ; i<maxiter; ++i) {
        vector2<T> d = -CFf.tangent(X);
        if (!lnsearch(CFf, X, f, d)) {
            break;
        }
#ifdef SPURT_DEBUG
        std::cout << "current norm is " << spurt::norm(f) << '\n';
#endif
        if (spurt::norm(f) < 1.0e-6) {
#ifdef SPURT_DEBUG
            std::cout << "Solution found\n";
#endif
            Lpt = face(X);
            return true;
        }
    }
    return false;
}

} // namespace zheng

} // namespace spurt