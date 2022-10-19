#ifndef __MLS_HPP__
#define __MLS_HPP__

#include <vector>
#include <set>
#include <list>
#include <math.h>
#include <iostream>
#include <sstream>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/timer.hpp>

// using LAPACK for SVD-based LS solution
#ifdef __INTEL_MKL__
#include "lapacke.h"
#else
#include "clapack.h"
#endif

#include <Eigen/Core>
#include <Eigen/SVD>

namespace MLS {
// number of degrees of freedom for requested polynomial precision

static size_t factorial(size_t n)
{
    if (n = 0) {
        return 1;
    }
    return n*factorial(n-1);
}

static size_t binomial(size_t n, size_t k)
{
    return factorial(n)/(factorial(n-k)*factorial(k));
}

static size_t dof(unsigned int dim, unsigned int prec)
{
    size_t n = 1;
    for (int k=1 ; k<=prec ; ++k) {
        n += binomial(dim, k);
    }
    return n;
}

inline unsigned int best_approx_order(int npts, int dim)
{
    unsigned int prec = 1;
    for (; prec<6 ; ++prec) {
        if (dof(dim, prec) > npts) {
            break;
        }
    }
    return prec-1;
}

// smooth (C^2) in dim 2 and 3
inline double wendland(double r, double radius)
{
    double __r = r / radius;
    if (__r >= 1.) {
        return 0;
    }
    return pow(1. - __r, 4.)*(4.*__r + 1);
}

inline double franke_nielson(double r, double radius)
{
    double d = radius - r;
    if (d <= 0) {
        return 0;
    }
    d /= radius * r;
}

template< typename V >
static inline void linear_2d(std::vector<double>& b, const V& p)
{
    b[0] = 1.;
    b[1] = p[0];
    b[2] = p[1];
}

template< typename V >
static inline void linear_3d(std::vector<double>& b, const V& p)
{
    b[0] = 1.;
    b[1] = p[0];
    b[2] = p[1];
    b[3] = p[2];
}

template< typename V >
static inline void quadratic_2d(std::vector<double>& b, const V& p)
{
    linear_2d(b, p);
    b[3] = p[0] * p[0];
    b[4] = p[0] * p[1];
    b[5] = p[1] * p[1];
}

template< typename V >
static inline void quadratic_3d(std::vector<double>& b, const V& p)
{
    linear_3d(b, p);
    b[4] = p[0] * p[0];
    b[5] = p[0] * p[1];
    b[6] = p[0] * p[2];
    b[7] = p[1] * p[1];
    b[8] = p[1] * p[2];
    b[9] = p[2] * p[2];
}

template< typename V >
inline void cubic_2d(std::vector<double>& b, const V& p)
{
    quadratic_2d(b, p);
    b[6] = p[0] * p[0] * p[0];
    b[7] = p[0] * p[0] * p[1];
    b[8] = p[0] * p[1] * p[1];
    b[9] = p[1] * p[1] * p[1];
}

template< typename V >
inline void cubic_3d(std::vector<double>& b, const V& p)
{
    quadratic_3d(b, p);
    b[10] = p[0] * p[0] * p[0];
    b[11] = p[0] * p[0] * p[1];
    b[12] = p[0] * p[0] * p[2];
    b[13] = p[0] * p[1] * p[1];
    b[14] = p[0] * p[1] * p[2];
    b[15] = p[0] * p[2] * p[2];
    b[16] = p[1] * p[1] * p[1];
    b[17] = p[1] * p[1] * p[2];
    b[18] = p[1] * p[2] * p[2];
    b[19] = p[2] * p[2] * p[2];
}

template< typename V >
inline void quartic_2d(std::vector<double>& b, const V& p)
{
    cubic_2d(b, p);
    b[10] = p[0] * p[0] * p[0] * p[0];
    b[11] = p[0] * p[0] * p[0] * p[1];
    b[12] = p[0] * p[0] * p[1] * p[1];
    b[13] = p[0] * p[1] * p[1] * p[1];
    b[14] = p[1] * p[1] * p[1] * p[1];
}

template< typename V >
inline void quartic_3d(std::vector<double>& b, const V& p)
{
    cubic_3d(b, p);
    b[20] = p[0] * p[0] * p[0] * p[0];
    b[21] = p[0] * p[0] * p[0] * p[1];
    b[22] = p[0] * p[0] * p[0] * p[2];
    b[23] = p[0] * p[0] * p[1] * p[1];
    b[24] = p[0] * p[0] * p[1] * p[2];
    b[25] = p[0] * p[0] * p[2] * p[2];
    b[26] = p[0] * p[1] * p[1] * p[1];
    b[27] = p[0] * p[1] * p[1] * p[2];
    b[28] = p[0] * p[1] * p[2] * p[2];
    b[29] = p[0] * p[2] * p[2] * p[2];
    b[30] = p[1] * p[1] * p[1] * p[1];
    b[31] = p[1] * p[1] * p[1] * p[2];
    b[32] = p[1] * p[1] * p[2] * p[2];
    b[33] = p[1] * p[2] * p[2] * p[2];
    b[34] = p[2] * p[2] * p[2] * p[2];
}

template< typename V >
inline void quintic_2d(std::vector<double>& b, const V& p)
{
    quartic_2d(b, p);
    b[15] = p[0] * p[0] * p[0] * p[0] * p[0];
    b[16] = p[0] * p[0] * p[0] * p[0] * p[1];
    b[17] = p[0] * p[0] * p[0] * p[1] * p[1];
    b[18] = p[0] * p[0] * p[1] * p[1] * p[1];
    b[19] = p[0] * p[1] * p[1] * p[1] * p[1];
    b[20] = p[1] * p[1] * p[1] * p[1] * p[1];
}

template< typename V >
inline void quintic_3d(std::vector<double>& b, const V& p)
{
    quartic_3d(b, p);
    b[35] = p[0] * p[0] * p[0] * p[0] * p[0];
    b[36] = p[0] * p[0] * p[0] * p[0] * p[1];
    b[37] = p[0] * p[0] * p[0] * p[0] * p[2];
    b[38] = p[0] * p[0] * p[0] * p[1] * p[1];
    b[39] = p[0] * p[0] * p[0] * p[1] * p[2];
    b[40] = p[0] * p[0] * p[0] * p[2] * p[2];
    b[41] = p[0] * p[0] * p[1] * p[1] * p[1];
    b[42] = p[0] * p[0] * p[1] * p[1] * p[2];
    b[43] = p[0] * p[0] * p[1] * p[2] * p[2];
    b[44] = p[0] * p[0] * p[2] * p[2] * p[2];
    b[45] = p[0] * p[1] * p[1] * p[1] * p[1];
    b[46] = p[0] * p[1] * p[1] * p[1] * p[2];
    b[47] = p[0] * p[1] * p[1] * p[2] * p[2];
    b[48] = p[0] * p[1] * p[2] * p[2] * p[2];
    b[49] = p[0] * p[2] * p[2] * p[2] * p[2];
    b[50] = p[1] * p[1] * p[1] * p[1] * p[1];
    b[51] = p[1] * p[1] * p[1] * p[1] * p[2];
    b[52] = p[1] * p[1] * p[1] * p[2] * p[2];
    b[53] = p[1] * p[1] * p[2] * p[2] * p[2];
    b[54] = p[1] * p[2] * p[2] * p[2] * p[2];
    b[55] = p[2] * p[2] * p[2] * p[2] * p[2];
}

template<typename V>
inline void set_basis(std::vector<double>& b, const V& p, int dim, int prec)
{
    if (prec == 0) {
        b[0] = 1.;
    } else if (dim == 2) {
        switch (prec) {
            case 1:
                linear_2d(b, p);
                break;
            case 2:
                quadratic_2d(b, p);
                break;
            case 3:
                cubic_2d(b, p);
                break;
            case 4:
                quartic_2d(b, p);
                break;
            case 5:
                quintic_2d(b, p);
                break;
            default:
                assert(false);
        }
    } else {
        switch (prec) {
            case 1:
                linear_3d(b, p);
                break;
            case 2:
                quadratic_3d(b, p);
                break;
            case 3:
                cubic_3d(b, p);
                break;
            case 4:
                quartic_3d(b, p);
                break;
            case 5:
                quintic_3d(b, p);
                break;
            default:
                assert(false);
        }
    }
}

#ifndef __INTEL_MKL__
typedef __CLPK_integer      integer_type;
typedef __CLPK_doublereal   float_type;
#else
typedef lapack_int          integer_type;
typedef double              float_type;
#endif

struct CLAPACK_helper {
    float_type* A, *b, *work, *s, rcond;
    integer_type ncols, nrows, nrhs, lwork, info, rank;
    bool svd;
    char trans;
    unsigned int dim, prec;
    int verbose;
    
    CLAPACK_helper(unsigned int max_nvals,
                   unsigned int _dim,
                   unsigned int _prec, unsigned int _nrhs,
                   bool use_svd = false) : verbose(0) {
        dim = _dim;
        prec = _prec;
        nrows = max_nvals;
        ncols = MLS::dof(dim,prec);
        nrhs = _nrhs;
        lwork = -1;
        trans = 'N';
        svd = use_svd;
        A = (float_type*)malloc(nrows * ncols * sizeof(float_type));
        b = (float_type*)malloc(nrows * nrhs * sizeof(float_type));
        if (use_svd) {
            s = (float_type*)malloc(ncols * sizeof(float_type));
            rcond = 1.0e-4;
        }
        work = (float_type*)malloc(10 * sizeof(float_type));
        if (!use_svd) {
            dgels_(&trans, &nrows, &ncols, &nrhs, A, &nrows, b, &nrows, work, &lwork, &info);
        } else
            dgelss_(&nrows, &ncols, &nrhs, A, &nrows, b, &nrows, s, &rcond, &rank,
                    work, &lwork, &info);
        lwork = (integer_type)work[0];
        std::cout << "optimal work space is " << lwork << std::endl;
        delete[] work;
        work = (float_type*)malloc(lwork * sizeof(float_type));
    }
    
    ~CLAPACK_helper() {
        delete[] A;
        delete[] b;
        delete[] work;
        delete[] s;
    }
};

template<typename T>
struct distance_traits {};

template<typename T, nvis::size_type N>
struct distance_traits< nvis::fixed_vector<T, N> > {
    typedef double                      value_type;
    typedef nvis::fixed_vector<T, N>    vec_type;
    
    static value_type dist(const vec_type& v0, const vec_type& v1) {
        return nvis::norm(v1-v0);
    }
};

template<typename T, unsigned int N>
struct distance_traits< Eigen::Matrix<T, N, 1> > {
    typedef double                      value_type;
    typedef Eigen::Matrix<T,N,1>        vec_type;
    
    static value_type dist(const vec_type& v0, const vec_type& v1) {
        return (v1-v0).norm();
    }
};

template<typename ValType, typename PosType,
         typename DistanceTraits = distance_traits<PosType> >
class weighted_least_squares {
public:
    weighted_least_squares(int dimension, int precision, int nrhs)
        : _dim(dimension), _prec(precision), _nrhs(nrhs) {}
        
    int operator()(Eigen::MatrixXd& fitted_coef,
                   const std::vector<PosType>& points,
                   const std::vector<ValType>& values,
                   const PosType& x0, double radius) const {
                   
        assert(_prec <= 5);
        
        typedef ValType         value_type;
        typedef PosType         pos_type;
        typedef DistanceTraits  dtraits;
        
        unsigned int npts = points.size();
        unsigned int prec = std::min((unsigned int)_prec, (unsigned int)best_approx_order(npts, _dim));
        unsigned int nbasisfn = dof(_dim, prec);
        
        if (!npts) {
            return -1;
        }
        
        std::vector<double> basis(nbasisfn);
        fitted_coef.resize(nbasisfn, _nrhs);
        
        Eigen::MatrixXd A(npts, nbasisfn);
        Eigen::MatrixXd rhs(npts, _nrhs);
        
        distance_traits<nvis::fixed_vector<double, 3> > traits_instance;
        
        for (int i=0 ; i<npts ; ++i) {
            double r = traits_instance.dist(x0, points[i]);
            double w = sqrt(MLS::wendland(r, radius));
            MLS::set_basis(basis, points[i], _dim, prec);
            for (int c=0 ; c<nbasisfn ; ++c) {
                A(i, c) = w * basis[c];
            }
            for (int c=0 ; c<_nrhs ; ++c) {
                rhs(i, c) = w * values[i][c];
            }
        }
        
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        fitted_coef = svd.solve(rhs);
        
        return prec;
    }
    
protected:
    int _dim, _prec, _nrhs;
};

template<typename ValType, typename PosType>
class least_squares {
public:
    least_squares(int dimension, int precision, int nrhs)
        : _dim(dimension), _prec(precision), _nrhs(nrhs) {}
        
    int operator()(Eigen::MatrixXd& fitted_coef,
                   const std::vector<PosType>& points,
                   const std::vector<ValType>& values) const {
                   
        assert(_prec <= 5);
        
        typedef ValType         value_type;
        typedef PosType         pos_type;
        
        unsigned int npts = points.size();
        unsigned int prec = std::min((unsigned int)_prec, (unsigned int)best_approx_order(npts, _dim));
        unsigned int nbasisfn = dof(_dim, prec);
        
        if (!npts) {
            return -1;
        }
        
        std::vector<double> basis(nbasisfn);
        fitted_coef.resize(nbasisfn, _nrhs);
        
        Eigen::MatrixXd A(npts, nbasisfn);
        Eigen::MatrixXd rhs(npts, _nrhs);
        
        distance_traits<nvis::fixed_vector<double, 3> > traits_instance;
        
        for (int i=0 ; i<npts ; ++i) {
            MLS::set_basis(basis, points[i], _dim, prec);
            for (int c=0 ; c<nbasisfn ; ++c) {
                A(i, c) = basis[c];
            }
            for (int c=0 ; c<_nrhs ; ++c) {
                rhs(i, c) = values[i][c];
            }
        }
        
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        fitted_coef = svd.solve(rhs);
        
        return prec;
    }
    
protected:
    int _dim, _prec, _nrhs;
};

}



#endif