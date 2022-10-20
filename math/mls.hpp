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

// Boost
#include <boost/math/special_functions/binomial.hpp>

// using Eigen SVD-based LS solution
#include <Eigen/Core>
#include <Eigen/SVD>

#include "rbf_basis.hpp"
#include <misc/meta_utils.hpp>

namespace spurt { namespace MLS {

// number of degrees of freedom for requested polynomial precision

static inline size_t dof(size_t dim, size_t prec)
{
    size_t n = 1;
    for (size_t k=1 ; k<=(prec > dim ? dim : prec) ; ++k) {
        n += boost::math::binomial_coefficient<double>(dim, k);
    }
    return n;
}

inline size_t best_approx_order(size_t npts, size_t dim)
{
    size_t prec = 1;
    for (; prec<6 ; ++prec) {
        if (dof(dim, prec) > npts) {
            break;
        }
    }
    return prec-1;
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
    linear_2d(b, p);
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

template<typename ValueType, typename PositionType,
         typename WeightType = spurt::RBF::wendland_function<double>,
         typename DistanceTraits = distance_traits<PositionType> >
class weighted_least_squares {
public:
    typedef ValueType            value_type;
    typedef PositionType         pos_type;
    typedef WeightType           weight_type;
    typedef DistanceTraits       traits_type;

    weighted_least_squares(int dimension, int precision, int nrhs)
        : m_dim(dimension), m_prec(precision), m_nrhs(nrhs) {}

    inline double weight(const pos_type& p, const pos_type& x0, double radius) const {
        traits_type traits;
        weight_type w(radius);
        double r = traits.dist(x0, p);
        return sqrt(w(r));
    }

    int operator()(Eigen::MatrixXd& fitted_coef,
                   const std::vector<pos_type>& points,
                   const std::vector<value_type>& values,
                   const pos_type& x0, double radius) const {

        assert(m_prec <= 5);

        unsigned int npts = points.size();
        unsigned int prec = std::min((unsigned int)m_prec, (unsigned int)best_approx_order(npts, m_dim));
        unsigned int nbasisfn = dof(m_dim, prec);

        if (!npts) {
            return -1;
        }

        std::vector<double> basis(nbasisfn);
        fitted_coef.resize(nbasisfn, m_nrhs);

        Eigen::MatrixXd A(npts, nbasisfn);
        Eigen::MatrixXd rhs(npts, m_nrhs);
        spurt::data_traits<value_type> wrapper;

        for (int i=0 ; i<npts ; ++i) {
            double w = weight(points[i], x0, radius);
            MLS::set_basis(basis, points[i], m_dim, prec);
            for (int c=0 ; c<nbasisfn ; ++c) {
                A(i, c) = w * basis[c];
            }
            for (int c=0 ; c<m_nrhs ; ++c) {
                rhs(i, c) = w * wrapper.value(values[i], c);
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        fitted_coef = svd.solve(rhs);

        return prec;
    }

protected:
    int m_dim, m_prec, m_nrhs;
};

template<typename ValType, typename PosType>
class least_squares {
public:
    typedef ValType value_type;
    typedef PosType pos_type;

    least_squares(int dimension, int precision, int nrhs)
        : m_dim(dimension), m_prec(precision), m_nrhs(nrhs) {}

    int operator()(Eigen::MatrixXd& fitted_coef,
                   const std::vector<pos_type>& points,
                   const std::vector<value_type>& values) const {

        assert(m_prec <= 5);

        unsigned int npts = points.size();
        unsigned int prec = std::min((unsigned int)m_prec, (unsigned int)best_approx_order(npts, m_dim));
        unsigned int nbasisfn = dof(m_dim, prec);

        if (!npts) {
            return -1;
        }

        std::vector<double> basis(nbasisfn);
        fitted_coef.resize(nbasisfn, m_nrhs);

        Eigen::MatrixXd A(npts, nbasisfn);
        Eigen::MatrixXd rhs(npts, m_nrhs);

        distance_traits<nvis::fixed_vector<double, 3> > traits;
        spurt::data_traits<value_type> wrapper;

        for (int i=0 ; i<npts ; ++i) {
            MLS::set_basis(basis, points[i], m_dim, prec);
            for (int c=0 ; c<nbasisfn ; ++c) {
                A(i, c) = basis[c];
            }
            for (int c=0 ; c<m_nrhs ; ++c) {
                rhs(i, c) = wrapper.value(values[i], c);
            }
        }

        Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);
        fitted_coef = svd.solve(rhs);

        return prec;
    }

protected:
    int m_dim, m_prec, m_nrhs;
};

} // MLS
} // spurt



#endif
