#ifndef __XAVIER_MATH_HPP__
#define __XAVIER_MATH_HPP__

#include <complex>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>

namespace spurt {
template< typename T >
bool LtIndexed(const T& v1, const T& v2)
{
    return (v1.first < v2.first);
}

template< typename T >
bool GtIndexed(const T& v1, const T& v2)
{
    return (v1.first > v2.first);
}

template<typename Index, typename ArrayLike, typename Comparison=std::less<typename ArrayLike::value_type>>
inline void argsort(std::vector<Index>& ids, 
                    const ArrayLike& _array,
                    Comparison cmp=Comparison())
{
    ids.resize(std::distance(_array.begin(), _array.end()));
    std::iota(ids.begin(), ids.end(), 0);
    std::sort(ids.begin(), ids.end(), [&](Index a, Index b)
    {
        return cmp(_array[a], _array[b]);
    });
}

template< typename T >
inline void sort_ids(std::vector< unsigned int >& sorted,
                     const std::vector< T >& values, bool increasing = true)
{
    typedef std::pair< T, unsigned int > indexed_type;
    
    std::vector< indexed_type > tmp(values.size());
    for (unsigned int i = 0 ; i < values.size() ; i++) {
        tmp[i] = indexed_type(values[i], i);
    }
    
    if (increasing) {
        std::sort(tmp.begin(), tmp.end(), LtIndexed< indexed_type >);
    } else {
        std::sort(tmp.begin(), tmp.end(), GtIndexed< indexed_type >);
    }
    sorted.resize(tmp.size());
    for (unsigned int i = 0 ; i < sorted.size() ; i++) {
        sorted[i] = tmp[i].second;
    }
}

template<typename I>
inline I ipow(I n, I k)
{
    I result = 1;
    while (k) {
        if (k & 1) {
            result *= n;
        }
        k >>= 1;
        n *= n;
    }
    
    return result;
}

template< typename T >
inline unsigned int min_id(const std::vector< T >& values)
{
    return std::distance(values.begin(), std::min_element(values.begin(), values.end()));
}

template< typename T >
inline unsigned int max_id(const std::vector< T >& values)
{
    return std::distance(values.begin(), std::max_element(values.begin(), values.end()));
}

template<typename Vec>
inline bool bary_coords(Vec& b, const Vec& p,
                        const Vec& p0, const Vec& p1, const Vec& p2)
{
    // bring problem back to 2D
    Vec q, e1(p1), e2(p2);
    double x, y, x1, x2, y2, /*det,*/ delta;
    
    e1 -= p0;
    e2 -= p0;
    x1 = norm(e1);
    e1 *= 1 / x1; // e1 normalized
    x2 = inner(e2, e1);
    e2 -= x2 * e1; // e2 orthogonal to e1
    y2 = norm(e2);
    e2 *= 1 / y2;
    
    q = p - p0;
    x = inner(q, e1);
    y = inner(q, e2);
    
    delta = x1 * y2; // NB: y1==0
    b[1] = (x * y2 - y * x2) / delta;
    b[2] = x1 * y / delta;
    b[0] = 1.0 - b[1] - b[2];
    
#if 0
    if (b[0] >= 0 && b[1] >= 0 && b[2] >= 0) {
        std::cout
                << "p = " << p << std::endl
                << "p0 = " << p0 << std::endl
                << "p1 = " << p1 << std::endl
                << "p2 = " << p2 << std::endl
                << "e1 = " << e1 << std::endl
                << "e2 = " << e2 << std::endl
                << "b0 = " << b[0] << ", b1 = " << b[1] << ", b2 = " << b[2]
                << std::endl;
    }
#endif
    
    return (b[0] >= 0 && b[1] >= 0 && b[2] >= 0);
}


inline int quadratic_equation(double a, double b, double c,
                              std::complex<double> x[2])
// PAR: a*x^2+b*x+c=0.0, x[0..1] gets the complex solutions
// POST: x[0], x[1] contain the solution
// RETURN: number of different solutions
//         (-1,0,1,2) -1 stands for infinite many!
{
    // mario
    if (fabs(a) < 1e-9) {
        // linear equation
        if (fabs(b) < 1e-9) {
            return 0;
        }
        // FIXME: there may be infinite, how to deal with that?
        x[0] = std::complex<double>(-c / b, 0.);
        return 1;
    }
    // end mario
    
    double d = b * b - 4 * a * c;
    if (d > 0.0) {
        double q = b >= 0.0 ? -0.5 * (b + sqrt(d)) : -0.5 * (b - sqrt(d));
        x[0] = (q / a); // x[0]=q/a+i*0
        x[1] = (c / q); // x[1]=c/q+i*0
        return 2;
    }
    if (d == 0.0) {
        if (a != 0.0) {
            x[0] = std::complex<double>(-0.5 * b / a);
            //although there is only one solution
            //set x[1] for convenience to the same value as x[0]
            x[1] = x[0];
            return 1;
        } else if (c == 0.0) {
            return -1;
        } else {
            return 0;
        }
    }
    // ASSERT: d<0.0
    std::complex<double> dd = d;
    dd = sqrt(dd);
    dd = (real(b * dd) >= 0) ? dd : -dd;
    std::complex<double> q = -0.5 * (b + dd);
    x[0] = (q / a);
    x[1] = (c / q);
    return 2;
}

inline int cubic_equation(double a3, double a2, double a1, double a0,
                          std::complex<double> x[3])
// PAR: a3*x^3+a2*x^2+a1*x+a0=0.0, x[0..2] gets the solution
// POST: x[0], x[1], x[2] contain the solutions
// RETURN: number of different solutions
//         (-1,0,1,2,3) -1 stands for infinite many!
// REMARK:  ideas taken from "Numerical Recipes in C", p.184/185
{
    if (a3 == 0.0) {
        int n = quadratic_equation(a2, a1, a0, x);
        x[2] = x[1];
        return n;
    }
    a2 /= a3;
    a1 /= a3;
    a0 /= a3;
    double Q = (a2 * a2 - 3.0 * a1) / 9.0;
    double R = (2.0 * a2 * a2 * a2 - 9.0 * a2 * a1 + 27.0 * a0) / 54.0;
    double comp = R * R - Q * Q * Q;
    if (comp >= 0.0) {
        double sgn_R = (R >= 0.0) ? 1.0 : -1.0;
        double A = fabs(R) + sqrt(comp);
        A = pow(A, 1.0 / 3.0);
        A *= (-sgn_R);
        double B = (A != 0.0) ? Q / A : 0.0;
        x[0] = std::complex<double>((A + B) - a2 / 3.0);
        x[1] = std::complex<double>(-0.5 * (A + B) - a2 / 3.0, 0.5 * sqrt(3.0) * (A - B));
        x[2] = std::complex<double>(-0.5 * (A + B) - a2 / 3.0, -0.5 * sqrt(3.0) * (A - B));
        return 3;
    }
    double theta = acos(R / sqrt(Q * Q * Q));
    x[0] = std::complex<double>(-2.0 * sqrt(Q) * cos(theta / 3.0) - a2 / 3);
    x[1] = std::complex<double>(-2.0 * sqrt(Q) * cos((theta + 2 * M_PI) / 3.0) - a2 / 3);
    x[2] = std::complex<double>(-2.0 * sqrt(Q) * cos((theta - 2 * M_PI) / 3.0) - a2 / 3);
    return 3;
}

}

#endif








