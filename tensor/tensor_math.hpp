#ifndef __TENSOR_MATH_HPP__
#define __TENSOR_MATH_HPP__

#include <math/types.hpp>

namespace spurt {

template<typename T>
spurt::mat2 to_matrix(const spurt::small_vector<T, 3>& t)
{
    spurt::mat2 m;
    m(0, 0) = t[0];
    m(1, 0) = m(0, 1) = t[1];
    m(1, 1) = t[2];
    return m;
}

template<typename T>
spurt::mat3 to_matrix(const spurt::small_vector<T, 7>& t)
{
    spurt::mat3 m;
    m(0, 0) = t[1];
    m(1, 0) = m(0, 1) = t[2];
    m(2, 0) = m(0, 2) = t[3];
    m(1, 1) = t[4];
    m(2, 1) = m(1, 2) = t[5];
    m(2, 2) = t[6];
    return m;
}

template<typename T>
spurt::mat3 to_matrix(const spurt::small_vector<T, 6>& t)
{
    spurt::mat3 m;
    m(0, 0) = t[0];
    m(1, 0) = m(0, 1) = t[1];
    m(2, 0) = m(0, 2) = t[2];
    m(1, 1) = t[3];
    m(2, 1) = m(1, 2) = t[4];
    m(2, 2) = t[5];
    return m;
}

inline spurt::vec3 PQR(const spurt::mat3& m)
{
    spurt::vec3 r;
    r[0] = trace(m); // P
    r[1] = m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1) +
           m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1) +
           m(2, 2) * m(0, 0) - m(0, 2) * m(2, 0); // Q
    r[2] = determinant(m); // R;
    return r;
}

inline double discriminant(const spurt::mat3& m)
{
    spurt::vec3 r = PQR(m);
    double P = r[0];
    double Q = r[1];
    double R = r[2];
    return Q*Q*P*P - 4.*R*P*P*P - 4.*Q*Q*Q + 18.*P*Q*R - 27.*R*R;
}


}


#endif

