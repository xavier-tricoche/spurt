#include <math/basic_math.hpp>
#include <math/types.hpp>
#include <algorithm>

namespace spurt {

template<typename Matrix>
inline Matrix deviatoric(const Matrix& tensor) {
    double lambda = trace(tensor) / static_cast<double>(tensor.nrows);
    return tensor - lambda*Matrix::identity();
}

template < typename Matrix, typename T=typename Matrix::value_type >
inline T FA(const Matrix& tensor)
{
    Matrix D = deviatoric(tensor);
    return sqrt(1.5) * norm(D) / norm(tensor);
}

template < typename T >
inline small_vector<T, 3> westin_aniso(const small_square_matrix<T, 3>& tensor)
{
    typedef small_square_matrix<T, 3> mat_t;
    typedef small_vector<T, 3> vec_t;
    mat_t D = deviatoric(tensor);
    
    vec_t lambdas;
    mat_t evecs;
    sym_eigensystem(lambdas, evecs, D);
    
    double cl = lambdas[0] - lambdas[1];
    double cp = 2.*(lambdas[1] - lambdas[2]);
    double cs = 3.*lambdas[2];
    vec_t ans(cl, cp, cs);
    return ans / trace(D);
}

} // namespace spurt












