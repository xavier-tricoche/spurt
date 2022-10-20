#include <math/math.hpp>
#include <teem/ell.h>
#include <math/fixed_vector.hpp>

namespace spurt {
template<typename V, typename T>
inline T outer(const V& v1, const V& v2)
{
    int N = v1.size();
    T out;
    for (int i = 0 ; i < N ; ++i) {
        for (int j = 0 ; j < N ; ++j) {
            out(i, j) = v1[i] * v2[j];
        }
    }
    return out;
}

template < typename T >
inline double frob_norm(const T& tensor)
{
    double n = 0;
    int N = tensor.size();
    for (int i = 0 ; i < N ; ++i) {
        for (int j = 0 ; j < N ; ++j) {
            n += tensor(i, j) * tensor(i, j);
        }
    }
    return sqrt(n);
}

template < typename T >
inline double trace(const T& tensor)
{
    double tr = 0;
    int N = tensor.size();
    for (int i = 0 ; i < N ; ++i) {
        tr += tensor(i, i);
    }
    return tr;
}

template<typename T>
inline T deviatoric(const T& tensor)
{
    T D(tensor);
    double lambda = trace(tensor) / (float)tensor.size();
    for (int i = 0 ; i < 3 ; ++i) {
        D(i, i) -= lambda;
    }
    return D;
}

template < typename T >
inline double FA(const T& tensor)
{
    typedef T   tensor_type;
    
    tensor_type D = deviatoric(tensor);
    return sqrt(1.5)*frob_norm(D) / frob_norm(tensor);
}

template < typename T >
inline nvis::vec3 westin_aniso(const T& tensor)
{
    double tr = trace(tensor);
    T D = deviatoric(tensor);
    std::vector<double> eval;
    eigenvalues(eval, D);
    
    double cl = eval[0] - eval[1];
    double cp = 2.*(eval[1] - eval[2]);
    double cs = 3.*eval[2];
    nvis::vec3 ans(cl, cp, cs);
    return ans / tr;
}

}












