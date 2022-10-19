#ifndef __TEEM_MATRIX_HPP__
#define __TEEM_MATRIX_HPP__

#include <teem/ell.h>
#include <math/fixed_vector.hpp>
#include <ostream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <assert.h>

namespace xavier {
extern bool display_eigen_stuff;

// ultra-basic matrix type
template< typename T >
struct matrix {
    matrix() : data(0) {}
    
    matrix(unsigned int m, unsigned int n)
        : _m(m), _n(n) {
        data = (T*)calloc(_m * _n, sizeof(T));
    }
    
    matrix(T* begin, unsigned int m, unsigned int n)
        : _m(m), _n(n) {
        data = (T*)calloc(_m * _n, sizeof(T));
        memcpy(data, begin, _m*_n*sizeof(T));
    }
    
    matrix(const matrix<T>& a)
        : _m(a._m), _n(a._n) {
        data = (T*)calloc(_m * _n, sizeof(T));
        memcpy(data, a.data, _m*_n*sizeof(T));
    }
    
    ~matrix() {
        delete[] data;
    }
    
    int size() const {
        if (_m == _n) {
            return _m;
        }
        throw std::runtime_error("xavier::matrix::size(): matrix is not symmetric");
    }
    
    void resize(unsigned int m, unsigned int n) {
        if (m == _m && n == _n) {
            return;
        }
        delete[] data;
        _m = m;
        _n = n;
        data = (T*)calloc(m * n, sizeof(T));
    }
    
    T& operator()(unsigned int i, unsigned int j) {
        assert(i < _m && j < _n);
        return data[j*_m+i];
    }
    
    const T& operator()(unsigned int i, unsigned int j) const {
        assert(i < _m && j < _n);
        return data[j*_m+i];
    }
    
    matrix<T>& operator=(const matrix<T>& A) {
        if (A._m != _m && A._n != _n) {
            resize(A._m, A._n);
        }
        for (unsigned int i = 0 ; i < _m*_n ; ++i) {
            data[i] = A.data[i];
        }
        return *this;
    }
    
    matrix<T> operator*(const matrix<T>& A) const {
        assert(A._m == _n);
        matrix<T> res(_m, A._n);
        for (unsigned int i = 0 ; i < _m ; ++i) {
            for (unsigned int j = 0 ; j < A._n ; ++j) {
                for (unsigned int k = 0 ; k < _n ; ++k) {
                    res(i, j) += (*this)(i, k) * A(k, j);
                }
            }
        }
        return res;
    }
    
    matrix<T>& operator+=(const matrix<T>& A) {
        assert(A._m == _m && A._n == _n);
        for (unsigned int i = 0 ; i < _m*_n ; ++i) {
            data[i] += A.data[i];
        }
        return *this;
    }
    
    unsigned int _m, _n;
    T* data;
};

// Simple helper class for teem/ell (3x3) matrix class
class mat3 {
public:
    mat3() {
        std::fill(__mat, __mat + 9, 0);
    }
    
    mat3(const mat3& M) {
        std::copy(M.__mat, M.__mat + 9, __mat);
    }
    
    ~mat3() {
    }
    
    static int size() {
        return 3;
    }
    
    const double& operator()(unsigned int i, unsigned int j) const {
        // row first
        return __mat[j+i*3];
    }
    
    double& operator()(unsigned int i, unsigned int j) {
        // row first
        return __mat[j+i*3];
    }
    
    mat3 operator*(const mat3& M) const {
        mat3 _m;
        ell_3m_mul_d(_m.__mat, __mat, M.__mat);
        return _m;
    }
    
    nvis::vec3 operator*(const nvis::vec3& v) const {
        nvis::vec3 _w;
        ell_3mv_mul_d(_w.begin(), __mat, v.begin());
        return _w;
    }
    
    mat3& operator*=(mat3& M) {
        ell_3m_mul_d(__mat, __mat, M.__mat);
        return *this;
    }
    
    mat3& operator+=(const mat3& M) {
        for (unsigned int i = 0 ; i < 9 ; i++) {
            __mat[i] += M.__mat[i];
        }
        return *this;
    }
    
    mat3 operator+(const mat3& M) const {
        mat3 _m(*this);
        _m += M;
        return _m;
    }
    
    mat3& operator-=(const mat3& M) {
        for (unsigned int i = 0 ; i < 9 ; i++) {
            __mat[i] += M.__mat[i];
        }
        return *this;
    }
    
    mat3 operator-(const mat3& M) const {
        mat3 _m(*this);
        _m -= M;
        return _m;
    }
    
    void eigensystem(std::vector< double >& evals,
                     std::vector< nvis::vec3 >& evecs) const;
                     
    friend double det(const mat3& M);
    friend mat3 invert(const mat3& M);
    friend mat3 transpose(const mat3& M);
    friend int eigenvalues(std::vector< double >& evals, const mat3& M);
    friend int eigensystem(std::vector<nvis::vec3>& evecs,
                           std::vector<double>& evals,
                           const mat3& M);
    friend void eigen(std::vector< double >& evals,
                      std::vector< nvis::vec3 >& evecs, const mat3& M);
    friend nvis::vec3 eigenvector(const mat3& M, const double lambda);
    
private:
    double __mat[9];
};

nvis::vec3 eigenvector(const mat3& M, const double lambda);

inline double norm(const mat3& M)
{
    double frob = 0;
    for (unsigned int i = 0 ; i < 3 ; i++)
        for (unsigned int j = 0 ; j < 3 ; j++) {
            frob += M(i, j) * M(i, j);
        }
    return sqrt(frob);
}

inline double det(const mat3& M)
{
    return ell_3m_det_d(const_cast<double *>(M.__mat));
}

inline mat3 invert(const mat3& M)
{
    mat3 _m;
    ell_3m_inv_d(_m.__mat, M.__mat);
    return _m;
}

inline mat3 transpose(const mat3& M)
{
    mat3 _m;
    for (unsigned int i = 0 ; i < 3 ; i++) {
        for (unsigned int j = i + 1 ; j < 3 ; j++) {
            _m(i, j) = M(j, i);
            _m(j, i) = M(i, j);
        }
        _m(i, i) = M(i, i);
    }
    return _m;
}

inline int eigenvalues(std::vector< double >& evals,
                       const mat3& M)
{
    double e[3];
    int ret = ell_3m_eigenvalues_d(e, M.__mat, AIR_TRUE);
    evals.resize(3);
    for (int i = 0 ; i < 3 ; ++i) {
        evals[i] = e[i];
    }
    return ret;
}

inline int eigensystem(std::vector<nvis::vec3>& evecs,
                       std::vector<double>& evals,
                       const mat3& M)
{
    evecs.resize(3);
    evals.resize(3);
    double _eval[3];
    double _evec[9];
    int roots = ell_3m_eigensolve_d(_eval, _evec, M.__mat, 0);
    for (int i = 0 ; i < 3 ; ++i) {
        evals[i] = _eval[i];
        evecs[i] = nvis::vec3(_evec[3*i], _evec[3*i+1], _evec[3*i+2]);
    }
    return roots;
}

inline void eigen(std::vector< double >& evals,
                  std::vector< nvis::vec3 >& evecs,
                  const mat3& M)
{
    double _val[3];
    double _vec[9];
    int roots = ell_3m_eigensolve_d(_val, _vec, M.__mat, 1);
    std::vector< unsigned int > valid;
    switch (roots) {
        case ell_cubic_root_unknown:
            break;
        case ell_cubic_root_single:
            valid.push_back(0);
            break;
        case ell_cubic_root_single_double:
            if (_val[0] > _val[1]) {
                valid.push_back(0);
            } else {
                valid.push_back(1);
            }
            break;
        case ell_cubic_root_three:
            valid.push_back(0);
            valid.push_back(1);
            valid.push_back(2);
            break;
        default:
            break;
    }
    
    evals.clear();
    evecs.clear();
    for (unsigned int j = 0 ; j < valid.size() ; j++) {
        // check results
        unsigned int i = valid[j];
        nvis::vec3 ev(_vec[3*i], _vec[3*i+1], _vec[3*i+2]);
        double lambda = _val[i];
        double dot = nvis::inner(M * ev, ev);
        if ((lambda == 0 && fabs(dot) > 1.0e-6) ||
                (lambda != 0 && fabs((dot - lambda) / lambda) > 1.0e-6)) {
            continue;
        }
        
        evals.push_back(lambda);
        evecs.push_back(ev);
    }
}

inline std::ostream& operator<<(std::ostream& out, const mat3& M)
{
    out << std::endl
        << std::setprecision(30)
        << "(\t" << M(0, 0) << ",\t" << M(0, 1) << ",\t" << M(0, 2) << "\t)" << std::endl
        << "(\t" << M(1, 0) << ",\t" << M(1, 1) << ",\t" << M(1, 2) << "\t)" << std::endl
        << "(\t" << M(2, 0) << ",\t" << M(2, 1) << ",\t" << M(2, 2) << "\t)" << std::endl;
        
    return out;
}

void check_eigen(const mat3& M);
};

#endif
























