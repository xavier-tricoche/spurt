#pragma once

#include <teem/ell.h>
#include <math/types.hpp>
#include <ostream>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <assert.h>

namespace spurt {

template<typename T, int N=3>
inline double det(const fixed_matrix<T, N>& M)
{
    return M.det();
}

template<typename T, int N=3>
inline fixed_matrix<T, N> invert(const fixed_matrix<T, N>& M)
{
    return M.invert();
}

template< typename T, int N=3>
inline fixed_matrix<T, N> transpose(const fixed_matrix<T, N>& M)
{
    return M.transpose();
}

template< typename T, int N=3 >
inline int eigenvalues(std::vector< double >& evals,
                       const fixed_matrix<T, N>& M)
{
    typedef fixed_matrix<T, N> matrix_t;
    Eigen::EigenSolver<matrix_t> solver(M);
    auto _evals = solver.eigenvalues();
    evals.clear();
    for (int i=0; i<N; ++i) {
        if (_evals[i].imag() == 0) {
            evals.push_back(_evals[i].real());
        }
    }
    std::sort(evals.begin(), evals.end(), [&](double a, double b) {
        return b < a;
    });
    return evals.size();
}

template <typename T, int N = 3>
inline int eigensystem(std::vector<T> &evals,
                       std::vector<fixed_vector<T, 3>> &evecs,
                       const fixed_matrix<T, N> &M)
{
    typedef fixed_matrix<T, N> matrix_t;
    typedef fixed_vector<T, N> vec_t;
    Eigen::EigenSolver<matrix_t> solver(M);

    evecs.clear();
    evals.clear();
    auto _evals = solver.eigenvalues();
    auto _evecs = solver.eigenvectors();
    for (int i=0; i<N; ++i) 
    {
        if (evals[i].imag() == 0) {
            evals[i] = _evals[0].real();
        }
    }

    std::vector<int> ids(evals.size());
    std::iota(ids.begin(), ids.end(), 0);
    if (evals.size() > 1) {
        std::sort(ids.begin(), ids.end(), [&](int a, int b) {
            return evals[b] < evals[a];
        });
    }
    std::vector<double> sorted_evals;
    for (int i=0; i<ids.size(); ++i) {
        sorted_evals.push_back(evals[ids[i]]);
        evecs.push_back(_evecs.col(ids[i]));
    }
    evals.swap(sorted_evals);
    return evals.size();
}

template<typename T, int N=3>
inline void eigen(std::vector< double >& evals,
                  std::vector< fixed_vector<T, N> >& evecs,
                  const fixed_matrix<T, N>& M)
{
    auto blah = eigensystem(evecs, evals, M);
}

} // namespace spurt
