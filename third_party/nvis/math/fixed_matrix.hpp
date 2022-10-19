#ifndef __fixed_matrix_hpp
#define __fixed_matrix_hpp

#include <algorithm>
#include <iosfwd>

#include <boost/operators.hpp>
#include <boost/call_traits.hpp>

#include <math/fixed_vector.hpp>
#include <stdexcept>

#include <iomanip>

namespace nvis
{

template < typename T, size_type Nrows, size_type Ncols = Nrows >
class fixed_matrix
{
public:

    typedef fixed_matrix<T, Nrows, Ncols>				self_type;

    typedef fixed_vector<T, Ncols>						row_type;
    typedef fixed_vector<T, Nrows>						col_type;

    typedef typename boost::call_traits<T>::param_type	value_arg;

    typedef T											value_type;

    // --- constructors ---

    fixed_matrix() {
    }

    fixed_matrix(const self_type& m) : __data(m.__data) {
    }

    fixed_matrix(value_arg v) : __data(row_type(v)) {
    }

    fixed_matrix(const fixed_vector<row_type, Nrows>& vv) : __data(vv) {
    }

    row_type& operator[](size_type i) {
        return __data[i];
    }

    const row_type& operator[](size_type i) const {
        return __data[i];
    }

    T& operator()(size_type i, size_type j) {
        return __data[i][j];
    }

    const T& operator()(size_type i, size_type j) const {
        return __data[i][j];
    }

    // --- matrix addition / subtraction ---

    self_type operator+(const self_type& A) const {
        self_type r;
        for (size_type i = 0 ; i < Nrows ; ++i) {
            r[i] = (*this)[i] + A[i];
        }
        return r;
    }

    self_type& operator+=(const self_type& A) {
        for (size_type i = 0 ; i < Nrows ; ++i) {
            (*this)[i] += A[i];
        }
        return *this;
    }

    self_type operator-(const self_type& A) const {
        self_type r;
        for (size_type i = 0 ; i < Nrows ; ++i) {
            r[i] = (*this)[i] - A[i];
        }
        return r;
    }

    self_type& operator-=(const self_type& A) {
        for (size_type i = 0 ; i < Nrows ; ++i) {
            (*this)[i] -= A[i];
        }
        return *this;
    }

    // --- matrix-matrix multiplication ---

    template<size_type N>
    fixed_matrix<T, Nrows, N> operator*(const fixed_matrix<T, Ncols, N>& A) const {
        fixed_matrix<T, Nrows, N> r;
        for (size_type i = 0 ; i < Nrows ; ++i) {
            for (size_type j = 0 ; j < N ; ++j) {
                r[i][j] = 0;
                for (size_type k = 0 ; k < Ncols ; ++k) {
                    r[i][j] += (*this)[i][k] * A[k][j];
                }
            }
        }
        return r;
    }

    // --- matrix-vector multiplication ---

    col_type operator*(const row_type& v) const {
        col_type r;

        for (size_type n = 0; n < Nrows; ++n)
            r[n] = inner((*this)[n], v);

        return r;
    }

    // --- matrix-scalar multiplication ---

    self_type operator*(const T& val) const {
        self_type r(*this);
        r.__data *= val;
    }

    self_type& operator*=(const T& val) {
        __data *= val;
        return *this;
    }

    self_type& operator/=(const T& val) {
        __data /= val;
        return *this;
    }

    // --- identity ---

    static fixed_matrix identity() {
        fixed_matrix id;

        for (unsigned int r = 0; r < Nrows; ++r)
            for (unsigned int c = 0; c < Ncols; ++c)
                id[r][c] = (r == c) ? 1 : 0;

        return id;
    }

private:
    nvis::fixed_vector<row_type, Nrows> __data;
};

// transpose
template<typename T, size_type Nrows, size_type Ncols>
inline fixed_matrix<T, Ncols, Nrows> transpose(const fixed_matrix<T, Nrows, Ncols>& _m) {
    fixed_matrix<T, Ncols, Nrows> tm;
    for (int r = 0; r < Nrows; ++r)
        for (unsigned int c = 0; c < Ncols; ++c)
            tm[c][r] = _m[r][c];
    return tm;
}

// solve linear equation system
// Gauss elimination, with pivoting

template<typename T, size_type N>
fixed_vector<T, N> solve(const fixed_matrix<T, N>& _m,
                         const fixed_vector<T, N>& _v)
{
    fixed_matrix<T, N, N> m(_m);
    fixed_vector<T, N>   v(_v);

    // outer loop on rows
    for (int i = 0; i < N; i++) {
        double pvt = m[i][i];

        if (!pvt) {
            int j;

            for (j = i + 1; j < N; j++) {
                if ((pvt = m[j][i]) != 0.0)
                    break;
            }

            if (!pvt)
                throw std::runtime_error("singular matrix inversion attempted");

            // swap matrix rows and result vector entries
            std::swap(m[j], m[i]);
            std::swap(v[j], v[i]);
        }

        for (int k = i + 1; k < N; k++) {
            double tmp = m[k][i] / pvt;

            for (int j = i + 1; j < N; j++)
                m[k][j] -= tmp * m[i][j];

            v[k] -= tmp * v[i];
        }
    }

    // back substitution
    for (int i = N - 1; i >= 0; i--) {
        v[i] = v[i];

        for (int j = N - 1; j > i; j--)
            v[i] -= m[i][j] * v[j];

        v[i] /= m[i][i];
    }

    return v;
}

// trace

template<typename T, size_type Nrows, size_type Ncols>
T trace(const fixed_matrix<T, Nrows, Ncols>& m)
{
    T tr = 0;
    size_type n = std::min(Nrows, Ncols);
    for (size_type i = 0 ; i < n ; ++i) {
        tr += m[i][i];
    }
    return tr;
}

// determinant by Gauss elimination

template<typename T, size_type N>
typename fixed_matrix<T, N>::value_type det(const fixed_matrix<T, N>& _m)
{
    fixed_matrix<T, N, N> m(_m);

    double det = 1.0;

    // outer loop on rows
    for (int i = 0; i < N; i++) {
        double pvt = m[i][i];

        for (int k = i + 1; k < N; k++) {
            double tmp = m[k][i] / pvt;

            for (int j = i + 1; j < N; j++)
                m[k][j] -= tmp * m[i][j];
        }

        det *= m[i][i];
    }

    return det;
}

// -----------------------------------------------------------------

template<typename T, size_type Nrows, size_type Ncols>
std::ostream& operator<<(std::ostream& out, const fixed_matrix<T, Nrows, Ncols>& m)
{
    out << "\n";
    for (size_type i = 0 ; i < Nrows ; ++i) {
        out << "(";
        for (size_type j = 0 ; j < Ncols ; ++j) {
            out << std::setw(10) << m[i][j] << '\t';
        }
        out << ")\n";
    }

    return out;
}

// -----------------------------------------------------------------

template<typename T, size_type Nrows, size_type Ncols>
fixed_matrix<T, Nrows, Ncols> operator*(const T& val, const fixed_matrix<T, Nrows, Ncols>& m)
{
    fixed_matrix<T, Nrows, Ncols> r(m);
    r *= val;
    return r;
}

// --- Frobenius norm ---

template<typename T, size_type Nrows, size_type Ncols>
T norm(const fixed_matrix<T, Nrows, Ncols>& m)
{
    T n = 0.;
    for (size_type i = 0 ; i < Nrows ; ++i) {
        n += inner(m[i], m[i]);
    }
    return sqrt(n);
}

// --- outer product of vectors ---
template< typename T, size_type Nrows, size_type Ncols >
fixed_matrix< T, Nrows, Ncols > outer(const fixed_vector<T, Nrows>& v1,
                                      const fixed_vector<T, Ncols>& v2)
{
    fixed_matrix<T, Nrows, Ncols> m;
    for (size_type i = 0 ; i < Nrows ; ++i)
        for (size_type j = 0 ; j < Ncols ; ++j)
            m[i][j] = v1[i] * v2[j];

    return m;
}

typedef fixed_matrix<double, 2>  mat2;
typedef fixed_matrix<double, 3>  mat3;
typedef fixed_matrix<double, 4>  mat4;

}  // namespace nvis

#endif // __fixed_matrix_hpp

















