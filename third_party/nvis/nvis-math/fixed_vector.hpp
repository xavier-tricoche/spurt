#ifndef __fixed_vector_hpp
#define __fixed_vector_hpp

#include <algorithm>
#include <array>
#include <functional>
#include <numeric>
#include <iosfwd>
#include <iostream>
#include <cmath>
#include <stdexcept>

#include <type_traits>
#include <nvis-math/base_types.hpp>

namespace nvis
{

template<typename T, size_type N>
class fixed_vector
{
public:

    typedef T                           value_type;
    typedef T*                          pointer;
    typedef T&                          reference;
    typedef const T&                    const_reference;
    typedef T*                          iterator;
    typedef const T*                    const_iterator;
    typedef fixed_vector<T, N>          self_type;

protected:

    value_type data_[N];

public:

    // --- constructors ---

    fixed_vector() {
    }

    fixed_vector(const T& v) {
        std::fill(begin(), end(), v);
    }

	fixed_vector(const std::array<T, N>& a) {
		std::copy(a.begin(), a.end(), begin());
	}

    // --- multi-argument initialization constructors

    fixed_vector(const T& v0, const T& v1) {
        static_assert(N >= 2, "Invalid dimension (<2) in fixed_vector constructor");

        iterator i = begin();

        *(i++) = v0;
        *(i++) = v1;

        std::fill(i, end(), T());
    }

    fixed_vector(const T& v0, const T& v1, const T& v2) {
        static_assert(N >= 3, "Invalid dimension (<3) in fixed_vector constructor");

        iterator i = begin();

        *(i++) = v0;
        *(i++) = v1;
        *(i++) = v2;

        std::fill(i, end(), T());
    }

    fixed_vector(const T& v0, const T& v1, const T& v2, const T& v3) {
        static_assert(N >= 4, "Invalid dimension (<4) in fixed_vector constructor");

        iterator i = begin();

        *(i++) = v0;
        *(i++) = v1;
        *(i++) = v2;
        *(i++) = v3;

        std::fill(i, end(), T());
    }

    fixed_vector(const T& v0, const T& v1, const T& v2, const T& v3, const T& v4) {
        static_assert(N >= 5, "Invalid dimension (<5) in fixed_vector constructor");

        iterator i = begin();

        *(i++) = v0;
        *(i++) = v1;
        *(i++) = v2;
        *(i++) = v3;
        *(i++) = v4;

        std::fill(i, end(), T());
    }

    // --- conversion/copy constructor

    template<typename S>
    fixed_vector(const fixed_vector<S, N>& rhs) {
        std::copy(rhs.begin(), rhs.end(), begin());
    }

    // --- element access

    reference operator[](size_type n) {
        return *(begin() + n);
    }

    const_reference operator[](size_type n) const {
        return *(begin() + n);
    }

    // --- iterators ---

    iterator begin() {
        return iterator(data_);
    }
    iterator end()   {
        return iterator(data_ + N);
    }

    const_iterator begin() const {
        return const_iterator(data_);
    }
    const_iterator end()   const {
        return const_iterator(data_ + N);
    }

    // --- size ---

    static size_type size() {
        return N;
    }

    // --- vector-vector operators ---

    self_type& operator+=(const self_type& rhs) {
        const_iterator ri = rhs.begin();

        for (iterator di = begin(); di != end(); ++di, ++ri)
            *di += *ri;

        return *this;
    }

    self_type& operator-=(const self_type& rhs) {
        const_iterator ri = rhs.begin();

        for (iterator di = begin(); di != end(); ++di, ++ri)
            *di -= *ri;

        return *this;
    }

    self_type& operator*=(const self_type& rhs) {
        const_iterator ri = rhs.begin();

        for (iterator di = begin(); di != end(); ++di, ++ri)
            *di *= *ri;

        return *this;
    }

    self_type& operator/=(const self_type& rhs) {
        const_iterator ri = rhs.begin();

        for (iterator di = begin(); di != end(); ++di, ++ri)
            *di /= *ri;

        return *this;
    }

    // --- vector-scalar operators ---

    template<typename S>
    self_type& operator*=(const S& s) {
        for (iterator i = begin(); i != end(); ++i)
            (*i) *= s;

        return *this;
    }

    template<typename S>
    self_type& operator/=(const S& s) {
        for (iterator i = begin(); i != end(); ++i)
            (*i) /= s;

        return *this;
    }
};

// --- operators ---------------------------------------------------

// +

template<typename S, typename T, size_type N> inline
fixed_vector<T, N> operator+(const S& s, const fixed_vector<T, N>& v)
{
    return fixed_vector<T, N>(v) += s;
}

template<typename S, typename T, size_type N> inline
fixed_vector<T, N> operator+(const fixed_vector<T, N>& v, const S& s)
{
    return fixed_vector<T, N>(v) += s;
}

template<typename T, size_type N> inline
fixed_vector<T, N> operator+(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    return fixed_vector<T, N>(v1) += v2;
}

// -

template<typename S, typename T, size_type N> inline
fixed_vector<T, N> operator-(const fixed_vector<T, N>& v, const S& s)
{
    return fixed_vector<T, N>(v) -= s;
}

template<typename S, typename T, size_type N> inline
fixed_vector<T, N> operator-(const S& s, const fixed_vector<T, N>& v)
{
    return fixed_vector<T, N>(s) -= v;
}

template<typename T, size_type N> inline
fixed_vector<T, N> operator-(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    return fixed_vector<T, N>(v1) -= v2;
}

// *

template<typename S, typename T, size_type N> inline
fixed_vector<T, N> operator*(const S& s, const fixed_vector<T, N>& v)
{
    return fixed_vector<T, N>(v) *= s;
}

template<typename S, typename T, size_type N> inline
fixed_vector<T, N> operator*(const fixed_vector<T, N>& v, const S& s)
{
    return fixed_vector<T, N>(v) *= s;
}

template<typename T, size_type N> inline
fixed_vector<T, N> operator*(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    return fixed_vector<T, N>(v1) *= v2;
}

// /

template<typename S, typename T, size_type N> inline
fixed_vector<T, N> operator/(const fixed_vector<T, N>& v, const S& s)
{
    return fixed_vector<T, N>(v) /= s;
}

template<typename S, typename T, size_type N> inline
fixed_vector<T, N> operator/(const S& s, const fixed_vector<T, N>& v)
{
    return fixed_vector<T, N>(s) /= v;
}

template<typename T, size_type N> inline
fixed_vector<T, N> operator/(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    return fixed_vector<T, N>(v1) /= v2;
}

// --- accumulators ------------------------------------------------

template<typename T, size_type N> inline
typename fixed_vector<T, N>::value_type sum(const fixed_vector<T, N>& v)
{
    typename fixed_vector<T, N>::value_type s =
        typename fixed_vector<T, N>::value_type(0.0);

    typename fixed_vector<T, N>::const_iterator i;

    for (i = v.begin(); i < v.end(); ++i)
        s += *i;

    return s;
}

template<typename T, size_type N> inline
typename fixed_vector<T, N>::value_type prod(const fixed_vector<T, N>& v)
{
    typename fixed_vector<T, N>::value_type p =
        typename fixed_vector<T, N>::value_type(1.0);

    typename fixed_vector<T, N>::const_iterator i;

    for (i = v.begin(); i < v.end(); ++i)
        p *= *i;

    return p;
}

template<typename T, size_type N> inline
typename fixed_vector<T, N>::value_type inner(const fixed_vector<T, N>& v1,
        const fixed_vector<T, N>& v2)
{
    typename fixed_vector<T, N>::value_type s =
        typename fixed_vector<T, N>::value_type(0.0);

    typename fixed_vector<T, N>::const_iterator i1, i2;

    for (i1 = v1.begin(), i2 = v2.begin(); i1 < v1.end(); ++i1, ++i2)
        s += (*i1) * (*i2);

    return s;
}

// --- min/max -----------------------------------------------------

template<typename T, size_type N> inline
typename fixed_vector<T, N>::value_type min(const fixed_vector<T, N>& v)
{
    return *std::min_element(v.begin(), v.end());
}

template<typename T, size_type N> inline
fixed_vector<T, N> min(const fixed_vector<T, N>& v1,
                       const fixed_vector<T, N>& v2)
{
    fixed_vector<T, N> r;

    typename fixed_vector<T, N>::const_iterator i1, i2;
    typename fixed_vector<T, N>::iterator       i;

    for (i = r.begin(), i1 = v1.begin(), i2 = v2.begin(); i1 < v1.end(); ++i, ++i1, ++i2)
        *i = *i1 < *i2 ? *i1 : *i2;

    return r;
}

template<typename T, size_type N> inline
typename fixed_vector<T, N>::value_type max(const fixed_vector<T, N>& v)
{
    return *std::max_element(v.begin(), v.end());
}

template<typename T, size_type N> inline
fixed_vector<T, N> max(const fixed_vector<T, N>& v1,
                       const fixed_vector<T, N>& v2)
{
    fixed_vector<T, N> r;

    typename fixed_vector<T, N>::const_iterator i1, i2;
    typename fixed_vector<T, N>::iterator       i;

    for (i = r.begin(), i1 = v1.begin(), i2 = v2.begin(); i1 < v1.end(); ++i, ++i1, ++i2)
        *i = *i1 > *i2 ? *i1 : *i2;

    return r;
}

// --- cross product -----------------------------------------------

template<typename T> inline
T cross(const fixed_vector<T, 2>& v1, const fixed_vector<T, 2>& v2)
{
    return v1[0]*v2[1] - v2[0]*v1[1];
}

template<typename T> inline
fixed_vector<T, 3> cross(const fixed_vector<T, 3>& v1, const fixed_vector<T, 3>& v2)
{
    fixed_vector<T, 3> c;

    c[0] = v1[1] * v2[2] - v2[1] * v1[2];
    c[1] = v1[2] * v2[0] - v2[2] * v1[0];
    c[2] = v1[0] * v2[1] - v2[0] * v1[1];

    return c;
}

template<typename T, size_type N> inline
double norm(const fixed_vector<T, N>& v)
{
    return std::sqrt(inner(v, v));
}

template<typename T, size_type N> inline
double norm_inf(const fixed_vector<T, N>& v)
{
    size_type mi = 0;

    for (int i = 1; i < N; ++i)
        if (std::abs(v[i]) > std::abs(v[i-1]))
            mi = i;

    return std::abs(v[mi]);
}

// --- comparison --------------------------------------------------

template<typename T, size_type N> inline
fixed_vector<bool, N> operator<(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_type i = 0; i < N; ++i)
        r[i] = v1[i] < v2[i];

    return r;
}

template<typename T, size_type N> inline
fixed_vector<bool, N> operator<=(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_type i = 0; i < N; ++i)
        r[i] = v1[i] <= v2[i];

    return r;
}

template<typename T, size_type N> inline
fixed_vector<bool, N> operator>(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_type i = 0; i < N; ++i)
        r[i] = v1[i] > v2[i];

    return r;
}

template<typename T, size_type N> inline
fixed_vector<bool, N> operator>=(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_type i = 0; i < N; ++i)
        r[i] = v1[i] >= v2[i];

    return r;
}


template<typename T, size_type N> inline
fixed_vector<bool, N> operator==(const fixed_vector<T, N>& v1,
                                 const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_type i = 0; i < N; ++i)
        r[i] = v1[i] == v2[i];

    return r;
}

template<typename T, size_type N> inline
fixed_vector<bool, N> operator!=(const fixed_vector<T, N>& v1,
                                 const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_type i = 0; i < N; ++i)
        r[i] = v1[i] != v2[i];

    return r;
}

template<size_type N> inline
bool any(const fixed_vector<bool, N>& v)
{
    for (size_type i = 0; i < N; ++i)
        if (v[i])
            return true;

    return false;
}

template<size_type N> inline
bool all(const fixed_vector<bool, N>& v)
{
    for (size_type i = 0; i < N; ++i)
        if (!v[i])
            return false;

    return true;
}

// -----------------------------------------------------------------

template<typename T, size_type N> inline
fixed_vector<T, N> reverse(const fixed_vector<T, N>& v)
{
    fixed_vector<T, N> r;

    for (size_type i = 0; i < N; ++i)
        r[i] = v[N-i-1];

    return r;
}

template<typename T, size_type N> inline
fixed_vector<T, N> shift(const fixed_vector<T, N>& v)
{
    fixed_vector<T, N> r;

    r[N-1] = v[0];

    for (size_type i = 1; i < N; ++i)
        r[i-1] = v[i];

    return r;
}

// -----------------------------------------------------------------

template<size_type M1, size_type M2, typename T, size_type N> inline
fixed_vector<T, M2>& subv(fixed_vector<T, N>& v)
{
    static_assert(M1 + M2 <= N, "size mismatch in fixed_vector::subv");
    return *((fixed_vector<T, M2>*)(&(v[M1])));
}

template<size_type M1, size_type M2, typename T, size_type N> inline
const fixed_vector<T, M2>& subv(const fixed_vector<T, N>& v)
{
    static_assert(M1 + M2 <= N, "size mismatch in fixed_vector::subv");
    return *((const fixed_vector<T, M2>*)(&(v[M1])));
}

template<typename T, size_type N> inline
fixed_vector<T, N>& suba(T* array, size_t i)
{
    return *((fixed_vector<T, N>*)(&(array[N*i])));
}

template<typename T, size_type N> inline
const fixed_vector<T, N>& suba(const T* array, size_t i)
{
    return *((const fixed_vector<T, N>*)(&(array[N*i])));
}

template<typename T, size_type N, typename S> inline
fixed_vector < T, N + 1 > prepend(const S& s, const fixed_vector<T, N>& v)
{
    fixed_vector < T, N + 1 > r;
    r[0] = s;

    std::copy(v.begin(), v.end(), r.begin() + 1);
    return r;
}

// -----------------------------------------------------------------

template<typename T, size_type N>
std::ostream& operator<<(std::ostream& out, const fixed_vector<T, N>& v)
{
    out << "[ ";

    for (size_type i = 0; i < N - 1; ++i)
        out << v[i] << ", ";

    out << v[N-1] << " ]";

    return out;
}

template<typename T, size_type N>
std::istream& operator>>(std::istream& in, fixed_vector<T, N>& v)
{
    // several syntaxes are supported here:
    // - delimiters: "[...]" or "(...)" or none in case the vector is
    //   provided as a set of (possibly comma-separated) values
    // - separators: "," or " "
    // - with our without spaces in between

    char delim, c, separ;

    // determine which delimiter is used
    in >> c;
    if (c == '[') delim = ']';
    else if (c == '(') delim = ')';
    else {
        in.putback(c);
        delim = ' ';
    }

    // extract a sequence containing N values separated by commas or spaces
    // and ending with the proper delimiter character
    size_t i;
    for (i=0 ; !in.eof() && i<N ; ++i) {
        // formatted input
        in >> v[i]; // i-th value
        if (in.fail()) throw std::runtime_error
            ("Invalid format in input stream");
        // if this is not the last value, skip following comma or
        // space character(s)
        if (i < N-1) {
            // skip leading white spaces
            int skipped = 0;
            for (in >> c; !in.eof() && c==' ' ; in >> c, ++skipped) {}
            // if the last non-space character read is not a comma,
            // the separator used must be space and the last character
            // read is the first digit of the next entry: put
            // this character back into the input stream
            if (c != ',') {
                if (!skipped) throw std::runtime_error
                    ("Invalid format in input stream");
                else {
                    if (!i) separ = ' ';
                    if (separ == ' ') in.putback(c);
                    else throw std::runtime_error
                        ("Inconsistent separator used in input stream");
                }
            }
            else if (i && separ != ',') {
                throw std::runtime_error
                    ("Inconsistent separator used in input stream");
            }
            else if (!i) separ = ',';
        }
        else if (delim != ' '){
            // we have read all expected values. look for delimiter if
            // it is not a space character.
            // skip leading white spaces
            for (in >> c; !in.eof() && c==' '; in >> c) {}
            // check that the last non-space character is the expected
            // delimiter
            if (c != delim) throw std::runtime_error
                ("Invalid delimiter used in input stream");
        }
    }
    if (i < N) throw std::runtime_error("Invalid format in input stream");

    return in;
}

// -----------------------------------------------------------------

struct lexicographical_order {
    template<typename T, size_type N>
    bool operator()(const fixed_vector<T, N>& v1,
                    const fixed_vector<T, N>& v2) const {
        for (unsigned int i = 0; i < N; ++i) {
            if (v1[i] < v2[i])
                return true;
            if (v1[i] > v2[i])
                return false;
        }

        return false;
    }
};

// -----------------------------------------------------------------

struct eps_lexicographical_order {

    eps_lexicographical_order(double eps)
            : _eps(fabs(eps)) {}

    template<typename T, size_type N>
    bool operator()(const fixed_vector<T, N>& v1,
                    const fixed_vector<T, N>& v2) const {
        for (unsigned int i = 0; i < N; ++i) {
            if (v1[i] + _eps < v2[i])
                return true;
            if (v1[i] > v2[i] + _eps)
                return false;
        }

        return false;
    }

    double _eps;
};

// -----------------------------------------------------------------

template<typename T, size_type N>
fixed_vector<T, N> abs(const fixed_vector<T, N>& v)
{
    fixed_vector<T, N> r(0);
    for (unsigned int i = 0 ; i < N ; ++i) {
        r[i] = fabs(v[i]);
    }
    return r;
}

// -----------------------------------------------------------------

typedef fixed_vector<double, 1> vec1;
typedef fixed_vector<double, 2> vec2;
typedef fixed_vector<double, 3> vec3;
typedef fixed_vector<double, 4> vec4;
typedef fixed_vector<double, 5> vec5;

typedef fixed_vector<float, 1> fvec1;
typedef fixed_vector<float, 2> fvec2;
typedef fixed_vector<float, 3> fvec3;
typedef fixed_vector<float, 4> fvec4;
typedef fixed_vector<float, 5> fvec5;

typedef fixed_vector<int, 1> ivec1;
typedef fixed_vector<int, 2> ivec2;
typedef fixed_vector<int, 3> ivec3;
typedef fixed_vector<int, 4> ivec4;
typedef fixed_vector<int, 5> ivec5;

typedef fixed_vector<unsigned int, 1> uvec1;
typedef fixed_vector<unsigned int, 2> uvec2;
typedef fixed_vector<unsigned int, 3> uvec3;
typedef fixed_vector<unsigned int, 4> uvec4;
typedef fixed_vector<unsigned int, 5> uvec5;

typedef fixed_vector<size_type, 1> size1;
typedef fixed_vector<size_type, 2> size2;
typedef fixed_vector<size_type, 3> size3;
typedef fixed_vector<size_type, 4> size4;
typedef fixed_vector<size_type, 5> size5;

}  // namespace nvis

#endif // __fixed_vector_hpp
