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

namespace {
    template<typename T>
    using must_be_number = typename std::enable_if< std::is_arithmetic< T >::value >::type;
}

namespace spurt
{

// Convenience wrapper class for std::array that endowes it with some limited
// linear algebra semantic.
// This code is based on an implementation generously provided by Prof. Dr. Christoph Garth, 
// Univ. of Kaiserslautern, Germany.
template<typename T, size_t N>
class fixed_vector : public std::array<T, N>
{
public:
    static constexpr size_t dimension = N;
    
    static constexpr size_t size() { return dimension; }
    
    // static constexpr size_t size = N;
    
    typedef std::array<T, N>   base_type;
    
    typedef typename base_type::value_type      value_type;
    typedef typename base_type::pointer         pointer;
    typedef typename base_type::reference       reference;
    typedef typename base_type::reference       const_reference;
    typedef typename base_type::iterator        iterator;
    typedef typename base_type::const_iterator  const_iterator;
    typedef fixed_vector<T, N>                  self_type;

    // --- constructors ---

    fixed_vector() : base_type() {
    }

    fixed_vector(const T& v) : base_type() {
        std::fill(this->begin(), this->end(), v);
    }
    
    fixed_vector(const base_type& other) : base_type(other) {
    }
    
    fixed_vector(const self_type& other) : base_type(static_cast<base_type>(other)) {}
    
    template<typename Q, must_be_number<Q> >
    fixed_vector(const fixed_vector<Q, N>& other) {
        for (size_t i=0; i<dimension; ++i) 
            base_type::operator[](i) = static_cast<T>(other[i]);
    }
    
    template<typename Q1, size_t N1, typename Q2, typename = typename std::enable_if< N1<N >::type>
    fixed_vector(const fixed_vector<Q1, N1>& first, const fixed_vector<Q2, N-N1>& second) {
        for (size_t i=0; i<N1; ++i) 
            base_type::operator[](i) = static_cast<T>(first[i]);
        for (size_t i=N1; i<N; ++i) 
            base_type::operator[](i) = static_cast<T>(second[i-N1]);
    }
    
    template<typename Q1, typename Q2>
    fixed_vector(const fixed_vector<Q1, N-1>& first, Q2 last) {
        for (size_t i=0; i<N-1; ++i) 
            base_type::operator[](i) = static_cast<T>(first[i]);
        base_type::operator[](N-1) = static_cast<T>(last);
    }
    
    template<typename Q>
    fixed_vector(const Q* _carray) {
        for (size_t i=0; i<N-1; ++i)
            base_type::operator[](i) = static_cast<T>(_carray[i]);
    }
    
    template<typename Q, size_t M, typename = typename std::enable_if< M >= N >::type>
    fixed_vector(const std::array<Q, M>& _array) {
        for (size_t i=0; i<N-1; ++i)
            base_type::operator[](i) = static_cast<T>(_array[i]);
    } 
    
    // --- multi-argument initialization constructors

    fixed_vector(const T& v0, const T& v1) : base_type({v0, v1}) {}
    fixed_vector(const T& v0, const T& v1, const T& v2) : base_type({v0, v1, v2}) {}
    fixed_vector(const T& v0, const T& v1, const T& v2, const T& v3) : base_type({v0, v1, v2, v3}) {}
    
    // --- coefficient-wise vector-vector operators ---

    self_type& operator+=(const self_type& rhs) {
        for (auto i=0; i<base_type::size(); ++i) this->operator[](i) += rhs[i];
        return *this;
    }
    
    self_type operator+(const self_type& rhs) const {
        self_type r(*this);
        r+=rhs;
        return r;
    }

    self_type& operator-=(const self_type& rhs) {
        for (auto i=0; i<base_type::size(); ++i) this->operator[](i) -= rhs[i];
        return *this;
    }
    
    self_type operator-(const self_type& rhs) const {
        self_type r(*this);
        r-=rhs;
        return r;
    }

    self_type& operator*=(const self_type& rhs) {
        for (auto i=0; i<base_type::size(); ++i) this->operator[](i) *= rhs[i];
        return *this;
    }

    self_type& operator*=(const value_type a) {
        for (auto i=0; i<base_type::size(); ++i) this->operator[](i) *= a;
        return *this;
    }
    
    self_type operator*(const self_type& rhs) const {
        self_type r(*this);
        r*=rhs;
        return r;
    }
    
    self_type operator*(const value_type a) const {
        self_type r(*this);
        r*=a;
        return r;
    }

    self_type& operator/=(const self_type& rhs) {
        for (auto i=0; i<base_type::size(); ++i) this->operator[](i) /= rhs[i];
        return *this;
    }

    self_type& operator/=(const value_type a) {
        for (auto i=0; i<base_type::size(); ++i) this->operator[](i) /= a;
        return *this;
    }
    
    self_type operator/(const value_type a) const {
        self_type r(*this);
        r/=a;
        return r;
    }
};

// --- right hand side scalar vector product

template<typename T, size_t N> inline
fixed_vector<T, N> operator*(const typename fixed_vector<T, N>::value_type& a, const fixed_vector<T, N>& v) {
    fixed_vector<T, N> r(v);
    r*=a;
    return r;
}

template<typename T, size_t N> inline
fixed_vector<T, N> operator/(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2) {
    fixed_vector<T, N> r(v1);
    r/=v2;
    return r;
}

// --- accumulators ------------------------------------------------

template<typename T, size_t N> inline
typename fixed_vector<T, N>::value_type sum(const fixed_vector<T, N>& v)
{
    auto it=v.begin();
    T r = *it++;
    for (; it!=v.end(); r += *it++) {}
    std::for_each(it, v.end(), [&](const T& t) { r += t; });
    return r;
}

template<typename T, size_t N> inline
typename fixed_vector<T, N>::value_type prod(const fixed_vector<T, N>& v)
{
    auto it=v.begin();
    T r = *it++;
    for (; it!=v.end(); r *= *it++) {}
    return r;
}

template<typename T, size_t N> inline
typename fixed_vector<T, N>::value_type inner(const fixed_vector<T, N>& v1,
        const fixed_vector<T, N>& v2)
{
    auto it1 = v1.begin();
    auto it2 = v2.begin();
    T r=*it1++ * *it2++;
    for (; it1 != v1.end(); r += *it1++ * *it2++) {}
    
    return r;
}

// --- min/max -----------------------------------------------------

template<typename T, size_t N> inline
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

template<typename T, size_t N> inline
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
    return fixed_vector<T, 3>(v1[1]*v2[2]-v1[2]*v2[1],
                              v1[2]*v2[0]-v1[0]*v2[2],
                              v1[0]*v2[1]-v1[1]*v2[0]);
}

template<typename T, size_t N> inline
double norm(const fixed_vector<T, N>& v)
{
    auto it=v.begin();
    double r = *it * *it;
    for (it++; it!=v.end(); it++) { r += *it * *it; }
    return sqrt(r);
}

template<typename T, size_t N> inline
double norm_inf(const fixed_vector<T, N>& v)
{
    size_t mi = 0;

    for (int i = 1; i < N; ++i)
        if (std::abs(v[i]) > std::abs(v[mi]))
            mi = i;

    return std::abs(v[mi]);
}

// --- comparison --------------------------------------------------


template<typename T, size_t N> inline
fixed_vector<bool, N> operator<(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_t i = 0; i < N; ++i)
        r[i] = v1[i] < v2[i];

    return r;
}

template<typename T, size_t N> inline
fixed_vector<bool, N> operator<=(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_t i = 0; i < N; ++i)
        r[i] = v1[i] <= v2[i];

    return r;
}

template<typename T, size_t N> inline
fixed_vector<bool, N> operator>(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_t i = 0; i < N; ++i)
        r[i] = v1[i] > v2[i];

    return r;
}

template<typename T, size_t N> inline
fixed_vector<bool, N> operator>=(const fixed_vector<T, N>& v1, const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_t i = 0; i < N; ++i)
        r[i] = v1[i] >= v2[i];

    return r;
}


template<typename T, size_t N> inline
fixed_vector<bool, N> operator==(const fixed_vector<T, N>& v1,
                                 const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_t i = 0; i < N; ++i)
        r[i] = v1[i] == v2[i];

    return r;
}

template<typename T, size_t N> inline
fixed_vector<bool, N> operator!=(const fixed_vector<T, N>& v1,
                                 const fixed_vector<T, N>& v2)
{
    fixed_vector<bool, N> r;

    for (size_t i = 0; i < N; ++i)
        r[i] = v1[i] != v2[i];

    return r;
}


template<size_t N> inline
bool any(const fixed_vector<bool, N>& v)
{
    for (size_t i = 0; i < N; ++i)
        if (v[i])
            return true;

    return false;
}

template<size_t N> inline
bool all(const fixed_vector<bool, N>& v)
{
    for (size_t i = 0; i < N; ++i)
        if (!v[i])
            return false;

    return true;
}

// -----------------------------------------------------------------


template<typename T, size_t N> inline
fixed_vector<T, N> reverse(const fixed_vector<T, N>& v)
{
    fixed_vector<T, N> r;

    for (size_t i = 0; i < N; ++i)
        r[i] = v[N-i-1];

    return r;
}

template<typename T, size_t N> inline
fixed_vector<T, N> shift(const fixed_vector<T, N>& v)
{
    fixed_vector<T, N> r;

    r[N-1] = v[0];

    for (size_t i = 1; i < N; ++i)
        r[i-1] = v[i];

    return r;
}

// -----------------------------------------------------------------

template<size_t M1, size_t M2, typename T, size_t N> inline
fixed_vector<T, M2>& subv(fixed_vector<T, N>& v)
{
    static_assert(M1 + M2 <= N, "dimensions mismatch");
    return *((fixed_vector<T, M2>*)(&(v[M1])));
}

template<size_t M1, size_t M2, typename T, size_t N> inline
const fixed_vector<T, M2>& subv(const fixed_vector<T, N>& v)
{
    static_assert(M1 + M2 <= N, "dimensions mismatch");
    return *((const fixed_vector<T, M2>*)(&(v[M1])));
}

template<typename T, size_t N> inline
fixed_vector<T, N>& suba(T* array, size_t i)
{
    return *((fixed_vector<T, N>*)(&(array[N*i])));
}

template<typename T, size_t N> inline
const fixed_vector<T, N>& suba(const T* array, size_t i)
{
    return *((const fixed_vector<T, N>*)(&(array[N*i])));
}

template<typename T, size_t N, typename S> inline
fixed_vector < T, N + 1 > prepend(const S& s, const fixed_vector<T, N>& v)
{
    fixed_vector < T, N + 1 > r;
    r[0] = s;

    std::copy(v.begin(), v.end(), r.begin() + 1);
    return r;
}

// -----------------------------------------------------------------

template<typename T, size_t N>
std::ostream& operator<<(std::ostream& out, const fixed_vector<T, N>& v)
{
    out << "[ ";

    for (size_t i = 0; i < N - 1; ++i)
        out << v[i] << ", ";

    out << v[N-1] << " ]";

    return out;
}

template<typename T, size_t N>
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
    template<typename T, size_t N>
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

    template<typename T, size_t N>
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

template<typename T, size_t N>
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

typedef fixed_vector<char, 1> cvec1;
typedef fixed_vector<char, 2> cvec2;
typedef fixed_vector<char, 3> cvec3;
typedef fixed_vector<char, 4> cvec4;
typedef fixed_vector<char, 5> cvec5;

typedef fixed_vector<short, 1> svec1;
typedef fixed_vector<short, 2> svec2;
typedef fixed_vector<short, 3> svec3;
typedef fixed_vector<short, 4> svec4;
typedef fixed_vector<short, 5> svec5;

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

typedef fixed_vector<size_t, 1> size1;
typedef fixed_vector<size_t, 2> size2;
typedef fixed_vector<size_t, 3> size3;
typedef fixed_vector<size_t, 4> size4;
typedef fixed_vector<size_t, 5> size5;

}  // namespace spurt

#endif // __fixed_vector_hpp



