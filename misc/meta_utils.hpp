#ifndef __XAVIER_METAPROGRAMMING_UTILS_HPP__
#define __XAVIER_METAPROGRAMMING_UTILS_HPP__

#include <array>
#include <map>
#include <type_traits>
#include <complex>
#include <math/fixed_vector.hpp>
#include <Eigen/Core>

namespace spurt {

using namespace std;

// uniform treatment of data attribute as (possibly unary) array
template<typename T, typename Enable=void>
struct data_traits {};

template<typename T>
struct data_traits< T, typename std::enable_if<std::is_scalar<T>::value>::type > {
    typedef T value_type;
    typedef T data_type;

    constexpr static size_t size() { return 1; };
    constexpr static size_t mem_size() { return sizeof(value_type); }
    static const value_type& value(const data_type& v, size_t) { return v; }
    static value_type& value(data_type& v, size_t) { return v; }
    static data_type& assign(data_type& inout, value_type val) { inout=val; return inout; }
    static value_type norm(const data_type& v) { return std::labs(v); }
};

template<typename T, size_t N>
struct data_traits< nvis::fixed_vector<T, N> > {
    typedef T value_type;
    typedef nvis::fixed_vector<value_type, N> data_type;

    constexpr static size_t size() { return N; }
    constexpr static size_t mem_size() { return N*sizeof(value_type); }
    static const value_type& value(const data_type& v, size_t i) {
        return v[i];
    }
    static value_type& value(data_type& v, size_t i) { return v[i]; }
    static data_type& assign(data_type& inout, T val) {
        std::fill(inout.begin(), inout.end(), val);
        return inout;
    }
    static value_type norm(const data_type& v) { return nvis::norm(v); }
};

template<typename T, size_t N>
struct data_traits< std::array<T, N> > {
    typedef T value_type;
    typedef std::array<value_type, N> data_type;

    constexpr static size_t size() { return N; }
    constexpr static size_t mem_size() { return N*sizeof(value_type); }
    static const value_type& value(const data_type& v, size_t i) {
        return v[i];
    }
    static value_type& value(data_type& v, size_t i) { return v[i]; }
    static data_type& assign(data_type& inout, value_type val) {
        std::fill(inout.begin(), inout.end(), val);
        return inout;
    }
    static value_type norm(const data_type& v) {
        return std::inner_product(v.begin(), v.end(), v.begin(), 0);
    }
};

template<typename T, int NRows>
struct data_traits< Eigen::Matrix<T, NRows, 1> > {
    typedef T value_type;
    typedef Eigen::Matrix<T, NRows, 1> data_type;

    constexpr static int size() { return NRows; }
    constexpr static size_t mem_size() { return NRows*sizeof(value_type); }
    static const value_type& value(const data_type& v, size_t i) {
        return v[i];
    }
    static value_type& value(data_type& v, size_t i) { return v[i]; }
    static data_type& assign(data_type& inout, value_type val) {
        inout.setConstant(val);
        return inout;
    }
    static value_type norm(const data_type& v) {
        return v.norm();
    }
};

template<typename T>
struct data_wrapper {
    typedef data_traits<T> traits_type;
    typedef typename traits_type::value_type value_type;
    typedef T data_type;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;

    data_wrapper(const data_type& v) : m_v(v) {}
    int size() const { return traits_type::size(m_v); }
    constexpr static size_t mem_size() { return traits_type::mem_size(); }
    const value_type& operator[](size_t i) const { return traits_type::value(m_v, i); }
    value_type& operator[](size_t i) { return traits_type::value(m_v, i); }
    data_type& assign(value_type s) { traits_type::assign(m_v, s); return m_v; }
    data_type& set(value_type s) { return this->assign(s); }
    value_type norm() const { return traits_type::norm(m_v); }
    const_iterator begin() const { return &traits_type::value(m_v, 0); }
    const_iterator end() const { return &traits_type::value(m_v, 0) + this->size(); }

    data_type& m_v;
};

template<typename T>
struct is_array : public std::false_type {};

template<typename Value_, size_t N>
struct is_array< nvis::fixed_vector<Value_, N> > : public std::true_type {};

template<typename Value_, size_t N>
struct is_array< std::array<Value_, N> > : public std::true_type {};

template<typename Value, int NRows>
struct is_array< Eigen::Matrix<Value, NRows, 1> > : public std::true_type {};

// type to human-readable string conversion
template<typename T>
struct type2string {};

#define LIST_OF_STRING_MAPABLE_TYPES \
    X(bool); \
    X(unsigned char); \
    X(char); \
    X(unsigned short); \
    X(short); \
    X(unsigned int); \
    X(int); \
    X(long); \
    X(long long); \
    X(float); \
    X(double); \
    X(long double);

#define X(type) \
    template<> \
    struct type2string<type> { \
        typedef type data_type; \
        static const std::string type_name() { return #type ; } \
    }; \
    template<> \
    struct type2string< std::complex<type> >{ \
        typedef type data_type; \
        static const std::string type_name() { return "std::complex<" #type ">"; } \
    };
    LIST_OF_STRING_MAPABLE_TYPES
#undef X

// predicate to verify that object is a pair
template<typename T>
struct is_pair : public std::false_type {};

template<typename T1, typename T2>
struct is_pair< std::pair<T1,T2> > : public std::true_type {};

// predicate to verify that unsigned integral constant expression is
// strictly positive
template<size_t>
struct is_strictly_positive : public std::true_type {};

template<>
struct is_strictly_positive<0> : public std::false_type {};

// predicate to verify that multiple boolean constant expression are true
template<bool b1, bool b2,
         bool b3=true, bool b4=true, bool b5=true,
         bool b6=true, bool b7=true>
struct are_all_true : public std::false_type {};

template<>
struct are_all_true<true, true, true, true, true, true, true>
    : public std::true_type {};

// convert boolean predicate into type for SFINAE template signature
// verification
#define REQATTR(_attr,_type) \
    typename = typename enable_if<is_##_attr<_type>::value>::type
#define REQ2ATTR(_attr1,_type1,_attr2,_type2) \
    typename = typename enable_if< \
        are_all_true< \
            is_##_attr1<_type1>::value, \
            is_##_attr2<_type2>::value>::value>::type
#define REQ3ATTR(_attr1,_type1,_attr2,_type2,_attr3,_type3) \
    typename = typename enable_if< \
        are_all_true< \
            is_##_attr1<_type1>::value, \
            is_##_attr2<_type2>::value, \
            is_##_attr3<_type3>::value>::value>::type
#define REQSUBATTR(_attr,_type) \
    typename = typename enable_if<is_##_attr<typename _type>::value>::type
#define REQ2SUBATTR(_attr1,_type1,_attr2,_type2) \
    typename = typename enable_if< \
        are_all_true< \
            is_##_attr1<typename _type1>::value, \
            is_##_attr2<typename _type2>::value>::value>::type
#define REQ3SUBATTR(_attr1,_type1,_attr2,_type2,_attr3,_type3) \
    typename = typename enable_if< \
        are_all_true< \
            is_##_attr1<typename _type1>::value, \
            is_##_attr2<typename _type2>::value, \
            is_##_attr3<typename _type3>::value>::value>::type
#define REQPRED(_pred,_type) \
    typename = typename enable_if<_pred<_type>::value>::type
#define REQSUBPRED(_pred,_type) \
    typename = typename enable_if<_pred<typename _type>::value>::type
#define REQTYPE(_type1,_type2) \
    typename = typename enable_if<is_same<_type1, _type2>::value>::type
#define REQSUBTYPE(_type1,_type2) \
    typename = typename enable_if<is_same<typename _type1, _type2>::value>::type

} // spurt


#endif
