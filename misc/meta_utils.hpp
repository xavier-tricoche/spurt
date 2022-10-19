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
    static const T& value(const T& v, size_t) { return v; }
    static T& value(T& v, size_t) { return v; }
    static T& assign(T& inout, T val) { inout=val; return inout; }
};

template<typename T, size_t N>
struct data_traits< fixed_vector<T, N> > {
    typedef T value_type;
    typedef fixed_vector<value_type, N> data_type;
    constexpr static size_t size() { return N; }
    static const T& value(const fixed_vector<T, N>& v, size_t i) {
        return v[i];
    }
    static T& value(fixed_vector<T, N>& v, size_t i) { return v[i]; }
    static data_type& assign(data_type& inout, T val) {
        std::fill(inout.begin(), inout.end(), val);
        return inout;
    }
};

template<typename T, size_t N>
struct data_traits< std::array<T, N> > {
    typedef T value_type;
    typedef std::array<value_type, N> data_type;
    constexpr static size_t size() { return N; }
    static const T& value(const std::array<T, N>& v, size_t i) {
        return v[i];
    }
    static T& value(std::array<T, N>& v, size_t i) { return v[i]; }
    static data_type& assign(data_type& inout, T val) {
        std::fill(inout.begin(), inout.end(), val);
        return inout;
    }
};

template<typename T, int NRows>
struct data_traits< Eigen::Matrix<T, NRows, 1> > {
    typedef T value_type;
    typedef Eigen::Matrix<T, NRows, 1> data_type;
    constexpr static int size() { return NRows; }
    static const T& value(const data_type& v, size_t i) {
        return v[i];
    }
    static T& value(data_type& v, size_t i) { return v[i]; }
    static data_type& assign(data_type& inout, T val) {
        inout.setConstant(val);
        return inout;
    }
};

template<typename T>
struct is_array : public std::false_type {};

template<typename Value_, size_t N>
struct is_array< fixed_vector<Value_, N> > : public std::true_type {};

template<typename Value_, size_t N>
struct is_array< std::array<Value_, N> > : public std::true_type {};

template<typename Value, int NRows>
struct is_array< Eigen::Matrix<Value, NRows, 1> > : public std::true_type {};

// need to figure out how to use std::decltype to test existence of a given
// method in class definition (here size())
// template<T, typename = typename std::enable_if<std::decltype(T().size())>::type
// struct is_array<T> : public std::true_type {};

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
        static const std::string type_name; \
    }; \
    const std::string type2string<type>::type_name = #type;
    LIST_OF_STRING_MAPABLE_TYPES
#undef X

template<>
struct type2string< std::complex<float> > {
    typedef float data_type;
    static const std::string type_name;
};
template<>
struct type2string< std::complex<double> > {
    typedef double data_type;
    static const std::string type_name;
};
template<>
struct type2string< std::complex<long double> > {
    typedef long double data_type;
    static const std::string type_name;
};
const std::string type2string< std::complex<float> >::type_name = "std::complex<float>";
const std::string type2string< std::complex<double> >::type_name = "std::complex<double>";
const std::string type2string< std::complex<long double> >::type_name = "std::complex<long double>";


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

} // xavier


#endif
