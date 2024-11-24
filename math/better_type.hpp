#pragma once

#include <type_traits>

namespace spurt {

// better_type:
// deciding which type should be
// used to perform a binary operation on a pair of types
template<typename T1, typename T2, typename Enable=void>
struct better_type {};

// always pick floating point over integral
template <typename T1, typename T2>
struct better_type < T1, T2,
    typename std::enable_if<std::is_floating_point<T1>::value && 
                            std::is_integral<T2>::value >::type >
{
    typedef T1 type;
};
template <typename T1, typename T2>
struct better_type<T1, T2,
    typename std::enable_if<std::is_floating_point<T2>::value &&
                            std::is_integral<T1>::value >::type >
{
    typedef T2 type;
};

// all integral types handled below:
// bool, 
// unsigned char, char, 
// unsigned short, short, 
// unsigned int, int, 
// unsigned long, long, 
// unsigned long long, long long

#define INTEGRAL_TYPES_MATCHES \
    X(bool,               bool,               bool); \
    X(bool,               unsigned char,      char); \
    X(bool,               char,               char); \
    X(bool,               unsigned short,     short); \
    X(bool,               short,              short); \
    X(bool,               unsigned int,       int); \
    X(bool,               int,                int); \
    X(bool,               unsigned long,      long); \
    X(bool,               long,               long); \
    X(bool,               unsigned long long, long long); \
    X(bool,               long long,          long long); \
    X(unsigned char,      bool,               char); \
    X(unsigned char,      unsigned char,      char); \
    X(unsigned char,      char,               char); \
    X(unsigned char,      unsigned short,     short); \
    X(unsigned char,      short,              short); \
    X(unsigned char,      unsigned int,       int); \
    X(unsigned char,      int,                int); \
    X(unsigned char,      unsigned long,      long); \
    X(unsigned char,      long,               long); \
    X(unsigned char,      unsigned long long, long long); \
    X(unsigned char,      long long,          long long); \
    X(char,               bool,               char); \
    X(char,               unsigned char,      char); \
    X(char,               char,               char); \
    X(char,               unsigned short,     short); \
    X(char,               short,              short); \
    X(char,               unsigned int,       int); \
    X(char,               int,                int); \
    X(char,               unsigned long,      long); \
    X(char,               long,               long); \
    X(char,               unsigned long long, long long); \
    X(char,               long long,          long long); \
    X(unsigned short,     bool,               short); \
    X(unsigned short,     unsigned char,      char); \
    X(unsigned short,     char,               char); \
    X(unsigned short,     unsigned short,     short); \
    X(unsigned short,     short,              short); \
    X(unsigned short,     unsigned int,       int); \
    X(unsigned short,     int,                int); \
    X(unsigned short,     unsigned long,      long); \
    X(unsigned short,     long,               long); \
    X(unsigned short,     unsigned long long, long long); \
    X(unsigned short,     long long,          long long); \
    X(short,              bool,               short); \
    X(short,              unsigned char,      short); \
    X(short,              char,               short); \
    X(short,              unsigned short,     short); \
    X(short,              short,              short); \
    X(short,              unsigned int,       int); \
    X(short,              int,                int); \
    X(short,              unsigned long,      long); \
    X(short,              long,               long); \
    X(short,              unsigned long long, long long); \
    X(short,              long long,          long long); \
    X(unsigned int,       bool,               int); \
    X(unsigned int,       unsigned char,      int); \
    X(unsigned int,       char,               int); \
    X(unsigned int,       unsigned short,     int); \
    X(unsigned int,       short,              unsigned int); \
    X(unsigned int,       unsigned int,       int); \
    X(unsigned int,       int,                int); \
    X(unsigned int,       unsigned long,      long); \
    X(unsigned int,       long,               long); \
    X(unsigned int,       unsigned long long, long long); \
    X(unsigned int,       long long,          long long); \
    X(int,                bool,               int); \
    X(int,                unsigned char,      int); \
    X(int,                char,               int); \
    X(int,                unsigned short,     int); \
    X(int,                short,              int); \
    X(int,                unsigned int,       int); \
    X(int,                int,                int); \
    X(int,                unsigned long,      long); \
    X(int,                long,               long); \
    X(int,                unsigned long long, long long); \
    X(int,                long long,          long long); \
    X(unsigned long,      bool,               long); \
    X(unsigned long,      unsigned char,      long); \
    X(unsigned long,      char,               long); \
    X(unsigned long,      unsigned short,     long); \
    X(unsigned long,      short,              long); \
    X(unsigned long,      unsigned int,       long); \
    X(unsigned long,      int,                long); \
    X(unsigned long,      unsigned long,      long); \
    X(unsigned long,      long,               long); \
    X(unsigned long,      unsigned long long, long long); \
    X(unsigned long,      long long,          long long); \
    X(long,               bool,               long); \
    X(long,               unsigned char,      long); \
    X(long,               char,               long); \
    X(long,               unsigned short,     long); \
    X(long,               short,              long); \
    X(long,               unsigned int,       long); \
    X(long,               int,                long); \
    X(long,               unsigned long,      long); \
    X(long,               long,               long); \
    X(long,               unsigned long long, long long); \
    X(long,               long long,          long long); \
    X(unsigned long long, bool,               long long); \
    X(unsigned long long, unsigned char,      long long); \
    X(unsigned long long, char,               long long); \
    X(unsigned long long, unsigned short,     long long); \
    X(unsigned long long, short,              long long); \
    X(unsigned long long, unsigned int,       long long); \
    X(unsigned long long, int,                long long); \
    X(unsigned long long, unsigned long,      long long); \
    X(unsigned long long, long,               long long); \
    X(unsigned long long, unsigned long long, long long); \
    X(unsigned long long, long long,          long long); \
    X(long long,          bool,               long long); \
    X(long long,          unsigned char,      long long); \
    X(long long,          char,               long long); \
    X(long long,          unsigned short,     long long); \
    X(long long,          short,              long long); \
    X(long long,          unsigned int,       long long); \
    X(long long,          int,                long long); \
    X(long long,          unsigned long,      long long); \
    X(long long,          long,               long long); \
    X(long long,          unsigned long long, long long); \
    X(long long,          long long,          long long); \
    X(float,              double,             double);   \
    X(float,              float,              float);    \
    X(double,             float,              double);   \
    X(double,             double,             double);

#define X(type1,type2,type3) \
    template <> \
    struct better_type<type1,type2,void> \
    { \
        typedef type3 type; \
    };
    INTEGRAL_TYPES_MATCHES
#undef X

} // namespace spurt