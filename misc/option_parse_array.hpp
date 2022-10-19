#ifndef __XAVIER_OPTION_PARSE_ARRAY_HPP__
#define __XAVIER_OPTION_PARSE_ARRAY_HPP__

#include <boost/lexical_cast.hpp>
#include <boost/config.hpp>
#include <boost/array.hpp>
#include <boost/program_options/errors.hpp>
#include <boost/shared_ptr.hpp>
#define __CPP11_ARRAY__
#include <array>

#include <math/fixed_vector.hpp>

namespace spurt { namespace command_line {

namespace po = boost::program_options;

// Add support for fixed-size array types in boost::program options
// Since validation occurs through a template function, explicit
// specialization is required. Macro meta-programming to the rescue...

template<typename arrayT, size_t N>
struct array_validator {
    static void validate(boost::any& v, 
                         const std::vector<std::string>& s) {                             
        typedef typename arrayT::value_type value_type;
           
        // check that number of arguments matches array size
        if (s.size() < N) 
            throw po::invalid_syntax(po::invalid_syntax::missing_parameter);
        else if (s.size() > N) 
            throw po::invalid_syntax(po::invalid_syntax::extra_parameter);
        
        if (v.empty()) {
            v = boost::any(arrayT());
        }
        arrayT* tv = boost::any_cast<arrayT>(&v);
        assert(NULL != tv);
        for (size_t i = 0; i < N; ++i)
        {
            try {
                // Use type-specific validator, if available, to validate
                // individual array entries
                boost::any a;
                std::vector<std::string> cv;
                cv.push_back(s[i]);
                po::validate(a, cv, (value_type*)0, 0);                
                (*tv)[i] = boost::any_cast<value_type>(a);
            }
            catch(const boost::bad_lexical_cast& /*e*/) {
                boost::throw_exception(po::invalid_option_value(s[i]));
            }
        }
    }
};

template<typename T>
struct vector_validator {
    static void validate(boost::any& v, 
                         const std::vector<std::string>& s) {                             
        typedef T val_t;
        typedef std::vector< T > vec_t;
           
        if (v.empty()) {
            v = boost::any(vec_t());
        }
        vec_t* tv = boost::any_cast<vec_t>(&v);
        assert(NULL != tv);
        for (size_t i = 0; i < s.size(); ++i)
        {
            try {
                // Use type-specific validator, if available, to validate
                // individual array entries
                boost::any a;
                std::vector<std::string> cv;
                cv.push_back(s[i]);
                po::validate(a, cv, (val_t*)0, 0);                
                (*tv)[i] = boost::any_cast<val_t>(a);
            }
            catch(const boost::bad_lexical_cast& /*e*/) {
                boost::throw_exception(po::invalid_option_value(s[i]));
            }
        }
    }
};

} // commandline
} // xavier

#define XAVIER_COMMAND_LINE_VECTOR_TYPES \
    X(bool); \
    X(char); \
    X(short); \
    X(int); \
    X(unsigned); \
    X(size_t); \
    X(float); \
    X(double); \
        
// Add fixed-size array validators to boost::program_options namespace
namespace boost { namespace program_options {
    
#define X(_Type) \
    template<> \
    void validate<>(boost::any& v, const std::vector<std::string>& s, \
                    std::vector<_Type>*, long) { \
         using namespace spurt::command_line; \
         std::cout << "Parsing bounded sequence" << std::endl; \
         vector_validator< std::vector<_Type> >::validate(v, s); \
    }
XAVIER_COMMAND_LINE_VECTOR_TYPES
#undef X
    
} }
            

#define XAVIER_COMMAND_LINE_ARRAY_TYPES_AND_SIZES \
    X(bool,     1); \
    X(bool,     2); \
    X(bool,     3); \
    X(bool,     4); \
    X(bool,     5); \
    X(bool,     6); \
    X(bool,     7); \
    X(bool,     8); \
    X(bool,     9); \
    X(int,      1); \
    X(int,      2); \
    X(int,      3); \
    X(int,      4); \
    X(int,      5); \
    X(int,      6); \
    X(int,      7); \
    X(int,      8); \
    X(int,      9); \
    X(unsigned, 1); \
    X(unsigned, 2); \
    X(unsigned, 3); \
    X(unsigned, 4); \
    X(unsigned, 5); \
    X(unsigned, 6); \
    X(unsigned, 7); \
    X(unsigned, 8); \
    X(unsigned, 9); \
    X(size_t,   1); \
    X(size_t,   2); \
    X(size_t,   3); \
    X(size_t,   4); \
    X(size_t,   5); \
    X(size_t,   6); \
    X(size_t,   7); \
    X(size_t,   8); \
    X(size_t,   9); \
    X(float,    1); \
    X(float,    2); \
    X(float,    3); \
    X(float,    4); \
    X(float,    5); \
    X(float,    6); \
    X(float,    7); \
    X(float,    8); \
    X(float,    9); \
    X(double,   1); \
    X(double,   2); \
    X(double,   3); \
    X(double,   4); \
    X(double,   5); \
    X(double,   6); \
    X(double,   7); \
    X(double,   8); \
    X(double,   9); \

// Add fixed-size array validators to boost::program_options namespace
namespace boost { namespace program_options {

#define X(_Type, _Size) \
    template<> \
    void validate<>(boost::any& v, const std::vector<std::string>& s, \
                    spurt::fixed_vector<_Type, _Size>*, long) { \
        using namespace spurt::command_line; \
        array_validator<spurt::fixed_vector<_Type, _Size>, _Size>::validate(v, s); \
    } \
    template<> \
    void validate<>(boost::any& v, const std::vector<std::string>& s, \
                    boost::array<_Type, _Size>*, long) { \
        using namespace spurt::command_line; \
        array_validator<boost::array<_Type, _Size>, _Size>::validate(v, s); \
    } 
XAVIER_COMMAND_LINE_ARRAY_TYPES_AND_SIZES
#undef X

#ifdef __CPP11_ARRAY__
#define X(_Type, _Size) \
    template<> \
    void validate<>(boost::any& v, const std::vector<std::string>& s, \
                    std::array<_Type, _Size>*, long) { \
        using namespace spurt::command_line; \
        array_validator<std::array<_Type, _Size>, _Size>::validate(v, s); \
    } 
XAVIER_COMMAND_LINE_ARRAY_TYPES_AND_SIZES
#undef X
#endif

} // program_options
} // boost


namespace spurt { namespace command_line {

#define STRINGIFY(x) #x
#define TO_STRING(x) STRINGIFY(x)
#define XAVIER_COMMAND_LINE_TYPES_TO_LONG_SYMBOLS \
    X(char, char); \
    X(unsigned char, uchar); \
    X(short, short); \
    X(unsigned short, ushort); \
    X(int, int); \
    X(unsigned int, uint); \
    X(long, long); \
    X(unsigned long, ulong); \
    X(long long, llong); \
    X(unsigned long long, ullong); \
    X(float, float); \
    X(double, double); \
    X(std::string, string); \
    X(char*, string);

template<typename T>
inline std::string type_as_long_string(bool use_brackets=true)
{
    return use_brackets ? "<arg>" : "arg";
}

#define X(type, typestr) \
    template<> \
    inline std::string type_as_long_string<type>(bool use_brackets) { \
        return use_brackets ? \
            std::string("<" TO_STRING(typestr) ">") : \
            std::string( TO_STRING(typestr) ); \
    }; \
    template<> \
    inline std::string type_as_long_string<std::vector<type> >(bool use_brackets) { \
        return use_brackets ? \
            std::string("[<" TO_STRING(typestr) ">]*") : \
            std::string("[" TO_STRING(typestr) "]*"); \
    }
XAVIER_COMMAND_LINE_TYPES_TO_LONG_SYMBOLS
#undef X

#define XAVIER_COMMAND_LINE_TYPES_TO_SHORT_SYMBOLS \
    X(char, c); \
    X(unsigned char, uc); \
    X(short, s); \
    X(unsigned short, us); \
    X(int, i); \
    X(unsigned int, ui); \
    X(long, l); \
    X(unsigned long, ul); \
    X(long long, ll); \
    X(unsigned long long, ull); \
    X(float, f); \
    X(double, d); \
    X(std::string, str); \
    X(char*, str);

template<typename T>
inline std::string type_as_short_string(bool use_brackets)
{
    return use_brackets ? "<?>" : "?";
}

#define X(type, typestr) \
    template<> \
    inline std::string type_as_short_string<type>(bool use_brackets) { \
        return use_brackets ? \
            std::string("<" TO_STRING(typestr) ">") : \
            std::string( TO_STRING(typestr) ); \
    }; \
    template<> \
    inline std::string type_as_short_string<std::vector<type> >(bool use_brackets) { \
        return use_brackets ? \
            std::string("[<" TO_STRING(typestr) ">]*") : \
            std::string("[" TO_STRING(typestr) "]*"); \
    }
XAVIER_COMMAND_LINE_TYPES_TO_SHORT_SYMBOLS
#undef X

template<typename T>
inline std::string type_as_string(bool use_brackets=true, bool use_short=false) {
    return use_short ? 
        type_as_short_string<T>(use_brackets) :
        type_as_long_string<T>(use_brackets);
}

} // commandline
} // xavier

#endif // __XAVIER_OPTION_PARSE_ARRAY_HPP__
