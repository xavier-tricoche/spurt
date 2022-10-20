// stdlib
#include <exception>
#include <sstream>
#include <string>
#include <vector>
// boost
#include <boost/program_options.hpp>
#include <boost/tuple.hpp>

namespace {
#define LIST_OF_AUTOMATICALLY_STRING_MAPPED_TYPES \
    X(char, "char"); \
    X(unsigned char, "uchar"); \
    X(short, "short"); \
    X(unsigned short, "ushort"); \
    X(int, "int"); \
    X(unsigned int, "uint"); \
    X(long, "long"); \
    X(unsigned long, "ulong"); \
    X(long long, "longlong"); \
    X(unsigned long long, "ulonglong"); \
    X(float, "float"); \
    X(double, "double"); \
    X(std::string, "string"); \
    X(char*, "string");

template<typename T>
std::string type_as_string() { return "X"; }

#define X(type, typestr) \
    template<> \
    std::string type_as_string<type> { \
        return std::string(typestr); \
    };
    LIST_OF_AUTOMATICALLY_STRING_MAPPED_TYPES
#undef X

template<typename T, unsigned N>
std::string array_as_string() { return "X ..."; }

#define X(type, typestr) \
    template<unsigned N> \
    std::string array_as_string<type, N> { \
        if (N == 1) return typestr; \
        else if (N == 2) return std::string(typestr) + " " + std::string(typestr); \
        else if (N == 3) return std::string(typestr) + " " + std::string(typestr) + " " + std::string(typestr); \
        else if (N > 3) { \
            std::ostringstream os; \
            os << typestr << " x" << N; \
            return os.str(); \
        } \
    };
    LIST_OF_AUTOMATICALLY_STRING_MAPPED_TYPES
#undef X
}

namespace spurt {

class cmdline_parser {
    namespace po = boost::program_options;
    
    bpo::options_description __desc;
    
    // function needed to turn multitoken items into
    // fixed or bounded length arrays since program_options
    // does not offer that possibility.
    template<typename Array, typename T, unsigned N>
    Array array_parser(const std::vector<T>& v,
                       const std::string& opt_name,
                       bool strict=true) {
        Array r;
        if (v.size() < N) {
            std::ostringstream os;
            os << "not enough values (" << v.size() << "<" << N
               << ") given for option " << opt_name;
            throw std::runtime_error(os.str());
        }
        for (size_t i=0 ; i<N ; ++i) {
            r[i] = v[i];
        }
        if (strict && v.size() > N) {
            std::ostringstream os;
            os << "too many values (" << v.size() << ">" << N
               << ") given for option " << opt_name;
            throw std::runtime_error(os.str());
        }
        return r;
    };
    
    std::map<std::string, std::pair<unsigned, unsigned> > array_bounds;
    std::map<std::string, boost::any> array_storage;
    
public:
    cmdline_parser(const std::string& caption) : __desc(caption) {}
    
    // add option corresponding to istream'able type
    template<typename T>
    void add(
        const std::string& opt_name,        // opt_name="canonical,short" enables --canonical and -short
        T& storage,                         // user supplied storage for the option value
        const std::string& description="",  // plain text description of the option
        bool required=true,                 // is the option required?
        const std::string& token_display="" // how should the value token be displayed?
    );
    
    // same as above with a fixed size random access container ("array")
    template<typename Array, typename T, unsigned N>
    void add(
        const std::string& opt_name,
        Array& storage,
        const std::string& description="", 
        unsigned min_token=N,               // least number of values needed to initialize the array
        bool required=true,
        const std::string& token_display=""
    );
    
    // same as above with C-style array
    template<typename T, unsigned N>
    void add(
        const std::string& opt_name,
        T* storage,
        const std::string& description="", 
        unsigned min_token=N,
        bool required=true,
        const std::string& token_display=""
    );
    
    // same as above with std::vector (there is special support for it in boost::program_options)
    template<typename T, unsigned N=-1>
    void add(
        const std::string& opt_name,
        std::vector<T>& storage,
        const std::string& description="", 
        unsigned min_token=N,
        bool required=true,
        const std::string& token_display=""
    );
};

template<typename T>
void spurt::cmdline_parser::add(const std::string& opt_name,
                                 T& storage, const std::string& description, 
                                 bool required, const std::string& token_display) {
                                     
    bo::typed_value<T>* semantic = bo::value<T>(&storage);
    if (required) semantic->required();
    if (token_display.empty()) semantic->value_name(type_as_string<T>());
    else semantic->value_name(token_display);
    _desc.add_options(opt_name, semantic, description);
}

template<typename Array, typename T, unsigned N>
void spurt::cmdline_parser::add(const std::string& opt_name,
                                 Array& storage, const std::string& description,
                                 unsigned min_token, bool required, 
                                 const std::string& token_display) {
                                     
    bo::typed_value<T>* semantic = bo::value<T>();
    semantic->multi_token();
    if (required) semantic->required();
    if (token_display.empty()) semantic->value_name(array_as_string<T, N>());
    else semantic->value_name(token_display);
    _desc.add_options(opt_name, semantic, description);
}






} // namespace spurt
