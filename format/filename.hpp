#ifndef __XAVIER_FILENAME_HPP__
#define __XAVIER_FILENAME_HPP__

#include <string>

namespace {
    std::string __before_char(const std::string& str, char c) {
        size_t n = str.find_last_of(c);
        if (n == std::string::npos) return str;
        else return str.substr(0, n);
    }
    std::string __after_char(const std::string& str, char c) {
        size_t n = str.find_last_of(c);
        if (n == std::string::npos) return "";
        else return str.substr(n+1);
    }
    void __around_char(std::string& before, std::string& after, const std::string& str, char c) {
        size_t n = str.find_last_of(c);
        if (n == std::string::npos) {
            before = str;
            after = "";
        }
        else {
            before = str.substr(0, n);
            after = str.substr(n+1);
        }
    }
}

// basic helper functions to manipulate file names
namespace spurt { namespace filename {
    
inline std::string parent_path(const std::string& name) {
    return __before_char(name, '/');
}

inline std::string extension(const std::string& name)
{
    return __after_char(name, '.');
}

inline std::string filename(const std::string& name) {
    return __after_char(name, '/');
}

inline std::string remove_extension(const std::string& name)
{
    return __before_char(name, '.');
}

inline std::string replace_extension(const std::string& name, const std::string& ext)
{
    return remove_extension(name) + '.' + ext;
}

} // filename 
} // spurt

#endif
