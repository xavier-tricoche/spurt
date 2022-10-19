#ifndef __FUNCTORS_HPP__
#define __FUNCTORS_HPP__

#include <string>
#include <vector>
#include <sstream>
#include <math/fixed_vector.hpp>
#include "maps_lib/definitions.hpp"

namespace map_analysis {

// measure p-map vector within period tolerance interval
struct p_step_func {
    p_step_func(unsigned int period, double eps)
        : __period(period) {
        for (int i = 1 ; i <= period ; ++i) {
            double q = (double)period / (double)i;
            __qs.push_back(interval_type(q - eps, q + eps));
        }
    }
    
    int order() const {
        return 1;
    }
    
    std::string name() const {
        std::ostringstream os;
        os << __period << "-vector";
        return os.str();
    }
    
    bool is_valid(const point_data& d) const {
        bool is_inside = false;
        for (int i = 0 ; i < __qs.size() && !is_inside ; ++i) {
            is_inside = __qs[i].inside(d.period());
        }
        if (!is_inside) {
            return false;
        }
        
        nvis::vec2 v = vector_value(d, __period, static_data::metric);
        return (v[0] != std::numeric_limits<double>::max());
    }
    
    std::string value_string(const point_data& d) const {
        bool is_inside = false;
        for (int i = 0 ; i < __qs.size() && !is_inside ; ++i) {
            is_inside = __qs[i].inside(d.period());
        }
        if (!is_inside) {
            return "0 0 0";
        }
        
        nvis::vec2 v = vector_value(d, __period, static_data::metric);
        std::ostringstream os;
        os << v[0] << " " << v[1] << " 0";
        return os.str();
    }
    
    nvis::vec2 value(const point_data& d) const {
        bool is_inside = false;
        for (int i = 0 ; i < __qs.size() && !is_inside ; ++i) {
            is_inside = __qs[i].inside(d.period());
        }
        if (!is_inside) {
            return nvis::vec2(0, 0);
        } else {
            return vector_value(d, __period, static_data::metric);
        }
    }
    
    unsigned int __period;
    std::vector<interval_type> __qs;
};

inline nvis::vec2 sanitize(const nvis::vec2& v)
{
    nvis::vec2 w(v);
    for (int i = 0 ; i < 2 ; ++i) {
        if (std::isnan(w[i]) || std::isinf(w[i])) {
            w[i] = 0.;
        }
    }
    return w;
}

// measure p-error vector
struct error_func {
    error_func(double maxnorm = 10.) : __maxnorm(maxnorm) {}
    
    int order() const {
        return 1;
    }
    
    std::string name() const {
        std::ostringstream os;
        os << "error_function";
        return os.str();
    }
    
    bool is_valid(const point_data& d) const {
        const nvis::vec2& err = sanitize(d.error());
        return (err[0] != std::numeric_limits<double>::max());
    }
    
    std::string value_string(const point_data& d) const {
        nvis::vec2 e = sanitize(d.error());
        std::ostringstream os;
        if (e[0] == std::numeric_limits<double>::max() || nvis::norm(e) > __maxnorm) {
            os << "0 0 0";
        } else {
            os << e[0] << " " << e[1] << " 0";
        }
        return os.str();
    }
    
    nvis::vec2 value(const point_data& d) const {
        nvis::vec2 e = sanitize(d.error());
        if (e[0] == std::numeric_limits<double>::max() || nvis::norm(e) > __maxnorm) {
            return nvis::vec2(0, 0);
        } else {
            return e;
        }
    }
    
    double __maxnorm;
};

}


#endif





