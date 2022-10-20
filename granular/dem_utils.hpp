#ifndef __XAVIER_GRANULAR_DEM_UTILS_HPP__
#define __XAVIER_GRANULAR_DEM_UTILS_HPP__

#include <math.h>
#include <algorithm>

namespace spurt {
namespace granular {
namespace dem {
namespace utils {
    
    template<typename T, typename Less>
    struct less_than_pair {
        typedef std::pair<T, T> pair_t;
        bool operator()(const pair_t& p0, const pair_t& p1) const {
            if (_less(p0.first, p1.first)) return true;
            else if (_less(p1.first, p0.first)) return false;
            return _less(p0.second, p1.second);
        }
        Less _less;
    };
    
    template<typename vector_t, 
             typename value_t = typename vector_t::value_type>
    inline value_t inner_product(const vector_t& v0, const vector_t& v1) {
        return v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
    }
    
    template<typename vector_t, 
             typename value_t = typename vector_t::value_type>
    inline value_t norm_sq(const vector_t& v) {
        return inner_product<vector_t, value_t>(v, v);
    }
    template<typename vector_t, 
             typename value_t = typename vector_t::value_type>
    inline value_t norm(const vector_t& v) {
        return sqrt(norm_sq<vector_t, value_t>(v));
    }
    
    template<typename vector_t, 
             typename value_t = typename vector_t::value_type>
    inline value_t distance_linf(const vector_t& v0, const vector_t& v1) {
        value_t d = fabs(v0[0]-v1[0]);
        d = std::max(d, fabs(v0[1]-v1[1]));
        d = std::max(d, fabs(v0[2]-v1[2]));
        return d;
    }
    
    template<typename vector_t, 
             typename value_t = typename vector_t::value_type>
    inline value_t distance_sq(const vector_t& v0, const vector_t& v1) {
        return norm_sq<vector_t, value_t>(v0-v1);
    }
    
    template<typename vector_t, 
             typename value_t = typename vector_t::value_type>
    inline value_t distance(const vector_t& v0, const vector_t& v1) {
        return sqrt(distance_sq<vector_t, value_t>(v0, v1));
    }
    
    template<typename vector_t>
    inline vector_t outer_product(const vector_t& u, const vector_t& v) {
        vector_t w;
        w[0] = u[1]*v[2] - u[2]*v[1];
        w[1] = u[2]*v[0] - u[0]*v[2];
        w[2] = u[0]*v[1] - u[1]*v[0];
        return w;
    }
    
    template<typename vector_t>
    inline vector_t project_on_vector(const vector_t& u, 
                                      const vector_t& n) {
        return inner_product<vector_t>(u, n)*n;
    }
    
    template<typename vector_t>
    inline vector_t project_on_plane(const vector_t& u,
                                     const vector_t& normal) {
        vector_t v(u);
        v -= project_on_vector<vector_t>(u, normal);
        return v;
    }
    
} // utils
} // dem
} // granular
} // spurt