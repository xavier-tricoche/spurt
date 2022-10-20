#ifndef __METRIC_HPP__
#define __METRIC_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace spurt {
template<typename T, int N>
class metric {
    static T __modulo(T a, T b) {
        T r = fmod(a, b);
        if (!r) {
            return 0;
        } else {
            return a >= 0 ? r : b + r;
        }
    }
    
    static T __min_norm(T a, T b) {
        if (2*a > b) {
            return a - b;
        } else if (2*a + b < 0) {
            return a + b;
        } else {
            return a;
        }
    }
    
public:
    typedef T                                   scalar_type;
    typedef nvis::fixed_vector<scalar_type, N>  vec_type;
    typedef nvis::fixed_vector<bool, N>         bvec_type;
    typedef nvis::bounding_box<vec_type>        bounds_type;
    
    metric() :
        __bounds(vec_type(0), vec_type(0)), __periodic(false) {}
    metric(const bounds_type& bounds) :
        __bounds(bounds), __periodic(false) {}
    metric(const bounds_type& bounds, const bvec_type& periodic) :
        __bounds(bounds), __periodic(periodic) {}
        
    const bounds_type& bounds() const {
        return __bounds;
    }
    
    bounds_type& bounds() {
        return __bounds;
    }
    
    const bvec_type& periodic() const {
        return __periodic;
    }
    
    bvec_type& periodic() {
        return __periodic;
    }
    
    vec_type size() const {
        return __bounds.size();
    }
    
    vec_type modulo(const vec_type& x) const {
        vec_type y = x - __bounds.min();
        const vec_type& sz = size();
        for (int i = 0 ; i < N ; ++i) {
            if (__periodic[i]) {
                y[i] = __modulo(y[i], sz[i]);
            }
        }
        return y + __bounds.min();
    }
    
    vec_type displacement(const vec_type& a, const vec_type& b) const {
        vec_type d = modulo(b) - modulo(a);
        const vec_type& sz = size();
        for (int i = 0 ; i < N ; ++i) {
            if (__periodic[i]) {
                d[i] = __min_norm(d[i], sz[i]);
            }
        }
        return d;
    }
    
    scalar_type distance(const vec_type& a, const vec_type& b) const {
        return nvis::norm(displacement(a, b));
    }
    
private:
    bounds_type __bounds;
    bvec_type   __periodic;
};
} // namespace spurt

#endif




