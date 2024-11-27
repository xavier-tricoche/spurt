#ifndef m_METRIC_HPPm_
#define m_METRIC_HPPm_

#include <math/types.hpp>
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
    typedef T                                    scalar_type;
    typedef spurt::small_vector<scalar_type, N>  vec_type;
    typedef spurt::small_vector<bool, N>         bvec_type;
    typedef spurt::bounding_box<vec_type>        bounds_type;
    
    metric() :
        m_bounds(vec_type(0), vec_type(0)), m_periodic(false) {}
    metric(const bounds_type& bounds) :
        m_bounds(bounds), m_periodic(false) {}
    metric(const bounds_type& bounds, const bvec_type& periodic) :
        m_bounds(bounds), m_periodic(periodic) {}
        
    const bounds_type& bounds() const {
        return m_bounds;
    }
    
    bounds_type& bounds() {
        return m_bounds;
    }
    
    const bvec_type& periodic() const {
        return m_periodic;
    }
    
    bvec_type& periodic() {
        return m_periodic;
    }
    
    vec_type size() const {
        return m_bounds.size();
    }
    
    vec_type modulo(const vec_type& x) const {
        vec_type y = x - m_bounds.min();
        const vec_type& sz = size();
        for (int i = 0 ; i < N ; ++i) {
            if (m_periodic[i]) {
                y[i] = __modulo(y[i], sz[i]);
            }
        }
        return y + m_bounds.min();
    }
    
    vec_type displacement(const vec_type& a, const vec_type& b) const {
        vec_type d = modulo(b) - modulo(a);
        const vec_type& sz = size();
        for (int i = 0 ; i < N ; ++i) {
            if (m_periodic[i]) {
                d[i] = __min_norm(d[i], sz[i]);
            }
        }
        return d;
    }
    
    scalar_type distance(const vec_type& a, const vec_type& b) const {
        return spurt::norm(displacement(a, b));
    }
    
private:
    bounds_type m_bounds;
    bvec_type   m_periodic;
};
} // namespace spurt

#endif




