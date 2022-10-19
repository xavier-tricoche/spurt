#ifndef __MAP_FIELD_WRAPPER_HPP__
#define __MAP_FIELD_WRAPPER_HPP__

#include <math/fixed_vector.hpp>
#include <data/raster.hpp>

namespace xavier {
namespace map {
typedef raster_grid<3, double>     grid_type;

template<typename F>
class wrapper {
    static double __mod(double a, double b) {
        return a >= 0 ? fmod(a, b) : b + fmod(a, b);
    }
    
public:
    typedef F   field_type;
    
    wrapper(const field_type& field) :
        __field(field), __periodic(false) {}
        
    wrapper(const wrapper& w) :
        __field(w.__field), __periodic(w.__periodic) {}
        
    void periodic(bool p) {
        __periodic = p;
    }
    
    nvis::vec3 unproject(const nvis::vec2& v) const {
        return nvis::vec3(v[1], v[0], 0.0);
    }
    
    nvis::vec2 project(const nvis::vec3& v) const {
        return (__periodic ?
                nvis::vec2(__mod(v[1], __field.get_grid().dimensions()[1] - 1),
                           __mod(v[0], __field.get_grid().dimensions()[0] - 1)) :
                nvis::vec2(v[1], v[0]));
    }
    
    nvis::mat2 project(const nvis::mat3& m) const {
        nvis::mat2 _m;
        // exchange role of x and y
        _m[0][0] = m[1][1];
        _m[0][1] = m[1][0];
        _m[1][0] = m[0][1];
        _m[1][1] = m[0][0];
        return _m;
    }
    
    nvis::vec2 modulo(const nvis::vec2& x) const {
        return project(__field.get_grid().dmodulo(unproject(x)));
    }
    
    double distance(const nvis::vec2& x, const nvis::vec2& y) const {
        return __field.get_grid().get_metric().distance(unproject(x), unproject(y));
    }
    
    bool intersect(const nvis::vec3& p0, const nvis::vec3& p1) const {
        int z0 = (int)trunc(p0[2] / (__field.get_grid().dimensions()[2] - 1));
        int z1 = (int)trunc(p1[2] / (__field.get_grid().dimensions()[2] - 1));
        
        return z0 != z1;
    }
    
    nvis::vec3 operator()(const double&, const nvis::vec3& p, bool verbose = false) const {
        try {
            if (verbose) {
                // __field.verbose(true);
            }
            nvis::vec3 f = __field.interpolate(p);
            if (verbose) {
                // __field.verbose(false);
            }
            return f;
        } catch (std::runtime_error& e) {
            // __field.verbose(false);
            throw;
        }
    }
    
    nvis::mat3 derivative(const double&, const nvis::vec3& p, bool verbose = false) const {
        try {
            if (verbose) {
                __field.verbose(true);
            }
            nvis::mat3 J = nvis::transpose(nvis::mat3(__field.derivative(p)));
            if (verbose) {
                __field.verbose(false);
            }
            return J;
        } catch (std::runtime_error& e) {
            __field.verbose(false);
            throw;
        }
    }
    
    const grid_type& grid() const {
        return __field.get_grid();
    }
    
private:
    const field_type&   __field;
    bool                __periodic;
};

} // map
} // xavier

#endif

