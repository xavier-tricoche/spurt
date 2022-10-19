#ifndef __FIELD_WRAPPER_HPP__
#define __FIELD_WRAPPER_HPP__

#include <stdexcept>
#include <iostream>
#include <vector>
#include <teem/nrrd.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <image/nrrd_wrapper.hpp>

namespace xavier {
template<typename _Format> class nrrd_data_traits {};

template<>
class nrrd_data_traits<Nrrd*> {
public:
    typedef double                              scalar_type;
    typedef nvis::fixed_vector<scalar_type, 3>  vector_type;
    typedef nvis::vec3                          point_type;
    
    nrrd_data_traits(const Nrrd* nin) {
        if (nin->dim < 3) {
            throw std::runtime_error("nrrd_data_traits<Nrrd*>: invalid dimension");
        }
        for (int i = 0 ; i < 3 ; ++i) {
            __size[i] = nin->axis[nin->dim-3+i].size;
        }
        __bounds = xavier::nrrd_utils::get_bounds<3>(nin);
//      std::cerr << "bounds are " << __bounds << "\n";
        __step = __bounds.size() / point_type(__size - nvis::ivec3(1, 1, 1));
        __scalar_values = xavier::nrrd_utils::to_array<scalar_type>(nin);
        __offset[0] = 0;
        __offset[1] = 1;
        __offset[2] = 1+__size[0];
        __offset[3] = __size[0];
        for (int i=0 ; i<4 ; ++i) {
            __offset[4+i] = __offset[i] + __size[1];
        }
    }
    
    ~nrrd_data_traits() {
        delete[] __scalar_values;
    }
    
    size_t size() const {
        size_t s = 0;
        s += 3*__size[0]*__size[1]*__size[2]*sizeof(scalar_type);
        s += sizeof(__bounds);
        s += sizeof(__offset);
        std::cerr << __size[0] << " x " << __size[1] << " x " << __size[2] << " x 3 x " << sizeof(scalar_type) << std::endl;
        return s;
    }
    
    bool get_value(const point_type& x, scalar_type& f) const {
        point_type y;
        if (!g2l(y, x)) {
            return false;
        }
        nvis::ivec3 id(floor(y[0]), floor(y[1]), floor(y[2]));
        point_type z = y - point_type(id);
        scalar_type u = z[0], v = z[1], w = z[2];
        scalar_type U = 1. - u, V = 1. - v, W = 1. - w;
        scalar_type weight[8] = {
            U* V*W, u* V*W, u* v*W, U* v*W,
            U* V*w, u* V*w, u* v*w, U* v* w
        };
        int n = index(id);
        f = 0;
        for (int i=0 ; i<8 ; ++i) {
            f += weight[i] * __scalar_values[n+__offset[i]];
        }
        
        return true;
    }
    
    bool get_value(const point_type& x, vector_type& f) const {
        point_type y;
        if (!g2l(y, x)) {
            return false;
        }
        nvis::ivec3 id(floor(y[0]), floor(y[1]), floor(y[2]));
        point_type z = y - point_type(id);
        scalar_type u = z[0], v = z[1], w = z[2];
        scalar_type U = 1. - u, V = 1. - v, W = 1. - w;
        scalar_type weight[8] = {
            U* V*W, u* V*W, u* v*W, U* v*W,
            U* V*w, u* V*w, u* v*w, U* v* w
        };
        int n = index(id);
        f = vector_type(0);
        for (int i=0 ; i<8 ; ++i) {
            f += weight[i] * nvis::suba<scalar_type, 3>(__scalar_values, n+__offset[i]);
        }
        return true;
    }
    
    const nvis::bbox3& bounds() const {
        return __bounds;
    }
    
private:

    bool g2l(point_type& l, const point_type& g) const {
        if (!__bounds.inside(g)) {
            return false;
        }
        l = (g - __bounds.min()) / __step;
        return true;
    }
    
    int index(const nvis::ivec3& id) const {
        return id[0] + __size[0]*(id[1] + __size[1]*id[2]);
    }
    
    nvis::ivec3                             __size;
    nvis::bbox3                             __bounds;
    point_type                              __step;
    scalar_type*                            __scalar_values;
    int                                     __offset[8];
};

template<typename _nrrd_data_traits>
struct rhs_wrapper {
    typedef _nrrd_data_traits                        nrrd_data_traits;
    typedef typename nrrd_data_traits::vector_type   vector_type;
    typedef typename nrrd_data_traits::point_type    point_type;
    
    rhs_wrapper(const nrrd_data_traits& traits) : __traits(traits) {}
    
    bool operator()(double t, const point_type& x, vector_type& f) const {
        return __traits.get_value(x, f);
    }
    
    const nrrd_data_traits& __traits;
};

}



#endif








