#ifndef __REGULAR_MESH_HPP__
#define __REGULAR_MESH_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

namespace spurt {

typedef nvis::vec2              point_type;
typedef nvis::bbox2             bounds_type;

struct regular_mesh {
    regular_mesh(const bounds_type& bounds, size_t res_x, size_t res_y)
        : __bounds(bounds), __res_x(res_x), __res_y(res_y) {
        point_type diagonal = __bounds.size();
        __dx = diagonal[0] / (double)__res_x;
        __dy = diagonal[1] / (double)__res_y;
    }
    
    int idx(const point_type& x) const {
        if (!__bounds.inside(x)) {
            return -1;
        }
        point_type p = x - __bounds.min();
        int i = (int)floor(p[0] / __dx);
        int j = (int)floor(p[1] / __dy);
        return i + j*__res_x;
    }
    
    point_type center(int n) const {
        int i = n % __res_x;
        int j = n / __res_x;
        return __bounds.min() + point_type(((double)i + 0.5) * __dx,
                                           ((double)j + 0.5) * __dy);
    }
    
    point_type random_point(int n) const {
        int i = n % __res_x;
        int j = n / __res_x;
        return __bounds.min() + point_type(((double)i + drand48()) * __dx,
                                           ((double)j + drand48()) * __dy);
    }
    
    bounds_type     __bounds;
    size_t          __res_x, __res_y;
    double          __dx, __dy;
};


}


#endif











