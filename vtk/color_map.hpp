#ifndef __XAVIER_COLORMAP_HPP__
#define __XAVIER_COLORMAP_HPP__

#include <map>
#include <vector>
#include <math/fixed_vector.hpp>

namespace spurt {
struct color_map_bwr {
    typedef nvis::fvec3             color_type;
    
    color_map_bwr(const std::vector<float>& values, float gamma = 1) : __gamma(gamma) {
        std::vector<float> __vals(values);
        std::sort(__vals.begin(), __vals.end());
        float min = __vals.front();
        float max = __vals.back();
        
        __width = std::max(fabs(min), fabs(max));
    }
    
    color_map(float min, float max, float gamma = 1) : __gamma(gamma) {
        __width = std::max(fabs(min), fabs(max));
    }
    
    color_type operator()(float v) const {
        color_type hue;
        if (v < 0) {
            hue = color_type(0, 0, 1);
        } else {
            hue = color_type(1, 0, 0);
        }
        
        float u = fabs(v) / __width;
        float x = pow(u, __gamma);
        return (1. - x)*color_type(1, 1, 1) + x*hue;
    }
    
    const float& width() const {
        return __width;
    }
    
    std::map<float, color_type>     __map;
    float                           __width, __gamma;
};

}


#endif



