#ifndef __COLORS_HPP__
#define __COLORS_HPP__

#include <math/fixed_vector.hpp>
#include <vector>
#include <map>
#include <cassert>

namespace spurt {
template<typename T>
inline double luminosity(const nvis::fixed_vector<T,3>& color, double gamma=2.2) {
    return 0.2126 * pow(color[0], gamma) + 
           0.7152 * pow(color[1], gamma) + 
           0.0722 * pow(color[2], gamma);
}

inline nvis::vec3 hsv2rgb(double hue, double sat, double val)
{
    double chroma = val * sat;
    double h_ = hue / 60.;
    double x = chroma * (1.-fabs(fmod(h_, 2.)-1));
    
    nvis::vec3 rgb_;
    if(val == 0) {
        rgb_ = nvis::vec3(0, 0, 0);
    } else if(0 <= h_ && h_ < 1) {
        rgb_ = nvis::vec3(chroma, x, 0);
    } else if(1 <= h_ && h_ < 2) {
        rgb_ = nvis::vec3(x, chroma, 0);
    } else if(2 <= h_ && h_ < 3) {
        rgb_ = nvis::vec3(0, chroma, x);
    } else if(3 <= h_ && h_ < 4) {
        rgb_ = nvis::vec3(0, x, chroma);
    } else if(4 <= h_ && h_ < 5) {
        rgb_ = nvis::vec3(x, 0, chroma);
    } else if(5 <= h_ && h_ < 6) {
        rgb_ = nvis::vec3(chroma, 0, x);
    }
    
    double m = val-chroma;
    return rgb_+  nvis::vec3(m, m, m);
}

const nvis::fvec3 red(1, 0, 0);
const nvis::fvec3 green(0, 1, 0);
const nvis::fvec3 blue(0, 0, 1);
const nvis::fvec3 yellow(1, 1, 0);
const nvis::fvec3 cyan(0, 1, 1);
const nvis::fvec3 white(1, 1, 1);
const nvis::fvec3 black(0, 0, 0);
const nvis::fvec3 orange(1, 0.5, 0);
const nvis::fvec3 magenta(1, 0, 1);
// some funky color names after Apple's color editor
const nvis::fvec3 cayenne(0.5, 0, 0);
const nvis::fvec3 midnight(0, 0, 0.5);
const nvis::fvec3 aqua(0, 0.5, 1);
const nvis::fvec3 sea_green(0, 1, 0.5);
const nvis::fvec3 lime(0.5, 1, 0);
const nvis::fvec3 framboise(1, 0, 0.5);
const nvis::fvec3 carnation(1, 0.5, 1);

const nvis::fvec3 rainbow[] = {
    black, midnight, blue, aqua,
    cyan, sea_green, green, lime,
    yellow, orange, red, framboise,
    magenta, carnation, white
}; // 15 colors spanning the rainbow

static void spiral_scale(std::vector<nvis::fvec3>& colors, int n, 
                         double minval, int r=1, 
                         double initsat=1, double endsat=1,
                         double order=1) {
    // double dhue = (double)r*360./(double)n;
    // double dval = (1-minval)/(double)n;
    // double dsat = (endsat-initsat)/(double)n;
    colors.resize(n);
    // double hue = 0;
    // double val = minval;
    // double sat = initsat;
    double hue, val, sat;
    for (int i=0 ; i<n ; ++i) {
        double u = pow((double)i/(double)(n-1), order);
        hue = u*r*360.;
        if (hue > 360) {
            hue -= 360;
        }
        val = (1-u)*minval + u*1;
        sat = (1-u)*initsat + u*endsat;
        colors[i] = hsv2rgb(hue, sat, val);
    }
}

// 0------1------2------3------4------5
// 0-----------1----------2-----------3

template<typename T>
struct adaptive_color_map {
    struct More {
        bool operator()(const T& v1, const T& v2) {
            return v2 > v1;
        }
    };
	
	nvis::fvec3 interpolate(double u, const std::vector<nvis::fvec3>& scale) {
		assert(!scale.empty());
		if (u<=0 || scale.size()==1) return scale[0];
		else if (u>=1) return scale.back();
		else {
			double du = 1./static_cast<double>(scale.size()-1);
			int id = static_cast<int>(std::floor(u / du));
			u -= id*du;
			return (1.-u)*scale[id] + u*scale[id+1];
		}
	}
	
    adaptive_color_map(const std::vector<T>& vals, 
                       const std::vector<nvis::fvec3>& scale,
                       bool _ascending=true,
					   int ncolors=-1)
        : t(scale.size()), colors(scale), ascending(_ascending) {
		assert(scale.size() > 0);
		if (ncolors > 0) {
			t.resize(ncolors);
			colors.resize(ncolors);
			// consider a few trivial special cases
			if (ncolors == 1) colors[0] = scale[scale.size()/2];
			else {
				colors.front() = scale.front();
				colors.back() = scale.back();
				double du = 1./static_cast<double>(ncolors-1);
				for (int i=1; i<colors.size()-1; ++i) {
					colors[i] = interpolate(i*du, scale);
				}
			}
		}	
		
        std::vector<T> tmp(vals);
        std::sort(tmp.begin(), tmp.end());
        unsigned int step = tmp.size() / (colors.size()-1);
        for (int i = 0 ; i < t.size() ; ++i) {
            t[i] = tmp[i*step];
        }
        t.back() = tmp.back();
    }
    
    nvis::fvec3 operator()(const T& val) const {
        unsigned int bin = 
            std::distance(t.begin(), 
                          std::lower_bound(t.begin(), t.end(), val));
        if (bin > t.size()-2) bin = t.size()-2;
        T u = (val-t[bin]) / (t[bin+1]-t[bin]);
        if (ascending) {
            return (1.-u)*colors[bin]+  u*colors[bin+1];
        }
        else {
            return (1.-u)*colors[colors.size()-1-bin]+u*colors[colors.size()-bin-2];
        }
    }

    std::vector<T>           t;
    std::vector<nvis::fvec3> colors;
    bool ascending;
};

template<typename T>
struct discrete_color_map {

    discrete_color_map(const std::vector<T>& cps, 
                       const std::vector<nvis::fvec3>& scale)
        : colors(scale) {
        assert(cps.size() == colors.size());
        for (int i=0 ; i<cps.size() ; ++i) {
            lookup[cps[i]] = i;
        }
    }
    
    nvis::fvec3 operator()(const T& val) const {
        typename std::map<T, int>::const_iterator it = lookup.find(val);
        if (it == lookup.end()) {
            return nvis::fvec3(0,0,0);
        } else {
            return colors[it->second];
        }
    }

    std::map<T, int>            lookup;
    std::vector<nvis::fvec3>     colors;
};

template<typename T>
struct fixed_color_map {

    fixed_color_map(const std::vector<T>& cps, 
                    const std::vector<nvis::fvec3>& scale)
        : t(cps), colors(scale) {
        std::sort(t.begin(), t.end());
        assert(t.size() == colors.size());
    }
    
    fixed_color_map(T min, T max, const std::vector<nvis::fvec3>& scale)
        : t(scale.size()), colors(scale) {    
        int n = scale.size();
        t[0] = min;
        t.back() = max;
        T dv = (max-min)/(T)(n-1);
        T v = min;
        for (int i=1 ; i<t.size()-1 ; ++i) {
            v += dv;
            t[i] = v;
        }    
    }
    
    nvis::fvec3 operator()(const T& val) const {
        unsigned int bin = 
            std::distance(t.begin(),
                          std::lower_bound(t.begin(), t.end(), val));
        if (bin > t.size()-2) bin = t.size()-2;
        T u = (val-t[bin]) / (t[bin+1]-t[bin]);
        return (1.-u)*colors[bin]+  u*colors[bin+1];
    }

    std::vector<T>                t;
    std::vector<nvis::fvec3>     colors;
};

template<typename T>
struct band_color_map {
    band_color_map(const std::vector<T>& cps, const std::vector<nvis::fvec3>& scale)
        : t(cps), colors(scale) {
        std::sort(t.begin(), t.end());
        assert(t.size()+2 == colors.size());
    }
    
    nvis::fvec3 operator()(const T& val) const {
        unsigned int bin = 
            std::distance(t.begin(),
                          std::lower_bound(t.begin(), t.end(), val));
        return colors[bin];
    }

    std::vector<T>                t;
    std::vector<nvis::fvec3>     colors;
};

} // spurt

#endif
