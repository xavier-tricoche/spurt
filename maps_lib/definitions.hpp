#ifndef __DEFINITIONS_HPP__
#define __DEFINITIONS_HPP__

#include "maps_lib/misc.hpp"
#include "poincare/metric.hpp"
#include "maps_lib/interval.hpp"
#include <boost/rational.hpp>
#include "maps_lib/triangulation.hpp"
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <tokamak/map2d.hpp>
#include <data/kdtree.hpp>

using namespace spurt;

namespace map_analysis {

class orbit;
class point_data;

typedef boost::rational<int>                                    rational_type;
typedef spurt::interval<double>                                interval_type;
typedef spurt::point_locator<double, int, 2>                   locator_type;
typedef spurt::triangulation<point_data, locator_type>         mesh_type;
typedef spurt::default_metric_type                             metric_type;

namespace static_data {

extern std::vector<orbit>               central_map_orbits;
extern metric_type                      metric;

} // map_analysis::static_data

namespace newton {
template<typename T>
struct cached_map {
    typedef T   map_type;
    
    cached_map(const map_type& map, unsigned int p)
        : __period(p) {
        __map = map.clone();
        clear_cache();
    }
    
    void clear_cache() const {
        __last_f = nvis::vec2(std::numeric_limits<double>::min(),
                              std::numeric_limits<double>::min());
    }
    
    void update_cache(const nvis::vec2& x) const {
        if (nvis::all(x == __last_x)) {
            return;
        }
        __last_x = x;
        map2d::value_type result;
        __map->jmap_complete(x, result, __period);
        __last_return = result.x;
        __last_f = static_data::metric.displacement(x, result.x);
        __last_J = result.J - nvis::mat2::identity();
        __last_err = result.err;
    }
    
    nvis::vec2 operator()(const nvis::vec2& x) const {
        update_cache(x);
        return __last_return;
    }
    
    nvis::vec2 f(const nvis::vec2& x) const {
        update_cache(x);
        return __last_f;
    }
    
    nvis::mat2 J(const nvis::vec2& x) const {
        update_cache(x);
        return __last_J;
    }
    
    nvis::vec2 error(const nvis::vec2& x) const {
        update_cache(x);
        return __last_err;
    }
    
    nvis::vec2 one_step(const nvis::vec2& x) const {
        clear_cache();
        nvis::vec2 y;
        try {
            y = static_data::metric.modulo(__map->map(x, 1));
        } catch (...) {
            std::cerr << "caught exception in one_step\n";
        }
        return y;
    }
    
    void set_precision(double eps) const {
        __map->precision(eps);
        clear_cache();
    }
    
    double get_precision() const {
        return __map->get_precision();
    }
    
    int period() const {
        return __period;
    }
    
    const metric_type& metric() const {
        return static_data::metric;
    }
    
    mutable map_type*   __map;
    unsigned int        __period;
    mutable nvis::vec2  __last_x;
    mutable nvis::vec2  __last_f;
    mutable nvis::mat2  __last_J;
    mutable nvis::vec2  __last_err;
    mutable nvis::vec2  __last_return;
};


template<typename T>
struct rhs_map {
    typedef T   cached_map_type;
    
    rhs_map(const cached_map_type& map)
        : __map(map) {}
        
    nvis::vec2 operator()(const nvis::vec2& x) const {
        return __map.f(x);
    }
    
    nvis::vec2 one_step(const nvis::vec2& x) const {
        return __map.one_step(x);
    }
    
    nvis::vec2 p_step(const nvis::vec2& x) const {
        return __map(x);
    }
    
    nvis::vec2 error(const nvis::vec2& x) const {
        return __map.error(x);
    }
    
    void set_precision(double eps) const {
        __map.set_precision(eps);
    }
    
    double get_precision()  const {
        return __map.get_precision();
    }
    
    int period() const {
        return __map.period();
    }
    
    const metric_type& metric() const {
        return __map.metric();
    }
    
    const cached_map_type& __map;
};

template<typename T>
struct J_map {
    typedef T   cached_map_type;
    
    J_map(const cached_map_type& map)
        : __map(map) {}
        
    nvis::mat2 operator()(const nvis::vec2& x) const {
        nvis::mat2 j = __map.J(x);
        return j;
    }
    
    void set_precision(double eps) const {
        __map.set_precision(eps);
    }
    
    double get_precision()  const {
        return __map.get_precision();
    }
    
    const cached_map_type& __map;
};

} // newton

} // map_analysis


#endif





























