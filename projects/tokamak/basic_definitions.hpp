#ifndef __MAPS_BASIC_DEFINITIONS_HPP__
#define __MAPS_BASIC_DEFINITIONS_HPP__

#include <vector>
#include <list>
#include <math/types.hpp>
#include <boost/rational.hpp>
#include <poincare/metric.hpp>
#include "period.hpp"
#include "tokamak/map2d.hpp"

namespace spurt {
class orbit {
public:
    orbit() : __points(), __q(-1.) {}
    
    orbit(const std::vector<spurt::vec2>& points, double period = -1)
        : __points(points), __q(period) {}
        
    orbit(const std::vector<spurt::vec2>& points, const std::vector<spurt::vec2>& errors,
          double period = -1)
        : __points(points), __errors(errors), __q(period) {}
        
    orbit(const std::list<spurt::vec2>& points, double period = -1)
        : __points(points.begin(), points.end()), __q(period) {}
        
    orbit(const std::list<spurt::vec2>& points, const std::list<spurt::vec2>& errors, double period = -1)
        : __points(points.begin(), points.end()), __errors(errors.begin(), errors.end()), __q(period) {}
        
    size_t size() const {
        return __points.size();
    }
    
    const spurt::vec2& operator[](size_t i) const {
        assert(i < __points.size());
        return __points[i];
    }
    
    spurt::vec2& operator[](size_t i) {
        assert(i < __points.size());
        return __points[i];
    }
    
    const std::vector<spurt::vec2>& points() const {
        return __points;
    }
    
    const std::vector<spurt::vec2>& errors() const {
        return __errors;
    }
    
    double period() const {
        return __q;
    }
    
    double& period() {
        return __q;
    }
    
private:
    std::vector<spurt::vec2>     __points;
    std::vector<spurt::vec2>     __errors;
    double __q;
};

extern std::vector<orbit> __map_orbits;

class point_data {
public:
    point_data() : __orbit(-1), __index(-1) {}
    point_data(size_t which_orbit, size_t idx = -1)
        : __orbit(which_orbit), __index(idx) {}
        
    size_t orbit_id() const {
        assert(__orbit < __map_orbits.size());
        return __orbit;
    }
    
    void set_orbit(size_t which_orbit) {
        __orbit = which_orbit;
    }
    
    double period() const {
        assert(__orbit < __map_orbits.size());
        return __map_orbits[__orbit].period();
    }
    
    const spurt::vec2& pos() const {
        assert(__orbit < __map_orbits.size());
        return __map_orbits[__orbit][__index];
    }
    
    size_t index() const {
        return __index;
    }
    
    double scalar_value() const {
        return period();
    }
    
    const spurt::vec2& error() const {
        const std::vector<spurt::vec2>& errors = __map_orbits[__orbit].errors();
        return errors[__index];
    }
    
private:
    size_t  __orbit, __index;
};

template<typename map_metric>
spurt::vec2 vector_value(const point_data& pt, int period, const map_metric& metric)
{
    const orbit& chain = __map_orbits[pt.orbit_id()];
    size_t N = chain.size();
    if (period >= N) {
        std::cerr << "WARNING: attempting to determine " << period
                  << "-map on a chain of length " << N
                  << " seeded at " << __map_orbits[pt.orbit_id()][0]
                  << std::endl;
        double large = std::numeric_limits<double>::max();
        return spurt::vec2(large, large);
    }
    int i = pt.index();
    if (i >= period && i + period < N) {
        const spurt::vec2& x0 = chain[i-period];
        const spurt::vec2& x2 = chain[i+period];
        return 0.5*(metric.displacement(x0, x2));
    } else if (i + period < N) {
        const spurt::vec2& x1 = chain[i];
        const spurt::vec2& x2 = chain[i+period];
        return metric.displacement(x1, x2);
    } else if (i >= period) {
        const spurt::vec2& x0 = chain[i-period];
        const spurt::vec2& x1 = chain[i];
        return metric.displacement(x0, x1);
    } else {
        double large = std::numeric_limits<double>::max();
        return spurt::vec2(large, large);
        
        // std::cerr << "WARNING: lacking values for period " << period
        //           << " on a chain of length " << N
        //           << " seeded at " << __map_orbits[pt.orbit_id()][0]
        //           << std::endl;
        // return spurt::vec2(0, 0);
    }
}

template<typename map_metric>
std::pair<spurt::vec2, spurt::vec2> vector_and_error_value(const point_data& pt, int period, const map_metric& metric)
{
    const orbit& chain = __map_orbits[pt.orbit_id()];
    size_t N = chain.size();
    if (period >= N) {
        std::cerr << "WARNING: attempting to determine " << period
                  << "-map on a chain of length " << N
                  << " seeded at " << __map_orbits[pt.orbit_id()][0]
                  << std::endl;
        double large = std::numeric_limits<double>::max();
        return std::make_pair(spurt::vec2(large, large), large);
    }
    int i = pt.index();
    if (i >= period && i + period < N) {
        const spurt::vec2& x0 = chain[i-period];
        const spurt::vec2& x2 = chain[i+period];
        return std::make_pair(0.5*(metric.displacement(x0, x2)), chain.errors()[i-period] + chain.errors()[i+period]);
    } else if (i + period < N) {
        const spurt::vec2& x1 = chain[i];
        const spurt::vec2& x2 = chain[i+period];
        return std::make_pair(metric.displacement(x1, x2), chain.errors()[i] + chain.errors()[i+period]);
    } else if (i >= period) {
        const spurt::vec2& x0 = chain[i-period];
        const spurt::vec2& x1 = chain[i];
        return std::make_pair(metric.displacement(x0, x1), chain.errors()[i-period] + chain.errors()[i]);
    } else {
        double large = std::numeric_limits<double>::max();
        return std::make_pair(spurt::vec2(large, large), large);
        
        // std::cerr << "WARNING: lacking values for period " << period
        //           << " on a chain of length " << N
        //           << " seeded at " << __map_orbits[pt.orbit_id()][0]
        //           << std::endl;
        // return spurt::vec2(0, 0);
    }
}

template<typename map_metric>
spurt::vec2 vector_value(size_t orbit_id, size_t idx, size_t period, const map_metric& metric)
{
    point_data pd(orbit_id, idx);
    return vector_value(pd, period, metric);
}

template<typename map_metric>
std::pair<spurt::vec2, spurt::vec2> vector_and_error_value(size_t orbit_id, size_t idx, size_t period, const map_metric& metric)
{
    point_data pd(orbit_id, idx);
    return vector_and_error_value(pd, period, metric);
}

template<typename IntType>
inline double value(const boost::rational<IntType>& r)
{
    return (double)r.numerator() / (double)r.denominator();
}

template<typename T>
struct interval {
    interval() : __min(T(0)), __max(T(0)) {}
    interval(T min, T max) : __min(min), __max(max) {}
    interval(const interval<T>& _int) : __min(_int.__min), __max(_int.__max) {}
    
    bool inside(T val) const {
        return (val >= __min && val <= __max);
    }
    
    bool empty() const {
        return (__max <= __min);
    }
    
    T length() const {
        return std::max(0., __max - __min);
    }
    
    T __min, __max;
};

template<typename T>
inline interval<T> intersect(const interval<T>& i0, const interval<T>& i1)
{
    return interval<T>(std::max(i0.__min, i1.__min),
                       std::min(i0.__max, i1.__max));
}

template<typename I>
struct orbit_integrator {
    typedef I   integrator_type;
    
    orbit_integrator(integrator_type& integrator, size_t nsteps,
                     const spurt::map_metric& metric)
        : __nsteps(nsteps), __integ(integrator), __metric(metric) {}
        
    void operator()(const spurt::vec2& x0,
                    std::vector<spurt::vec2>& points,
                    std::vector<spurt::point_data>& data) const {
                    
        const integrator_type* pmap = __integ.clone();
        std::vector<typename map2d::value_type> returned_values;
        try {
            pmap->map_complete(x0, returned_values, __nsteps);
        } catch (...) {
            // std::cerr << "orbit_integrator: unable to integrate from " << x0 << std::endl;
            return;
        }
        if (returned_values.empty()) {
            // std::cerr << "orbit_integrator: unable to integrate from " << x0 << std::endl;
            return;
        }
        
        std::vector<spurt::vec2> errors(returned_values.size() + 1);
        points.resize(returned_values.size() + 1);
        data.resize(returned_values.size() + 1);
        points[0] = x0;
        errors[0] = 0.;
        for (int i = 1 ; i < points.size() ; ++i) {
            points[i] = returned_values[i-1].x;
            errors[i] = returned_values[i-1].err;
        }
        
        // period computation is done prior to applying congruence relation to coordinates
        double q = spurt::dist_based_x_period(points, __metric);
        
        size_t orbit_id = spurt::__map_orbits.size();
        for (int i = 0 ; i < points.size() ; ++i) {
            points[i] = __metric.modulo(points[i]);
            data[i] = spurt::point_data(orbit_id, i);
        }
        spurt::__map_orbits.push_back(spurt::orbit(points, errors, q));
    }
    
    void precision(double h) {
        __integ.precision(h);
    }
    
    size_t                  __nsteps;
    integrator_type&        __integ;
    spurt::map_metric      __metric;
};

template<typename T>
inline void push_front(const T& val, std::vector<T>& _in)
{
    if (!_in.size()) {
        _in.push_back(val);
    } else if (nvis::all(val != _in[0])) {
        std::vector<T> _out;
        _out.reserve(_in.size() + 1);
        _out.push_back(val);
        std::copy(_in.begin(), _in.end(), std::back_inserter(_out));
        std::swap(_in, _out);
    }
}

}




#endif






























