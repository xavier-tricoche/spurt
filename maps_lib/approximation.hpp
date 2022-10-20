#ifndef __APPROXIMATION_HPP__
#define __APPROXIMATION_HPP__

#include "maps_lib/definitions.hpp"
#include "maps_lib/orbits.hpp"
#include "poincare/metric.hpp"
#include "maps_lib/quality_control.hpp"
#include "maps_lib/misc.hpp"
#include <math/angle.hpp>

namespace map_analysis {

// functors for adaptive triangulation

// goal: quality of piecewise linear approximation
struct error_measure_linear {
    typedef point_data  data_type;
    
    error_measure_linear(unsigned int period, double threshold, const spurt::map_metric& metric)
        : __period(period), __eps(threshold), __metric(metric) {}
        
    double operator()(const data_type& d,
                      const nvis::vec3& beta,
                      const data_type data[3]) const {
        nvis::vec2 ref = vector_value(d, __period, __metric);
        nvis::vec2 v[3];
        for (int i = 0 ; i < 3 ; ++i) {
            v[i] = vector_value(data[i], __period, __metric);
        }
        
        nvis::vec2 approx = beta[0] * v[0] + beta[1] * v[1] + beta[2] * v[2];
        return nvis::norm(approx - ref) / nvis::norm(ref) - __eps;
    }
    
    unsigned int __period;
    double __eps;
    spurt::map_metric __metric;
};

// goal: piecewise uniform p-map norm
struct error_measure_norm {
    typedef point_data  data_type;
    
    error_measure_norm(unsigned int period, const metric_type& metric)
        : __period(period), __metric(metric) {}
        
    double operator()(const data_type& d,
                      const nvis::vec3& beta,
                      const data_type data[3]) const {
        nvis::vec2 ref = vector_value(d, __period, __metric);
        double dref = nvis::norm(ref);
        double dist[3];
        for (int i = 0 ; i < 3 ; ++i) {
            nvis::vec2 v = vector_value(data[i], __period, __metric);
            dist[i] = nvis::norm(v);
        }
        
        if (dref < *std::min_element(&dist[0], &dist[3])) {
            return 1. / dref;
        } else {
            return -1;
        }
    }
    
    unsigned int __period;
    metric_type __metric;
};

// goal: limited rotation of p-map direction
struct angular_variation_priority {
    typedef nvis::vec2                      point_type;
    typedef map_analysis::point_data        data_type;
    typedef map_analysis::interval_type     interval_type;
    typedef map_analysis::rational_type     rational_type;
    
    angular_variation_priority(unsigned int period, double eps,
                               const metric_type& metric, double dq)
        : _p(period), _eps(eps), _dq(dq), _metric(metric) {
        for (int i = 1 ; i <= period ; ++i) {
            rational_type r(period, i);
            if (r.numerator() != period) {
                continue;
            }
            double q = spurt::value(r);
            _valid_q.push_back(interval_type(q - dq, q + dq));
        }
    }
    
    double operator()(const point_type p[3], const data_type d[3]) const {
    
        // check the easy part: is the triangle degenerate?
        if (max_length(p) > 0.05*static_data::metric.diameter()) {
            return -1;
        }
        
        double dtheta = 0;
        std::list<double> qs;
        for (int i = 0 ; i < 3 ; ++i) {
            qs.push_back(d[i].period());
        }
        double minq = *std::min_element(qs.begin(), qs.end());
        double maxq = *std::max_element(qs.begin(), qs.end());
        interval_type span(minq - _dq, maxq + _dq);
        
        bool valid = false;
        for (int i = 0 ; i < _valid_q.size() && !valid ; ++i) {
            valid = !(intersect(span, _valid_q[i]).empty());
        }
        if (!valid) {
            return -1;
        }
        
        for (int i = 0 ; i < 3 ; ++i) {
            dtheta += fabs(spurt::signed_angle(vector_value(d[i], _p, _metric),
                                                vector_value(d[(i+1)%3], _p, _metric)));
        }
        return area(p)*(dtheta - _eps);
    }
    
    double _eps, _dq;
    unsigned int _p;
    std::vector<interval_type> _valid_q;
    const metric_type& _metric;
};

// goal: limited rotation of p-map rotation coupled with requirements
//       on map precision and approximated period
//
// NB:   does not seem to be working all that well...
struct conditional_angular_variation_priority {
    typedef nvis::vec2                  point_type;
    typedef point_data          data_type;
    
    conditional_angular_variation_priority(unsigned int period,
                                           const std::vector<double>& qs,
                                           double qeps, double eps,
                                           const metric_type& metric)
        : _p(period), _ints(qs.size()), _eps(eps), _metric(metric) {
        for (int i = 0 ; i < qs.size() ; ++i) {
            _ints[i] = interval_type(qs[i] - qeps, qs[i] + qeps);
        }
    }
    
    double operator()(const point_type p[3], const data_type d[3]) const {
    
        // bool display_stuff = (_p == 1) && (saddle_box.inside(p[0]) ||
        //                                    saddle_box.inside(p[1]) ||
        //                                    saddle_box.inside(p[2]));
        
        // double min, max;
        // min = max = d[0].period();
        // for (int i = 1 ; i < 3 ; ++i) {
        //  double q = d[i].period();
        //  if (q < min) min = q;
        //  else if (q > max) max = q;
        // }
        // interval_type cell_int(min, max);
        // bool relevant = false;
        // for (int i = 0 ; i < _ints.size() && !relevant ; ++i) {
        //  relevant = !(spurt::intersect(_ints[i], cell_int).empty());
        // }
        
        // if (!relevant) {
        //  if (display_stuff) {
        //      std::cerr << "triangle: " << p[0] << ", " << p[1] << ", " << p[2]
        //                << " has been deemed irrelevant because its associated periods are "
        //                << d[0].period() << ", " << d[1].period() << ", and " << d[2].period()
        //                << '\n';
        //  }
        //  return 0;
        // }
        
        double dtheta = 0;
        double maxerror = 0;
        for (int i = 0 ; i < 3 ; ++i) {
            std::pair<nvis::vec2, nvis::vec2> vecerr0 = vector_and_error_value(d[i], _p, _metric);
            std::pair<nvis::vec2, nvis::vec2> vecerr1 = vector_and_error_value(d[(i+1)%3], _p, _metric);
            dtheta += fabs(spurt::signed_angle(vecerr0.first, vecerr1.first));
            double err = nvis::norm(vecerr0.second);
            if (err > maxerror) {
                maxerror = err;
            }
        }
        
        // // if triangle area is smaller than half of square of side maxerror, bail
        // if (2*spurt::area(p) < maxerror*maxerror) return 0;
        
        // if (display_stuff) {
        //  std::cerr << "triangle: " << p[0] << ", " << p[1] << ", " << p[2]
        //            << " has a priority value of " << dtheta - _eps << '\n';
        // }
        return dtheta - _eps;
    }
    
    unsigned int _p;
    std::vector<interval_type> _ints;
    double _eps;
    const metric_type& _metric;
};

// checks if period span of triangle is compatible with provided period.
// test is based on period approximation tolerance dq
struct triangle_filter {
    typedef nvis::vec2  point_type;
    typedef point_data  data_type;
    
    triangle_filter(unsigned int period, double eps,
                    const metric_type& metric, double dq)
        : _p(period), _eps(eps), _dq(dq), _metric(metric) {
        for (int i = 1 ; i <= period ; ++i) {
            rational_type r(period, i);
            if (r.numerator() != period) {
                continue;
            }
            double q = spurt::value(r);
            _valid_q.push_back(interval_type(q - dq, q + dq));
        }
    }
    
    bool valid(const point_type p[3], const data_type d[3]) const {
        double minq, maxq;
        minq = maxq = d[0].period();
        for (int i = 1 ; i < 3 ; ++i) {
            if (d[i].period() > maxq) {
                maxq = d[i].period();
            } else if (d[i].period() < minq) {
                minq = d[i].period();
            }
        }
        interval_type span(minq - _dq, maxq + _dq);
        
        bool valid = false;
        for (int i = 0 ; i < _valid_q.size() && !valid ; ++i) {
            valid = !(intersect(span, _valid_q[i]).empty());
        }
        return valid;
    }
    
    unsigned int _p;
    double _eps, _dq;
    const metric_type& _metric;
    std::vector<interval_type> _valid_q;
};

}

#endif



