#ifndef __XAVIER_QUALITY_CONTROL_HPP__
#define __XAVIER_QUALITY_CONTROL_HPP__

#include <math/fixed_vector.hpp>
#include <vector>
#include <set>
#include <maps/basic_definitions.hpp>
#include <algorithm>
#include <math/bounding_box.hpp>


namespace xavier {

inline double area(const nvis::vec2 p[3])
{
    nvis::vec2 e0 = p[1] - p[0];
    nvis::vec2 e1 = p[2] - p[0];
    return 0.5*(e0[0] * e1[1] - e0[1] * e1[0]);
}

inline nvis::vec3 edge_lengths(const nvis::vec2 p[3])
{
    nvis::vec3 l;
    for (int i = 0 ; i < 3 ; ++i) {
        int j = (i + 1) % 3;
        l[i] = nvis::norm(p[i] - p[j]);
    }
    return l;
}

template<typename T>
struct value_stats {

    value_stats() : min(0), max(0), mean(0), var(0) {}
    value_stats(const std::vector<T>& vals) {
        if (vals.empty()) {
            return;
        }
        double n = vals.size();
        mean = min = max = vals[0];
        typename std::vector<T>::const_iterator it = vals.begin();
        for (++it ; it != vals.end() ; ++it) {
            mean += *it;
            if (*it < min) {
                min = *it;
            } else if (*it > max) {
                max = *it;
            }
        }
        mean /= n;
        var = 0;
        for (it = vals.begin() ; it != vals.end() ; ++it) {
            double d = (*it - mean);
            var += d * d;
        }
        var /= n;
    }
    
    double min, max, mean, var;
};

template<typename M>
inline void statistics(const M& mesh, value_stats<double>& area_stats,
                       value_stats<double>& edge_stats, const nvis::bbox2& bounds)
{
    int ntri = mesh.get_nb_triangles();
    std::vector<double> _area, _edge;
    _area.reserve(ntri);
    _edge.reserve(3*ntri);
    // we are assuming that the imposed bounding box excludes the boundary
    // and as such ensures that each edge is going to be counted exactly twice
    // thus removing the need for a uniqueness check
    nvis::ivec3 ids;
    nvis::vec2 p[3];
    nvis::vec3 l;
    for (int i = 0 ; i < ntri ; ++i) {
        ids = mesh.get_triangle_vertices(i);
        bool skip = false;
        for (int k = 0 ; k < 3 ; ++k) {
            p[k] = mesh.get_vertex(ids[k]);
            if (!bounds.inside(p[k])) {
                skip = true;
                break;
            }
        }
        if (skip) {
            continue;
        }
        _area.push_back(area(p));
        l = edge_lengths(p);
        for (int k = 0 ; k < 3 ; ++k) {
            _edge.push_back(l[k]);
        }
    }
    area_stats = value_stats<double>(_area);
    edge_stats = value_stats<double>(_edge);
}


inline double max_length(const nvis::vec2 p[3])
{
    nvis::vec3 l = edge_lengths(p);
    return *std::max_element(&l[0], &l[3]);
}

inline double aspect_ratio(const nvis::vec2 p[3])
{
    double l, min, max;
    min = max = nvis::norm(p[0] - p[1]);
    for (int i = 1 ; i < 3 ; ++i) {
        int j = (i + 1) % 3;
        l = nvis::norm(p[i] - p[j]);
        if (l < min) {
            min = l;
        } else if (l > max) {
            max = l;
        }
    }
    return max / min;
}

inline bool zero_vector(double b[3], nvis::vec2 v[3])
{
    b[1]    = (-v[2][0] * (v[1][1] - v[2][1])) - (-v[2][1] * (v[1][0] - v[2][0]));
    b[2]    = (-v[2][1] * (v[0][0] - v[2][0])) - (-v[2][0] * (v[0][1] - v[2][1]));
    double denom = ((v[0][0] - v[2][0]) * (v[1][1] - v[2][1])) -
                   ((v[0][1] - v[2][1]) * (v[1][0] - v[2][0]));
    if (denom != 0) {
        b[1] /= denom;
        b[2] /= denom;
        b[0] = 1. - b[1] - b[2];
        return b[0] >= 0 && b[1] >= 0 && b[2] >= 0;
    } else {
        return false;
    }
}

inline interval<double> span(const double v[3])
{
    return interval<double> (*std::min_element(&v[0], &v[3]),
                             *std::max_element(&v[0], &v[3]));
}

struct triangle_area_controller {
    typedef xavier::point_data  data_type;
    
    triangle_area_controller(double max) : __max_area(max) {}
    
    double operator()(const nvis::vec2 p[3], const data_type v[3]) const {
        return std::max(area(p) - __max_area, 0.);
    }
    
    double __max_area;
};

struct triangle_aspect_ratio_controller {
    typedef xavier::point_data  data_type;
    
    triangle_aspect_ratio_controller(double max_r) : __max_ratio(max_r) {}
    
    double operator()(const nvis::vec2 p[3], const data_type v[3]) const {
        return std::max(aspect_ratio(p) - __max_ratio, 0.);
    }
    
    double __max_ratio;
};

struct triangle_edge_length_controller {
    typedef xavier::point_data  data_type;
    
    triangle_edge_length_controller(double max) : __max_length(max) {
        std::cerr << "initializing edge_length_controller with threshold = " << max << std::endl;
    }
    
    double operator()(const nvis::vec2 p[3], const data_type v[3]) const {
        return std::max(max_length(p) - __max_length, 0.);
    }
    
    double __max_length;
};

struct triangle_shape_controller {
    typedef xavier::point_data  data_type;
    typedef nvis::bbox2         bounds_type;
    
    triangle_shape_controller(double max_a, double max_r)
        : __area_controller(max_a), __ratio_controller(max_r), __bounds(), __use_bounds(false) {}
        
    triangle_shape_controller(double max_a, double max_r, const bounds_type& bounds)
        : __area_controller(max_a), __ratio_controller(max_r), __bounds(bounds), __use_bounds(true) {}
        
    double operator()(const nvis::vec2 p[3], const data_type v[3]) const {
        if (__use_bounds) {
            for (int i = 0 ; i < 3 ; ++i) {
                if (!__bounds.inside(p[i])) {
                    return 0.;
                }
            }
        }
        return std::max(__area_controller(p, v), __ratio_controller(p, v));
    }
    
    triangle_area_controller            __area_controller;
    triangle_aspect_ratio_controller    __ratio_controller;
    bounds_type                         __bounds;
    bool                                __use_bounds;
};

struct triangle_multiple_orbits_controller {
    typedef xavier::point_data  data_type;
    
    triangle_multiple_orbits_controller(double min_a)
        : __min_area(min_a) {}
        
    double operator()(const nvis::vec2 p[3], const data_type v[3]) const {
        if (v[1].orbit_id() == v[0].orbit_id() ||
                v[2].orbit_id() == v[0].orbit_id() ||
                v[2].orbit_id() == v[1].orbit_id()) {
            return std::max(area(p) - __min_area, 0.);
        } else {
            return 0;
        }
    }
    
    double __min_area;
};

struct triangle_value_inclusion_controller {
    typedef xavier::point_data  data_type;
    typedef interval<double>    interval_type;
    
    triangle_value_inclusion_controller(double val, double min_a)
        : __value(val), __min_area(min_a) {}
        
    double operator()(const nvis::vec2 p[3], const data_type v[3]) const {
        double cell_area = area(p);
        if (cell_area <= __min_area) {
            return 0.;
        }
        double q[3] = {v[0].period(), v[1].period(), v[2].period()};
        interval<double> cell_span = span(q);
        if (!cell_span.inside(__value)) {
            return 0.;
        }
        return cell_span.length();
    }
    
    double __value;
    double __min_area;
};

struct triangle_interval_inclusion_controller {
    typedef xavier::point_data  data_type;
    typedef interval<double>    interval_type;
    
    triangle_interval_inclusion_controller(const interval_type& valid, double min_a)
        : __valid_range(valid), __min_area(min_a) {}
        
    double operator()(const nvis::vec2 p[3], const data_type v[3]) const {
        double cell_area = area(p);
        if (cell_area <= __min_area) {
            return 0.;
        }
        double q[3] = {v[0].period(), v[1].period(), v[2].period()};
        interval<double> cell_span = span(q);
        interval<double> intersection = intersect(cell_span, __valid_range);
        if (intersection.empty()) {
            return 0.;
        }
        return cell_span.length();
    }
    
    interval_type __valid_range;
    double __min_area;
};

struct triangle_value_span_controller {
    typedef xavier::point_data  data_type;
    typedef interval<double>    interval_type;
    
    triangle_value_span_controller(double max_s, double min_a)
        : __max_span(max_s), __min_area(min_a) {}
        
    double operator()(const nvis::vec2 p[3], const data_type v[3]) const {
        double cell_area = area(p);
        if (cell_area <= __min_area) {
            return 0.;
        }
        double q[3] = {v[0].period(), v[1].period(), v[2].period()};
        return std::max(0., span(q).length() - __max_span);
    }
    
    double __max_span;
    double __min_area;
};

struct triangle_contains_fixed_point_controller {
    typedef xavier::point_data  data_type;
    
    triangle_contains_fixed_point_controller(double min_a, int period,
            const map_metric& metric)
        : __min_area(min_a), __period(period), __metric(metric) {}
        
    double operator()(const nvis::vec2 p[3], const data_type d[3]) const {
        double cell_area = area(p);
        if (cell_area <= __min_area) {
            return 0.;
        }
        
        nvis::vec2 v[3];
        for (int i = 0 ; i < 3 ; ++i) {
            v[i] = vector_value(d[i], __period, __metric);
        }
        double b[3];
        bool zero = zero_vector(b, v);
        if (!zero) {
            return 0.;
        } else {
            return (cell_area - __min_area);
        }
    }
    
    double              __min_area;
    int                 __period;
    const map_metric&   __metric;
};

}


#endif


































































