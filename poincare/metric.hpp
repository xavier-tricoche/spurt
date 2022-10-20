#ifndef __XAVIER_METRIC_HPP__
#define __XAVIER_METRIC_HPP__

#include <vector>
#include <limits>
#include <map>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <utility>
#include <assert.h>
#include <misc/misc_helper.hpp>
#include <math/vector_manip.hpp>

namespace {

    template<typename T, typename Bounds>
    struct bounds_traits{};
    
    template<typename T>
    struct bounds_traits<T, nvis::bounding_box<T> > {
        typedef T value_type;
        typedef nvis::bounding_box<T> bounds_type;
        
        static T& min(bounds_type& b) { return b.min(); }
        static const T& min(const bounds_type& b) { return b.min(); }
        static T& max(bounds_type& b) { return b.max(); }
        static const T& max(const bounds_type& b) { return b.max(); }
        static bool inside(const bounds_type& b, const value_type& v) { return b.inside(v); }
    };
    
    template<typename T>
    struct bounds_traits<T, std::pair<T, T> > {
        typedef T value_type;
        typedef std::pair<T, T> bounds_type;
        
        static T& min(bounds_type& b) { return b.first(); }
        static const T& min(const bounds_type& b) { return b.first(); }
        static T& max(bounds_type& b) { return b.second(); }
        static const T& max(const bounds_type& b) { return b.second(); }
        static bool inside(const bounds_type& b, const value_type& v) {
            value_type v1 = v - b.first();
            if (spurt::vector::min(v1) < 0) return false;
            v1 = b.second() - v;
            if (spurt::vector::min(v1) < 0) return false;
            return true;
        }
    };

    template<typename T>
    struct bounds_traits<T, std::array<T, 2> > {
        typedef T value_type;
        typedef std::array<T, 2> bounds_type;
        
        static T& min(bounds_type& b) { return b[0]; }
        static const T& min(const bounds_type& b) { return b[0]; }
        static T& max(bounds_type& b) { return b[1]; }
        static const T& max(const bounds_type& b) { return b[1]; }
        static bool inside(const bounds_type& b, const value_type& v) {
            value_type v1 = v - b[0];
            if (spurt::vector::min(v1) < 0) return false;
            v1 = b[1] - v;
            if (spurt::vector::min(v1) < 0) return false;
            return true;
        }
    };

}


namespace spurt {
    
template<int N, typename Pos_ = nvis::fixed_vector<double, N>, 
         typename Bounds_ = nvis::bounding_box<nvis::fixed_vector<double, N> > >
class map_metric {
public:
    static constexpr int dimension = N;
    typedef Pos_     pos_type;
    typedef Bounds_  bounds_type;
    typedef bounds_traits<pos_type, bounds_type> traits_type;
    typedef std::array<bool, N> bvecN;

    map_metric()
            : m_bounds(pos_type(-std::numeric_limits<double>::max()),
                       pos_type(std::numeric_limits<double>::max())) {
        std::fill(m_periodic.begin(), m_periodic.end(), false);
    }

    map_metric(const bounds_type& bounds, bool periodic[2])
            : m_bounds(bounds) {
        static_assert(dimension == 2, "Invalid periodic array size");
        m_periodic[0] = periodic[0];
        m_periodic[1] = periodic[1];
    }
    
    map_metric(const bounds_type& bounds, const bvecN& periodic)
        : m_bounds(bounds), m_periodic(periodic) {}
    
    map_metric(const map_metric& metric)
            : m_bounds(metric.m_bounds) {
        m_periodic = metric.m_periodic;
    }

    const bounds_type& bounds() const {
        return m_bounds;
    }

    bounds_type& bounds() {
        return m_bounds;
    }
    
    pos_type size() const {
        return traits_type::max(m_bounds) - traits_type::min(m_bounds);
    }

    double width() const {
        return m_bounds.size()[0];
    }

    double height() const {
        return m_bounds.size()[1];
    }

    double diameter() const {
        return spurt::vector::norm(size());
    }

    bool periodic(size_t i) const {
        assert(i < dimension);
        return m_periodic[i];
    }

    bool& periodic(size_t i) {
        assert(i < dimension);
        return m_periodic[i];
    }

    pos_type displacement(const pos_type& a, const pos_type& b, const bvecN& per) const {
        pos_type dx = (b - a);
        pos_type p = dx;
        pos_type sz = size();
        
        for (int i=0; i<dimension; ++i) {
            if (per[i]) {
                p[i] = sign(dx[i]) * fmod(fabs(dx[i]), sz[i]);
                p[i] = mindist(p[i], sz[i]);
            }
        }
        return p;
    }

    pos_type displacement(const pos_type& a, const pos_type& b) const {
        return displacement(a, b, m_periodic);
    }

    double distance(const pos_type& a, const pos_type& b) const {
        return spurt::vector::norm(displacement(a, b));
    }

    double distance(const nvis::vec2& a, const nvis::vec2& b, const bvecN& per) const {
        return nvis::norm(displacement(a, b, per));
    }

    double angle(const pos_type& x0, const pos_type& x1, const pos_type& x2) const {
        pos_type v0 = displacement(x0, x1);
        pos_type v1 = displacement(x1, x2);
        v0 /= spurt::vector::norm(v0);
        v1 /= spurt::vector::norm(v1);
        
        pos_type v2 = v1;
        v2 -= spurt::vector::dot(v1, v0)*v0;
        v2 /= spurt::vector::norm(v2);
        
        double cos_alpha = spurt::vector::dot(v0, v1);
        double sin_alpha = spurt::vector::dot(v1, v2);
        double alpha = acos(cos_alpha);
        if (sin_alpha < 0) {
            alpha *= -1;
        }
        return alpha;
    }
    
    double modulo(double x, int dim) const {
        if (!m_periodic[dim]) return x;
        double y = x - traits_type::min(m_bounds)[dim];
        double span = traits_type::max(m_bounds)[dim] - traits_type::min(m_bounds)[dim];
        y = fmod(y, span);
        if (y<0) y += span;
        y += traits_type::min(m_bounds)[dim];
        return y;
    }

    pos_type modulo(const pos_type& x) const {
        pos_type p;
        for (int i=0; i<dimension; ++i) {
            p[i] = modulo(x[i], i);
        }
        return p;
    }

    void clip_segment(std::vector<std::pair<pos_type, pos_type> >& segments,
                      const pos_type& a, const pos_type& b) const {
        segments.clear();
        pos_type d = displacement(a, b);
        pos_type __a = modulo(a);
        if (!traits_type::inside(m_bounds, __a)) return;
        if (traits_type::inside(m_bounds, __a + d)) {
            segments.resize(1);
            segments[0] = std::make_pair(__a, __a + d);
        }
        else {
            std::array<double, dimension> t;
            std::fill(t.begin(), t.end(), std::numeric_limits<double>::max());
            pos_type __b = modulo(b);
            pos_type __x, __y;
            pos_type dir = d / spurt::vector::norm(d);
            // compute up/down and left/right intersection with boundary
            for (int i = 0 ; i < dimension ; ++i) {
                double w = (dir[i] > 0 ? traits_type::max(m_bounds)[i] : traits_type::min(m_bounds)[i]);
                t[i] = (w - __a[i]) / dir[i];
            }
            // which intersection is closest?
            int imin = spurt::vector::min_id(t);
            __x = __a + t[imin] * dir;
            // split segment accordingly
            if (m_periodic[imin]) {
                __y = __x;
                __y[imin] = (dir[imin] > 0 ? traits_type::min(m_bounds)[imin] : traits_type::max(m_bounds)[imin]);
                segments.resize(2);
                segments[0] = std::make_pair(__a, __x);
                segments[1] = std::make_pair(__y, __b);
            }
            // if not periodic along that direction, simply discard second half
            else {
                segments.resize(1);
                segments[0] = std::make_pair(__a, __x);
            }
        }
    }

private:
    bounds_type                 m_bounds;
    std::array<bool, dimension> m_periodic;
};

typedef map_metric<2, nvis::vec2, nvis::bbox2> default_metric_type;

static default_metric_type __default_metric;

static const default_metric_type euclidean_metric = default_metric_type();

}

#endif








































