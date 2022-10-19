#ifndef __XAVIER_COLLAB_GMIG_UTILS_HPP__
#define __XAVIER_COLLAB_GMIG_UTILS_HPP__

// stl
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
// nvis
#include <math/fixed_vector.hpp>

namespace xavier { namespace gmig {

inline double deg2rad(double deg) {
    static const double k = M_PI/180.;
    return deg*k;
}

inline double d_acos(double z) {
    return -1/sqrt(1-z*z);
}

// compute km distance between two points expressed in longitude / latitude
// formula found (by Kim) here: 
// http://bluemm.blogspot.com/2007/01/excel-formula-to-calculate-distance.html
// and successfully tested on Hui Huang's data
inline double distance_in_km(const nvis::vec2& lola0, const nvis::vec2& lola1) {
    static const double MEAN_EARTH_RADIUS = 6371;
    double theta0 = deg2rad(lola0[0]);
    double theta1 = deg2rad(lola1[0]);
    double phi0 = deg2rad(90-lola0[1]);
    double phi1 = deg2rad(90-lola1[1]);
    return acos(cos(phi0)*cos(phi1) + sin(phi0)*sin(phi1)*cos(theta0-theta1))*
           MEAN_EARTH_RADIUS;
}

inline nvis::vec2 d_distance_in_km(const nvis::vec2& lola_var, const nvis::vec2& lola_fix) {
    static const double MEAN_EARTH_RADIUS = 6371;
    double v = deg2rad(lola_var[0]);
    double v0 = deg2rad(lola_fix[0]);
    double u = deg2rad(90-lola_var[1]);
    double u0 = deg2rad(90-lola_fix[1]);
    
    double a = cos(u0);
    double b = sin(v0);
    double c = v0;
    
    nvis::vec2 r;
    r[0] = (-a*sin(u) + b*cos(u)*cos(v-c))*
            d_acos(cos(u)*a + sin(u)*b*cos(v-c))*M_PI/180.;
    r[1] = (b*sin(u)*sin(v-c))*
            d_acos(cos(u)*a + sin(u)*b*cos(v-c))*M_PI/180.;
    r *= MEAN_EARTH_RADIUS;
    return r;
}

inline double km_per_angle_ratio(const nvis::vec2& lola, const nvis::vec2& direction, double eps=0.0001) {
    // note: selecting a value of epsilon less than 1.0e-4 leads to 
    //       wrong results
    const nvis::vec2& from = lola;
    nvis::vec2 to = lola + eps*direction;
    double km = distance_in_km(from, to);
    return km/eps;
}

struct ProgressDisplay {
    bool active;
    float progress;
    int bar_width;
    int precision;
    size_t size;
    // std::ostream& os;
      
    ProgressDisplay(bool activate=true/*, std::ostream& _os=std::cout*/) 
        : active(activate), progress(0), bar_width(60), precision(3)
          /*, os(_os)*/{}

    void start(size_t _size) {
        size=_size;
        progress=0;
        if (active) display(); 
    }

    void end() { 
        if (active) {
            display(); 
            std::cout << std::endl;
        }
    }

    void update(size_t at) { 
        progress=(float)at/(float)size;
        if (active) display(); 
    }
    
    void display() {
        std::cout << "[";
        int pos=bar_width * progress;
        std::cout << std::string(pos, '=') << '>' 
                  << std::string(bar_width-pos-1, ' ')
                  << "] " 
                  << std::setprecision(precision)
                  << std::fixed
                  << progress*100.0
                  << " %\r"
                  << std::flush;
    }
};

class progress_message {
    std::string        _what;
    size_t             _size;
    
    mutable std::ostringstream _os;
    
public:
    progress_message(size_t size, std::string what = "")
    : _what(what), _size(size) {}
    
    std::string operator()(size_t n, double elapsed=0) const {
        _os.clear();
        _os.str("");
        _os << "\rCompleted " << 100.*(float)n/float(_size)
        << "\% of " << _size;
        if (_what.size()) {
            _os << " " << _what;
        }
        if (elapsed > 0) {
            _os << " in " << elapsed << " seconds (" << (float)n/elapsed << " Hz)";
        }
        _os << "                       \r";
        return _os.str();
    }
    
    void reset(const std::string& what = "") {
        if (what.size()) {
            _what = what;
        }
    }
};

template<typename Point_, typename Scalar_ = typename Point_::scalar_t>
struct default_distance {
    typedef Scalar_ scalar_t;
    typedef Point_  point_t;
    
    scalar_t operator()(const point_t& p1, const point_t& p2) const {
        return (p1-p2).norm();
    }
};

template<typename Point_, typename Distance_ = default_distance<Point_>,
         typename Scalar_ = typename Point_::scalar_t>
inline Scalar_ point_curve_distance(const Point_& p, const std::vector<Point_>& c) {
    typedef Scalar_   scalar_t;
    typedef Point_    point_t;
    typedef Distance_ distance_t;
      
    distance_t dm;
    scalar_t _min = std::numeric_limits<scalar_t>::max();
    std::for_each(c.begin(), c.end(), [&](const point_t& x) 
    {
        scalar_t d = dm(x,p);
        if (d<_min) _min=d;
    });  
    
    return _min;
}

template<typename Point_, typename Distance_ = default_distance<Point_>,
         typename Scalar_ = typename Point_::scalar_t>
inline Scalar_ hausdorff_distance(const std::vector<Point_>& c1, const std::vector<Point_>& c2, Scalar_& mean_dist, Scalar_& std_dev) {
    typedef Point_              point_t;
    typedef Distance_           distance_t;
    typedef Scalar_             scalar_t;
    typedef std::vector<Point_> curve_t;
    
    std::vector<scalar_t> p_dist;
    
    // Note: brute force implementation
    scalar_t _max12 = std::numeric_limits<scalar_t>::min();
    std::for_each(c1.begin(), c1.end(), [&](const point_t& p)
    {
        scalar_t d = 
            point_curve_distance<point_t, distance_t, scalar_t>(p, c2);
        if (d > _max12) _max12 = d;
        p_dist.push_back(d);
    });

    scalar_t _max21 = std::numeric_limits<scalar_t>::min();
    std::for_each(c2.begin(), c2.end(), [&](const point_t& p)
    {
        scalar_t d = 
            point_curve_distance<point_t, distance_t, scalar_t>(p, c1);
        if (d > _max21) _max21 = d;
        p_dist.push_back(d);
    });
    
    mean_dist = 
        std::accumulate(p_dist.begin(), p_dist.end(), 0) /
            static_cast<scalar_t>(p_dist.size());
    std_dev = sqrt(
        std::accumulate(p_dist.begin(), p_dist.end(), 0,
                        [&](scalar_t a, scalar_t b)
        {
            return a + (b-mean_dist)*(b-mean_dist);
        }) / static_cast<scalar_t>(p_dist.size()));
    
    return std::max(_max12, _max21);
}

template<typename Key_>
class Counter : private std::map<Key_, size_t> 
{
    typedef typename std::map<Key_, size_t> base_t;
public:
    typedef Key_ key_t;
    typedef typename base_t::value_type value_t;
    typedef typename base_t::iterator iterator;
    typedef typename base_t::const_iterator const_iterator;
private:
    typedef std::pair<iterator, bool> answer_t;
public:
    Counter() : base_t(), __total_size(0), 
    __this(static_cast<base_t*>(this)),
    __cthis(static_cast<const base_t*>(this)) {}
    
    size_t increment(const key_t& key) {
        ++__total_size;
        answer_t what=__this->insert(value_t(key, 1));
        if (!what.second) ++(what.first->second);
        return what.first->second;
    }
    
    size_t decrement(const key_t& key) {
        iterator it=__this->find(key);
        if (it == __this->end() || 
            it->second==0) return 0;
        --__total_size;
        return --(it->second);
    }
    
    size_t operator[](const key_t& key) const {
        const_iterator it=__cthis->find(key);
        if (it != __cthis->end()) 
            return it->second;
        else return 0;
    }
    
    std::pair<size_t, size_t> size() const {
        return std::make_pair(__cthis->size(), __total_size);
    }
    
    iterator begin() { return __this->begin(); }
    iterator end() { return __this->end(); }
    const_iterator begin() const { return __cthis->begin(); }
    const_iterator end() const { return __cthis->end(); }
    
private:
    size_t __total_size;
    base_t* __this;
    const base_t* __cthis;
};

template<typename Key_>
std::ostream& operator<<(std::ostream& os, const Counter<Key_>& count) {
    size_t total_size=count.size().second;
    os << "Count distribution: (" << total_size << " items)\n";
    std::for_each(count.begin(), count.end(), 
        [&](const typename Counter<Key_>::value_t& value) {
        os << value.first << ": " << value.second 
           << " (" << (float)value.second/(float)total_size*100.
           << "%)\n";
    });
    return os;
} 

template<typename Pos_, typename Point_>
struct LessDistFrom {
    LessDistFrom() : from(0,0,0) {}
    LessDistFrom(const Pos_& _from) : from(_from) {}
    LessDistFrom(const Point_& _point) : from(_point.coordinate()) {}
    
    bool operator()(const Point_& a, const Point_& b) const {
        return nvis::norm(a.coordinate()-from) < 
                          nvis::norm(b.coordinate()-from);
    }
    
    Pos_ from;
};

template<typename Scalar_, typename Pos_, typename Order_>
struct LessPosTolerance {
    LessPosTolerance() : eps(1.0e-3) {}
    LessPosTolerance(Scalar_ _eps) : eps(_eps) {}
    bool operator()(const Pos_& p1, const Pos_& p2) const {
        if (nvis::norm(p1-p2) < eps) return false;
        else return order(p1, p2);
    }
    Scalar_ eps;
    Order_ order;
};

} // namespace gmig
} // namespace xavier

#endif