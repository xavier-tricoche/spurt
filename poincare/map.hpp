#ifndef ____MAP_HPP__
#define ____MAP_HPP__

#include <iostream>
#include <vector>
#include <list>

#include <boost/bind.hpp>
#include <boost/thread.hpp>

#include <util/timer.hpp>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>

#include <math/angle.hpp>
#include <poincare/macros.hpp>
#include <poincare/basic_math.hpp>
#include <poincare/ls.hpp>
#include <poincare/metric.hpp>

#define PARALLEL

namespace spurt {

namespace map_debug {
extern int verbose_level;
extern std::vector< nvis::vec2 > jacobian_sample_pos;
extern std::vector< nvis::vec2 > jacobian_sample_vals;
extern double jacobian_timer;
}

/////////////////////////////////////////////////////////////////////
//
//                          miscellaneous
//
/////////////////////////////////////////////////////////////////////
inline nvis::vec4 operator*(const nvis::vec4& a, const nvis::vec4& b)
{
    nvis::vec4 c;
    c[0] = a[0] * b[0] + a[1] * b[2];
    c[1] = a[0] * b[1] + a[1] * b[3];
    c[2] = a[2] * b[0] + a[3] * b[2];
    c[3] = a[2] * b[1] + a[3] * b[3];
    return c;
}

/////////////////////////////////////////////////////////////////////
//
//                    useful map wrappers
//
/////////////////////////////////////////////////////////////////////
template< typename Map>
class map_wrapper {
public:
    map_wrapper(const Map& map, unsigned int p = 1)
        : _map(map), _p(p) {}
        
    nvis::vec2 operator()(const nvis::vec2& x, int n = 1) const {
        return _map.map(x, n*_p);
    }
    
    nvis::mat2 derivative(const nvis::vec2& x, int n = 1) const {
        nvis::mat2 r;
        nvis::vec2 hit;
        _map.jmap(x, hit, r, n*_p);
        return r;
    }
    
    nvis::vec2 single_step(const nvis::vec2& x) const {
        return _map.map(x, 1);
    }
    
    void reset_period(unsigned int p) {
        _p = p;
    }
    
private:
    const Map& _map;
    unsigned int _p;
};

template< typename Map >
class map_from_last_wrapper {
public:
    map_from_last_wrapper(const Map& map, unsigned int p = 1)
        : _map(map), _p(p) {}
        
    nvis::vec2 operator()(const nvis::vec2& x, int n = 1) const {
        return spurt::__default_metric.displacement(x, _map.map(x, n*_p));
    }
    
    nvis::mat2 derivative(const nvis::vec2& x, int n = 1) const {
        nvis::mat2 r;
        nvis::vec2 hit;
        _map.jmap(x, hit, r, n*_p);
        r[0][0] -= 1.;
        r[1][1] -= 1.;
        return r;
    }
    
    nvis::vec2 single_step(const nvis::vec2& x) const {
        return _map.map(x, 1);
    }
    
    void reset_period(unsigned int p) {
        _p = p;
    }
    
private:
    const Map& _map;
    unsigned int _p;
};

template< typename Map >
class map_from_fixed_location_wrapper {
public:
    map_from_fixed_location_wrapper(const Map& map, const nvis::vec2& x0, unsigned int p = 1)
        : _map(map), _x0(x0), _p(p) {}
        
    nvis::vec2 operator()(const nvis::vec2& x, int n = 1) const {
        return spurt::__default_metric.displacement(_x0, _map.map(x, n*_p));
    }
    
    nvis::mat2 derivative(const nvis::vec2& x, int n = 1) const {
        nvis::mat2 r;
        nvis::vec2 hit;
        _map.jmap(x, hit, r, n*_p);
        return r;
    }
    
    nvis::vec2 single_step(const nvis::vec2& x) const {
        return _map.map(x, 1);
    }
    
    void reset_period(unsigned int p) {
        _p = p;
    }
    
private:
    const Map& _map;
    const nvis::vec2& _x0;
    unsigned int _p;
};

/////////////////////////////////////////////////////////////////////
//
//                      Jacobian computation(s)
//
/////////////////////////////////////////////////////////////////////
template< typename Map >
class central_diff_jacobian {
public:
    central_diff_jacobian(const Map& map, double dx, double dy)
        : _map(map), _dx(dx), _dy(dy) {}
        
    nvis::mat2 operator()(const nvis::vec2& x) const {
        nvis::mat2 J;
        nvis::vec2 y[] = { spurt::__default_metric.modulo(x - nvis::vec2(_dx, 0)),
                           spurt::__default_metric.modulo(x + nvis::vec2(_dx, 0)),
                           spurt::__default_metric.modulo(x - nvis::vec2(0, _dy)),
                           spurt::__default_metric.modulo(x + nvis::vec2(0, _dy))
                         };
                         
        // std::cerr << "central difference Jacobian at " << x
        // << " with dx = " << _dx << " and dy = " << _dy << '\n';
        // std::cerr << "sampled points are "
        //           << y[0] << ", " << y[1] << ", " << y[2] << ", " << y[3]
        //           << '\n';
        
        nvis::vec2 d_dx = 0.5 * (_map(y[1]) - _map(y[0])) / _dx;
        nvis::vec2 d_dy = 0.5 * (_map(y[3]) - _map(y[2])) / _dy;
        
        // std::cerr << "d_dx = " << d_dx << '\n';
        // std::cerr << "d_dy = " << d_dy << '\n';
        
        J[0][0] = d_dx[0];
        J[1][0] = d_dx[1];
        J[0][1] = d_dy[0];
        J[1][1] = d_dy[1];
        return J;
    }
    
private:
    const Map& _map;
    double _dx, _dy;
};

template< typename Map >
class least_sq_jacobian {
public:
    least_sq_jacobian(const Map& map, double dx, double dy)
        : _map(map), _dx(dx), _dy(dy) {
        double dtheta = 2.*M_PI / 8.;
        double theta = 0.;
        for (unsigned int i = 0 ; i < 8 ; i++, theta += dtheta) {
            p[i] = nvis::vec2(_dx * cos(theta), _dy * sin(theta));
        }
    }
    
    nvis::mat2 operator()(const nvis::vec2& x) const {
        std::vector< nvis::vec2 > pos(8);
        std::vector< nvis::vec2 > vec(8);
        
        try {
            for (unsigned int i = 0 ; i < 8 ; ++i) {
                pos[i] = spurt::__default_metric.modulo(x + p[i]);
                vec[i] = _map(pos[i]);
            }
        } catch (...) {
            WARNING_MACRO(0, "invalid position in LS_jacobian: " << x << '\n');
            nvis::mat2 zero;
            zero[0][0] = zero[0][1] = zero[1][0] = zero[1][1] = 0.;
            return zero;
        }
        
        nvis::vec4 _j = ls_jacobian(pos, vec, x);
        nvis::mat2 J;
        J[0][0] = _j[0];
        J[0][1] = _j[1];
        J[1][0] = _j[2];
        J[1][1] = _j[3];
    }
    
private:
    const Map& _map;
    double _dx, _dy;
    nvis::vec2 p[8];
};

template< typename Map, typename Jacobian >
class cumulative_jacobian {
public:
    // NB: this constructor must be called with a Jacobian corresponding to
    // a single iteration of the map.
    cumulative_jacobian(const Map& map, const Jacobian& one_jacobian,
                        unsigned int p)
        : _map(map), _jacobian(one_jacobian), _p(p) {}
        
    nvis::mat2 operator()(const nvis::vec2& x) const {
        nvis::mat2 j[_p], J_p(nvis::mat2::identity());
        nvis::vec2 y = x;
        try {
            for (unsigned int i = 0 ; i < _p ; ++i) {
                j[i] = _jacobian(y);
                y = _map(y, 1);
                J_p = j[i] * J_p;
            }
        } catch (...) {
            nvis::mat2 zero;
            zero[0][0] = zero[0][1] = zero[1][0] = zero[1][1] = 0.;
            return zero;
        }
        return J_p;
    }
    
    void evaluate_Js(const nvis::vec2& x, std::vector<nvis::mat2>& J) const {
        nvis::vec2 y = x;
        J.resize(_p);
        for (unsigned int i = 0 ; i < _p ; ++i) {
            J[i] = _jacobian(y);
            std::cout << "jacobian(" << y << ") = " << J[i] << std::endl;
            y = _map(y, 1);
        }
    }
    
private:
    const Map& _map;
    const Jacobian& _jacobian;
    unsigned int _p;
};

template< typename Map >
class integral_jacobian {
public:
    integral_jacobian(const Map& map)
        : _map(map) {}
        
    nvis::mat2 operator()(const nvis::vec2& x) const {
        return _map.derivative(x);
    }
    
private:
    const Map& _map;
};

/////////////////////////////////////////////////////////////////////
//
//                      Periods measure
//
/////////////////////////////////////////////////////////////////////

template< typename Map >
unsigned int period(const Map& map, const nvis::vec2& x, unsigned int pmax,
                    double tolerance = 1.0e-6)
{
    if (map_debug::verbose_level > 1) {
        std::cout << "period called at " << x << std::endl;
    }
    
    unsigned int p;
    double dmin = std::numeric_limits<double>::max();
    unsigned int pmin = 0;
    
    std::vector< nvis::vec2 > hits;
    try {
        map.map(x, hits, pmax);
    } catch (...) {
        return 0;
    }
    for (p = 0 ; p < pmax ; ++p) {
        const nvis::vec2& y = hits[p];
        double d = spurt::__default_metric.distance(x, y);
        if (d < tolerance) {
            return p + 1;
        } else if (d < dmin) {
            dmin = d;
            pmin = p + 1;
        }
    }
    // we did not meet the prescribed tolerance - return the best we got
    if (map_debug::verbose_level > 1) {
        std::cout << "unable to determine period at " << x
                  << " - min dist was " << dmin
                  << " (>" << tolerance << ") for period " << pmin
                  << std::endl;
    }
    
    return pmin;
}

/////////////////////////////////////////////////////////////////////
//
//              Classical examples of analytical maps
//
/////////////////////////////////////////////////////////////////////

class LinearMap {
public:
    LinearMap(double dt) : _dt(dt), __metric() {
    }
    nvis::vec2 map(const nvis::vec2& p) const {
        static nvis::vec2 _p;
        nvis::vec4 J(0, -1, 1, 0);
        _p[0] = p[0] + _dt * J[0] * (p[0] - 0.5) + _dt * J[1] * (p[1] - 0.5);
        _p[1] = p[1] + _dt * J[2] * (p[0] - 0.5) + _dt * J[3] * (p[1] - 0.5);
        return _p;
    }
    
    const spurt::default_metric_type& metric() const {
        return __metric;
    }
    
private:
    double _dt;
    spurt::default_metric_type __metric;
};

class StandardMap
{
	static constexpr double twopi = 2.*M_PI;

	inline void forward(double& x, double& y) const {
		y -= _k / twopi * sin(twopi * x);
		x += y;
	}

	inline void backward(double& x, double& y) const {
		x -= y;
		y += _k / twopi * sin(twopi * x);
	}

	inline nvis::vec2 _map(const nvis::vec2& x, bool fwd = true) const {
		nvis::vec2 y = x;
		if (fwd)
			forward(y[0], y[1]);
		else
			backward(y[0], y[1]);
		if (is_periodic) {
			y = spurt::__default_metric.modulo(y);
		}
		return y;
	}


public:
    typedef spurt::default_metric_type metric_type;
    
    StandardMap(double k) : _k(k), __metric() {
        is_periodic = true;
        __metric.bounds().min() = nvis::vec2(0, 0);
        __metric.bounds().max() = nvis::vec2(1, 1);
        __metric.periodic(0) = true;
        __metric.periodic(1) = true;
    }
    
    double K() const {
        return _k;
    }
    
    nvis::vec2 map(const nvis::vec2& x, int n = 1) const {
        nvis::vec2 y = x;
        for (unsigned int i = 0 ; i < abs(n) ; ++i) {
            y = _map(y, (n > 0));
        }
        return y;
    }
    
    void map(const nvis::vec2& x, std::vector< nvis::vec2 >& hits, int n = 1) const {
        hits.resize(abs(n));
        nvis::vec2 y = x;
        for (unsigned int i = 0; i < abs(n) ; ++i) {
            y = _map(y, (n > 0));
            hits[i] = y;
        }
    }
    
    void set_to_periodic(bool periodic = true) {
        is_periodic = periodic;
        __metric.periodic(0) = periodic;
        __metric.periodic(1) = periodic;
    }
    
    const metric_type& metric() const {
        return __metric;
    }
    
private:
    double _k;
    bool is_periodic;
    metric_type __metric;
};

class Tokamap {
public:
    typedef spurt::default_metric_type metric_type;
    
    Tokamap(const double l) : _l(l), __metric() {
        __metric.periodic(0) = true;
        __metric.periodic(1) = true;
    }
    
    nvis::vec2 map(const nvis::vec2 X, unsigned int _p = 1) const {
        static nvis::vec2 Y;
        
        _psi = X[0] * X[0] + X[1] * X[1];
        _theta = atan(X[1] / X[0]);
        if (X[0] < 0) {
            _theta += M_PI;
        }
        
        for (unsigned int i = 0 ; i < _p ; i++) {
            _b = (1 - _psi + _l * sin(2 * M_PI * _theta));
            _delta = _b * _b + 4 * _psi;
            
            _psi = 0.5 * (-_b + sqrt(_delta));
            _theta = _theta + 1 / q_profile(_psi) -
                     0.5 * _l / M_PI / ((1 + _psi) * (1 + _psi)) * cos(2 * M_PI * _theta);
        }
        
        Y[0] = sqrt(_psi) * cos(_theta);
        Y[1] = sqrt(_psi) * sin(_theta);
        return Y;
    }
    
    const metric_type& metric() const {
        return __metric;
    }
    
private:
    double _l;
    mutable double _psi, _theta;
    mutable double _delta, _b;
    
    double q_profile(double psi) const {
        return 4 / ((2 - psi)*(2 - 2*psi + psi*psi));
    }
    metric_type __metric;
};

template< typename Map >
bool track_map_rotation(const Map& the_map,
                        const nvis::vec2& x, double theta0, double theta1, double radius,
                        std::vector< std::pair<double, nvis::vec2> >& _samples,
                        unsigned int max_samples, double dtheta,
                        nvis::vec2& minnormpos)
{
    typedef std::pair< double, nvis::vec2 >     Sample;
    typedef std::list< Sample >::iterator       Iterator;
    
    std::list< Sample > samples;
    _samples.clear();
    nvis::vec2 y0, y1, f0, f1, f;
    y0 = x + radius * nvis::vec2(cos(theta0), sin(theta0));
    f0 = the_map(y0);
    samples.push_back(Sample(theta0, f0));
    y1 = x + radius * nvis::vec2(cos(theta1), sin(theta1));
    f1 = the_map(y1);
    samples.push_back(Sample(theta1, f1));
    double theta = unsigned_angle<2>(f0 , f1);
    if (theta < dtheta) {
        // std::cout << "initial angle at " << x << " is " << theta << ", which is valid\n";
        if (nvis::norm(f0) > nvis::norm(f1)) {
            minnormpos = y1;
        } else {
            minnormpos = y0;
        }
        std::copy(samples.begin(), samples.end(), std::back_inserter(_samples));
        return true;
    }
    
    unsigned int sz = 2;
    double maxangle = theta;
    double nf0 = nvis::norm(f0);
    double nf1 = nvis::norm(f1);
    double minnorm;
    if (nf0 < nf1) {
        minnorm = nf0;
        minnormpos = y0;
    } else {
        minnorm = nf1;
        minnormpos = y1;
    }
    while (maxangle > dtheta && sz < max_samples) {
        maxangle = 0.;
        Iterator cur = samples.begin();
        Iterator next = cur;
        ++next;
        for (; next != samples.end() ; ++cur, ++next) {
            double alpha = unsigned_angle<2>(cur->second, next->second);
            if (alpha > maxangle) {
                maxangle = alpha;
            }
            if (alpha > dtheta) {
                double& __theta0 = cur->first;
                double& __theta1 = next->first;
                nvis::vec2 x0 = x + radius * nvis::vec2(cos(__theta0), sin(__theta0));
                nvis::vec2 x1 = x + radius * nvis::vec2(cos(__theta1), sin(__theta1));
                
                // std::cout << "TRACK: angle around " << x << " between vector " << cur->second << " at " << x0
                // << " (=" << 180*__theta0 / M_PI << " degrees) and vector " << next->second << " at "
                // << x1 << " (=" << 180*__theta1 / M_PI << " degrees) is too large: " << 180*alpha / M_PI << " > "
                // << (int)floor(degrees(dtheta)) << " degrees" << std::endl;
                
                double theta = 0.5 * (cur->first + next->first);
                nvis::vec2 y = x + radius * nvis::vec2(cos(theta), sin(theta));
                nvis::vec2 f = the_map(y);
                double nf = nvis::norm(f);
                if (nf < minnorm) {
                    minnorm = nf;
                    minnormpos = y;
                    // std::cout << "minimum norm so far is " << minnorm << " sfound at " << y
                    // << " with an angular gap = " << 180*alpha / M_PI << "\n";
                }
                cur = samples.insert(next, Sample(theta, f));
                map_debug::jacobian_sample_pos.push_back(y);
                map_debug::jacobian_sample_vals.push_back(f);
                ++sz;
            }
        }
    }
    
    std::copy(samples.begin(), samples.end(), std::back_inserter(_samples));
    return (maxangle < dtheta);
}


template< typename T >
bool sample_map_on_circle(std::vector< std::pair< nvis::vec2, nvis::vec2 > >& _samples,
                          const T& the_map, const nvis::vec2& x, double radius,
                          double dtheta, double eps, unsigned int max_nb_pts = 50)
{
    assert(dtheta > 0 && dtheta <= M_PI);
    
    typedef std::pair< double, nvis::vec2 >     Sample;
    typedef std::list< Sample >::iterator       Iterator;
    std::list< Sample > samples;
    
    // start with coarse and uniform sampling
    for (unsigned int i = 0 ; i < 8 ; ++i) {
        double theta = (double)i * M_PI / 4.;
        nvis::vec2 y = x + radius * nvis::vec2(cos(theta), sin(theta));
        nvis::vec2 f = the_map(y);
        samples.push_back(Sample(theta, f));
        map_debug::jacobian_sample_pos.push_back(y);
        map_debug::jacobian_sample_vals.push_back(f);
    }
    samples.push_back(Sample(2*M_PI, samples.front().second)); // cyclic list
    
    unsigned int sz = 8;
    double maxangle = M_PI;
    while (maxangle > dtheta && sz < max_nb_pts) {
        maxangle = 0.;
        Iterator cur = samples.begin();
        Iterator next = cur;
        ++next;
        for (; next != samples.end() ; ++cur, ++next) {
            double alpha = unsigned_angle<2>(cur->second, next->second);
            if (alpha > maxangle) {
                maxangle = alpha;
            }
            if (alpha > dtheta) {
                double& theta0 = cur->first;
                double& theta1 = next->first;
                nvis::vec2 x0 = x + radius * nvis::vec2(cos(theta0), sin(theta0));
                nvis::vec2 x1 = x + radius * nvis::vec2(cos(theta1), sin(theta1));
                
                // std::cout << "SAMPLE: angle around " << x << " between vector " << cur->second << " at " << x0
                // << " (=" << 180*theta0 / M_PI << " degrees) and vector " << next->second << " at "
                // << x1 << " (=" << 180*theta1 / M_PI << " degrees) is too large: " << 180*alpha / M_PI << " > "
                // << (int)floor(degrees(dtheta)) << " degrees" << std::endl;
                
                double theta = 0.5 * (cur->first + next->first);
                nvis::vec2 y = x + radius * nvis::vec2(cos(theta), sin(theta));
                nvis::vec2 f = the_map(y);
                double nf = nvis::norm(f);
                // std::cout << "norm at added point is " << nf << std::endl;
                //          if (nf < eps) {
                //              std::cout << "norm underflow at " << y << std::endl;
                //          }
                
                cur = samples.insert(next, Sample(theta, f));
                map_debug::jacobian_sample_pos.push_back(y);
                map_debug::jacobian_sample_vals.push_back(f);
                ++sz;
            }
        }
    }
    
    _samples.clear();
    for (Iterator it = samples.begin() ; it != samples.end() ; ++it) {
        double theta = it->first;
        _samples.push_back(std::pair<nvis::vec2, nvis::vec2>(
                               x + radius * nvis::vec2(cos(theta),
                                       sin(theta)),
                               it->second));
    }
    _samples.pop_back();
    
    if (maxangle > dtheta) {
        return false;
    } else {
        return true;
    }
}

template< typename T >
nvis::vec4 index_based_jacobian(const T& the_map, const nvis::vec2& x, double dx, double dy,
                                double dtheta = 0.6)
{
    std::vector< std::pair< nvis::vec2, nvis::vec2 > > samples;
    sample_map_on_circle(samples, the_map, x, std::min(dx, dy), dtheta);
    
    std::vector< nvis::vec2 > ps, fs;
    
    for (unsigned int i = 0 ; i < samples.size() ; ++i) {
        ps.push_back(samples[i].first);
        fs.push_back(samples[i].second);
    }
    
    return ls_jacobian(ps, fs, x);
}

template< typename Map >
int poincare_index(std::vector< std::pair< nvis::vec2, nvis::vec2 > >& samples,
                   const Map& map, const nvis::vec2& x,
                   double radius, double dtheta, unsigned int nb_pts = 50)
{
    samples.clear();
    bool ok = sample_map_on_circle(samples, map, x, radius, dtheta, nb_pts);
    if (!ok) {
        return 0;
    }
    samples.push_back(samples[0]); // close the loop
    
    double sigma_theta = 0.;
    for (unsigned int i = 0 ; i < samples.size() - 1 ; ++i) {
        sigma_theta += signed_angle(samples[i].second, samples[i+1].second);
    }
    samples.pop_back();
    
    sigma_theta /= 2.*M_PI;
    return (int)round(sigma_theta);
}


};

#endif
































































































































































































































