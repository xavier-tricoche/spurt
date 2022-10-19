#ifndef __new_universal_poincare_map_hpp
#define __new_universal_poincare_map_hpp

#include <vector>
#include <exception>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <math/dopri5.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <iomanip>

#include "mapNd.hpp"

template< typename RHS, int N, typename Section >
class new_universal_poincare_map : public mapNd<N> {
public:
    typedef RHS                             rhs_type;
    typedef Section                         section_type;
    typedef typename section_type::vec_type vec_type;

    new_universal_poincare_map(const rhs_type* rhs, const section_type* sec)
        : _rhs(rhs), _prec(1.0e-6), _sec(sec) {}

    new_universal_poincare_map(const new_universal_poincare_map& other)
    {
        _prec = other._prec;
        _rhs  = other._rhs->clone();
        _sec  = other._sec->clone();
    }

    void precision(const double& prec)
    {
        _prec = prec;
    }

    vec_type map(const vec_type& in, int niter) const
    {
        std::vector<vec_type> out;
        std::vector<double> ts;
        std::vector<vec_type> orbit;

        bool success;

        map_RKDP45(in, out, ts, niter, orbit, false);

        if (out.size() < std::abs(niter)) {
            typename mapNd<N>::map_undefined m_u;
            throw m_u;
        }
        return out.back();
    }

    const rhs_type* rhs() const
    {
        return _rhs;
    }

    void map(const vec_type& in, std::vector<vec_type>& out, int niter) const
    {
        std::vector< double > ts;
        std::vector< vec_type> orbit;
        map_RKDP45(in, out, ts, niter, orbit, false);
    }

    void map(const vec_type& in, std::vector<vec_type>& out, int niter,
             std::vector<vec_type>& orbit) const
    {
        std::vector< double > ts;
        map_RKDP45(in, out, ts, niter, orbit, true);
    }

    void map(const vec_type& seed, std::vector<vec_type>& steps,
             std::vector< double >& times, int niter) const
    {
        std::vector<vec_type> orbit;
        map_RKDP45(seed, steps, times, niter, orbit, false);
    }

    void map(const vec_type& seed, std::vector<vec_type>& steps,
             std::vector< double >& times, int niter,
             std::vector< vec_type>& orbit) const
    {
        map_RKDP45(seed, steps, times, niter, orbit, true);
    }

    nvis::bounding_box<vec_type> bounds() const
    {
        return _rhs->bounds();
    }

    new_universal_poincare_map* clone() const
    {
        return new new_universal_poincare_map<rhs_type, N, section_type>(*this);
    }

private:
    void map_RKDP45(const vec_type& seed, std::vector<vec_type>& steps,
                    std::vector< double >& times, int niter,
                    std::vector<vec_type>& orbit, bool do_orbit=false) const;

    const rhs_type*     _rhs;
    double              _prec;
    const section_type* _sec;
};

template<int N>
class planar_section {
public:
    typedef nvis::fixed_vector<double, N>             vec_type;
    typedef nvis::fixed_vector < double, N + 1 >      plane_type;

    planar_section(const plane_type& plane)
        : _normal(nvis::subv<0, N>(plane)), _offset(plane[N]), _sign(1.) {}

    planar_section(const vec_type& normal, double offset)
        : _normal(normal), _offset(offset), _sign(1.) {}

    void direction(bool dir) const
    {
        _sign = dir ? 1. : -1;
    }

    double dist(const vec_type& x) const
    {
        return _sign*nvis::inner(x, _normal) - _offset;
    }

    const vec_type& ddist(const vec_type&) const
    {
        return _normal;
    }

    planar_section* clone() const
    {
        return new planar_section(_normal, _offset);
    }

private:
    vec_type _normal;
    double _offset;
    mutable double _sign;
};

namespace {

template<typename Section>
class intersect_stop {
    typedef Section                         section_type;
    typedef typename section_type::vec_type vec_type;

    const section_type* _sec;
    bool                _stopped;
    vec_type            _inter;
    double              _inter_time;

public:

    intersect_stop(const section_type* sec) :
        _sec(sec), _stopped(false)
    {
    }

    static double secant_method(double v0, double t0, double v1, double t1)
    {
        return (v1*t0 - v0*t1)/(v1-v0);
    }

    template<typename STEP>
    bool operator()(const STEP& s)
    {
        _stopped = false;

        double t0 = std::min(s.t0(), s.t1());
        double t1 = std::max(s.t0(), s.t1());

        double f0 = _sec->dist(s.y(t0));
        double f1 = _sec->dist(s.y(t1));

		double ff0 = _sec->dist(s.y(s.t0()));
		double ff1 = _sec->dist(s.y(s.t1()));
		int sign_ff0 = ff0 > 0.0 ? 1 : -1;
		int sign_ff1 = ff1 > 0.0 ? 1 : -1;

        // std::cerr << "y(" << t0 << ")=" << f0 << " -- y("
        //           << t1 << ")=" << f1 << std::endl;

        if (f0 < 0 && f1 > 0) {
            // std::cerr << "intersection detected...\n";

            double dt, t, df;
            df = f1 - f0;
            t = secant_method(f0, t0, f1, t1);
            double oldf, f = _sec->dist(s.y(t));
            // std::cerr << "y(" << t << ")=" << f << std::endl;

            // apply secant method
            while (fabs(f) > 1e-8) {
                // std::cerr << "iterating..." << std::endl;
                if (f > 0) {
                    t1 = t;
                    f1 = f;
                } else {
                    t0 = t;
                    f0 = f;
                }
                t = secant_method(f0, t0, f1, t1);
                oldf = f;
                f = _sec->dist(s.y(t));
                // std::cerr << "y(" << t << ")=" << f << std::endl;
                if (oldf <= f) {
                    // std::cerr << "endless loop detected.\n";
                    break;
                }
            }
            // std::cerr << "intersection found.\n";

			int sign_f = f > 0.0 ? 1 : -1;
			if (sign_f == sign_ff1) {
				_inter = s.y(t);
				_inter_time = t;
			}
			else {
				int sign_f0 = f0 > 0.0 ? 1 : -1;
				int sign_f1 = f1 > 0.0 ? 1 : -1;
				if (sign_f0 == sign_ff1) {
					_inter = s.y(t0);
					_inter_time = t0;
				}
				else {
					_inter = s.y(t1);
					_inter_time = t1;
				}
			}
            _stopped = true;
        }

        return _stopped;
    }

    bool did_stop() const
    {
        return _stopped;
    }

    const vec_type& where() const
    {
        return _inter;
    }

    double when() const
    {
        return _inter_time;
    }
};
} // anonymous namespace

// --------------------------------------------------------------------------

template< typename RHS, int N, typename Section>
void new_universal_poincare_map<RHS, N, Section>::
map_RKDP45(const vec_type& seed, std::vector<vec_type>& steps,
           std::vector< double >& times, int niter,
           std::vector<vec_type>& orbit, bool do_orbit) const
{
    steps.clear();
    steps.reserve(std::abs(niter));
    times.clear();
    times.reserve(std::abs(niter));

    if (do_orbit) orbit.clear();

    nvis::dopri5<vec_type> intg;

    intg.t = 0;
    intg.y = seed;
    if (do_orbit) orbit.push_back(seed);

    intg.t_max = niter < 0 ? -1e8 : 1e8;
    intg.h = 1e-5;
    intg.abstol = _prec;
    intg.reltol = _prec;

    typename nvis::dopri5<vec_type>::result res;
    typename nvis::dopri5<vec_type>::step   step;

    const section_type* section = _sec->clone();
    section->direction(niter < 0 ? false : true);
    intersect_stop<section_type> stop(section);

    while (steps.size() < std::abs(niter)) {
        try {
            res = intg.do_step(*_rhs, step);
            if (do_orbit) orbit.push_back(step.y1());
        } catch (...) {
            std::cout << "caught exception in RKDP45" << std::endl;
            break;
        }

        if (res != nvis::dopri5<vec_type>::OK) {
            break;
        }

        stop(step);

        if (stop.did_stop()) {
            steps.push_back(stop.where());
            times.push_back(stop.when());
        }
    }

    delete section;
}




#endif
