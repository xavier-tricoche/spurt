#include "universal_poincare_map.hpp"

#include <math/dopri5.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

using namespace nvis;

// --------------------------------------------------------------------------

namespace {

template< int N >
class intersect_stop {
    typedef nvis::fixed_vector<double, N>   vecN;

    vecN        normal;
    double      offset;

    bool        stopped;
    vecN        intersection;

public:

    intersect_stop(const nvis::fixed_vector < double, N + 1 > & plane) :
        normal(subv<0, N>(plane)), offset(plane[N]), stopped(false) {
    }

    template<typename STEP>
    bool operator()(const STEP& s) {
        stopped = false;

        double tmin = std::min(s.t0(), s.t1());
        double tmax = std::max(s.t0(), s.t1());

        double f0 = inner(normal, s.y0()) - offset;
        double f1 = inner(normal, s.y1()) - offset;

        if (f0 < 0 && f1 > 0) {
            double h, t, df;

            // initial guess
            t  = 0.5 * (s.t0() + s.t1());
            df = -inner(normal, s.y(t)) + offset;

            // Newton iteration
            while (fabs(df) > 1e-8) {
                h = df / (inner(normal, s.dy(t)) - offset);
                t = std::min(tmax, std::max(tmin, t + h));

                df = -inner(normal, s.y(t)) + offset;
            }

            intersection = s.y(t);
            stopped = true;
        }

        return stopped;
    }

    bool did_stop() {
        return stopped;
    }

    const vecN& where() {
        return intersection;
    }
};

} // anonymous namespace

// --------------------------------------------------------------------------

template< typename RHS, int N >
void universal_poincare_map<RHS, N>::map_RK56(const vecN& in, std::vector<vecN>& out, int niter) const
{
    out.clear();
    out.reserve(std::abs(niter));

    nvis::dopri5<vecN> intg;

    intg.t = 0;
    intg.y = in;

    intg.t_max = niter < 0 ? -1e8 : 1e8;
    intg.h = 1e-5;
    intg.abstol = _prec;
    intg.reltol = _prec;

    typename nvis::dopri5<vecN>::result res;
    typename nvis::dopri5<vecN>::step   step;

    nvis::fixed_vector < double, N + 1 > plane = niter < 0 ? -1.0 * _rhs->plane() : _rhs->plane();

    intersect_stop<N> stop(plane);

    while (out.size() < std::abs(niter)) {
        try {
            res = intg.do_step(*_rhs, step);
        } catch (...) {
            break;
        }

        if (res != dopri5<vecN>::OK) {
            break;
        }

        stop(step);

        if (stop.did_stop()) {
            out.push_back(_rhs->project(stop.where()));
            std::cout << "Iteration #" << out.size() << " t=" << std::max(step.t0(), step.t1()) << '\n';
        }
    }
}
