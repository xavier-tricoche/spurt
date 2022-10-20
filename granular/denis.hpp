#ifndef __GRANULAR_DENIS_HPP__
#define __GRANULAR_DENIS_HPP__

#include <math/fixed_vector.hpp>
#include <poincare/metric.hpp>

namespace spurt { 
namespace denis {
    
// gamma = 0.5, 1, 1.5

class discrete_holmes {
    nvis::vec2 do_step(const nvis::vec2& x, bool fwd) const {
        const double& theta_0 = x[0];
        const double& v_0     = x[1];
        nvis::vec2 next;
        double& theta_1 = next[0];
        double& v_1     = next[1];
        
        if (fwd) {
            theta_1 = theta_0 + v_0;
            v_1     = _rho*v_0 + _gamma*W(theta_1);
        }
        else {
            v_1     = (v_0 - _gamma*W(theta_0))/_rho;
            theta_1 = theta_0 - v_1;
        }
        
        return next;
    }
    
public:
    // constant of gravity
    constexpr static double g = 9.80665;
    
    discrete_holmes(double T, double a, double rho, double gamma, size_t N)
        : _T(T), _a(a), _rho(rho), _gamma(gamma), _N(N) {
        // set value of dependent parameters
        _g_star = g/(double)_N;
        _omega  = sqrt(_gamma*_g_star/(2.*_a*(1.+_rho)));
        // set metric
        _metric.bounds().min()[0] = 0;
        _metric.bounds().max()[0] = _omega*_T;
        _metric.periodic(0) = true;
        _metric.periodic(1) = false;
        _period = _omega*_T;
    }
    
    const double& period() const     { return _period;   }
    const double& amplitude() const  { return _a;        }
    const double& gamma() const      { return _gamma;    }
    const double& rho() const        { return _rho;      }
    const double& omega() const      { return _omega;    }
    const map_metric& metric() const { return _metric;   }
    
    const size_t& number_of_particles() const { return _N; }
    
    double W(double t) const
    {   
        double s = _metric.modulo(t, 0);
        if (s <= M_PI) {
            return cos(s);
        } else {
            return 0.;
        }
    }
    
    nvis::vec2 operator()(const nvis::vec2& x, bool fwd=true, bool modulo=true) const 
    {   
        nvis::vec2 next = do_step(x, fwd);
        if (modulo) return _metric.modulo(next);
        else return next;
    }
    
    std::pair<nvis::vec2, nvis::vec2> step(const nvis::vec2& x, bool fwd=true) const
    {
        nvis::vec2 next = do_step(x, fwd);
        std::pair<nvis::vec2, nvis::vec2> r;
        r.first = _metric.displacement(x, next);
        r.second = _metric.modulo(next);
        return r;
    }
    
private:
    // Equation parameters
    double     _omega, _T, _gamma, _rho, _g_star, _a, _period;
    size_t     _N;
    map_metric _metric;
};
    

} // denis
} // spurt


#endif // __GRANULAR_DENIS_HPP__