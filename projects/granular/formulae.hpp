#ifndef __NJIT_FORMULAE_HPP__
#define __NJIT_FORMULAE_HPP__

#include <math.h>
#include <math/fixed_vector.hpp>

#define NEW_VERSION
namespace njit {

#define     g           9.80665     // gravity
#define     VERBOSE

struct tapping_parameters {
    int         N;      // number of particles in column
    double      d;      // particle diameter (m)
    double      f;      // tap frequency (Hz)
    double      rho;    // material density (kg/m^3)
    double      K1;     // loading stiffness (N/m^2)
    double      e;      // restitution coefficient
    double      a;      // tapping amplitude (m)
    double      T;      // total period (s)
    
    double omega() const {
        return 2.*M_PI*f;
    }
    
    double tau_star() const {
        return 0.005*d; // 1% of min particle radius
    }
    
    double mass() const {
        return rho/6.*M_PI*d*d*d; // rho*4/3*π*r^3
    }
    
    double total_mass() const {
        return N*mass();
    }
};

inline double y0(double t, const tapping_parameters& param)
{
    static const double omg = param.omega();
    double u = fmod(t, param.T);
    double v;
    if (u <= M_PI/omg) {
        v = param.a*sin(omg*u);
    } else {
        v = 0;
    }
    return v;
}

inline double dy0dt(double t, const tapping_parameters& param)
{
    static const double omg = param.omega();
    double u = fmod(t, param.T);
    double v;
    if (u <= M_PI/omg) {
        v = param.a*omg*cos(omg*u);
    } else {
        v = 0;
    }
    return v;
}

inline double sigma(double s, double alpha=1.0e3)
{
#ifdef VERBOSE
    std::cerr << "sigma(" << s << ", " << alpha << ")=" << tanh(alpha*s) << '\n';
#endif
    return tanh(alpha*s);
}

inline double chi(double s, double alpha=1.0e3)
{
#ifdef VERBOSE
    std::cerr << "chi(" << s << ", " << alpha << ")=" << 0.5*(1 + tanh(alpha*s)) << '\n';
#endif
    return 0.5*(1 + tanh(alpha*s));
}

inline double F(double tau, const tapping_parameters& param)
{
    static const double t_ = param.tau_star();
    double v;
    if (tau <= t_) {
        v = 1 + 1/sqrt(tau) - 1/sqrt(t_);
    } else {
        v = 1;
    }
#ifdef VERBOSE
    std::cerr << "F(" << tau << ")=" << v << '\n';
#endif
    return v;
}

inline double F(double tau, double nu, const tapping_parameters& param)
{
    static const double t_ = param.tau_star();
#ifdef VERBOSE
    std::cerr << "F(" << tau << ", " << nu << ")=" << 1.+chi(-nu)*chi(t_-tau)*(1.+1/sqrt(tau)-1/sqrt(t_)) << '\n';
#endif
    return 1.+chi(-nu)*chi(t_-tau)*(1.+1/sqrt(tau)-1/sqrt(t_));
}

inline double e(const tapping_parameters& param)
{
    const double& rho = param.e; // coefficient of restitution here
    double v = (1.-rho*rho)/(1.+rho*rho);
#ifdef VERBOSE
    std::cerr << "e=" << v << '\n';
#endif
    return v;
}

inline nvis::vec2 ode_BRTU_31(const nvis::vec2& x, double t, const tapping_parameters& param)
{
    static const double M   = param.total_mass(); // total mass
    static const double N   = param.N;            // number of particles
    static const double K   = param.K1;           // stiffness constant
    static const double r1  = param.d/2;          // radius of first particle
    static const double rho = param.e;            // coefficient of restitution
    
    // const double& xi  = x[0];
    const double& eta = x[1];
    const double  v   = (1.-rho*rho)/(1.+rho*rho);
    
    double y    = y0(t, param);
    double ydot = dy0dt(t, param);
    
#ifdef VERBOSE
    std::cerr << "v=" << v << '\n';
    std::cerr << "K/M=" << K/M << '\n';
    std::cerr << "y=" << y << '\n';
    std::cerr << "ydot=" << ydot << '\n';
    std::cerr << "r1=" << r1 << '\n';
#endif
    
    nvis::vec2 f;
    f[0] = eta; // xi dot
    
#ifndef NEW_VERSION
    f[1] = -g + K/M*((1-chi(2.*eta/N - ydot))*(r1-(2.*xi/N-y)) +
                     chi(2.*eta/N-ydot)*(1-chi(2.*xi/N-y-rho*rho*r1))*((2.*xi/N-y)/(rho*rho)-r1))*
           chi(r1-2.*xi/N+y)*F(fabs(2.*xi/N-y), param);
#else
    f[1] = -g + K/M*(1.-v*sigma(2.*eta/N-ydot))*(r1-(2.*eta/N-y))*
           chi(r1-(2.*eta/N-y))*F(fabs(1.*eta/N-y), 2.*eta/N-ydot, param);
#endif
           
#ifdef VERBOSE
    std::cerr << x[0] << "," << x[1] << "," << t << "," << f[0] << "," << f[1] << '\n';
#endif
    
    return f;
}

};



#endif