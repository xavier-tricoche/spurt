#ifndef ____NEWTON_HPP
#define ____NEWTON_HPP

#include "math/fixed_vector.hpp"
#include "math/fixed_matrix.hpp"
#include <vector>

#include <poincare/macros.hpp>
#include <poincare/metric.hpp>
#include <poincare/map.hpp>

namespace spurt {

namespace map_debug {
extern std::vector< std::vector< nvis::vec2 > > newton_steps;
extern std::vector< std::vector< nvis::vec2 > > newton_rhs;
extern bool save_newton_steps;
}

inline nvis::mat2 convert(const nvis::vec4& A)
{
    nvis::mat2 _A;
    _A(0, 0) = A[0];
    _A(1, 0) = A[2];
    _A(0, 1) = A[1];
    _A(1, 1) = A[3];
    
    return _A;
}

template<typename Map>
bool lnsearch(const Map& map, nvis::vec2& x, nvis::vec2& f,
              const nvis::vec2& dd)
{
    double lambda = 1.0;
    const double alpha = 1e-4;
    static double maxlength = spurt::__default_metric.diameter() * 0.05;
    
    nvis::vec2 xsave = x, fsave = f, d = norm(dd) > maxlength ? dd * maxlength / norm(dd) : dd;
    
    for (unsigned int i = 0; i < 7; ++i) {
        x = xsave + lambda * d;
        f = map(x);
        if (norm(f) < (1 - alpha*lambda)*norm(fsave)) {
            return true;
        }
        
        lambda *= 0.5;
    }
    
    return false;
}

template<typename Map, typename Jacobian>
bool newton(const Map& map, const Jacobian& jacobian,
            nvis::vec2& x, nvis::vec2& f,
            double eps, size_t maxiter)
{
    nvis::vec2 d; // d is destination (f is rhs)
    nvis::mat2 J;
    
    if (map_debug::save_newton_steps) {
        map_debug::newton_steps.push_back(std::vector< nvis::vec2 >());
        map_debug::newton_rhs.push_back(std::vector< nvis::vec2 >());
    }
    
    unsigned int k;
    try {
        f = map(x);
        double dinit = nvis::norm(f);
        for (k = 0; k < maxiter; k++) {
            if (map_debug::save_newton_steps) {
                map_debug::newton_rhs.back().push_back(f);
                map_debug::newton_steps.back().push_back(x);
            }
            
            WARNING_MACRO(2, "newton: k = " << k << ", norm(f) = " << norm(f) << " (eps=" << eps << ") at " << x << '\n');
            if (norm(f) < eps) {
                break;
            }
            
            // determine local search direction
            J = jacobian(x);
            d = nvis::solve(J, J * x - f);
            
            // do a relaxation linesearch
            // (updates x and f)
            lnsearch(map, x, f, d - x);
        }
        
        double dfin = nvis::norm(map(x));
        if (k == maxiter) {
            WARNING_MACRO(1, "\t\t initial distance = " << dinit
                          << ", final distance = " << dfin << ". failed.\n");
        }
    } catch (...) {
        WARNING_MACRO(0, "\texception caught in Newton. current position is " << x << std::endl);
        return false;
    }
    
    return !(k == maxiter);
}

template< typename Map, typename RHS, typename Jacobian >
bool compute_iterates(const Map& map, const RHS& rhs, const Jacobian& jacobian,
                      unsigned int period,
                      const nvis::vec2& x,
                      std::vector< nvis::vec2 >& chain,
                      double eps, unsigned int niter = 20)
{
    try {
        chain.clear();
        chain.reserve(period);
        nvis::vec2 y = x;
        nvis::vec2 f = rhs(x);
        for (unsigned int i = 0 ; i < period ; ++i) {
            bool found = spurt::newton(rhs, jacobian, y, f, eps, niter);
            if (!found) {
                std::cout << "unable to converge from iterated position" << std::endl;
                return false;
            }
            f = rhs(y);
            std::cout << "after Newton correction, iterate #"
                      << i << " of " << period
                      << " is now associated with rhs norm "
                      << nvis::norm(f)
                      << std::endl;
            chain.push_back(y);
            y = map.map(y, 1);
        }
        
        return true;
    } catch (...) {
        std::cout << "exception caught in compute_iterates" << std::endl;
        return false;
    }
}

}

#endif




































































