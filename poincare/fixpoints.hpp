#ifndef __XAVIER_FIXPOINT_HPP__
#define __XAVIER_FIXPOINT_HPP__

#include <poincare/basic_math.hpp>
#include <math/fixed_vector.hpp>
#include "poincare/macros.hpp"
#include <iostream>

namespace spurt {

struct fixpoint {
    fixpoint() : isolated(true) {}
    nvis::vec2 pos;
    bool saddle;
    nvis::mat2 J;
    nvis::vec2 evec[2];
    unsigned int K;
    bool isolated;
};

std::ostream& operator<<(std::ostream& os, const fixpoint& fp)
{
    os << "fixed point: \tpos=" << fp.pos
       << ", \tperiod=" << fp.K
       << ", \t" << (fp.saddle ? "saddle" : "center")
       << "\nJ=" << fp.J;
    if (fp.saddle) {
        os << ", \tev0=" << fp.evec[0] << ",\t ev1=" << fp.evec[1];
    }
    
    return os;
}

template< typename Jacobian >
void linear_analysis(const Jacobian& jac, unsigned int period,
                     const nvis::vec2& x, spurt::fixpoint& fp)
{
    nvis::mat2 J = jac(x);
    
    // std::cerr << "measured Jacobian at " << x << " is " << J << '\n';
    
    fp.pos = x;
    fp.K = period;
    fp.J = J;
    fp.saddle = eigen(fp.evec, fp.J);
    if (fp.saddle) {
        std::cout << "found a saddle point at " << fp.pos << "\n";
    } else {
        std::cout << "found a center at " << fp.pos << "\n";
    }
}

template< typename Jacobian >
bool linear_chain_analysis(const Jacobian& one_jacobian,
                           const std::vector< nvis::vec2 >& pos,
                           std::vector< spurt::fixpoint >& fps)
{
    unsigned int period = pos.size();
    fps.resize(period);
    std::vector< nvis::mat2 > J(period);
    
    try {
        for (int i = 0 ; i < period ; ++i) {
            J[i] = one_jacobian(pos[i]);
        }
    } catch (...) {
        WARNING_MACRO(0, "exception caught in linear_chain_analysis at "
                      << pos[0] << '\n');
        return false;
    }
    
    unsigned int nb_saddles = 0;
    for (unsigned int i = 0 ; i < period ; ++i) {
        nvis::mat2 J_p = nvis::mat2::identity();
        for (unsigned int j = 0 ; j < period ; ++j) {
            J_p = J[(i+j)%period] * J_p;
        }
        
        fps[i].pos = pos[i];
        fps[i].K = period;
        fps[i].J = J_p;
        fps[i].saddle = eigen(fps[i].evec, J_p);
        if (fps[i].saddle) {
            ++nb_saddles;
        }
        if (map_debug::verbose_level > 0) {
            if (fps[i].saddle) {
                std::cout << "found a saddle point at " << fps[i].pos
                          << ", jacobian = " << J_p << "\n";
            } else {
                std::cout << "found a center at " << fps[i].pos
                          << ", jacobian = " << J_p << "\n";
            }
        }
    }
    
    if (fps.size() > 1 && nb_saddles*(fps.size() - nb_saddles) > 0) {
        WARNING_MACRO(0, "inconsistent results: "
                      << nb_saddles << " saddles and "
                      << fps.size() - nb_saddles << " centers\n");
        return false;
    }
    
    // analysis successful
    return true;
}

}


#endif











