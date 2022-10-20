#ifndef __FIND_FIXED_POINTS_HPP__
#define __FIND_FIXED_POINTS_HPP__

#include <sstream>

#include "definitions.hpp"
#include <poincare/newton.hpp>
#include <poincare/topology.hpp>
#include <poincare/fixpoints.hpp>

namespace map_analysis {
template<typename Map, typename Jac>
inline bool robust_newton(const Map& rhs, const Jac& jacobian, nvis::vec2& x,
                          double tol, double epsilon)
{
    nvis::vec2 f = rhs(x);
    double norm = nvis::norm(f);
    
    double eps = tol;
    bool ok = false;
    double ratio;
    
    while (true) {
        nvis::vec2 y(x);
        bool found;
        try {
            found = spurt::newton(rhs, jacobian, y, f, eps, 20);
        } catch (...) {
            return false;
        }
        
        double new_norm = nvis::norm(f);
        
        if (new_norm < epsilon) {
            x = y;
            ok = true;
            break;
        }
        
        nvis::vec2 err = rhs.error(y);
        
        // evaluate ration between error and rhs
        ratio = nvis::norm(err) / new_norm;
        
        if (ratio >= 0.2 || (!found && rhs.get_precision() > 1.0e-10)) {
            rhs.set_precision(0.2*rhs.get_precision());
        }
        if (found) {
            ok = true;
            x = y;
            norm = new_norm;
            eps /= 5.;
        } else if (ratio > 0.2 || rhs.get_precision() < 1.0e-10) {
            // error is fine but Newton does not progress
            break;
        }
    }
    
    return ok;
}

template<typename Map, typename Jac>
inline bool find_fixed_points(std::list<spurt::fixpoint>& fps,
                              const Map& rhs, const Jac& jacobian, const nvis::vec2& seed,
                              double tol, double epsilon = 1.0e-6)
{
    fps.clear();
    nvis::vec2 x(seed);
    bool found = robust_newton(rhs, jacobian, x, tol, epsilon);
    if (!found) {
        return false;
    }
    
    nvis::vec2 f = rhs(x);
    std::ostringstream os;
    os << "fixed point found at " << x << " with norm "
       << std::setprecision(12) << nvis::norm(f) << '\n';
    std::cerr << os.str();
    
    double h = rhs.get_precision();
    if (h > 1.0e-10) {
        rhs.set_precision(0.1*h);
    }
    
    
    for (int i = 0 ; i < rhs.period() ; ++i) {
    
        spurt::fixpoint fpt;
        spurt::linear_analysis(jacobian, rhs.period(), x, fpt);
        fps.push_back(fpt);
        
        // move on to next one on the chain
        x = rhs.one_step(x);
    }
    
    return true;
}

inline double chain_distance(const std::list<spurt::fixpoint>& chain1,
                             const std::list<spurt::fixpoint>& chain2)
{
    typedef std::list<spurt::fixpoint>         chain_type;
    typedef chain_type::const_iterator          fp_iter_type;
    
    // quadratic complexity method. problem is tiny anyway
    double max = 0;
    for (fp_iter_type it1 = chain1.begin() ; it1 != chain1.end() ; ++it1) {
        std::list<double> dist;
        for (fp_iter_type it2 = chain2.begin() ; it2 != chain2.end() ; ++it2) {
            dist.push_back(static_data::metric.distance(it1->pos, it2->pos));
        }
        double d = *std::min_element(dist.begin(), dist.end());
        max = std::max(max, d);
    }
    
    return max;
}

template<typename Map>
inline double chain_norm(const std::list<spurt::fixpoint>& chain,
                         const Map& map)
{
    double l1norm = 0;
    for (std::list<spurt::fixpoint>::const_iterator it = chain.begin() ;
            it != chain.end() ; ++it) {
        l1norm += nvis::norm(map(it->pos));
    }
    
    return l1norm;
}

template<typename Map>
inline void uniquify_chains(std::list<std::list<spurt::fixpoint> >& all_chains,
                            const Map& map, double tolerance)
{
    typedef std::list<spurt::fixpoint>         chain_type;
    typedef chain_type::iterator                fp_iter_type;
    typedef std::list<chain_type>::iterator     chain_iter_type;
    
    std::cerr << "uniquifying " << all_chains.size() << " chains\n";
    
    // start by sorting the chains
    std::vector<chain_type> saddle_chains, center_chains, mixed_chains;
    for (chain_iter_type it_ch = all_chains.begin() ; it_ch != all_chains.end() ; ++it_ch) {
        if (it_ch->front().saddle) {
            bool all_saddle = true;
            for (fp_iter_type it_fp = it_ch->begin() ; it_fp != it_ch->end() ; ++it_fp) {
                if (!it_fp->saddle) {
                    all_saddle = false;
                    break;
                }
            }
            if (all_saddle) {
                saddle_chains.push_back(*it_ch);
            } else {
                mixed_chains.push_back(*it_ch);
            }
        } else {
            bool all_centers = true;
            for (fp_iter_type it_fp = it_ch->begin() ; it_fp != it_ch->end() ; ++it_fp) {
                if (it_fp->saddle) {
                    all_centers = false;
                    break;
                }
            }
            if (all_centers) {
                center_chains.push_back(*it_ch);
            } else {
                mixed_chains.push_back(*it_ch);
            }
        }
    }
    
    std::cerr << "we have identified " << saddle_chains.size() << " saddle chains, "
              << center_chains.size() << " center chains, and " << mixed_chains.size()
              << " mixed / broken chains\n";
              
    std::cerr << "tolerance is " << tolerance << '\n';
    
    std::vector<bool> discarded_saddles(saddle_chains.size(), false);
    std::vector<chain_type> unique_saddle_chains, unique_center_chains;
    for (unsigned int i = 0 ; i < saddle_chains.size() ; ++i) {
        if (discarded_saddles[i]) {
            continue;
        }
        double ref_norm = chain_norm(saddle_chains[i], map);
        bool discarded = false;
        for (int j = i + 1 ; j < saddle_chains.size() ; ++j) {
            double d = chain_distance(saddle_chains[i], saddle_chains[j]);
            if (d < tolerance) {
                if (chain_norm(saddle_chains[j], map) < ref_norm) {
                    discarded_saddles[i] = true;
                    break;
                } else {
                    discarded_saddles[j] = true;
                }
            }
        }
        if (!discarded_saddles[i]) {
            unique_saddle_chains.push_back(saddle_chains[i]);
        }
    }
    std::vector<bool> discarded_centers(center_chains.size(), false);
    for (unsigned int i = 0 ; i < center_chains.size() ; ++i) {
        if (discarded_centers[i]) {
            continue;
        }
        double ref_norm = chain_norm(center_chains[i], map);
        bool discarded = false;
        for (int j = i + 1 ; j < center_chains.size() ; ++j) {
            double d = chain_distance(center_chains[i], center_chains[j]);
            if (d < tolerance) {
                if (chain_norm(center_chains[j], map) < ref_norm) {
                    discarded_centers[i] = true;
                    break;
                } else {
                    discarded_centers[j] = true;
                }
            }
        }
        if (!discarded_centers[i]) {
            unique_center_chains.push_back(center_chains[i]);
        }
    }
    
    all_chains.clear();
    std::copy(unique_saddle_chains.begin(), unique_saddle_chains.end(),
              std::back_inserter(all_chains));
    std::copy(unique_center_chains.begin(), unique_center_chains.end(),
              std::back_inserter(all_chains));
              
    std::cerr << "after uniquification we have " << unique_saddle_chains.size()
              << " unique saddle chains and " << unique_center_chains.size()
              << " unique center chains\n";
}

}

#endif


















































