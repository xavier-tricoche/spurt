#ifndef __XAVIER__TOPOLOGY__HPP__
#define __XAVIER__TOPOLOGY__HPP__

#include <set>

#include <poincare/map.hpp>
#include <poincare/newton.hpp>
#include <poincare/macros.hpp>
#include <poincare/fixpoints.hpp>
#include <poincare/chains.hpp>
#include <tokamak/map2d.hpp>

namespace xavier {
struct MarkedPos {
    MarkedPos(double eps) : _eps(eps) {}
    bool known(nvis::vec2& x, unsigned int p) const {
        for (unsigned int i = 0 ; i < _pos.size() ; ++i) {
            double d = xavier::__default_metric.distance(x, _pos[i]);
            if (d < _eps && ((_per[i] % p == 0) || (p % _per[i] == 0))) {
                return true;
            }
        }
        return false;
    }
    
    void add(const nvis::vec2& x, unsigned int p) {
        _pos.push_back(x);
        _per.push_back(p);
    }
    
    std::vector< nvis::vec2 > _pos;
    std::vector< unsigned int > _per;
    double _eps;
};

template< typename Map, typename Jacobian >
bool compute_iterates(const Map& map, const Jacobian& jacobian,
                      unsigned int period, const nvis::vec2& x,
                      std::vector< nvis::vec2 >& chain,
                      double eps, unsigned int niter = 20)
{
    try {
        chain.clear();
        nvis::vec2 y(x);
        nvis::vec2 f;
        for (unsigned int i = 0 ; i < period ; ++i) {
            bool found = xavier::newton(map, jacobian, y, f, eps, niter);
            if (!found) {
                WARNING_MACRO(0, "unable to converge from iterated position" << std::endl);
                return false;
            }
            if (map_debug::verbose_level > 1) {
                WARNING_MACRO(0, "after Newton correction, iterate #"
                              << i + 1 << " of " << period
                              << " is now associated with rhs norm "
                              << nvis::norm(f)
                              << std::endl);
            }
            chain.push_back(y);
            y = map.single_step(y);
        }
        
        return true;
    } catch (...) {
        WARNING_MACRO(0, "exception caught in compute_iterates" << std::endl);
        return false;
    }
}

template<typename Map>
void find_fixed_points(const Map& map_coarse, const Map& map_fine,
                       const MarkedPos& all_marked_positions,
                       const std::vector< nvis::vec2 >& seeds,
                       unsigned int p, unsigned int pmax, double maxnorm,
                       double dx0, double dy0, double eps0, // coarse Newton
                       double dx1, double dy1, double eps1, // fine Newton
                       std::vector< std::vector< xavier::fixpoint > >& chains)
{
    chains.clear();
    // we keep a local list of marked positions to prevent conflicts
    // in case the extraction of fixed points of several periods is run in
    // parallel.
    MarkedPos loc_marked_pos(all_marked_positions._eps);
    std::cout << "processing " << seeds.size() << " candidates for period " << p << '\n';
    
    // wrappers needed for Newton search
    // coarse precision
    map_from_last_wrapper<Map> pmap_coarse(map_coarse, p);
    // least_sq_jacobian< map_from_last_wrapper<Map> > jacobian_coarse(pmap_coarse, dx0, dy0);
    integral_jacobian< map_from_last_wrapper<Map> > jacobian_coarse(pmap_coarse);
    
    // fine precision
    map_from_last_wrapper<Map> pmap_fine(map_fine, p);
    // least_sq_jacobian< map_from_last_wrapper<Map> > jacobian_fine(pmap_fine, dx1, dy1);
    integral_jacobian< map_from_last_wrapper<Map> > jacobian_fine(pmap_fine);
    
    // temporary storage of fixed points prior to uniqueness filtering
    std::set< xavier::fixpoint > fixpoints;
    nvis::vec2 f;
    
    for (unsigned int n = 0 ; n < seeds.size() ; ++n) {
        nvis::vec2 x(seeds[n]);
        
        // check for redundancy
        if (loc_marked_pos.known(x, p) || all_marked_positions.known(x, p)) {
            WARNING_MACRO(0, "seed point #" << n + 1
                          << " lies in direct vicinity of previously discovered location\n");
        }
        
        // check for maximum norm violation
        try {
            f = pmap_coarse(x);
        } catch (...) {
            WARNING_MACRO(0, "could not compute map at " << x << '\n');
            continue;
        }
        double normf = nvis::norm(f);
        if (normf > maxnorm) {
            WARNING_MACRO(0, "candidate #" << n + 1 << " of " << seeds.size()
                          << " associated with norm " << normf
                          << " which is larger than threshold (=" << maxnorm << ")\n");
            break; // list of seeds assumed to be sorted by norm
        }
        
        // coarse numerical search
        WARNING_MACRO(0, "\n**candidate #" << n + 1 << " / " << seeds.size() << " located at " << x
                      << " has period = " << p << " and value " << normf
                      << std::endl);
        if (map_debug::verbose_level > 0) {
            WARNING_MACRO(1, "RHS norm BEFORE Newton is " << normf << std::endl);
        }
        bool found;
        if (normf < eps0) {
            found = true;
        } else {
            found = xavier::newton(pmap_coarse, jacobian_coarse, x, f, eps0, 15);
        }
        if (!found) {
            WARNING_MACRO(0, "Coarse Newton search did not converge at " << x << '\n');
            continue;
        }
        WARNING_MACRO(0, "Coarse Newton search converged at " << x << std::endl);
        normf = nvis::norm(f);
        WARNING_MACRO(1, "RHS norm at " << x << " AFTER coarse Newton is " << normf << std::endl);
        
        // intermediate redundancy check
        if (loc_marked_pos.known(x, p) || all_marked_positions.known(x, p)) {
            WARNING_MACRO(0, "point " << x << " found after coarse search for period "
                          << p << "is already known. skipping...\n");
            continue;
        }
        
        // verify that estimated period remains valid at newly reached location
        unsigned int q = xavier::period(map_coarse, x, pmax, eps0);
        if (q == 0) {
            WARNING_MACRO(0, "unable to determine actual period at " << x << '\n');
            continue;
        }
        if (q != p) {
            WARNING_MACRO(0, "estimated period was contradicted after coarse Newton search\n");
            continue;
        }
        
        // fine Newton search
        // did we already converge by any chance...?
        if (normf < eps1) {
            found = true;
        } else {
            found = xavier::newton(pmap_fine, jacobian_fine, x, f, eps1, 50);
        }
        if (!found) {
            WARNING_MACRO(0, "WARNING: initially found fixpoint was lost after period correction!\n");
            continue;
        } else if (loc_marked_pos.known(x, p) || all_marked_positions.known(x, p)) {
            WARNING_MACRO(0, "fixpoint " << x << " for period " << p
                          << "found after refinement is already known. skipping...\n");
            continue;
        } else {
            WARNING_MACRO(0, "position " << x << " passed uniqueness test\n");
            WARNING_MACRO(0, "final norm of rhs at " << x << " after refined Newton step is "
                          << nvis::norm(f)
                          << std::endl);
        }
        
        // chain extraction and linear analysis
        std::vector< nvis::vec2 > iterates;
        WARNING_MACRO(0, "computing iterates from " << x << '\n');
        chains.push_back(std::vector< xavier::fixpoint >());
        std::vector< xavier::fixpoint >& chain = chains.back();
        if (!xavier::compute_iterates(pmap_fine, jacobian_fine, p, x, iterates, eps1, 50)) {
            WARNING_MACRO(0, "unable to recover complete " << p << "-chain. skipping\n");
            chains.pop_back();
            continue;
        }
        map_wrapper<Map> one_map(map_fine, 1);
        integral_jacobian< map_wrapper<Map> > jacobian_one_fine(one_map);
        if (!xavier::linear_chain_analysis(jacobian_one_fine, iterates, chain)) {
            chain.clear();
            WARNING_MACRO(0, "unable to determine consistent linear type\n");
            chains.pop_back();
            continue;
        }
        WARNING_MACRO(0, "Numerical search was successful.\n");
        
        central_diff_jacobian< map_wrapper<Map> > cdj(one_map, 0.05, 0.05);
        for (unsigned int n = 0 ; n < iterates.size() ; ++n) {
            nvis::mat2 jint = jacobian_one_fine(iterates[n]);
            nvis::mat2 japp = cdj(iterates[n]);
            std::cout << "at " << iterates[n] << ":\n"
                      << "integrated Jacobian   = " << jint << '\n'
                      << "approximated Jacobian = " << japp << '\n';
        }
        
        // add found positions to the list
        for (unsigned int i = 0 ; i < chain.size() ; ++i) {
            loc_marked_pos.add(chain[i].pos, p);
        }
    }
    
    WARNING_MACRO(0, "\n\nfiltering chains of period " << p << std::endl);
    filter_duplicates(chains, maxnorm);
    WARNING_MACRO(0, "after filtering there are " << chains.size()
                  << " valid chains of period " << p << '\n');
}


template<typename Map>
void find_fixed_points(const Map& map,
                       const MarkedPos& all_marked_positions,
                       const std::vector< nvis::vec2 >& seeds,
                       unsigned int p, unsigned int pmax, double maxnorm,
                       double eps,
                       std::vector< std::vector< xavier::fixpoint > >& chains,
                       unsigned int verbose_level = 0)
{
    xavier::map_debug::verbose_level = verbose_level;
    
    chains.clear();
    // we keep a local list of marked positions to prevent conflicts
    // in case the extraction of fixed points of several periods is run in
    // parallel.
    MarkedPos loc_marked_pos(all_marked_positions._eps);
    std::cout << "processing " << seeds.size() << " candidates for period " << p << '\n';
    
    // wrappers needed for Newton search
    map_from_last_wrapper<Map> pmap(map, p);
    integral_jacobian< map_from_last_wrapper<Map> > jacobian(pmap);
    
    // temporary storage of fixed points prior to uniqueness filtering
    std::set< xavier::fixpoint > fixpoints;
    nvis::vec2 f;
    
    for (unsigned int n = 0 ; n < seeds.size() ; ++n) {
        nvis::vec2 x(seeds[n]);
        
        // check for redundancy
        if (loc_marked_pos.known(x, p) || all_marked_positions.known(x, p)) {
            WARNING_MACRO(0, "seed point #" << n + 1
                          << " lies in direct vicinity of previously discovered location\n");
        }
        
        // check for maximum norm violation
        try {
            f = pmap(x);
        } catch (...) {
            WARNING_MACRO(0, "could not compute map at " << x << '\n');
            continue;
        }
        double normf = nvis::norm(f);
        if (normf > maxnorm) {
            WARNING_MACRO(0, "candidate #" << n + 1 << "/" << seeds.size() << " at " << x
                          << " associated with norm " << normf
                          << " which is larger than threshold (=" << maxnorm << "). skipping\n");
            continue;
            // break; // list of seeds assumed to be sorted by norm
        }
        
        // compute average error estimate in the neighborhood of this point
        {
            nvis::vec2 avg_ee(0,0);
            double avg_norm = 0;
            bool failed = false;
            nvis::vec2 __x;
            map2d::value_type __v;
            for (int i = 0 ; i < 10 ; ++i) {
                __x[0] = x[0] + normf * (-1. + 2.*drand48());
                __x[1] = x[1] + normf * (-1. + 2.*drand48());
                try {
                    __v = map.map_complete(__x, p);
                } catch (...) {
                    WARNING_MACRO(0, "could not compute map at " << x << '\n');
                    failed = true;
                    break;
                }
                
                avg_norm += nvis::norm(__v.err);
                avg_ee += __v.err;
            }
            if (failed) {
                continue;
            }
            
            avg_norm /= 10.;
            avg_ee /= 10.;
            std::cout << "average error estimate is " << avg_norm << " (" << avg_ee << "), requested accuracy is " << eps << std::endl;
            eps = avg_norm;
        }
        
        // coarse numerical search
        
        std::cerr << "blahblahblah" << std::endl;
        
        
        WARNING_MACRO(0, "\n**candidate #" << n + 1 << " / " << seeds.size() << " located at " << x
                      << " has period = " << p << " and value " << normf
                      << std::endl);
        if (map_debug::verbose_level > 0) {
            WARNING_MACRO(1, "RHS norm BEFORE Newton is " << normf << std::endl);
        }
        bool found;
        if (normf < eps) {
            found = true;
        } else {
            found = xavier::newton(pmap, jacobian, x, f, eps, 50);
        }
        if (!found) {
            WARNING_MACRO(0, "Newton search did not converge at " << x << '\n');
            continue;
        }
        WARNING_MACRO(0, "Newton search converged at " << x << std::endl);
        normf = nvis::norm(f);
        WARNING_MACRO(1, "RHS norm at " << x << " AFTER Newton is " << normf << std::endl);
        
        // intermediate redundancy check
        if (loc_marked_pos.known(x, p) || all_marked_positions.known(x, p)) {
            WARNING_MACRO(0, "point " << x << " found after coarse search for period "
                          << p << "is already known. skipping...\n");
            continue;
        }
        
        // verify that estimated period remains valid at newly reached location
        unsigned int q = xavier::period(map, x, pmax, eps);
        if (q == 0) {
            WARNING_MACRO(0, "unable to determine actual period at " << x << '\n');
            continue;
        }
        if (q != p) {
            WARNING_MACRO(0, "estimated period was contradicted after coarse Newton search\n");
            continue;
        }
        
        // chain extraction and linear analysis
        std::vector< nvis::vec2 > iterates;
        WARNING_MACRO(0, "computing iterates from " << x << '\n');
        chains.push_back(std::vector< xavier::fixpoint >());
        std::vector< xavier::fixpoint >& chain = chains.back();
        if (!xavier::compute_iterates(pmap, jacobian, p, x, iterates, eps, 50)) {
            WARNING_MACRO(0, "unable to recover complete " << p << "-chain. skipping\n");
            chains.pop_back();
            continue;
        }
        map_wrapper<Map> one_map(map, 1);
        integral_jacobian< map_wrapper<Map> > jacobian_one(one_map);
        if (!xavier::linear_chain_analysis(jacobian_one, iterates, chain)) {
            chain.clear();
            WARNING_MACRO(0, "unable to determine consistent linear type\n");
            chains.pop_back();
            continue;
        }
        WARNING_MACRO(0, "Numerical search was successful.\n");
        
        // for debugging purpose only...
        central_diff_jacobian< map_wrapper<Map> > cdj(one_map, 0.05, 0.05);
        for (unsigned int n = 0 ; n < iterates.size() ; ++n) {
            nvis::mat2 jint = jacobian_one(iterates[n]);
            nvis::mat2 japp = cdj(iterates[n]);
            std::cout << "at " << iterates[n] << ":\n"
                      << "integrated Jacobian   = " << jint << '\n'
                      << "approximated Jacobian = " << japp << '\n';
        }
        // ...
        
        // add found positions to the list
        for (unsigned int i = 0 ; i < chain.size() ; ++i) {
            loc_marked_pos.add(chain[i].pos, p);
        }
    }
    
    WARNING_MACRO(0, "\n\nfiltering chains of period " << p << std::endl);
    filter_duplicates(chains, maxnorm);
    WARNING_MACRO(0, "after filtering there are " << chains.size()
                  << " valid chains of period " << p << '\n');
}

}





#endif













































































