#ifndef __XAVIER_CHAINS_HPP__
#define __XAVIER_CHAINS_HPP__

#include <math/fixed_vector.hpp>
#include <vector>

#include <poincare/fixpoints.hpp>
#include <poincare/metric.hpp>
#include <misc/sort.hpp>

namespace xavier {
double chain_distance(const std::vector< xavier::fixpoint >& c1,
                      const std::vector< xavier::fixpoint>& c2)
{
    unsigned int p = c1.size();
    
    if (p == 1) {
        return nvis::norm(c1[0].pos - c2[0].pos);
    }
    
    // find best initial match
    nvis::vec2 x0 = c1[0].pos;
    std::vector< double > dist(p, 0);
    for (unsigned int i = 0 ; i < p ; ++i) {
        dist[i] = xavier::__default_metric.distance(x0, c2[i].pos);
    }
    std::vector<unsigned int > sorted;
    xavier::sort(dist, sorted);
    
    unsigned int i0 = sorted[0];
    double d = 0;
    for (unsigned int i = 0 ; i < p ; ++i) {
        d += xavier::__default_metric.distance(c1[i].pos, c2[(i+i0)%p].pos);
    }
    return d / (double)p;
}

void filter_duplicates(std::vector< std::vector< xavier::fixpoint > >& all_p_chains,
                       double mindist)
{
    if (all_p_chains.size() < 2) {
        return;
    }
    
    std::vector< bool > removed(all_p_chains.size(), false);
    
    for (unsigned int i = 0 ; i < all_p_chains.size() ; ++i) {
        const std::vector< xavier::fixpoint >& c1 = all_p_chains[i];
        for (unsigned int j = i + 1 ; j < all_p_chains.size() ; ++j) {
            if (removed[j]) {
                continue;
            }
            double d = chain_distance(c1, all_p_chains[j]);
            if (d < mindist) {
                removed[j] = true;
            }
        }
    }
    
    std::vector< std::vector< xavier::fixpoint > > filtered;
    for (unsigned int i = 0 ; i < all_p_chains.size() ; ++i) {
        if (!removed[i]) {
            filtered.push_back(std::vector< xavier::fixpoint>());
            std::copy(all_p_chains[i].begin(), all_p_chains[i].end(), std::back_inserter(filtered.back()));
        }
    }
    
    std::cout << "out of " << all_p_chains.size() << " chains " << filtered.size() << " were unique" << std::endl;
    std::swap(all_p_chains, filtered);
}
}


#endif



