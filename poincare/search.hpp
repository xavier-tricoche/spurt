#ifndef __SEARCH_HPP__
#define __SEARCH_HPP__

#include <vector>
#include <math/fixed_vector.hpp>
#include <iostream>
#include <poincare/newton.hpp>
#include <poincare/image.hpp>
#include <poincare/fixpoints.hpp>
#include <math/math.hpp>
#include <util/timer.hpp>

namespace search {
extern  std::vector< std::vector< xavier::fixpoint > > all_p_chains;

template< typename T >
void find_extrema(const Image<T>& image, const SubImage& subimage, int radius,
                  std::vector< unsigned int >& extrema, bool max = true)
{
    unsigned int N = image.m * image.n;
    std::vector< T > __copy(&image.data[image.shift], &image.data[image.shift+image.m*image.n]);
    std::vector< unsigned int > __ids(N);
    std::vector< bool > __checked(N, false);
    
    extrema.clear();
    
    nvis::timer __timer;
    
    xavier::sort_ids(__ids, __copy, !max);
    
    for (unsigned int l = 0 ; l < N ; ++l) {
        unsigned int n = __ids[l];
        
        bool skip = false;
        if (__checked[n]) {
            skip = true;
        }
        
        int i = n % image.m;
        int j = n / image.m;
        
        if (!subimage.inside(i, j) ||
                i < radius || i >= image.m - radius ||
                j < radius || j >= image.n - radius) {
            skip = true;
        }
        
        // remove its neighborhood from further considerations
        int mini = i - radius, maxi = i + radius, minj = j - radius, maxj = j + radius;
        if (map_metric::periodic[0]) {
            mini = std::max(mini, 0);
            maxi = std::min(maxi, (int)image.m - 1);
        }
        if (map_metric::periodic[1]) {
            minj = std::max(minj, 0);
            maxj = std::min(maxj, (int)image.n - 1);
        }
        for (int u = mini ; u <= maxi ; ++u) {
            int _u = (u + image.m) % image.m;
            for (int v = minj ; v <= maxj ; ++v) {
                int _v = (v + image.n) % image.n;
                __checked[_u+_v*image.m] = true;
            }
        }
        
        if (skip) {
            continue;
        }
        
        // we found a local extremum
        extrema.push_back(n);
    }
    
    std::cout << "found " << extrema.size() << " " << (max ? "maxima" : "minima")
              << " in " << __timer.elapsed() << "s." << std::endl;
}

}

#endif




















