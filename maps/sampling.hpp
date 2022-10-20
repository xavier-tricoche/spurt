#ifndef __SAMPLING_HPP__
#define __SAMPLING_HPP__

#include <maps/basic_definitions.hpp>
#include <math/bounding_box.hpp>
#include <tokamak/poincare_map.hpp>

namespace spurt {
void sample_along_line(const nvis::bbox2& bounds, const poincare_map& pmap,
                       double x0, int niter, int nprobes, double dx = 0)
{
    poincare_map __map(pmap);
    std::vector<nvis::vec2> steps;
    srand48(time(0));
    double min = bounds.min()[1];
    double max = bounds.max()[1];
    double length = max - min;
    double spacing = length / (double)(nprobes - 1);
    std::vector<spurt::orbit> orbits;
    
    orbits.reserve(nprobes);
    for (int i = 0 ; i < nprobes ; ++i) {
        double dy = 0.1 * spacing * (-1. + 2 * drand48());
        double y = (double)i * spacing + dy;
        y = std::max(min, std::min(y, max));
        double x = x0;
        if (dx > 0) {
            x += 0.5 * dx * (-1. + 2 * drand48());
            x = std::max(bounds.min()[0], std::min(bounds.max()[0], x));
        }
        nvis::vec2 seed(x, y);
        try {
            __map.map(seed, steps, niter);
        } catch (...) {
            continue;
        }
        
        if (steps.size() < niter) {
            continue;
        }
        steps.push_back(seed);
        orbits.push_back(spurt::orbit(steps, 0));
    }
    
    std::copy(orbits.begin(), orbits.end(), std::back_inserter(spurt::__map_orbits));
    std::cout << "1D sampling at x=" << x0 << " resulted in " << orbits.size() << " new orbits\n";
}

void sample_on_raster(const nvis::bbox2& bounds, const poincare_map& pmap,
                      int niter, int nprobes[2], nvis::vec2 dp = nvis::vec2(0, 0))
{
    poincare_map __map(pmap);
    std::vector<nvis::vec2> steps;
    srand48(time(0));
    
    double du = dp[0];
    double dv = dp[1];
    
    nvis::vec2 d = bounds.size();
    double spacing[] = { d[0] / (double)(nprobes[0] - 1),
                         d[1] / (double)(nprobes[1] - 1)
                       };
                       
    du = std::min(du, spacing[0]);
    dv = std::min(dv, spacing[1]);
    
    std::vector<spurt::orbit> orbits;
    orbits.reserve(nprobes[0]*nprobes[1]);
    for (int i = 0 ; i < nprobes[0] ; ++i) {
        for (int j = 0 ; j < nprobes[1] ; ++j) {
            nvis::vec2 seed(bounds.min());
            double dx = 0.5 * du * (-1. + 2.*drand48());
            double dy = 0.5 * dv * (-1. + 2.*drand48());
            seed[0] += i * spacing[0] + dx;
            seed[1] += j * spacing[1] + dy;
            seed[0] = std::max(bounds.min()[0], std::max(bounds.min()[0], seed[0]));
            seed[1] = std::max(bounds.min()[1], std::max(bounds.min()[1], seed[1]));
            
            try {
                __map.map(seed, steps, niter);
            } catch (...) {
                continue;
            }
            
            if (steps.size() < niter) {
                continue;
            }
            steps.push_back(seed);
            orbits.push_back(spurt::orbit(steps, 0));
        }
    }
    
    std::copy(orbits.begin(), orbits.end(), std::back_inserter(spurt::__map_orbits));
    std::cout << "2D sampling in " << bounds << " resulted in " << orbits.size() << " new orbits\n";
}

}

#endif




















