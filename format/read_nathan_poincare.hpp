#ifndef __read_nathan_poincare_hpp__
#define __read_nathan_poincare_hpp__

#include <math/fixed_vector.hpp>
#include <vector>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <math/bounding_box.hpp>

namespace xavier {

inline nvis::bbox3 read_nathan_poincare(std::vector<std::vector<nvis::vec3> >& orbits,
                                        const std::string& name)
{
    typedef std::vector<nvis::vec3> orbit_type;
    
    FILE* file = fopen(name.c_str(), "rb");
    
    unsigned int N = 0;
    fread(&N, sizeof(unsigned int), 1, file);
    
    if (N > 0) {
        orbits.resize(N);
    }
    
    nvis::bbox3 bounds;
    
    for (int i = 0; i < N ; i++) {
        unsigned int size = 0;
        fread(&size, sizeof(unsigned int), 1, file);
        if (!size) {
            continue;
        }
        orbit_type& orbit = orbits[i];
        orbit.resize(size);
        
        std::vector<nvis::vec4> tmp(size);
        
        fread(&tmp[0], sizeof(double), 4*size, file);
        
        for (int j = 0 ; j < size ; ++j) {
            orbit[j] = nvis::vec3(tmp[j][0], tmp[j][2], tmp[j][3]);
            // std::cerr << orbit[j] << '\n';
            bounds.add(orbit[j]);
        }
        // std::cerr << '\n';
    }
    
    fclose(file);
    return bounds;
}

} // namespace xavier

#endif





