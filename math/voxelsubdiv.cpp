#include <iostream>
#include <vector>
#include <math/fixed_vector.hpp>
#include <data/grid.hpp>

int main(int argc, char* argv[])
{
    using namespace xavier;
    
    nvis::vec3 v[8] = { nvis::vec3(0, 0, 0), nvis::vec3(1, 0, 0),
                        nvis::vec3(1, 1, 0), nvis::vec3(0, 1, 0),
                        nvis::vec3(0, 0, 1), nvis::vec3(1, 0, 1),
                        nvis::vec3(1, 1, 1), nvis::vec3(0, 1, 1),
                      };
                      
    nvis::vec3 f[6] = { nvis::vec3(0.5, 0.5, 0), nvis::vec3(0.5, 0.5, 1),
                        nvis::vec3(0.5, 0, 0.5), nvis::vec3(1, 0.5, 0.5),
                        nvis::vec3(0.5, 1, 0.5), nvis::vec3(0, 0.5, 0.5)
                      };
                      
    nvis::vec3 c(0.5, 0.5, 0.5);
    
    nvis::vec3 tets[24][4];
    for (int face = 0 ; face < 6 ; ++face) {
        const nvis::ivec4 aface = xavier::grid::voxel_faces[face];
        for (int i = 0 ; i < 4 ; ++i) {
            tets[4*face+i][0] = v[aface[i]];
            tets[4*face+i][1] = v[aface[(i+1)%4]];
            tets[4*face+i][2] = f[face];
            tets[4*face+i][3] = c;
        }
    }
    
    std::vector<int> id2tet(64, -1);
    for (int i = 0 ; i < 24 ; ++i) {
        nvis::vec3 center = 0.25 * (tets[i][0] + tets[i][1] + tets[i][2] + tets[i][3]);
        double u = center[0], v = center[1], w = center[2];
        unsigned int code = 0;
        if (u <= v) {
            code |= 1;
        }
        if (u <= 1 - v) {
            code |= 2;
        }
        if (v <= w) {
            code |= 4;
        }
        if (v <= 1 - w) {
            code |= 8;
        }
        if (w <= u) {
            code |= 16;
        }
        if (w <= 1 - u) {
            code |= 32;
        }
        
        std::cerr << "tet " << i << "\tface " << i / 4
                  << "\tvertices: 0: " << tets[i][0] << "\t1: " << tets[i][1] << "\t2: " << tets[i][2] << "\t3: " << tets[i][3]
                  << "\tx = " << center
                  << "\t(u<=v) " << (u <= v ? "TRUE" : "FALSE")
                  << "\t(u<=1-v) " << (u <= 1 - v ? "TRUE" : "FALSE")
                  << "\t(v<=w) " << (v <= w ? "TRUE" : "FALSE")
                  << "\t(v<=1-w) " << (v <= 1 - w ? "TRUE" : "FALSE")
                  << "\t(w<=u) " << (w <= u ? "TRUE" : "FALSE")
                  << "\t(w<=1-u) " << (w <= 1 - u ? "TRUE" : "FALSE")
                  << "\t code " << code << '\n';
        id2tet[code] = i;
    }
    
    for (int i=0 ; i<64 ; ++i) {
        std::cerr << id2tet[i] << ", ";
    }
    std::cerr << '\n';
    
    return 0;
}











