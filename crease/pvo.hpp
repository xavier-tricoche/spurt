#ifndef __PVO_HPP__
#define __PVO_HPP__

#include <vector>
#include <math/fixed_vector.hpp>
#include <math/matrix.hpp>
#include "crease.hpp"

namespace spurt {
namespace crease {

// linear parallel vector operator on a face
bool par_vec_op(std::vector< vec3 >& beta, const spurt::mat3& V,
                const spurt::mat3& W);
                
bool linear_parallel_operator(std::vector< vec3 >& b,
                              const vec3 v[3],
                              const vec3 w[3]);
                              
// Stetten's averaging method for eigenvectors
vec3 average(vec3 dirs[4]);


struct search_face_PVO {
    int operator()(std::vector< vec3 >& xing,
                   const vec3& p0, const vec3& p1,
                   const vec3& p2, const vec3& p3,
                   unsigned int maxdepth, bool& something_found) const;
};

}
}


#endif













