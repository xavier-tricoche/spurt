#ifndef __PVO_HPP__
#define __PVO_HPP__

#include <vector>
#include <math/fixed_vector.hpp>
#include <math/matrix.hpp>
#include "crease.hpp"

namespace spurt {
namespace crease {

// linear parallel vector operator on a face
bool par_vec_op(std::vector< nvis::vec3 >& beta, const spurt::mat3& V,
                const spurt::mat3& W);
                
bool linear_parallel_operator(std::vector< nvis::vec3 >& b,
                              const nvis::vec3 v[3],
                              const nvis::vec3 w[3]);
                              
// Stetten's averaging method for eigenvectors
nvis::vec3 average(nvis::vec3 dirs[4]);


struct search_face_PVO {
    int operator()(std::vector< nvis::vec3 >& xing,
                   const nvis::vec3& p0, const nvis::vec3& p1,
                   const nvis::vec3& p2, const nvis::vec3& p3,
                   unsigned int maxdepth, bool& something_found) const;
};

}
}


#endif













