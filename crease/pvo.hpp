#ifndef __PVO_HPP__
#define __PVO_HPP__

#include <vector>
#include <array>
#include <math/types.hpp>
#include "crease.hpp"

namespace spurt {
namespace crease
{
    // linear parallel vector operator on a face
    bool par_vec_op(std::vector<vec3> &beta, const mat3 &V,
                    const mat3 &W);

    bool linear_parallel_operator(std::vector<vec3> &b,
                                  const std::array<vec3, 3>& v,
                                  const std::array<vec3, 3>& w);

    // Stetten's averaging method for eigenvectors
    vec3 average(vec3 dirs[4]);

    struct search_face_PVO
    {
        int operator()(std::vector<vec3> &xing,
                       const vec3 &p0, const vec3 &p1,
                       const vec3 &p2, const vec3 &p3,
                       unsigned int maxdepth, bool &something_found) const;
    };

}
}


#endif













