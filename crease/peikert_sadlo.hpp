#ifndef __PEIKERT_SADLO_HPP__
#define __PEIKERT_SADLO_HPP__

#include <math/types.hpp>

namespace spurt {

namespace peikert_sadlo_criterion {

double surface_criterion(const vec3& g, const mat3& H) {
    mat3 M;
    vec3 Hg = H*g;
    vec3 HHg = H*Hg;
    for (int r=0; r<3; ++r) {
        M(r,0) = g[r];
        M(r,1) = Hg[r];
        M(r,2) = HHg[r];
    }
    return M.determinant();
}

void compute_surface_criterion(std::vector<double>& vals,
                               const std::vector<vec3>& grads,
                               const std::vector<mat3>& hess) {
    for (size_t i=0; i<grads.size(); ++i) {
        const vec3& g = grads[i];
        const mat3& H = hess[i];
        vals[i] = surface_criterion(g, H);
    }
}

} // peikert_sadlo_criterion

} // spurt


#endif // __PEIKERT_SADLO_HPP__
