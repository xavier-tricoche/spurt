#ifndef __PEIKERT_SADLO_HPP__
#define __PEIKERT_SADLO_HPP__

#include <Eigen/Eigen>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>

namespace xavier {

namespace peikert_sadlo_criterion {

double surface_criterion(const nvis::vec3& g, const nvis::mat3& H) {
    Eigen::Mat3 M;
    nvis::vec3 Hg = H*g;
    nvis::vec3 HHg = H*Hg;
    for (int r=0; r<3; ++r) {
        M(r,0) = g[r];
        M(r,1) = Hg[r];
        M(r,2) = HHg[r];
    }
    return M.determinant();
}

void compute_surface_criterion(std::vector<double>& vals,
                               const std::vector<nvis::vec3>& grads,
                               const std::vector<nvis::mat3>& hess) {
    for (size_t i=0; i<grads.size(); ++i) {
        const nvis::vec3& g = grads[i];
        const nvis::mat3& H = hess[i];
        vals[i] = surface_criterion(g, H);
    }
}

} // peikert_sadlo_criterion

} // xavier


#endif // __PEIKERT_SADLO_HPP__
