#ifndef ____XAVIER_LS_HPP__
#define ____XAVIER_LS_HPP__

#include <math/fixed_vector.hpp>
#include <math/matrix.hpp>

namespace xavier {

typedef xavier::matrix<double>  mat;

bool NR_svddcmp(const mat& A, std::vector<double> w, mat& U, mat& V);

inline nvis::vec4 ls_jacobian(const std::vector< nvis::vec2 >& pos,
                              const std::vector< nvis::vec2 >& vec,
                              const nvis::vec2& x0 = nvis::vec2(0, 0))
{
    assert(pos.size() == vec.size());
    
    // NOTE TO SELF: This is a very naive one. A SVD-based computation would require
    // to solve two equations independently.
    
    // d/dx, d/dy, K
    xavier::mat3 M;
    nvis::vec3 rhsx(0, 0, 0), rhsy(0, 0, 0);
    
    for (unsigned int i = 0 ; i < pos.size() ; ++i) {
        nvis::vec2 p = pos[i] - x0;
        const nvis::vec2& v = vec[i];
        
        M(0, 0) += p[0] * p[0];
        M(0, 1) += p[0] * p[1];
        M(0, 2) += p[0];
        rhsx[0] += p[0] * v[0];
        rhsy[0] += p[0] * v[1];
        
        M(1, 1) += p[1] * p[1];
        M(1, 2) += p[1];
        rhsx[1] += p[1] * v[0];
        rhsy[1] += p[1] * v[1];
        
        M(2, 2) += 1;
        rhsx[2] += v[0];
        rhsy[2] += v[1];
    }
    
    // LS matrix is symmetric
    M(1, 0) = M(0, 1);
    M(2, 0) = M(0, 2);
    M(2, 1) = M(1, 2);
    
    xavier::mat3 Inv = invert(M);
    nvis::vec4 J;
    nvis::vec3 tmp = Inv * rhsx;
    J[0] = tmp[0];
    J[1] = tmp[1];
    tmp = Inv * rhsy;
    J[2] = tmp[0];
    J[3] = tmp[1];
    
    return J;
}

inline nvis::vec2 ls_gradient(const std::vector< nvis::vec2 >& pos, const std::vector< double >& val)
{
    // d/dx, d/dy, K
    xavier::mat3 M;
    nvis::vec3 rhs(0, 0, 0);
    
    for (unsigned int i = 0 ; i < pos.size() ; ++i) {
        nvis::vec2 p = pos[i] - pos[0];
        
        M(0, 0) += p[0] * p[0];
        M(0, 1) += p[0] * p[1];
        M(0, 2) += p[0];
        rhs[0] += p[0] * val[i];
        
        M(1, 1) += p[1] * p[1];
        M(1, 2) += p[1];
        rhs[1] += p[1] * val[i];
        
        M(2, 2) += 1;
        rhs[2] += val[i];
    }
    
    // LS matrix is symmetric
    M(1, 0) = M(0, 1);
    M(2, 0) = M(0, 2);
    M(2, 1) = M(1, 2);
    
    xavier::mat3 Inv = invert(M);
    nvis::vec3 tmp = Inv * rhs;
    return nvis::vec2(tmp[0], tmp[1]);
}

}

#endif




















