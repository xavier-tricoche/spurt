#ifndef __MULTIDIMENSIONAL_SCALING_HPP__
#define __MULTIDIMENSIONAL_SCALING_HPP__

#include <memory>

// Eigen
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace xavier {

template <typename T>
class MultidimensionalScaling {
public:
    typedef T scalar_t;
    typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> matrix_t;
    typedef Eigen::Matrix<scalar_t, Eigen::Dynamic, 1> col_vector_t;
    typedef col_vector_t vector_t;
    typedef Eigen::Matrix<scalar_t, 1, Eigen::Dynamic> row_vector_t;
    typedef Eigen::SelfAdjointEigenSolver<matrix_t> eigensolver_t;
    
    MultidimensionalScaling() {}
    
    void embed(const matrix_t& distances) {
        Tau = tau_func(distances);
        solver=std::shared_ptr<eigensolver_t>(new eigensolver_t(Tau));
    }
    
    matrix_t eigenvectors() const {
        return solver->eigenvectors();
    }
    
    col_vector_t eigenvector(int k) const {
        return solver->eigenvectors().col(k);
    }
    
    col_vector_t eigenvalues() const {
        return solver->eigenvalues();
    }
    
    scalar_t eigenvalue(int k) const {
        return solver->eigenvalues()(k);
    }
    
    matrix_t coordinates(int ndim=2) const {
        return solver->eigenvectors().rightCols(ndim)*
            matrix_t(solver->eigenvalues().tail(ndim).array().sqrt().
                matrix().asDiagonal());
    }
    
    matrix_t tau_func(const matrix_t& D) const {
        matrix_t S, H;
        const size_t N=D.rows();
        S=D.array()*D.array();        
        H=matrix_t::Identity(N, N).array()-1./scalar_t(N);
        return -0.5*H*S*H;
    }
    
    scalar_t l2_error(int ndim=2) const {
        matrix_t Y = coordinates(ndim);
        const size_t N=Y.rows();
        // compute pairwise distances
        matrix_t Dy = matrix_t::Zero(N, N);
        for (int i=0; i<N; ++i) {
            const row_vector_t& yi = Y.row(i);
            for (int j=0; j<i ; ++j) {
                const row_vector_t& yj = Y.row(j);
                Dy(i,j) = Dy(j,i) = (yi-yj).norm();
            }
        }
        return (Tau-tau_func(Dy)).norm();
    }
    
protected:
    std::shared_ptr<eigensolver_t> solver;
    matrix_t Tau;
};

} // namespace xavier

#endif