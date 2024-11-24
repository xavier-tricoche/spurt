#include "math/svd.hpp"
#include <teem/hest.h>
#include <teem/ell.h>
#include <iostream>
#include <misc/progress.hpp>
#include <math/types.hpp>
#include <vector>


int main(int argc, char* argv[])
{

    if (argc != 2) {
        std::cerr << "USAGE: " << argv[0] << " <# tests>\n";
        exit(0);
    }
    
    std::cerr << "machine epsilon is " << DBL_EPSILON << '\n';
    
    spurt::mat3 A;
    Eigen::Matrix3d A_, Id;
    srand48(time(0));
    double err = 0;
    
    for (int i=0 ; i<3 ; ++i) {
        Id(i,i) = 1;
    }
    
    for (int n=0 ; n<atoi(argv[1]) ; ++n) {
    
        // random matrix A
        // for (int i=0 ; i<3 ; ++i) {
        //     for (int j=0 ; j<3 ; ++j) {
        //         A(i,j) = -1. + 2.*drand48();
        //     }
        // }
        A(0,0) = 2;
        A(0,1) = 15;
        A(0,2) = -4;
        A(1,0) = 15;
        A(1,1) = 127;
        A(1,2) = 0.65;
        A(2,0) = -4;
        A(2,1) = 0.65;
        A(2,2) = 127;
        
        A_(0,0) = 2;
        A_(0,1) = 15;
        A_(0,2) = -4;
        A_(1,0) = 15;
        A_(1,1) = 127;
        A_(1,2) = 0.65;
        A_(2,0) = -4;
        A_(2,1) = 0.65;
        A_(2,2) = 127;
        
        //         // make A degenerate
        //         double e[3];
        //         int ret = ell_3m_eigenvalues_d(e, A[0].begin(), AIR_TRUE);
        //         A -= e[0]*spurt::mat3::identity();
        // A_ -= e[0]*Id;
        
        // check SVD
        spurt::mat3 U,V;
        spurt::vec3 w;
        spurt::timer t;
        nr_svd::svdcmp<double, 3,3>(A, w, U, V, DBL_EPSILON);
        std::cerr << "NR SVD took " << t.elapsed() << '\n';
        spurt::mat3 Approx = U * nr_svd::to_mat<double, 3>(w) * spurt::transpose(V);
        std::cerr << "NR SVD error = " << spurt::norm(A-Approx)/spurt::norm(A) << std::endl;
        std::cerr << "NR w = " << w << std::endl;
        std::cerr << "NR U = " << U << std::endl;
        std::cerr << "NR V = " << V << std::endl;
        
        Eigen::Matrix3d U_, V_, w_;
        t.restart();
        eigen_svd::svdcmp(A_, U_, w_, V_);
        std::cerr << "EIGEN SVD took " << t.elapsed() << '\n';
        Eigen::Matrix3d A_approx = U_ * w_ * V_.transpose();
        std::cerr << "Eigen SVD error = " << (A_ - A_approx).norm()/A_.norm() << std::endl;
        std::cerr << "Eigen w = " << '\n' << w_ << std::endl;
        std::cerr << "Eigen U = " << '\n' << U_ << std::endl;
        std::cerr << "Eigen V = " << '\n' << V_ << std::endl;
        
        // compute pseudoinverse
        spurt::mat3 B = nr_svd::pseudoinv<double, 3, 3>(A, DBL_EPSILON);
        
        // ad-hoc pseudo inverse of diagonal matrix
        std::vector<double> svs;
        for (int i=0 ; i< 3 ; ++i) {
            svs.push_back(w_(i,i));
        }
        double max = *std::max_element(svs.begin(), svs.end());
        for (int i=0 ; i< 3 ; ++i) {
            if (w_(i,i)/max > 1.0e-9) {
                w_(i,i) = 1./w_(i,i);
            } else {
                w_(i,i) = 0;
            }
        }
        Eigen::Matrix3d B_ = V_ * w_ * U_.transpose();
        
        std::cerr << "NR pseudoinv = " << B << std::endl;
        std::cerr << "Eigen pseudoinv = " << B << std::endl;
        
        // check pseudoinverse
        spurt::mat3 a = A*B*A;
        spurt::mat3 b = B*A*B;
        double err1 = spurt::norm(a-A)/spurt::norm(A);
        double err2 = spurt::norm(b-B)/spurt::norm(B);
        std::cerr << "NR: " << err1 << '\t' << err2 << std::endl;
        err += err1 + err2;
        
        Eigen::Matrix3d a_ = A_*B_*A_;
        Eigen::Matrix3d b_ = B_*A_*B_;
        err1 = (a_ - A_).norm() / A_.norm();
        err2 = (b_ - B_).norm() / B_.norm();
        std::cerr << "Eigen: " << err1 << '\t' << err2 << std::endl;
        
    }
    std::cerr << "Average error: " << 0.5*err/(double)atoi(argv[1]) << '\n';
    return 0;
}
