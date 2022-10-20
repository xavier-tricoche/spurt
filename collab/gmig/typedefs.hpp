#ifndef __XAVIER_COLLAB_GMIG_TYPEDEFS__
#define __XAVIER_COLLAB_GMIG_TYPEDEFS__

// STL
#include <vector>
// Eigen
#include <Eigen/Dense>
// nvis
#include <math/fixed_vector.hpp>
// spurt
#include <math/RBF.hpp>
#include <math/RBFbasis.hpp>

namespace spurt { namespace gmig {
    
static constexpr double invalid_double = std::numeric_limits<double>::max();
    
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  matrix_type;
typedef Eigen::HouseholderQR<matrix_type>                      solver_type;

// solvers - from fast to slow
typedef Eigen::HouseholderQR<matrix_type>          fast_solver_type;
typedef Eigen::ColPivHouseholderQR<matrix_type>    medium_solver_type;
typedef Eigen::FullPivHouseholderQR<matrix_type>   slow_solver_type;

// kernels - both local and global support
typedef RBF::linear_function<double>             linear_type;
typedef RBF::cubic_function<double>              cubic_type;
typedef RBF::quintic_function<double>            quintic_type;
typedef RBF::gaussian_function<double>           gaussian_type;
typedef RBF::wendland_function<double>           wendland_type;
typedef RBF::truncated_gaussian_function<double> truncated_gaussian_type;

// RBF interpolators with infinite support
// linear
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, linear_type, fast_solver_type>     fast_linear_rbf_type;
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, linear_type, medium_solver_type>   medium_linear_rbf_type;
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, linear_type, slow_solver_type>     slow_linear_rbf_type;
// cubic
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, cubic_type, fast_solver_type>      fast_cubic_rbf_type;
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, cubic_type, medium_solver_type>    medium_cubic_rbf_type;
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, cubic_type, slow_solver_type>      slow_cubic_rbf_type;
// quintic
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, quintic_type, fast_solver_type>    fast_quintic_rbf_type;
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, quintic_type, medium_solver_type>  medium_quintic_rbf_type;
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, quintic_type, slow_solver_type>    slow_quintic_rbf_type;
// gaussian
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, gaussian_type, fast_solver_type>   fast_gaussian_rbf_type;
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, gaussian_type, medium_solver_type> medium_gaussian_rbf_type;
typedef RBF::InfiniteSupportRBFInterpolator<double, double, 2, gaussian_type, slow_solver_type>   slow_gaussian_rbf_type;

// RBF interpolators with compact support
// Wendland C2
typedef RBF::CompactSupportRBFInterpolator<double, double, 2, wendland_type>            wendland_rbf_type;
// truncated gaussian
typedef RBF::CompactSupportRBFInterpolator<double, double, 2, truncated_gaussian_type>  trunc_gaussian_rbf_type;

} // gmig
} // spurt


#endif