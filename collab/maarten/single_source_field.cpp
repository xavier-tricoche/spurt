#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#ifdef _OPENMP
#include <omp.h>
#endif
// nvis
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/timer.hpp>
// teem
#include <teem/nrrd.h>
// xavier
#include <math/RBF.hpp>
#include <math/RBFbasis.hpp>
#include <image/nrrd_wrapper.hpp>
#include <format/format.hpp>
#include "maarten_utils.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/bind.hpp>
#include <boost/ref.hpp>

const double TWO_PI = 2.*M_PI;

// global variables
nvis::ivec2 resolution(800, 800);
nvis::bbox2 bounds;
bool verbose = false;
bool do_reconstruct = false;
bool do_export_weights = false;
std::string kernel_name = "r3";
std::string accuracy = "low";
nvis::vec2 source;
unsigned poly_order = 0;

class progress_message {
    std::string        _what;
    size_t             _size;
    
    mutable std::ostringstream _os;
    
public:
    progress_message(size_t size, std::string what = "")
        : _size(size), _what(what) {}
        
    std::string operator()(size_t n, double elapsed=0) const {
        _os.clear();
        _os.str("");
        _os << "\rCompleted " << 100.*(float)n/float(_size)
            << "\% of " << _size;
        if (_what.size()) {
            _os << " " << _what;
        }
        if (elapsed > 0) {
            _os << " in " << elapsed << " seconds (" << (float)n/elapsed << " Hz)";
        }
        _os << "                       \r";
        return _os.str();
    }
    
    void reset(const std::string& what = "") {
        if (what.size()) {
            _what = what;
        }
    }
};

namespace xrbf = xavier::RBF;

typedef nvis::vec1 data_type;
typedef nvis::vec2 point_type;
const size_t N = 2;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>     matrix_type;

// solvers - from fast to slow
typedef Eigen::HouseholderQR<matrix_type>          fast_solver_type;
typedef Eigen::ColPivHouseholderQR<matrix_type>    medium_solver_type;
typedef Eigen::FullPivHouseholderQR<matrix_type>   slow_solver_type;

// kernels - both local and global support
typedef xrbf::linear_function<double>             linear_type;
typedef xrbf::cubic_function<double>              cubic_type;
typedef xrbf::quintic_function<double>            quintic_type;
typedef xrbf::gaussian_function<double>           gaussian_type;
typedef xrbf::wendland_function<double>           wendland_type;
typedef xrbf::truncated_gaussian_function<double> truncated_gaussian_type;

// RBF interpolators with infinite support

// linear
// 0-th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, fast_solver_type>        fast_linear_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, medium_solver_type>      medium_linear_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, slow_solver_type>        slow_linear_rbf_type;
// 1st order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, fast_solver_type, 1>     fast_linear_pp1_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, medium_solver_type, 1>   medium_linear_pp1_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, slow_solver_type, 1>     slow_linear_pp1_rbf_type;
// 2nd order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, fast_solver_type, 2>     fast_linear_pp2_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, medium_solver_type, 2>   medium_linear_pp2_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, slow_solver_type, 2>     slow_linear_pp2_rbf_type;
// 3rd order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, fast_solver_type, 3>     fast_linear_pp3_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, medium_solver_type, 3>   medium_linear_pp3_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, slow_solver_type, 3>     slow_linear_pp3_rbf_type;
// 4th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, fast_solver_type, 4>     fast_linear_pp4_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, medium_solver_type, 4>   medium_linear_pp4_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, slow_solver_type, 4>     slow_linear_pp4_rbf_type;
// 5th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, fast_solver_type, 5>     fast_linear_pp5_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, medium_solver_type, 5>   medium_linear_pp5_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, slow_solver_type, 5>     slow_linear_pp5_rbf_type;
// 6th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, fast_solver_type, 6>     fast_linear_pp6_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, medium_solver_type, 6>   medium_linear_pp6_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, linear_type, slow_solver_type, 6>     slow_linear_pp6_rbf_type;

// cubic
// 0th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, fast_solver_type>         fast_cubic_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, medium_solver_type>       medium_cubic_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, slow_solver_type>         slow_cubic_rbf_type;
// 1st order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, fast_solver_type, 1>      fast_cubic_pp1_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, medium_solver_type, 1>    medium_cubic_pp1_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, slow_solver_type, 1>      slow_cubic_pp1_rbf_type;
// 2nd order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, fast_solver_type, 2>      fast_cubic_pp2_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, medium_solver_type, 2>    medium_cubic_pp2_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, slow_solver_type, 2>      slow_cubic_pp2_rbf_type;
// 3rd order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, fast_solver_type, 3>      fast_cubic_pp3_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, medium_solver_type, 3>    medium_cubic_pp3_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, slow_solver_type, 3>      slow_cubic_pp3_rbf_type;
// 4th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, fast_solver_type, 4>      fast_cubic_pp4_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, medium_solver_type, 4>    medium_cubic_pp4_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, slow_solver_type, 4>      slow_cubic_pp4_rbf_type;
// 5th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, fast_solver_type, 5>      fast_cubic_pp5_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, medium_solver_type, 5>    medium_cubic_pp5_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, slow_solver_type, 5>      slow_cubic_pp5_rbf_type;
// 6th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, fast_solver_type, 6>      fast_cubic_pp6_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, medium_solver_type, 6>    medium_cubic_pp6_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, cubic_type, slow_solver_type, 6>      slow_cubic_pp6_rbf_type;

// quintic
// 0th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, fast_solver_type>       fast_quintic_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, medium_solver_type>     medium_quintic_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, slow_solver_type>       slow_quintic_rbf_type;
// 1st order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, fast_solver_type, 1>    fast_quintic_pp1_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, medium_solver_type, 1>  medium_quintic_pp1_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, slow_solver_type, 1>    slow_quintic_pp1_rbf_type;
// 2nd order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, fast_solver_type, 2>    fast_quintic_pp2_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, medium_solver_type, 2>  medium_quintic_pp2_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, slow_solver_type, 2>    slow_quintic_pp2_rbf_type;
// 3rd order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, fast_solver_type, 3>    fast_quintic_pp3_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, medium_solver_type, 3>  medium_quintic_pp3_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, slow_solver_type, 3>    slow_quintic_pp3_rbf_type;
// 4th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, fast_solver_type, 4>    fast_quintic_pp4_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, medium_solver_type, 4>  medium_quintic_pp4_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, slow_solver_type, 4>    slow_quintic_pp4_rbf_type;
// 5th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, fast_solver_type, 5>    fast_quintic_pp5_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, medium_solver_type, 5>  medium_quintic_pp5_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, slow_solver_type, 5>    slow_quintic_pp5_rbf_type;
// 6th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, fast_solver_type, 6>    fast_quintic_pp6_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, medium_solver_type, 6>  medium_quintic_pp6_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, quintic_type, slow_solver_type, 6>    slow_quintic_pp6_rbf_type;

// gaussian
// 0th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, fast_solver_type>      fast_gaussian_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, medium_solver_type>    medium_gaussian_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, slow_solver_type>      slow_gaussian_rbf_type;
// 1st order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, fast_solver_type, 1>   fast_gaussian_pp1_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, medium_solver_type, 1> medium_gaussian_pp1_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, slow_solver_type, 1>   slow_gaussian_pp1_rbf_type;
// 2nd order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, fast_solver_type, 2>   fast_gaussian_pp2_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, medium_solver_type, 2> medium_gaussian_pp2_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, slow_solver_type, 2>   slow_gaussian_pp2_rbf_type;
// 3rd order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, fast_solver_type, 3>   fast_gaussian_pp3_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, medium_solver_type, 3> medium_gaussian_pp3_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, slow_solver_type, 3>   slow_gaussian_pp3_rbf_type;
// 4th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, fast_solver_type, 4>   fast_gaussian_pp4_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, medium_solver_type, 4> medium_gaussian_pp4_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, slow_solver_type, 4>   slow_gaussian_pp4_rbf_type;
// 5th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, fast_solver_type, 5>   fast_gaussian_pp5_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, medium_solver_type, 5> medium_gaussian_pp5_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, slow_solver_type, 5>   slow_gaussian_pp5_rbf_type;
// 6th order
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, fast_solver_type, 6>   fast_gaussian_pp6_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, medium_solver_type, 6> medium_gaussian_pp6_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<data_type, double, N, gaussian_type, slow_solver_type, 6>   slow_gaussian_pp6_rbf_type;

// RBF interpolators with compact support

// Wendland C2
// 0th-order
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, wendland_type>     wendland_rbf_type;
// 1st-order                                                                         
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, wendland_type, 1>  wendland_pp1_rbf_type;
// 2nd-order                                                                         
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, wendland_type, 2>  wendland_pp2_rbf_type;
// 3rd-order                                                                         
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, wendland_type, 3>  wendland_pp3_rbf_type;
// 4th-order                                                                         
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, wendland_type, 4>  wendland_pp4_rbf_type;
// 5th-order                                                                         
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, wendland_type, 5>  wendland_pp5_rbf_type;
// 6th-order                                                                         
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, wendland_type, 6>  wendland_pp6_rbf_type;

// truncated gaussian
// 0th order
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, truncated_gaussian_type>     trunc_gaussian_rbf_type;
// 1st order
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, truncated_gaussian_type, 1>  trunc_gaussian_pp1_rbf_type;
// 2nd order
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, truncated_gaussian_type, 2>  trunc_gaussian_pp2_rbf_type;
// 3rd order
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, truncated_gaussian_type, 3>  trunc_gaussian_pp3_rbf_type;
// 4th order
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, truncated_gaussian_type, 4>  trunc_gaussian_pp4_rbf_type;
// 5th order
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, truncated_gaussian_type, 5>  trunc_gaussian_pp5_rbf_type;
// 6th order
typedef xrbf::CompactSupportRBFInterpolator<data_type, double, N, truncated_gaussian_type, 6>  trunc_gaussian_pp6_rbf_type;


template<typename _Interpolator>
float* reconstruct(const _Interpolator& interpolator, const nvis::ivec2& res)
{
    typedef typename _Interpolator::data_type       value_type;
    typedef typename _Interpolator::derivative_type derivative_type;
    
    point_type spacing = bounds.size() / point_type(res - nvis::ivec2(1, 1));
    size_t number_of_samples = res[0]*res[1];
    float* result = (float*)calloc(3*number_of_samples, sizeof(float));
    
    size_t number_of_threads = 1;
#ifdef _OPENMP
    number_of_threads = omp_get_max_threads();
#endif
    size_t counter = 0;
    nvis::timer _timer;
    
    progress_message msg(number_of_samples, "interpolations");
    
    const nvis::vec2& source = interpolator.points()[0];
    
    double LARGE = std::numeric_limits<double>::max();
    std::vector<double> min_val(number_of_threads, LARGE);
    std::vector<nvis::vec2> minval_pos(number_of_threads);
    std::vector<double> min_gm(number_of_threads, LARGE);
    std::vector<nvis::vec2> mingm_pos(number_of_threads);
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (size_t n=0 ; n<number_of_samples ; ++n) {
            int i = n%res[0];
            int j = n/res[0];
            point_type x = bounds.min() + point_type(i, j)*spacing;
            value_type t = interpolator(x);
            derivative_type dt = interpolator.derivative(x);
            result[3*n  ] = t[0];
            result[3*n+1] = dt[0][0];
            result[3*n+2] = dt[1][0];
            
            int thread = 0;
#if _OPENMP
                thread = omp_get_thread_num();
#endif
            
            if (t[0] < min_val[thread]) {
                min_val[thread] = t[0];
                minval_pos[thread] = x;
            }
            
            double gm = nvis::norm(nvis::vec2(dt[0][0], dt[1][0]));
            if (gm < min_gm[thread]) {
                min_gm[thread] = gm;
                mingm_pos[thread] = x;
            }
            
            if (verbose) {
                if (!thread) {
                    double elapsed = _timer.elapsed();
                    std::cout << msg(n, elapsed) << std::flush;
                }
            }
        }
    }
    
    std::cout << "\nRBF reconstruction completed in "
              << _timer.elapsed() << " seconds ("
              << (float)number_of_samples/_timer.elapsed() << " Hz)\n";
              
    int dist = std::distance(min_val.begin(), std::min_element(min_val.begin(), min_val.end()));
    std::cout << "min sampled value = " << min_val[dist] << " at " << minval_pos[dist] << '\n';
    dist = std::distance(min_gm.begin(), std::min_element(min_gm.begin(), min_gm.end()));
    std::cout << "min grad mag = " << min_gm[dist] << " at " << mingm_pos[dist] << '\n';
       
    return result;
}

nvis::bbox2
read_nrrd(std::vector<point_type>& points,
          std::vector<data_type>& times,
          std::vector<data_type>& weights,
          const std::string& filename)
{
    const size_t psizes [] = { 1, 3, 6, 10, 15, 21, 28 };
    size_t npolys = 0;
    if (poly_order) {
        npolys = psizes[poly_order];
    }
    
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, filename.c_str(), NULL)) {
        std::cerr << "read_nrrd: " << biffGetDone(NRRD) << std::endl;
        throw;
    }
    std::vector<double> data;
    xavier::to_vector(data, nin);
    size_t N = nin->axis[0].size;
    bool has_weights = (N == 4);
    size_t nb_rows = data.size()/N;
    size_t nb_pts = nb_rows - npolys;
    points.resize(nb_pts);
    times.resize(nb_pts);
    if (has_weights) {
        weights.resize(nb_rows);
    }
    if (verbose && has_weights) {
        std::cout << "input data contains precomputed weights\n";
    }
    nvis::bbox2 _bounds;
    source = nvis::vec2(data[0], data[1]); // source is 1st position
    for (size_t i=0 ; i<nb_pts ; ++i) {
        points[i][0] = data[N*i  ];
        points[i][1] = data[N*i+1];
        times[i][0]  = data[N*i+2];
        if (has_weights) {
            weights[i][0] = data[N*i+3];
        }
        _bounds.add(points[i]);
    }
    if (has_weights) {
        for (size_t i=nb_pts ; i<nb_rows ; ++i) weights[i][0] = data[N*i+3];
    }
    nrrdNuke(nin);
    return _bounds;
}

// save full information about RBF reconstruction
void save_rbf(const std::vector<point_type>& points,
              const std::vector<data_type>& values,
              const std::vector<data_type>& weights,
              const std::string& kernel_name,
              const std::string& filename)
{
    size_t npts = points.size();
    size_t nweights = weights.size();
    double* data = (double*)calloc(nweights*4, sizeof(double));
    for (size_t i=0 ; i<npts ; ++i) {
        data[4*i  ] = points[i][0];
        data[4*i+1] = points[i][1];
        data[4*i+2] = values[i][0];
        data[4*i+3] = weights[i][0];
    }
    for (size_t i=npts ; i<nweights ; ++i) {
        data[4*i  ] = 0;
        data[4*i+1] = 0;
        data[4*i+2] = 0;
        data[4*i+3] = weights[i][0];
    }
    xavier::nrrd_params<double, 2> params;
    params.sizes()[0] = 4;
    params.sizes()[1] = nweights;
    params.description() = "RBF reconstruction data for kernel " + kernel_name;
    params.labels()[0] = "x_rec;y_rec,time,weight";
    params.labels()[1] = "RBF centers";
    xavier::writeNrrd(data, filename, params);
    std::cout << filename << " has been exported\n";
}

std::string me;
void usage(const std::string& message="")
{
    if (!message.empty()) {
        std::cerr << "ERROR: " << message << '\n';
    }
    std::cout << '\n'
              << "DESCRIPTION: Compute smooth (value + derivative) reconstruction\n"
              << "coefficient of discrete tomography travel time data associated\n"
              << "with a single source and/or a high-resolution field\n"
              << "reconstruction of the data using these coefficients.\n"
              << '\n'
              << "USAGE: " << me << " [parameters] [options]\n"
              << '\n'
              << "PARAMETERS:\n"
              << " -i | --input <string>      Input file name\n"
              << "                            (precomputed RBF weights accepted)\n"
              << '\n'                         
              << "OPTIONS:\n"                 
              << " -h | --help                Print this information\n"
              << " -b | --bounds <float>x4    Area to reconstruct:\n"
              << "                              lon_min lon_max lat_min lat_max\n"
              << " -p | --path <string>       Path to prepend to input file names\n"
              << " -w | --weights <string>    Export computed RBF weights to this file\n"
              << " -o | --output <string>     Export reconstructed field to this file\n"
              << " -k | --kernel <string>     Reconstruction kernel (\"r\", \"r3\", \"r5\",\n"
              << "                              \"gaussian:\"<sigma>,\n"
              << "                              \"wendland:\"<radius>,\n"
              << "                              \"tgaussian:\"<sigma>:<radius>)\n"
              << " -n | --order <int>         Polynomial precision order\n"
              << " -a | --accuracy <string>   Accuracy of linear algebra solver\n"
              << "                              (\"low\", \"medium\", \"high\")\n"
              << " -r | --resolution <int>x2  Resampling resolution: width height\n"
              << " -v | --verbose             Turn on verbose mode\n"
              << std::endl;
    exit(!message.empty());
}

void point_distribution(const std::vector<nvis::vec2>& points) {
    using namespace boost::accumulators;
    const std::vector<nvis::vec2>::size_type N = points.size();
    std::vector<double> dist(N*(N-1)/2);
    int c=0;
    
    for (int i=0 ; i<N ; ++i) {
        for (int j=i+1 ; j<N ; ++j, ++c) {
            dist[c] = nvis::norm(points[i]-points[j]);
        }
    }
    
    // compute statistics
    accumulator_set<double, 
                    features<tag::min, 
                             tag::max, 
                             tag::mean, 
                             tag::variance > > acc;
    std::for_each(dist.begin(), dist.end(), boost::bind<void>(boost::ref(acc), _1));
    std::cout << "statistics of point distribution:\n"
              << "\tdiameter: " << max(acc) << '\n'
              << "\tmin:      " << min(acc) << '\n'
              << "\taverage:  " << mean(acc) << '\n'
              << "\tvariance: " << variance(acc) << '\n';
}

template< typename Func >
std::pair<nvis::vec2, double> 
find_minimum(const nvis::vec2& ic, const Func& f, double h=0.01, double eps=1.0e-6) {
    typedef typename Func::derivative_type derivative_type;

    double init_val = f(ic)[0];
    derivative_type d = f.derivative(ic);
    nvis::vec2 g(d[0][0], d[1][0]);
    double init_norm = nvis::norm(g);

    if (verbose) {
        std::cout << "find_minimum: starting at " << ic << ", value=" << init_val << ", norm=" << init_norm << '\n';
        std::cout << "d = " << d << ", g = " << g << '\n';
    }

    double norm = init_norm;
    nvis::vec2 x = ic;
    double val = init_val;
    for (int n=0 ; n<100 && norm > eps ; ++n) {
        x -= h*g;
        val = f(x)[0];
        d = f.derivative(x);
        g = nvis::vec2(d[0][0], d[1][0]);
        norm = nvis::norm(g);
        if (verbose) {
            std::cout << "\tat " << x << ", f=" << val << ", g=" << g << " (" << norm << ")\n";
        }
    }
    return std::make_pair(x, val);
}

template<typename Int_>
float* reconstruct(const std::vector<point_type>& points,
                   const std::vector<data_type>& times,
                   const std::vector<data_type>& weights,
                   const std::string& file_name,
                   const typename Int_::function_type& fun =
                       typename Int_::function_type() )
{
    typedef typename Int_::derivative_type      derivative_type;
    typedef std::pair<size_t, derivative_type>  pair_type;
    
    bool solved = !weights.empty();
    bool save = !file_name.empty();
    
    typedef boost::shared_ptr<Int_> Int_ptr;
    Int_ptr int_;
    
    if (solved) {
        int_.reset(new Int_(points, times, weights, fun, verbose));
    } else {
        int_.reset(new Int_(points, times, fun, verbose));
        if (do_export_weights) 
            save_rbf(points, times, int_->weights(), 
                     kernel_name, file_name +  "-" + kernel_name + ".nrrd");
    }

    derivative_type d = int_->derivative(points[0]);
    if (verbose) std::cout << "gradient norm at source: " << nvis::norm(nvis::vec2(d[0][0], d[1][0])) << '\n';
    
    if (do_reconstruct) return reconstruct(*int_, resolution);
    
    return (float*)0;
}

int main(int argc, char* argv[])
{
    std::string input_name="", output_name="", weights_name="", path="";
    me = argv[0];
    bounds.reset();
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage();
        } else if (arg == "-i" || arg == "--input") {
            if (i == argc-1) {
                usage("missing input filename");
            }
            input_name = argv[++i];
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                usage("missing output filename");
            }
            output_name = argv[++i];
            do_reconstruct = true;
        } else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-4) {
                usage("missing reconstruction bounds");
            }
            bounds.min()[0] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
        } else if (arg == "-w" || arg == "--weights") {
            if (i == argc-1) {
                usage("missing weights filename");
            }
            weights_name = argv[++i];
            do_export_weights = true;
        } else if (arg == "-p" || arg == "--path") {
            if (i == argc-1) {
                usage("missing data path");
            }
            path = argv[++i];
        } else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) {
                usage("missing sampling resolution");
            }
            resolution[0] = atoi(argv[++i]);
            resolution[1] = atoi(argv[++i]);
        } else if (arg == "-k" || arg == "--kernel") {
            if (i == argc-1) {
                usage("missing kernel name");
            }
            kernel_name = argv[++i];
        } else if (arg == "-n" || arg == "--order") {
            if (i == argc-1) {
                usage("missing polynomial order");
            }
            poly_order = atoi(argv[++i]);
        } else if (arg == "-a" || arg == "--accuracy") {
            if (i == argc-1) {
                usage("missing accuracy level");
            }
            accuracy = argv[++i];
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else {
            usage("unrecognized argument: " + arg);
        }
    }
    if (input_name.empty()) {
        usage("Missing input filename");
    }
    if (!do_export_weights && !do_reconstruct) {
        usage("Neither RBF weights nor actual reconstruction requested: why bother?");
    }
    
    if (path != "") {
        if (*path.rbegin() != '/') {
            path.push_back('/');
        }
        input_name = path + input_name;
        if (do_reconstruct)
            output_name = path + output_name;
        if (do_export_weights)
            weights_name = path + weights_name;
    }
    if (verbose) {
        std::cout << "input:  " << input_name << '\n';
        if (do_reconstruct) std::cout << "output: " << output_name << '\n';
        if (do_export_weights) std::cout << "weights: " << weights_name << '\n';
    }
    weights_name = xavier::get_basename(weights_name);
    
    size_t number_of_samples = resolution[0]*resolution[1];
    
    std::vector<point_type> points;
    std::vector<data_type> values, weights;
    bool import_weights = false;
    // determine type of file to import
    std::string ext = xavier::get_extension(input_name);
    nvis::bbox2 in_bounds;
    if (ext == "nrrd" || ext == ".nhdr") {
        in_bounds = read_nrrd(points, values, weights, input_name);
        import_weights = !weights.empty();
    }
    else {
        std::vector<double> dummy;
        in_bounds = xavier::maarten::read_text(points, values, dummy, input_name);
    }    
    source = points[0];
    if (verbose) {
        point_distribution(points);
    }
    
    // if no bounds were specified, use input bounds for reconstruction
    // by default
    if (nvis::any(bounds.min() >= bounds.max()) && do_reconstruct) {
        std::cout << "Warning: No bounds were provided for the reconstruction.\n"
                  << "         Using receivers' bounds by default.\n";
        bounds = in_bounds;
    }
    
    if (verbose) {
        std::cout << "source located at " << source << '\n';
        if (do_reconstruct) std::cout << "receivers' bounds: " << in_bounds << '\n';
    }
    
    if (import_weights && do_export_weights) {
        std::cout << "WARNING: odd request to import and export RBF weights\n";
        std::cout << "         export request will be ignored\n";
        do_export_weights = false;
    }
    
    if (verbose) {
        std::cout << "imported " << points.size() << " data points\n";
        if (do_reconstruct) std::cout << "reconstruction bounding box: " << bounds << '\n';
        if (import_weights) std::cout << "weights have been imported\n";
    }
    
    // sanity check for accuracy setting
    if (accuracy == "low" || accuracy == "fast") {
        accuracy = "low";
    } else if (accuracy == "medium" || accuracy == "regular" || accuracy == "average") {
        accuracy = "medium";
    } else if (accuracy == "high" || accuracy == "slow") {
        accuracy == "high";
    } else {
        usage("unrecognized accuracy setting:" + accuracy);
    }
    
    point_type domain_size = bounds.size();
    
    float* result;
    if (kernel_name == "r") {
        if (accuracy == "low") {
            switch (poly_order) {
            case 0: result = reconstruct<fast_linear_rbf_type>(points, values, weights, weights_name);
            break;
            case 1: result = reconstruct<fast_linear_pp1_rbf_type>(points, values, weights, weights_name);
            break;
            case 2: result = reconstruct<fast_linear_pp2_rbf_type>(points, values, weights, weights_name);
            break;
            case 3: result = reconstruct<fast_linear_pp3_rbf_type>(points, values, weights, weights_name);
            break;
            case 4: result = reconstruct<fast_linear_pp4_rbf_type>(points, values, weights, weights_name);
            break;
            case 5: result = reconstruct<fast_linear_pp5_rbf_type>(points, values, weights, weights_name);
            break;
            case 6: result = reconstruct<fast_linear_pp6_rbf_type>(points, values, weights, weights_name);
            break;
            default: usage("unsupported polynomial order");
            }
        } else if (accuracy == "medium") {
            switch (poly_order) {
            case 0: result = reconstruct<medium_linear_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 1: result = reconstruct<medium_linear_pp1_rbf_type>(points, values, weights, weights_name);
            break;                      
            case 2: result = reconstruct<medium_linear_pp2_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 3: result = reconstruct<medium_linear_pp3_rbf_type>(points, values, weights, weights_name);
            break;                      
            case 4: result = reconstruct<medium_linear_pp4_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 5: result = reconstruct<medium_linear_pp5_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 6: result = reconstruct<medium_linear_pp6_rbf_type>(points, values, weights, weights_name);
            break;
            default: usage("unsupported polynomial order");
            }
        } else {
            switch (poly_order) {
            case 0: result = reconstruct<slow_linear_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 1: result = reconstruct<slow_linear_pp1_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 2: result = reconstruct<slow_linear_pp2_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 3: result = reconstruct<slow_linear_pp3_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 4: result = reconstruct<slow_linear_pp4_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 5: result = reconstruct<slow_linear_pp5_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 6: result = reconstruct<slow_linear_pp6_rbf_type>(points, values, weights, weights_name);
            break;
            default: usage("unsupported polynomial order");
            }
        }
    } else if (kernel_name == "r3") {
        if (accuracy == "low") {
            switch (poly_order) {
            case 0: result = reconstruct<fast_cubic_rbf_type>(points, values, weights, weights_name);
            break;
            case 1: result = reconstruct<fast_cubic_pp1_rbf_type>(points, values, weights, weights_name);
            break;
            case 2: result = reconstruct<fast_cubic_pp2_rbf_type>(points, values, weights, weights_name);
            break;
            case 3: result = reconstruct<fast_cubic_pp3_rbf_type>(points, values, weights, weights_name);
            break;
            case 4: result = reconstruct<fast_cubic_pp4_rbf_type>(points, values, weights, weights_name);
            break;
            case 5: result = reconstruct<fast_cubic_pp5_rbf_type>(points, values, weights, weights_name);
            break;
            case 6: result = reconstruct<fast_cubic_pp6_rbf_type>(points, values, weights, weights_name);
            break;
            default: usage("unsupported polynomial order");
            }
        } else if (accuracy == "medium") {
            switch (poly_order) {
            case 0: result = reconstruct<medium_cubic_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 1: result = reconstruct<medium_cubic_pp1_rbf_type>(points, values, weights, weights_name);
            break;                      
            case 2: result = reconstruct<medium_cubic_pp2_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 3: result = reconstruct<medium_cubic_pp3_rbf_type>(points, values, weights, weights_name);
            break;                      
            case 4: result = reconstruct<medium_cubic_pp4_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 5: result = reconstruct<medium_cubic_pp5_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 6: result = reconstruct<medium_cubic_pp6_rbf_type>(points, values, weights, weights_name);
            break;
            default: usage("unsupported polynomial order");
            }
        } else {
            switch (poly_order) {
            case 0: result = reconstruct<slow_cubic_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 1: result = reconstruct<slow_cubic_pp1_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 2: result = reconstruct<slow_cubic_pp2_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 3: result = reconstruct<slow_cubic_pp3_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 4: result = reconstruct<slow_cubic_pp4_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 5: result = reconstruct<slow_cubic_pp5_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 6: result = reconstruct<slow_cubic_pp6_rbf_type>(points, values, weights, weights_name);
            break;
            default: usage("unsupported polynomial order");
            }
        }
    } else if (kernel_name == "r5" ) {
        if (accuracy == "low") {
            switch (poly_order) {
            case 0: result = reconstruct<fast_quintic_rbf_type>(points, values, weights, weights_name);
            break;
            case 1: result = reconstruct<fast_quintic_pp1_rbf_type>(points, values, weights, weights_name);
            break;
            case 2: result = reconstruct<fast_quintic_pp2_rbf_type>(points, values, weights, weights_name);
            break;
            case 3: result = reconstruct<fast_quintic_pp3_rbf_type>(points, values, weights, weights_name);
            break;
            case 4: result = reconstruct<fast_quintic_pp4_rbf_type>(points, values, weights, weights_name);
            break;
            case 5: result = reconstruct<fast_quintic_pp5_rbf_type>(points, values, weights, weights_name);
            break;
            case 6: result = reconstruct<fast_quintic_pp6_rbf_type>(points, values, weights, weights_name);
            break;
            default: usage("unsupported polynomial order");
            }
        } else if (accuracy == "medium") {
            switch (poly_order) {
            case 0: result = reconstruct<medium_quintic_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 1: result = reconstruct<medium_quintic_pp1_rbf_type>(points, values, weights, weights_name);
            break;                      
            case 2: result = reconstruct<medium_quintic_pp2_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 3: result = reconstruct<medium_quintic_pp3_rbf_type>(points, values, weights, weights_name);
            break;                      
            case 4: result = reconstruct<medium_quintic_pp4_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 5: result = reconstruct<medium_quintic_pp5_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 6: result = reconstruct<medium_quintic_pp6_rbf_type>(points, values, weights, weights_name);
            break;
            default: usage("unsupported polynomial order");
            }
        } else {
            switch (poly_order) {
            case 0: result = reconstruct<slow_quintic_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 1: result = reconstruct<slow_quintic_pp1_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 2: result = reconstruct<slow_quintic_pp2_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 3: result = reconstruct<slow_quintic_pp3_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 4: result = reconstruct<slow_quintic_pp4_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 5: result = reconstruct<slow_quintic_pp5_rbf_type>(points, values, weights, weights_name);
            break;                       
            case 6: result = reconstruct<slow_quintic_pp6_rbf_type>(points, values, weights, weights_name);
            break;
            default: usage("unsupported polynomial order");
            }
        }
    } else if (kernel_name.substr(0, 8) == "gaussian") {
        double sigma;
        if (kernel_name[8] == ':') {
            sigma = atof(kernel_name.substr(9).c_str());
        } else {
            usage("Syntax error in kernel definition: " + kernel_name);
        }
        gaussian_type gauss(1./(2.*sigma*sigma));
	if (verbose) std::cout << "parsed gaussian with sigma=" << sigma << '\n';
        if (accuracy == "low") {
            switch (poly_order) {
            case 0: result = reconstruct<fast_gaussian_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 1: result = reconstruct<fast_gaussian_pp1_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 2: result = reconstruct<fast_gaussian_pp2_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 3: result = reconstruct<fast_gaussian_pp3_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 4: result = reconstruct<fast_gaussian_pp4_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 5: result = reconstruct<fast_gaussian_pp5_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 6: result = reconstruct<fast_gaussian_pp6_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            default: usage("unsupported polynomial order");
            }
        } else if (accuracy == "medium") {
            switch (poly_order) {
            case 0: result = reconstruct<medium_gaussian_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 1: result = reconstruct<medium_gaussian_pp1_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 2: result = reconstruct<medium_gaussian_pp2_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 3: result = reconstruct<medium_gaussian_pp3_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 4: result = reconstruct<medium_gaussian_pp4_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 5: result = reconstruct<medium_gaussian_pp5_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 6: result = reconstruct<medium_gaussian_pp6_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            default: usage("unsupported polynomial order");
            }
        } else {
            switch (poly_order) {
            case 0: result = reconstruct<slow_gaussian_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 1: result = reconstruct<slow_gaussian_pp1_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 2: result = reconstruct<slow_gaussian_pp2_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 3: result = reconstruct<slow_gaussian_pp3_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 4: result = reconstruct<slow_gaussian_pp4_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 5: result = reconstruct<slow_gaussian_pp5_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            case 6: result = reconstruct<slow_gaussian_pp6_rbf_type>(points, values, weights, weights_name, gauss);
            break;
            default: usage("unsupported polynomial order");
            }
        }
    } else if (kernel_name.substr(0, 8) == "wendland") {
        double radius = 1;
        if (kernel_name[8] == ':') {
            radius = atof(kernel_name.substr(9).c_str());
        } else {
            usage("Syntax error in kernel definition: " + kernel_name);
        }
        if (verbose) {
            std::cout << "kernel radius = " << radius << '\n';
        }
        switch (poly_order) {
        case 0: result = reconstruct<wendland_rbf_type>(points, values, weights, weights_name, wendland_type(radius));
        break;
        case 1: result = reconstruct<wendland_pp1_rbf_type>(points, values, weights, weights_name, wendland_type(radius));
        break;
        case 2: result = reconstruct<wendland_pp2_rbf_type>(points, values, weights, weights_name, wendland_type(radius));
        break;
        case 3: result = reconstruct<wendland_pp3_rbf_type>(points, values, weights, weights_name, wendland_type(radius));
        break;
        case 4: result = reconstruct<wendland_pp4_rbf_type>(points, values, weights, weights_name, wendland_type(radius));
        break;
        case 5: result = reconstruct<wendland_pp5_rbf_type>(points, values, weights, weights_name, wendland_type(radius));
        break;
        case 6: result = reconstruct<wendland_pp6_rbf_type>(points, values, weights, weights_name, wendland_type(radius));
        break;
        default: usage("unsupported polynomial order");
        }
    } else if (kernel_name.substr(0, 9) == "tgaussian") {
        double sigma;
        double radius;
        if (kernel_name[9] != ':') {
            usage("Syntax error in kernel definition: " + kernel_name);
        }
        size_t next_colon = kernel_name.find(':', 10);
        if (next_colon == std::string::npos) {
            usage("Syntax error in kernel definition: " + kernel_name);
        }
        sigma = atof(kernel_name.substr(10, next_colon-10).c_str());
        radius = sigma*atof(kernel_name.substr(next_colon+1).c_str());
        if (verbose) {
            std::cout << "Selected kernel: truncated gaussian: sigma^2="
                      << sigma
                      << ", radius = " << radius << '\n';
        }
        truncated_gaussian_type tgauss(1./(2.*sigma*sigma), radius);
        switch (poly_order) {
        case 0: result = reconstruct<trunc_gaussian_rbf_type>(points, values, weights, weights_name, tgauss);
        break;
        case 1: result = reconstruct<trunc_gaussian_pp1_rbf_type>(points, values, weights, weights_name, tgauss);
        break;
        case 2: result = reconstruct<trunc_gaussian_pp2_rbf_type>(points, values, weights, weights_name, tgauss);
        break;
        case 3: result = reconstruct<trunc_gaussian_pp3_rbf_type>(points, values, weights, weights_name, tgauss);
        break;
        case 4: result = reconstruct<trunc_gaussian_pp4_rbf_type>(points, values, weights, weights_name, tgauss);
        break;
        case 5: result = reconstruct<trunc_gaussian_pp5_rbf_type>(points, values, weights, weights_name, tgauss);
        break;
        case 6: result = reconstruct<trunc_gaussian_pp6_rbf_type>(points, values, weights, weights_name, tgauss);
        break;
        default: usage("unsupported polynomial order");
        }
    } else {
        usage("Unrecognized kernel type: " + kernel_name);
    }
    
    if (do_reconstruct) { 
        xavier::nrrd_params<float, 3> params;
        params.sizes()[0] = 3;
        params.sizes()[1] = resolution[0];
        params.sizes()[2] = resolution[1];
        params.mins()[0] = 0; // airNaN() is broken
        params.mins()[1] = bounds.min()[0];
        params.mins()[2] = bounds.min()[1];
        point_type spacing = bounds.size() / point_type(resolution - nvis::ivec2(1,1));
        params.spacings()[0] = 1;
        params.spacings()[1] = spacing[0];
        params.spacings()[2] = spacing[1];
        xavier::writeNrrd(result, output_name, params);
        std::cout << output_name << " has been exported\n";
    }
    return 0;
}
