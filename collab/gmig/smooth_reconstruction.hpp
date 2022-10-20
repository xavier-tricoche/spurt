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
// spurt
#include <math/RBF.hpp>
#include <math/RBFbasis.hpp>
#include <image/nrrd_wrapper.hpp>
#define BPO_WRAPPER_IS_BROKEN
#ifndef BPO_WRAPPER_IS_BROKEN
#   include <misc/option_parse.hpp>
#endif
#include <format/format.hpp>
#include <boost/shared_ptr.hpp>
#include "maarten_utils.hpp"
#include "smooth_reconstruction.hpp"

// global variables
nvis::ivec2 resolution(800, 800);
nvis::bbox2 bounds;
bool verbose = false;
std::string kernel_name = "r3";

namespace xrbf = spurt::RBF;

template<typename _Interpolator>
void spurt::gmig::check_solution(const _Interpolator& f) {
    typedef typename _Interpolator::point_type   point_type;
    typedef typename _Interpolator::data_type    data_type;
    
    const std::vector<point_type>& points = f.points();
    const std::vector<data_type>& times = f.data();
    
    double mean=0, max=0;
    for (size_t n=0 ; n<points.size() ; ++n) {
        double err = fabs(times[n][0]-f(points[n])[0]);
        mean += err;
        max = std::max(max, err);
    }
    mean /= points.size();
    std::cout << "mean error = " << mean << ", max error = " << max << '\n';
}

template<typename _Interpolator>
float* _reconstruct(const _Interpolator& interpolator,
                    bool do_gradient=false) {
    nvis::vec2 spacing = bounds.size() / nvis::vec2(resolution - nvis::ivec2(1, 1));
    double eps = 1.0e-6*nvis::norm(spacing);
    nvis::vec2 dx(eps, 0);
    nvis::vec2 dy(0, eps);
    size_t number_of_samples = resolution[0]*resolution[1];
    size_t val_per_sample = 1;
    if (do_gradient) val_per_sample += 2;
    float* result = (float*)calloc(val_per_sample*number_of_samples, sizeof(float));
    
    size_t number_of_threads = 1;
    #ifdef _OPENMP
    number_of_threads = omp_get_max_threads();
    #endif
    
    if (verbose)
        std::cout << "there are " << number_of_threads << " threads available\n";
    
    size_t counter = 0;
    size_t nb_jobs = number_of_samples;
    nvis::timer _timer;
    
    progress_message msg(nb_jobs, "interpolations");
    std::vector<double> g_err(number_of_threads, 0);
    std::vector<double> g_rel(number_of_threads, 0);
        
    typedef typename _Interpolator::derivative_type derivative_type;
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (size_t n=0 ; n<number_of_samples ; ++n) {
            int thread = 0;
            #if _OPENMP
            thread = omp_get_thread_num();
            #endif
            
            size_t base_id = n*val_per_sample;
            
            int i = n%resolution[0];
            int j = n/resolution[0];
            nvis::vec2 x = bounds.min() + nvis::vec2(i, j)*spacing;
            result[base_id] = interpolator(x)[0];
            if (do_gradient) {
                derivative_type g = interpolator.derivative(x);
                result[base_id+1] = g[0][0];
                result[base_id+2] = g[1][0];
                nvis::vec2 rbf_g(g[0][0], g[1][0]);
                nvis::vec2 cd_g;
                cd_g[0] = interpolator(x+dx)[0] - interpolator(x-dx)[0];
                cd_g[1] = interpolator(x+dy)[0] - interpolator(x-dy)[0];
                cd_g /= 2.*eps;
                g_err[thread] += nvis::norm(cd_g-rbf_g);
                g_rel[thread] += nvis::norm(cd_g-rbf_g)/nvis::norm(cd_g);
            }
            
            if (verbose && !thread) {
                double elapsed = _timer.elapsed();
                std::cout << msg(n, elapsed) << std::flush;
            }
        }
    }
    std::cout << "\nRBF reconstruction completed in " << _timer.elapsed() << " seconds ("
    << (float)number_of_samples/_timer.elapsed() << " Hz)\n";
    
    if (verbose) {
         check_solution(interpolator);
         double total_err = 0;
         double total_rel = 0;
         for (size_t i=0 ; i<number_of_threads ; ++i) {
             total_err += g_err[i];
             total_rel += g_rel[i];
         }
         total_err /= (double)number_of_samples;
         total_rel /= (double)number_of_samples;
         std::cout << "average gradient error = " << total_err << '\n';
         std::cout << "average relative error = " << total_rel << '\n';
    }
    
    return result;
}

template<typename Int_>
float* spurt::gmig::
reconstruct(const std::vector<nvis::vec2>& points,
            const std::vector<nvis::vec1>& times,
            const std::vector<nvis::vec1>& weights,
            const std::string& file_name,
            bool do_gradient = false,
            const typename Int_::function_type& fun = 
                  typename Int_::function_type() ) {
    bool solved = !weights.empty();
    bool save = !file_name.empty();
    
    typedef boost::shared_ptr<Int_> Int_ptr;
    Int_ptr int_;
    
    if (solved)
        int_.reset(new Int_(points, times, weights, fun, verbose));
    else {
        int_.reset(new Int_(points, times, fun, verbose));
        if (save) 
            save_rbf(points, times, int_->weights(), 
                     kernel_name, file_name + "-" + kernel_name + ".nrrd");
    }
    return _reconstruct<Int_>(*int_, do_gradient);
}