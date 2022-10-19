#include <flow/lavd.hpp>

#include <sstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <cctype>
#include <stdexcept>
#include <locale>
#include <iomanip>
#include <random>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/filesystem.hpp>

#include <data/field_wrapper.hpp>
#include <data/raster.hpp>
#include <format/filename.hpp>
#include <image/nrrd_wrapper.hpp>
#include <image/probe.hpp>
#include <misc/option_parse.hpp>
#include <misc/time_helper.hpp>
#include <misc/log_helper.hpp>
#include <vtk/vtk_utils.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>

// #include <vtkPoints.h>
// #include <vtkCellArray.h>
// #include <vtkPolyLine.h>
// #include <vtkPolyData.h>
// #include <vtkDataSetWriter.h>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace spurt::lavd;

typedef std::mt19937_64 random_generator_t;

std::ofstream log_file;

spurt::log::dual_ostream spurt::lavd::_log_(log_file, std::cout, 1, 0, 0, true);

std::string name_in, name_out;
std::string me;
size_t nsamples=1000000;
size_t nb_threads;
spurt::ivec3 up(2,2,2);
size_t n_used = 10;
value_t t_between_files=3*spurt::lavd::HOUR;
bbox_t domain, region;
size_t support_radius;

inline vec3 append(const vec2& v, value_t t) {
    return vec3(v[0], v[1], t);
}

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;
        
    xcl::option_traits 
            required_group(true, false, "Required Options"), 
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compare teem interpolation to trilinear interpolation over upsampled grid");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input info filename", required_group);
        parser.add_value("output", name_out, "Output filename", optional_group);
        parser.add_value("n", nsamples, nsamples, "Number of samples", optional_group);
        parser.add_tuple<3>("x", up, up, "Upsampling factor", optional_group);
        
        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

// global static data storage for parallel computation
Nrrd* current_velocity_volume=0;
Nrrd* current_vorticity_volume=0;

typedef spurt::image< vec3, 3, value_t, size_t > vector_image_t;
typedef spurt::image< value_t, 3, value_t, size_t > scalar_image_t;

std::shared_ptr< vector_image_t > vec_img;
std::shared_ptr< scalar_image_t > scl_img;
value_t current_t_min, current_t_max;

inline void nuke_current_volumes() {
    if (current_velocity_volume != NULL) {
        nrrdNuke(current_velocity_volume);
    }
    if (current_vorticity_volume != NULL) {
        nrrdNuke(current_vorticity_volume);
    }
}

void clean_exit(int n) {
    nuke_current_volumes();
    exit(n);
}

void import_data(const std::vector< std::string >& vel_filenames,
                 const std::vector< std::string >& vor_filenames)
{    
    // initialize velocity and vorticity volumes
    std::vector< Nrrd* > vel_tsteps(n_used);
    // std::vector< Nrrd* > vor_tsteps(n_used);
    for (size_t i=0; i<n_used; ++i) {
        vel_tsteps[i]=spurt::readNrrd(vel_filenames[i]);
        std::cout << "Imported " << vel_filenames[i] << '\n';
        // vor_tsteps[i]=spurt::readNrrd(vor_filenames[i]);
        // std::cout << "Imported " << vor_filenames[i] << '\n';
    }
    size_t last_time_step = n_used - 1;
    current_velocity_volume = spurt::lavd::create_nrrd_volume(vel_tsteps, 0, t_between_files);
    // current_vorticity_volume = spurt::lavd::create_nrrd_volume(vor_tsteps, 0, t_between_files);
    current_t_min = 0;
    current_t_max = last_time_step*t_between_files;
    
    vec_img = std::shared_ptr< vector_image_t >(upsample_vector(current_velocity_volume, up));
    // scl_img = upsample_scalar(current_vorticity_volume, lvec3(up));
    
    // vec_img->save_as_nrrd("upsampled_" + name_out );
    
    for (int i=0; i<vel_tsteps.size(); ++i) {
        nrrdNuke(vel_tsteps[i]);
        // nrrdNuke(vor_tsteps[i]);
    }
}

void check_nrrd_interpolation() {
    std::pair<double, double> minmax;
    spurt::bbox3 bounds;
    spurt::vec3 spacing;
    spurt::fixed_vector<size_t, 3> size;
    for (int i=0; i<3; ++i) {
        minmax = axis_bounds(current_velocity_volume->axis[i+1]);
        bounds.min()[i] = minmax.first;
        bounds.max()[i] = minmax.second;
        spacing[i] = current_velocity_volume->axis[i+1].spacing;
        size[i] = current_velocity_volume->axis[i+1].size;
    }
    size_t N = size[0]*size[1]*size[2];
    
    float* err = (float*)calloc(N, sizeof(float));
    
    NrrdVectorField velocity(current_velocity_volume);
    spurt::progress_display progress(false);
    progress.start(N, "nrrd interpolation check");
    for (size_t n=0; n<N; ++n) {
        spurt::fixed_vector<size_t, 3> coord=spurt::index_to_coord(n, size);
        spurt::vec3 x = bounds.min() + spurt::vec3(coord)*spacing;
        spurt::vec3 v0, v1;
        progress.update(n);
        try {
            velocity(x, v0);
            v1 = nrrd_value< spurt::vec3 >(current_velocity_volume, n);
            err[n] = spurt::norm(v0-v1);
        }
        catch(...) {
            std::cout << "unable to interpolate at " << x << '\n';
        }
    }       
    progress.stop(); 
    
    spurt::writeNrrdFromContainers(err, "nrrd_err_"+ name_out, size, spacing, bounds.min()); 
}

inline double unif01(std::uint_fast64_t n) {
    static const std::uint_fast64_t _max = std::numeric_limits<size_t>::max();
    return static_cast<double>(n)/static_cast<double>(_max);
}

int main(int argc, const char* argv[])
{
    using namespace spurt;
    
    me=argv[0];
    
    nb_threads=1;
    
#if _OPENMP
    nb_threads = omp_get_max_threads();
#endif
    
    initialize(argc, argv);
    
    std::vector<std::string> velocity_filenames;
    std::vector<std::string> vorticity_filenames;
    
    std::fstream info_file(name_in, std::ios::in);
    if (!info_file) {
        exit(1);
    }

    boost::filesystem::path p(name_in);
    std::string parent_dir=p.parent_path().string();
    while (info_file.good()) {
        std::string velname, vortname;
        value_t avg;
        info_file >> velname >> vortname >> avg;
        velocity_filenames.push_back(parent_dir + '/' + velname);
        vorticity_filenames.push_back(parent_dir + '/' + vortname);
    }
    info_file.close();
    assert(!velocity_filenames.empty());
            
    // compute bounds of entire domain
    spurt::vec2 input_spc;
    get_spatial_info(domain, input_spc, velocity_filenames[0], 1);
    region.min() = domain.min();
    region.max() = domain.max();
            
    support_radius = spurt::lavd::compute_support_radius();
    
    // initialize velocity and vorticity volumes
    import_data(velocity_filenames, vorticity_filenames);
    nrrdSave(("velvol_" + name_out).c_str(), current_velocity_volume, NULL);
    
    std::cout << "grid has bounds:\n" << vec_img->grid().bounds() << '\n';
    
    NrrdVectorField velocity(current_velocity_volume);
    // NrrdScalarField<3> vorticity(current_vorticity_volume);
    
    std::vector< shared_ptr< NrrdVectorField > > vf_copies(nb_threads);
    for (int i=0; i<nb_threads; ++i) {
        vf_copies[i] = shared_ptr< NrrdVectorField >(new NrrdVectorField(current_velocity_volume, std::string("velocity volume for thread #") + std::to_string(i), false));
    }
    
    std::cout << "region=\n" << region << '\n';
    std::cout << "current_t_max = " << current_t_max << '\n';
    
    spurt::progress_display progress(true);
    
    auto timer_start = std::chrono::high_resolution_clock::now();
    
    std::cout << "grid properties:\n";
    std::cout << "spacing: " << vec_img->grid().spacing() << '\n';
    std::cout << "bounds: " << vec_img->grid().bounds().min() << " " << vec_img->grid().bounds().max() << '\n';
    std::cout << "resolution: " << vec_img->grid().resolution() << '\n';
        
    
    std::uint_fast64_t seed = 314159265358979;
        
    if (name_out.empty()) {
        std::vector<value_t> err(nsamples, 0);
        
        random_generator_t randgen;
        
        progress.start(nsamples, "sampling velocity nrrd");

        randgen.seed(seed);
        #pragma omp parallel
        {
        #pragma omp for schedule(dynamic,1) 
        for (size_t n = 0 ; n < nsamples ; ++n) {
            #if _OPENMP
            const int thread=omp_get_thread_num();
            #else
            const int thread=0;
            #endif
            
            if (!thread) progress.update(n);
            
            spurt::vec3 p, v;
            p[0] = region.min()[0] + std::generate_canonical<double, 64>(randgen)*region.size()[0];
            p[1] = region.min()[1] + std::generate_canonical<double, 64>(randgen)*region.size()[1];
            p[2] = std::generate_canonical<double, 64>(randgen)*current_t_max;
                      
            try {
                velocity(p, v);
            } 
            catch (std::runtime_error& e) {
                std::cout << "Caught runtime_error while interpolating in velocity at " << p << '\n';
                std::cout << "error message: " << e.what() << '\n';
                exit(1);
            }
        }
        }
        progress.stop();
        
        progress.start(nsamples, "sampling velocity raster");
        
        randgen.seed(seed);
        #pragma omp parallel
        {
        #pragma omp for schedule(dynamic,1) 
        for (size_t n = 0 ; n < nsamples ; ++n) {
            #if _OPENMP
            const int thread=omp_get_thread_num();
            #else
            const int thread=0;
            #endif
            
            if (!thread) progress.update(n);
            
            spurt::vec3 p, v;
            p[0] = region.min()[0] + std::generate_canonical<double, 64>(randgen)*region.size()[0];
            p[1] = region.min()[1] + std::generate_canonical<double, 64>(randgen)*region.size()[1];
            p[2] = std::generate_canonical<double, 64>(randgen)*current_t_max;
            
            try {
                v = vec_img->value(p);
            } 
            catch (std::runtime_error& e) {
                std::cout << "Caught runtime_error while interpolating in v_img at " << p << '\n';
                std::cout << "error message: " << e.what() << '\n';
                exit(1);
            }
        }
        }
        progress.stop();
        
        progress.start(nsamples, "Compute relative error");
        randgen.seed(seed);
        size_t nvalid=0;
        for (size_t n=0; n<nsamples; ++n) {
            progress.update(n);
            
            spurt::vec3 p, v0, v1;
            p[0] = region.min()[0] + std::generate_canonical<double, 64>(randgen)*region.size()[0];
            p[1] = region.min()[1] + std::generate_canonical<double, 64>(randgen)*region.size()[1];
            p[2] = std::generate_canonical<double, 64>(randgen)*current_t_max;
            
            try {
                velocity(p, v0);
                v1 = vec_img->value(p);
                double norm = spurt::norm(v0);
                if (norm != 0) {
                    ++nvalid;
                    err[n] = spurt::norm(v1-v0)/norm;
                }
            } 
            catch (std::runtime_error& e) {
                std::cout << "Caught runtime_error while interpolating in v_img at " << p << '\n';
                std::cout << "error message: " << e.what() << '\n';
                exit(1);
            }
        }
        progress.stop();
    
        value_t max_err = *std::max_element(err.begin(), err.end());
        value_t avg_err = std::accumulate(err.begin(), err.end(), 0);
        avg_err /= static_cast<value_t>(nvalid);
        std::cout << "max error = " << max_err << '\n';
        std::cout << "mean error = " << avg_err << '\n';
        std::cout << "sorting " << nsamples << " values...\n";
        std::sort(err.begin(), err.end());
        size_t offset = nsamples - nvalid;
        std::cout << "error percentiles:\n";
        for (int i=0; i<10; ++i) {
            std::cout << spurt::number_to_rank(i) << " percentile: " << err[offset+i*nvalid/10] << '\n';
        }
        std::cout << "100th percentile: " << err.back() << '\n';
    }
    else {
        size_t N = nsamples*nsamples*nsamples;
        double* nrrd_data = (double*)std::calloc(3*N, sizeof(double));
        double* tril_data = (double*)std::calloc(3*N, sizeof(double));
        spurt::vec3* nrrd_vec = reinterpret_cast<spurt::vec3 *>(nrrd_data);
        spurt::vec3* tril_vec = reinterpret_cast<spurt::vec3 *>(tril_data);
        
        std::vector<double> rel_err(N, 0);
        
        
        spurt::vec2 orig = region.min();
        spurt::vec3 orig3d(orig[0], orig[1], 0);
        double dx = region.size()[0]/static_cast<double>(nsamples-1);
        double dy = region.size()[1]/static_cast<double>(nsamples-1);
        double dt = current_t_max/static_cast<double>(nsamples-1);
        spurt::vec3 step(dx, dy, dt);
        
        std::cout << "dx=" << dx << "\ndy=" << dy << "\ndt=" << dt << '\n';
        std::cout << "Sampling bounds: " 
            << "(" << orig[0] << ", " << orig[0]+(nsamples-1)*dx << ") x "
            << "(" << orig[1] << ", " << orig[1]+(nsamples-1)*dy << ") x "
            << "(0, " << (nsamples-1)*dt << ")\n";
                    
        spurt::fixed_vector<size_t, 3> res(nsamples);
        std::vector<spurt::vec3> pos(N);
        
        progress.start(N, "regular sampling nrrd");
        #pragma omp parallel
        {
        #pragma omp for schedule(dynamic,1)      
        for (size_t n = 0 ; n < N ; ++n) {            
            #if _OPENMP
            const int thread=omp_get_thread_num();
            #else
            const int thread=0;
            #endif
            
            if (!thread) progress.update(n);
            
            spurt::fixed_vector<size_t, 3> coord = spurt::index_to_coord(n, res);
            pos[n] = orig3d + spurt::vec3(coord)*step;
            try {
                (*vf_copies[thread])(pos[n], nrrd_vec[n]);
                nrrd_vec[n][2] = 0;
            }
            catch (std::runtime_error& e) {
                std::cout << "Caught runtime_error while interpolating in velocity at " << p << '\n';
                std::cout << "error message: " << e.what() << '\n';
                exit(1);
            }
        }
        }
        progress.stop();
        
        progress.start(nsamples*nsamples*nsamples, "regular sampling trilinear");
        #pragma omp parallel
        {
        #pragma omp for schedule(dynamic,1)      
        for (size_t n = 0 ; n < N ; ++n) {
            #if _OPENMP
            const int thread=omp_get_thread_num();
            #else
            const int thread=0;
            #endif
            
            if (!thread) progress.update(n);
                                    
            try {
                tril_vec[n] = vec_img->value(pos[n]);
                tril_vec[n][2] = 0;
            } catch (std::runtime_error& e) {
                std::cout << "Caught runtime_error while interpolating in v_img at " << p << '\n';
                std::cout << "error message: " << e.what() << '\n';
                exit(1);
            }
            
            double norm = spurt::norm(nrrd_vec[n]);
            
            if (norm != 0) {
                rel_err[n] = spurt::norm(nrrd_vec[n]-tril_vec[n])/norm;
            }
        }
        }
        progress.stop();
        
        std::sort(rel_err.begin(), rel_err.end());
        auto _begin = std::upper_bound(rel_err.begin(), rel_err.end(), 0);
        size_t nzero = std::distance(rel_err.begin(), _begin);
        double err_sum = std::accumulate(rel_err.begin(), rel_err.end(), 0);
        
        std::cout << std::setprecision(16);
        std::cout << "max error=" << rel_err.back() << "\n";
        std::cout << "mean error=" << err_sum/static_cast<double>(N-nzero) << "\n";
        std::cout << "median error=" << rel_err[(nzero+N)/2] << "\n";
                
        std::vector<size_t> __res(4);
        __res[0] = 3;
        __res[1] = nsamples;
        __res[2] = nsamples;
        __res[3] = nsamples;
    
        std::vector<double> __mins(4);
        __mins[0] = AIR_NAN;
        __mins[1] = region.min()[0];
        __mins[2] = region.min()[1];
        __mins[3] = 0;
    
        std::vector<double> __spc(4);
        __spc[0] = AIR_NAN;
        __spc[1] = region.size()[0]/static_cast<double>(nsamples-1);
        __spc[2] = region.size()[1]/static_cast<double>(nsamples-1);
        __spc[3] = current_t_max/static_cast<double>(nsamples-1);
    
        std::vector<int> __ctr(4);
        __ctr[0] = nrrdCenterUnknown;
        __ctr[1] = nrrdCenterNode;
        __ctr[2] = nrrdCenterNode;
        __ctr[3] = nrrdCenterNode;
    
        spurt::writeNrrdFromContainers(nrrd_data, "nrrd_" + name_out, 
                          __res, __spc, __mins, __ctr); 
        spurt::writeNrrdFromContainers(tril_data, "tril_" + name_out,
                          __res, __spc, __mins, __ctr);
                          
        std::cout << "relative error distribution:\n";
        size_t nvalid = N-nzero;
        for (int i=0 ; i<10; ++i) {
            std::cout << spurt::number_to_rank(i*10) << " percentile: " << rel_err[nzero+i*nvalid/10] << '\n';
        }
        std::cout << "100th percentile: " << rel_err.back() << '\n';
        
        if (false) check_nrrd_interpolation();
        
        nrrd_vec = 0;
        tril_vec = 0;
        std::free(nrrd_data);
        std::free(tril_data);
    }

    return 0;
}
