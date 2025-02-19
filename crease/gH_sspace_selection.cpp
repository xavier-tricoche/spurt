#include <vector>
#include <algorithm>
#include <string>
#include <iostream>
#include <tuple>
#include <regex>
#include <mutex>
#include <thread>

// Boost filesystem (c++17 compiler flaky with std::filesystem )
#include <boost/filesystem.hpp>
// Teem
#include <teem/nrrd.h>
// spurt
#include <math/types.hpp>
#include <math/stat.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <data/image.hpp>

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
std::atomic<size_t> tbb_progress_counter;
std::mutex update_progress_mutex;

using namespace spurt;
typedef double scalar_type;
typedef long size_type;
typedef small_vector<scalar_type, 3> pos_type;
typedef small_vector<size_type, 3> coord_type;
typedef small_vector<scalar_type, 3> vec_type;
typedef small_matrix<scalar_type, 3, 3> mat_type;

typedef raster_data<size_type, scalar_type, 3, scalar_type> scalar_raster;
typedef raster_data<size_type, scalar_type, 3, vec_type> vector_raster;
typedef raster_data<size_type, scalar_type, 3, mat_type> matrix_raster;

typedef image<size_type, scalar_type, 3, scalar_type, kernels::MitchellNetravaliBC> smooth_image;

Nrrd* import_nrrd(const std::string& filename)
{
    Nrrd* nin = nrrdNew();
    std::cout << "importing " << filename << '\n';
    if (nrrdLoad(nin, filename.c_str(), NULL))
    {
        char* err = biffGetDone(NRRD);
        std::cerr << "Thomas Schultz's ridge method: " << err << std::endl;
        exit(-1);
    }
    return nin;
}

scalar_raster import_nrrd_as_raster(const std::string& filename) 
{
    Nrrd *nin = import_nrrd(filename);
    return nrrd_utils::to_raster<size_type, scalar_type, 3, scalar_type>(nin, true);
}

scalar_type default_lmin_treshold(const scalar_raster& values) {
    scalar_type mode = 
        spurt::mode(values.begin(), values.end(), 1024, 
                    std::numeric_limits<scalar_type>::min(), scalar_type(0));
    std::cout << "mode is " << mode << '\n';
    return mode;
}

// to compute best scale:
// 1. compute scale normalized Hessian across a range of scales
// 2. compute the corresponding minor eigenvalues
// 3. determine at each point which scale yields the smallest minor eigenvalue
// 4. compute a gradient value at each point that matches the picked scale 

/* 
 scale space theory
 ------------------
 t: time of isotropic diffusion (described by a Gaussian)
 s: spatial scale  
 s = sqrt(t), t = s^2
 one-dimensional Gaussian: g(x,t) = 1/(sqrt(2pi*t))*exp(-x^2/(2t))
 in 3D: g(x,y,z,t) = g(x,t)*g(y,t)*g(z,t) = 
    1/sqrt(2pi*t)^3 * ( exp(-x^2/(2t)) * exp(-y^2/(2t)) * exp(-z^2/(2t)) = 
    1/sqrt(2pi*t)^3 * ( exp(-(x^2 + y^2 + z^2)/(2t) ) =
    1/(2pi*t)^1.5 * exp(-||x||^2/(2t))
 using spatial scale instead:
 g(x,y,z,t) = 1/(2pi)^{3/2}/s^3*exp(-||x||^2/{2*s^2})
*/

void ridge_strength(matrix_raster& hessian, 
                    scalar_raster& lmin,
                    scalar_raster& input, 
                    scalar_type scale) 
{
    hessian.initialize(0);

    // We'll use B-spline kernels to compute derivatives, like vprobe
    smooth_image volume(input, kernels::MitchellNetravaliBC(1,0)); 
    scalar_type scsq = scale*scale; // scale normalization
    if (scale < 1) scsq = 1;
    scalar_raster h00 = image_convolution(input, kernels::MitchellNetravaliBC(1,0), coord_type(2,0,0));
    size_type npts = input.grid().size();
    coord_type res = input.grid().resolution();
    ProgressDisplay progress;
    progress.begin(npts, "Computing ridge strength", 500);
    tbb_progress_counter = 0;
    tbb::parallel_for(tbb::blocked_range<size_type>(0, npts),
                       [&](tbb::blocked_range<size_type> r) {
        for (size_type n=r.begin(); n!=r.end(); ++n)
        {
            coord_type ijk = input.grid().coordinates(n);
            coord_type u(0);
            for (int d=0; d<3; ++d) { 
                if (ijk[d]==res[d]-1) {
                    --ijk[d];
                    u[d] = 1;
                }
            }
            hessian[n] = volume.second_derivative_in_voxel(ijk, u);
            hessian[n] *= scsq;
            vec_type lambdas;
            real_eigenvalues(lambdas, hessian[n]);
            lmin[n] = lambdas[2]; // eigenvalues are sorted in decreasing order
            std::unique_lock<std::mutex> lock(update_progress_mutex, std::defer_lock);
            if (lock.try_lock())
            {
                progress.update(tbb_progress_counter);
            }

            ++tbb_progress_counter;
        }
    });
    progress.end();
}

void sample_at_fixed_scale(scalar_raster& value, 
                           vector_raster& gradient,
                           scalar_raster& blurred,
                           const scalar_raster& best_scale,
                           const scalar_type& scale) 
{
    smooth_image volume(blurred, kernels::MitchellNetravaliBC(1,0));
    size_type npts = blurred.grid().size();
    coord_type res = blurred.grid().resolution();
    ProgressDisplay progress;
    progress.begin(npts, "Sampling image at scale " + std::to_string(scale), 500);
    tbb_progress_counter = 0;
    tbb::parallel_for(tbb::blocked_range<size_type>(0, npts),
                       [&](tbb::blocked_range<size_type> r) {
        for (size_type n=r.begin(); n!=r.end(); ++n)
        { 
            std::unique_lock<std::mutex> lock(update_progress_mutex, std::defer_lock);
            if (lock.try_lock())
            {
                progress.update(tbb_progress_counter);
            }
            ++tbb_progress_counter;

            if (best_scale[n] != scale) continue;
            coord_type ijk = blurred.grid().coordinates(n);
            pos_type u = 0;
            for (int d=0; d<3; ++d) {
                if (ijk[d] == res[d]-1) {
                    --ijk[d];
                    u[d] = 1;
                }
            }
            value[n] = volume.value_in_voxel(ijk, u);
            gradient[n] = volume.derivative_in_voxel(ijk, u);
            if (scale > 1) {
                gradient[n] *= scale;
            }
        }
    });
}

void sample_scale_space(scalar_raster& value, 
                        vector_raster& gradient,
                        const scalar_raster& best_scale,
                        const std::vector<std::pair<scalar_type, std::string>>& blur_info)
{
    for ( auto info : blur_info ) {
        scalar_type s = info.first;
        scalar_raster blurred = import_nrrd_as_raster(info.second);
        std::cout << "\nSampling image at scale " << s << '\n';
        sample_at_fixed_scale(value, gradient, blurred, best_scale, s);
    }
}

void img_min(matrix_raster& inout_hessian, 
             scalar_raster& inout_lmin, 
             scalar_raster& inout_scale,
             const matrix_raster& other_hessian,
             const scalar_raster& other_lmin,
             scalar_type other_scale) 
{
    size_type npts = inout_hessian.grid().size();

    ProgressDisplay progress;
    progress.begin(npts, "Updating best scale", 500);
    tbb_progress_counter = 0;
    tbb::parallel_for(tbb::blocked_range<size_type>(0, npts),
                       [&](tbb::blocked_range<size_type> r) {
        for (size_type n=r.begin(); n!=r.end(); ++n)
        {
            if (other_lmin[n] < inout_lmin[n]) {
                inout_lmin[n] = other_lmin[n];
                inout_hessian[n] = other_hessian[n];
                inout_scale[n] = other_scale;
            }
            
            std::unique_lock<std::mutex> lock(update_progress_mutex, std::defer_lock);
            if (lock.try_lock())
            {
                progress.update(tbb_progress_counter);
            }
            ++tbb_progress_counter;
        }
    });
    progress.end();
}

void select_scale(matrix_raster& hessian,
                  scalar_raster& lmin,
                  scalar_raster& best_scale,
                  const std::vector<std::pair<scalar_type, std::string>>& blur_info)
{
    matrix_raster other_hessian(hessian.grid());
    scalar_raster other_lmin(hessian.grid());
    best_scale.initialize(-1);
    lmin.initialize(1); // a positive value of lmin is not a ridge
    for ( auto info : blur_info ) {
        scalar_type s = info.first;
        scalar_raster blurred = import_nrrd_as_raster(info.second);
        std::cout << "processing scale " << s << '\n';
        spurt::timer t;
        ridge_strength(other_hessian, other_lmin, blurred, s);
        img_min(hessian, lmin, best_scale, other_hessian, other_lmin, s);
        std::cout << "Processing scale " << s << " took ";
        t.print_self();
        std::cout << '\n';
    }
}

std::vector<std::string> 
regex_match_files(const std::string& basename) {
    namespace fs = boost::filesystem;

    fs::path apath(basename);
    fs::path directory = apath.parent_path();
    fs::path search_path = apath.filename().string();
    std::string search_str = search_path.string();

    std::vector<std::string> matching_files;
    std::regex pattern(search_str);

    for (const auto& entry : fs::directory_iterator(directory)) {
        if (fs::is_regular_file(entry)) {
            std::string filename = entry.path().filename().string();
            if (std::regex_match(filename, pattern)) {
                matching_files.push_back(entry.path().string());
            }
        }
    }
    return matching_files;
}

void import_scales(std::vector<std::pair<scalar_type, std::string> >& blurred, 
                 const std::string& basename) {
    std::cout << "basename is " << basename << '\n';
    typedef std::pair<scalar_type, std::string> scale_t;
    blurred.clear();
    std::vector<std::string> all_filenames = regex_match_files(basename);
    for ( const std::string& name : all_filenames ) {
        Nrrd *nin = import_nrrd(name);
        char* scale_as_str = nrrdKeyValueGet(nin, "scale");
        blurred.push_back(std::make_pair(std::stof(scale_as_str), name));
        nrrdNuke(nin);
    }

    std::sort(blurred.begin(), blurred.end(), [&](const scale_t& p, const scale_t& q) {
        return p.first < q.first;
    });
}

int main(int argc, char* argv[]) 
{
    std::string input_name, blurred_basename, output_basename;
    vec3 _scales;
    std::vector<scalar_type> scales;
    scalar_type cutoff;
    bool verbose = false;

    namespace cl = spurt::command_line;
    cl::option_traits
        required_group(true, false, "Required Options"),
        positional_group(true, true, "Positional Group"),
        optional_group(false, false, "Optional Group");

    cl::option_parser parser(argv[0],
                             "Automatically determine best scale for ridge surface extraction and measure value, gradient, Hessian, and ridge strength at that scale for further processing");
    try 
    {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("input", blurred_basename, "Blurred input volume basename (regex)", required_group);
        parser.add_value("output", output_basename, "Output basename", required_group);
        parser.add_flag("verbose", verbose, "Verbose output", optional_group);
        parser.parse(argc, const_cast<const char**>(argv));
    }
    catch (std::runtime_error &e)
    {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    catch (std::exception &e)
    {
        std::cerr << "ERROR(2): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options enteredso far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }

    std::vector<std::pair<scalar_type, std::string>> blurred_info;
    import_scales(blurred_info, blurred_basename);

    assert(blurred_info.size() > 0);
    scalar_raster first = import_nrrd_as_raster(blurred_info[0].second);
    matrix_raster best_hessian(first.grid());
    vector_raster best_gradient(first.grid());
    scalar_raster best_scale(first.grid());
    scalar_raster best_value(first.grid());
    scalar_raster best_strength(first.grid());

    std::cout << "blurred input volumes are:\n";
    for ( auto bi : blurred_info ) { 
        std::cout << "sigma: " << bi.first << ", name: " << bi.second << '\n';
    }
    
    select_scale(best_hessian, best_strength, best_scale, blurred_info);

    sample_scale_space(best_value, best_gradient,
                       best_scale, blurred_info);

    std::cout << "Exporting data\n";

    save_as_nrrd(output_basename + "_scsp_value.nrrd", best_value);
    save_as_nrrd(output_basename + "_scsp_strength.nrrd", best_strength);
    save_as_nrrd(output_basename + "_scsp_scale.nrrd", best_scale);
    save_as_nrrd(output_basename + "_scsp_gradient.nrrd", best_gradient);
    save_as_nrrd(output_basename + "_scsp_hessian.nrrd", best_hessian, small_vector<scalar_type, 9>(0));

    return 0;
}