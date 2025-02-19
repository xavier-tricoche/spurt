#include <iostream>
#include <mutex>
#include <set>
#include <string>
#include <vector>

#include <data/image.hpp>
#include <math/types.hpp>
#include <misc/progress.hpp>
#include <misc/option_parse.hpp>
#include <utils/functions.hpp>
#include <format/filename.hpp>

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
std::atomic<size_t> tbb_progress_counter;
std::atomic<size_t> tbb_total_time;

using namespace spurt;

typedef long size_type;
typedef double scalar_type;
typedef lvec3 coord_type;
typedef vec3 pos_type;
typedef vec3 vector_type;
typedef mat3 jacobian_type;
typedef vec3 gradient_type;
typedef mat3 hessian_type;
typedef bbox3 bounds_type;

bool sequential = false;

typedef spurt::raster_grid<size_type, scalar_type, 3, coord_type, pos_type> grid_type;

template<typename T, typename Enable=void>
struct _is_valid {};

template<typename T>
struct _is_valid<T, typename std::enable_if<std::is_scalar<T>::value>::type> {
    bool operator()(const T& val) { return !( std::isnan(val) || std::isinf(val) ); }
};

template<typename T>
struct _is_valid<T, typename std::enable_if<!std::is_scalar<T>::value>::type> {
    bool operator()(const T& val) {
        return !spurt::any(spurt::isinvalid(val));
    }
};

template<typename Value_>
using dataset_type = spurt::raster_data<size_type, scalar_type, 3, Value_, coord_type, pos_type>;

template<typename Value_, typename Kernel_>
using image_type = spurt::image<size_type, scalar_type, 3, Value_, Kernel_, coord_type, pos_type>;

template<typename Value_>
void fill_image_impl(dataset_type<Value_>& ds, const std::function<Value_(const pos_type&)>& op)
{
    std::cout << "Filling an image with " << ds.size() << " entries\n";
    spurt::ProgressDisplay progress;
    
    std::mutex mtx;
    tbb_progress_counter = 0;
    
    if (!sequential) {
        progress.begin(ds.size(), "Sampling");
        tbb::parallel_for(tbb::blocked_range<int>(0, ds.size()),
                           [&](tbb::blocked_range<int> r) {
            for (int i=r.begin(); i!=r.end(); ++i)
            {
                ds[i] = op(ds.grid()[i]);
                ++tbb_progress_counter;
                std::unique_lock<std::mutex> mlock(mtx, std::try_to_lock);
                if (mlock) {
                    progress.update(tbb_progress_counter);
                }
            }
        });
        progress.end();
    }
    else {
        for (int i=0; i!=ds.size(); ++i)
        {
            ds[i] = op(ds.grid()[i]);
        }
    }
}

template<typename Value_>
void fill_image(dataset_type<Value_>& ds, const std::string& fname) {}

template<>
void fill_image<scalar_type>(dataset_type<scalar_type>& ds, const std::string& fname) 
{
    if (fname == "marshner" || fname == "marshner_lobb")
    {
        functions3d::marshner_lobb ml;
        fill_image_impl<scalar_type>(ds, [&](const pos_type& p) -> scalar_type { return ml(p); });
    }
    else if (fname == "ball")
    {
        fill_image_impl<scalar_type>(ds, [&](const pos_type& p) -> scalar_type { return functions3d::ball(p); });
    }
    else if (fname =="sphere")
    {
        fill_image_impl<scalar_type>(ds, [&](const pos_type& p) -> scalar_type { return functions3d::sphere(p); });
    }
    else if (fname == "ellipsoid")
    {
        fill_image_impl<scalar_type>(ds, [&](const pos_type& p) -> scalar_type { return functions3d::ellipsoid(p); });
    }
    else if (fname =="torus")
    {
        fill_image_impl<scalar_type>(ds, [&](const pos_type& p) -> scalar_type { return functions3d::torus(p); });
    }
    else if (fname == "helix")
    {
        fill_image_impl<scalar_type>(ds, [&](const pos_type& p) -> scalar_type { return functions3d::helix(p); });
    }
    else
    {
        std::cerr << "Unrecognized function name \"" << fname << "\"\n";
        exit(1);
    }
}

template<>
void fill_image<vec3>(dataset_type<vector_type>& ds, const std::string& fname) 
{
    if (fname == "abc" || fname == "ABC") {
        spurt::functions3d::ABC_field abc;
        fill_image_impl<vec3>(ds, [&](const pos_type& p) -> vec3 { return abc.evaluate(p); });
    }
    else
    {
        std::cerr << "unrecognized function name for vector valued function: \"" << fname << "\"\n";
        exit(1);
    }
}


template<typename Kernel>
void test_image(const coord_type& res, const bounds_type& bounds, 
                const std::string& fname, const std::string& outname, float upres=2) 
{
    grid_type ingrid(res, bounds);
    grid_type outgrid(coord_type(floor(upres*res)), bounds);
    
    typedef image_type<vector_type, Kernel> vec_image_type;
    typedef typename vec_image_type::first_derivative_type vec_derivative_type;
    typedef typename vec_image_type::second_derivative_type vec_second_derivative_type;
    
    typedef image_type<scalar_type, Kernel> scl_image_type;
    typedef typename scl_image_type::first_derivative_type scl_derivative_type;
    typedef typename scl_image_type::second_derivative_type scl_second_derivative_type;
    
    std::string basename = spurt::filename::remove_extension(outname);
    std::string origname = basename + "_original.nrrd";
    std::ostringstream oss;
    oss << basename << "_x" << upres << "_";
    std::string valname  =  oss.str() + "f.nrrd";
    std::string dername  =  oss.str() + "gv.nrrd";
    std::string ddername =  oss.str() + "H.nrrd";
    
    std::cout << "original grid resolution: " << ingrid.resolution() << '\n';
    std::cout << "output grid resolution: " << outgrid.resolution() << '\n';
    if (fname == "ABC" || fname == "abc")
    {
        dataset_type<vector_type> indata(ingrid, vector_type(0));
        dataset_type<vector_type> outvalues(outgrid, vector_type(0));
        dataset_type<vec_derivative_type> outderiv(outgrid);
        fill_image<vector_type>(indata, "ABC");
        vec_image_type img(indata);
        fill_image_impl<vector_type>(outvalues, [&](const pos_type& p) { return img.value(p); });
        fill_image_impl<vec_derivative_type>(outderiv, [&](const pos_type& p) { return img.derivative(p); });
        spurt::save_as_nrrd(valname, outvalues);
        spurt::save_as_nrrd(origname, indata);
        spurt::save_as_nrrd<size_type, scalar_type, 3, vec_derivative_type, coord_type, pos_type, jacobian_type>(dername, outderiv);
    }
    else {
        dataset_type<scalar_type> indata(ingrid, scalar_type(0));
        dataset_type<scalar_type> outvalues(outgrid, scalar_type(0));
        dataset_type<scl_derivative_type> outderiv(outgrid, scl_derivative_type(0));
        dataset_type<scl_second_derivative_type> out2ndderiv(outgrid, scl_second_derivative_type(0));
        fill_image<scalar_type>(indata, fname);
        scl_image_type img(indata);
        fill_image_impl<scalar_type>(outvalues, [&](const pos_type& p) { return img.value(p); });
        fill_image_impl<scl_derivative_type>(outderiv, [&](const pos_type& p) { return img.derivative(p); });
        fill_image_impl<scl_second_derivative_type>(out2ndderiv, [&](const pos_type& p) { return img.second_derivative(p); });
        // fill_image_impl<scl_second_derivative_type>(out2ndderiv, [&](const pos_type& p) { return mat3(p[0]); });
        spurt::save_as_nrrd(valname, outvalues);
        spurt::save_as_nrrd(origname, indata);
        spurt::save_as_nrrd<size_type, scalar_type, 3, scl_derivative_type, coord_type, pos_type, gradient_type>(dername, outderiv);
        spurt::save_as_nrrd<size_type, scalar_type, 3, scl_second_derivative_type, coord_type, pos_type, hessian_type>(ddername, out2ndderiv);
    }
}

int main(int argc, const char* argv[])
{
    std::string dataname, outname, kernel;
    svec3 res;
    float upsample;
    vec3 bmin, bmax;

    namespace cl = command_line;

    cl::option_traits
        required_group(true, false, "Required Options"),
        positional_group(true, true, "Positional Group"),
        optional_group(false, false, "Optional Group");

    cl::option_parser parser(argv[0],
        "Test image class and its interpolation kernels");

    try {
        // parser.use_default_symbol();
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("data", dataname, "marshner", "Name of function to reconstruct", optional_group);
        parser.add_value("res", res, {64, 64, 64}, "Resolution", optional_group);
        parser.add_value("kernel", kernel, "linear", "Name of interpolation kernel", optional_group);
        parser.add_value("upsample", upsample, 2, "Upsampling coefficient", optional_group);
        parser.add_value("output", outname, "Output filename", required_group);
        parser.add_value("min", bmin, {-1, -1, -1}, "Min bounds", optional_group);
        parser.add_value("min", bmax, {1, 1, 1}, "Max bounds", optional_group);

        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    catch(std::exception& e) {
        std::cerr << "ERROR(2): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    
    bounds_type bounds(bmin, bmax);
    std::cout << "Bounds are " << bounds << '\n';
    
    if (kernel == "linear")
        test_image<spurt::kernels::Linear>(res, bounds, dataname, outname, upsample);
    else if (kernel == "mnbc")
        test_image<spurt::kernels::MitchellNetravaliBC>(res, bounds, dataname, outname, upsample);
    else if (kernel == "gauss" || kernel == "gaussian")
        test_image<spurt::kernels::Gaussian>(res, bounds, dataname, outname, upsample);
    else {
        std::cerr << "Unrecognized kernel name \"" << kernel << "\"\n";
        exit(1);
    }
    
    return 0;
}