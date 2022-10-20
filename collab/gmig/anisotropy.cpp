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
#include <data/raster.hpp>
#include <math/RBF.hpp>
#include <math/RBFbasis.hpp>
#include <image/nrrd_wrapper.hpp>
#include <format/format.hpp>
#include "data_IO.hpp"
#include "reconstruct.hpp"
#include "utils.hpp"
// boost
#include <boost/shared_ptr.hpp>

const double TWO_PI = 2.*M_PI;

// global variables
nvis::ivec2 __resolution(200, 200);
nvis::bbox2 __bounds;
bool        __verbose = false;
std::string __kernel_name = "r3";
size_t      __number_of_bins = 36;

using namespace spurt::gmig;

typedef std::vector<nvis::vec2>      value_type;
typedef spurt::raster2d<value_type> dataset_type;
typedef dataset_type::grid_type      grid_type;

typedef spurt::RBF::gaussian_function<double> gauss_t;
typedef spurt::RBF::wendland_function<double> wendland_t;
typedef spurt::RBF::truncated_gaussian_function<double> trgauss_t;

// map [-pi, pi] angles to [0, 2pi]
inline double convert_angle(double a) {
    if (a >= 0) return a;
    else return a + TWO_PI;
}

void vec2bins(double* out, const value_type& v) {
    static const double dt = TWO_PI/__number_of_bins;
    std::vector<size_t> counter(__number_of_bins, 0);
    for (size_t i=0 ; i<v.size() ; ++i) {
        const nvis::vec2& g = v[i];
        double t = convert_angle(atan2(g[1], g[0]));
        int j = floor(t/dt);
        out[j] += nvis::norm(g);
        counter[j]++;
    }
    for (size_t i=0 ; i<__number_of_bins ; ++i) {
        if (counter[i] > 1)
            out[i] /= static_cast<double>(counter[i]);
    }
}

template<typename _Int>
void anisotropy(const _Int& interpolator, 
                dataset_type& raster,
                const nvis::vec2& source)
{
    typedef typename _Int::derivative_type  derivative_type;
    typedef dataset_type::grid_type         grid_type;
    typedef grid_type::coord_type           coord_type;
    
    const grid_type& grid = raster.grid();
    
    const size_t number_of_samples = grid.size();
        
    size_t number_of_threads = 1;
#ifdef _OPENMP
    number_of_threads = omp_get_max_threads();
#endif
    size_t counter = 0;
    size_t nb_jobs = number_of_samples;
    nvis::timer _timer;
    
    progress_message msg(nb_jobs, "interpolations");
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (int n=0 ; n<number_of_samples ; ++n) {
            coord_type coord = grid.coordinates(n);
            
            nvis::vec2 x = grid(coord);            
            derivative_type nabla = interpolator.derivative(x);
            nvis::vec2 g(nabla[0], nabla[1]); // dt/drx, dt/dry
            g /= nvis::inner(g, g);
            raster(coord).push_back(g);
            
            if (__verbose) {
                int thread = 0;
#if _OPENMP
                thread = omp_get_thread_num();
#endif
                if (thread == 0) {
                    size_t my_counter = counter;
                    double elapsed = _timer.elapsed();
                    std::cout << msg(my_counter, elapsed) << std::flush;
                }
            }
        }
    }
    
    std::cout << "\nRBF reconstruction completed in " << _timer.elapsed() 
              << " seconds ("
              << (float)nb_jobs/_timer.elapsed() 
              << " Hz)\n";
}

void usage(const std::string&);

std::string me;
void usage(const std::string& message="")
{
    if (!message.empty()) {
        std::cerr << "ERROR: " << message << '\n';
    }
    std::cout << '\n'
              << "DESCRIPTION: Compute anisotropy at vertices of a raster grid\n"
              << "based on the smooth reconstruction of the mapping from\n"
              << "(x_source, y_source) x (x_receiver, y_receiver) to the\n"
              << "associated travel time.\n"
              << '\n'
              << "USAGE: " << me << " [parameters] [options]\n"
              << '\n'
              << "PARAMETERS:\n"
              << " -i | --input <string>      File containing the names of all RBF reconstruction info files\n"
              << " -o | --output <string>     Output file name\n"
              << " -a | --area <float>x4      Area to process: lon_min, lon_max, lat_min, lat_max\n"
              << '\n'                           
              << "OPTIONS:\n"                   
              << " -h | --help                Print this information\n"
              << " -p | --path <string>       Path to prepend to all file names\n"
              << " -k | --kernel <string>     Reconstruction kernel (\"r3\", \"gaussian:\"<sigma>,\n"
              << "                            \"wendland:\"<radius>, \"truncgaussian:\"<sigma>:<radius>)\n"
              << " -r | --resolution <int>x2  Resampling resolution\n"
              << " -b | --bins <int>          Number of bins to discretize angular range\n"
              << " -v | --verbose             Turn on verbose mode\n"
              << std::endl;
    exit(!message.empty());
}

template<typename Int_>
void reconstruct(const std::vector<std::string>& filenames, 
                 dataset_type& raster,
                 const typename Int_::function_type& fun =
                       typename Int_::function_type() )
{
    typedef travel_time_interpolator<Int_> tt_intp_type;
    
    for (int i=0 ; i<filenames.size() ; ++i) {
        tt_intp_type intp(filenames[i], fun, __verbose);
        std::cout << "Processing source #" << i+1 << "/" 
                  << filenames.size() << " at " << intp.data().source << '\n';
        anisotropy(intp, raster, intp.data().source);
    }
}

int main(int argc, char* argv[])
{
    std::string input_name="", output_name="", path="";
    me = argv[0];
    __bounds.reset();
    bool relative_distance = false;
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
        } else if (arg == "-a" || arg == "--area") {
            if (i >= argc-4) usage("missing area information");
            __bounds.min()[0] = atof(argv[++i]);
            __bounds.max()[0] = atof(argv[++i]);
            __bounds.min()[1] = atof(argv[++i]);
            __bounds.max()[1] = atof(argv[++i]);
        } else if (arg == "-p" || arg == "--path") {
            if (i == argc-1) {
                usage("missing data path");
            }
            path = argv[++i];
        } else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) {
                usage("missing sampling resolution");
            }
            __resolution[0] = atoi(argv[++i]);
            __resolution[1] = atoi(argv[++i]);
        } else if (arg == "-k" || arg == "--kernel") {
            if (i == argc-1) {
                usage("missing kernel name");
            }
            __kernel_name = argv[++i];
        } else if (arg == "-b" || arg == "--bins") {
            if (i == argc-1) {
                usage("missing bin number");
            }
            __number_of_bins = atoi(argv[++i]);
        } else if (arg == "-v" || arg == "--verbose") {
            __verbose = true;
        } else {
            usage("unrecognized argument: " + arg);
        }
    }
    if (input_name.empty()) {
        usage("Missing input filename");
    }
    if (output_name.empty()) {
        usage("Missing output filename");
    }
    if (nvis::any(__bounds.min() > __bounds.max())) {
        usage("Missing area information");
    }
    
    if (path != "") {
        if (*path.rbegin() != '/') {
            path.push_back('/');
        }
        input_name = path + input_name;
        output_name = path + output_name;
    }
    if (__verbose) {
        std::cout << "input:  " << input_name << '\n';
        std::cout << "output: " << output_name << '\n';
        std::cout << "bounds: " << __bounds << '\n';
    }
    
    grid_type::coord_type size(__resolution[0], __resolution[1]);
    grid_type grid(size, __bounds);
    dataset_type raster(grid);
    
    nvis::vec2 domain_size = __bounds.size();
    
    // import file names
    std::vector<std::string> filenames;
    std::fstream in(input_name.c_str(), std::ios::in);
    while (!in.eof()) {
        std::string name = "invalid name";
        in >> name;
        if (name == "invalid name") break;
        filenames.push_back(name);
    }
    in.close();
    
    if (__kernel_name == "r") {
        reconstruct<fast_linear_rbf_type>(filenames, raster);
    } else if (__kernel_name == "r3") {
        reconstruct<fast_cubic_rbf_type>(filenames, raster);
    } else if (__kernel_name == "r5") {
        reconstruct<fast_quintic_rbf_type>(filenames, raster);
    } else if (__kernel_name.substr(0, 8) == "gaussian") {
        double sigma;
        if (__kernel_name[8] == ':') {
            sigma = atof(__kernel_name.substr(9).c_str());
        } else {
            usage("Syntax error in kernel definition: " + __kernel_name);
        }
        reconstruct<fast_gaussian_rbf_type>(filenames, raster, 
                                            gauss_t(0.5/(sigma*sigma)));
    } else if (__kernel_name.substr(0, 8) == "wendland") {
        double radius = 1;
        if (__kernel_name[8] == ':') {
            radius = atof(__kernel_name.substr(9).c_str());
        } else {
            usage("Syntax error in kernel definition: " + __kernel_name);
        }
        if (__verbose) {
            std::cout << "kernel radius = " << radius << '\n';
        }
        reconstruct<wendland_rbf_type>(filenames, raster, 
                                       wendland_t(radius));
    } else if (__kernel_name.substr(0, 9) == "tgaussian") {
        double sigma;
        double radius;
        if (__kernel_name[9] != ':') {
            usage("Syntax error in kernel definition: " + __kernel_name);
        }
        size_t next_colon = __kernel_name.find(':', 10);
        if (next_colon == std::string::npos) {
            usage("Syntax error in kernel definition: " + __kernel_name);
        }
        sigma = atof(__kernel_name.substr(10, next_colon-10).c_str());
        radius = sigma*atof(__kernel_name.substr(next_colon+1).c_str());
        if (__verbose) {
            std::cout << "Selected kernel: truncated gaussian: sigma^2="
                      << sigma
                      << ", radius = " << radius << '\n';
        }
        truncated_gaussian_type tgauss(1./(2.*sigma*sigma), radius);
        reconstruct<trunc_gaussian_rbf_type>(filenames, raster, 
                                             trgauss_t(0.5/(sigma*sigma), 
                                                       radius));
    } else {
        usage("Unrecognized kernel type: " + __kernel_name);
    }
    
    spurt::nrrd_params<float, 3> params;
    params.sizes()[0] = __number_of_bins;
    params.sizes()[1] = __resolution[0];
    params.sizes()[2] = __resolution[1];
    params.mins()[0] = 0;
    params.mins()[1] = __bounds.min()[0];
    params.mins()[2] = __bounds.min()[1];
    nvis::vec2 spacing = __bounds.size() / 
                         nvis::vec2(__resolution - nvis::ivec2(1,1));
    params.spacings()[0] = TWO_PI/(double)__number_of_bins;
    params.spacings()[1] = spacing[0];
    params.spacings()[2] = spacing[1];
    
    double* data = static_cast<double*>(calloc(__number_of_bins*
                                               raster.size(), sizeof(double)));
    for (size_t i=0 ; i<__number_of_bins*raster.size() ; ++i) {
        vec2bins(&data[__number_of_bins*i], raster[i]);
    }
    spurt::writeNrrd(data, output_name, params);
    std::cout << output_name << " has been exported\n";
    
    return 0;
}
