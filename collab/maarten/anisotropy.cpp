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
#include <data/raster.hpp>
#include <math/RBF.hpp>
#include <math/RBFbasis.hpp>
#include <image/nrrd_wrapper.hpp>
#include <format/format.hpp>
#include "maarten_utils.hpp"
#include <boost/shared_ptr.hpp>

const double TWO_PI = 2.*M_PI;

// global variables
nvis::ivec2 __resolution(200, 200);
nvis::bbox2 __bounds;
bool        __verbose = false;
std::string __kernel_name = "r3";
size_t      __number_of_bins = 36;

using namespace xavier::maarten;

typedef std::vector<double> bins;

class raster_helper : private xavier::raster2d<bins>
{
public:
    typedef bins                               value_type;
    typedef typename bins::value_type          scalar_type;
    typedef xavier::raster2d<bins>             base_type;
    typedef typename base_type::grid_type      grid_type;
    typedef typename grid_type::size_type      size_type;
    typedef typename grid_type::vec_type       vec_type;
    typedef typename grid_type::bounds_type    bounds_type;
    typedef typename base_type::const_iterator const_iterator;
    typedef typename base_type::iterator       iterator;
    
    raster_helper(const bounds_type& bounds, const nvis::ivec3& size)
        : _size(size), base_type(grid_type(nvis::subv<1,2>(size), bounds))
    {
        _dtheta = 2*M_PI/(double)_size[0];
    }
    
    ~raster_helper() {
    }
    
    const bounds_type& bounds() const { return base_type::grid().bounds(); }
    const vec_type& spacing() const { return base_type::grid().spacing(); }
    const nvis::ivec3& resolution() const { return _size; }
    
    size_type nb_positions() const {
        return base_type::grid().size();
    }
    
    vec_type pos(size_type n) const {
        return base_type::grid()(n%_size[1], n/_size[1]);
    }
    
    scalar_type value(size_type n, size_type bin) const {
        size_type i = n%_size[1];
        size_type j = n/_size[1];
        return base_type::operator()(i,j)[bin];
    }
    
    scalar_type& value(size_type n, size_type bin) {
        size_type i = n%_size[1];
        size_type j = n/_size[1];
        return base_type::operator()(i,j)[bin];
    }
    
    scalar_type& bin(size_type n, const vec_type& v) {
        double theta = atan2(v[1], v[0]);
        if (theta < 0) theta += 2*M_PI;
        size_type b = floor(theta / _dtheta);
        return value(n, b);
    }

    vec_type pos(size_type i, size_type j) const {
        return base_type::grid()(i,j);
    }
    
    const_iterator begin() const { return base_type::begin(); }
    const_iterator end() const { return base_type::end(); }
    
    iterator begin() { return base_type::begin(); }
    iterator end() { return base_type::end(); }

private:
    nvis::ivec3 _size;
    double      _dtheta;
};

template<typename _Int>
void anisotropy(const _Int& interpolator, 
                raster_helper& raster,
                const nvis::vec2& source)
{
    typedef typename _Int::derivative_type derivative_type;
    int number_of_samples = raster.nb_positions();
        
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
            nvis::vec2 x = raster.pos(n);            
            derivative_type nabla = interpolator.derivative(x);
            nvis::vec2 g(nabla[0][0], nabla[1][0]); // dt/drx, dt/dry
            raster.bin(n, g) = 1/nvis::norm(g);
            
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
    
    std::cout << "\nRBF reconstruction completed in " << _timer.elapsed() << " seconds ("
              << (float)nb_jobs/_timer.elapsed() << " Hz)\n";
}

void usage(const std::string&);

nvis::vec2 read_nrrd(std::vector<nvis::vec2>& points,
                     std::vector<nvis::vec1>& times,
                     std::vector<nvis::vec1>& weights,
                     const std::string& filename)
{
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, filename.c_str(), NULL)) {
        std::cerr << "read_nrrd: " << biffGetDone(NRRD) << std::endl;
        throw;
    }
    std::vector<double> data;
    xavier::to_vector(data, nin);
    size_t N = nin->axis[0].size;
    bool has_weights = (N == 4);
    if (!has_weights) {
        usage("Input file does not contain RBF weight information");
    }
    size_t nb_pts = data.size()/N;
    points.resize(nb_pts);
    times.resize(nb_pts);
    weights.resize(nb_pts);
    
    nvis::vec2 source(data[0], data[1]);
    for (size_t i=0 ; i<nb_pts ; ++i) {
        points[i][0]  = data[N*i  ];
        points[i][1]  = data[N*i+1];
        times[i][0]   = data[N*i+2];
        weights[i][0] = data[N*i+3];
    }
    nrrdNuke(nin);
    return source;
}

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
                 raster_helper& raster,
                 const typename Int_::function_type& fun =
                       typename Int_::function_type() )
{
    for (int i=0 ; i<filenames.size() ; ++i) {
        std::vector<nvis::vec2> points;
        std::vector<nvis::vec1> times;
        std::vector<nvis::vec1> distances;
        nvis::bbox2 bounds = read_text(points, times, distances, 
                                       filenames[i], __verbose);
                                       
        
        result = reconstruct<Int_>(__resolution, bounds, interpolator);
                                       
        Int_ int_(points, times, weights, fun, __verbose);
        std::cout << "Processing source #" << i+1 << "/" 
                  << filenames.size() << " at " << source << '\n';
        anisotropy(int_, raster, source);
    }
}

int main(int argc, char* argv[])
{
    std::string input_name="", output_name="", path="";
    me = argv[0];
    bounds.reset();
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
            bounds.min()[0] = atof(argv[++i]);
            bounds.max()[0] = atof(argv[++i]);
            bounds.min()[1] = atof(argv[++i]);
            bounds.max()[1] = atof(argv[++i]);
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
            kernel_name = argv[++i];
        } else if (arg == "-b" || arg == "--bins") {
            if (i == argc-1) {
                usage("missing bin number");
            }
            number_of_bins = atoi(argv[++i]);
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
    if (nvis::any(bounds.min() > bounds.max())) {
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
        std::cout << "bounds: " << bounds << '\n';
    }
    
    nvis::ivec3 size(number_of_bins, __resolution[0], __resolution[1]);
    raster_helper raster(bounds, size);
    
    nvis::vec2 domain_size = bounds.size();
    if (min_dist_to_source == 0) {
        min_dist_to_source = 0.1*nvis::norm(domain_size);
    }
    else if (relative_distance) {
        min_dist_to_source *= nvis::norm(domain_size);
    }
    if (__verbose) {
        std::cout << "Minimum distance to source: " << min_dist_to_source << '\n';
    }
    
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
    
    if (kernel_name == "r") {
        reconstruct<fast_linear_rbf_type>(filenames, raster);
    } else if (kernel_name == "r3") {
        reconstruct<fast_cubic_rbf_type>(filenames, raster);
    } else if (kernel_name == "r5") {
        reconstruct<fast_quintic_rbf_type>(filenames, raster);
    } else if (kernel_name.substr(0, 8) == "gaussian") {
        double sigma;
        if (kernel_name[8] == ':') {
            sigma = atof(kernel_name.substr(9).c_str());
        } else {
            usage("Syntax error in kernel definition: " + kernel_name);
        }
        gaussian_type gauss(1./(2.*sigma*sigma));
        reconstruct<fast_gaussian_rbf_type>(filenames, raster, gauss);
    } else if (kernel_name.substr(0, 8) == "wendland") {
        double radius = 1;
        if (kernel_name[8] == ':') {
            radius = atof(kernel_name.substr(9).c_str());
        } else {
            usage("Syntax error in kernel definition: " + kernel_name);
        }
        if (__verbose) {
            std::cout << "kernel radius = " << radius << '\n';
        }
        reconstruct<wendland_rbf_type>(filenames, raster, wendland_type(radius));
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
        if (__verbose) {
            std::cout << "Selected kernel: truncated gaussian: sigma^2="
                      << sigma
                      << ", radius = " << radius << '\n';
        }
        truncated_gaussian_type tgauss(1./(2.*sigma*sigma), radius);
        reconstruct<trunc_gaussian_rbf_type>(filenames, raster, tgauss);
    } else {
        usage("Unrecognized kernel type: " + kernel_name);
    }
    
    xavier::nrrd_params<float, 3> params;
    params.sizes()[0] = number_of_bins;
    params.sizes()[1] = __resolution[0];
    params.sizes()[2] = __resolution[1];
    params.mins()[0] = 0;
    params.mins()[1] = bounds.min()[0];
    params.mins()[2] = bounds.min()[1];
    nvis::vec2 spacing = bounds.size() / nvis::vec2(__resolution - nvis::ivec2(1,1));
    params.spacings()[0] = TWO_PI/(double)number_of_bins;
    params.spacings()[1] = spacing[0];
    params.spacings()[2] = spacing[1];
    xavier::writeNrrd(raster.get_data(), output_name, params);
    std::cout << output_name << " has been exported\n";
    
    return 0;
}
