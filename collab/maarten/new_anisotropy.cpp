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
#include <format/format.hpp>
#include "maarten_utils.hpp"
#include <boost/shared_ptr.hpp>

const double TWO_PI = 2.*M_PI;

// global variables
nvis::ivec2 resolution(200, 200);
nvis::bbox2 bounds;
bool verbose = false;
std::string kernel_name = "r3";
size_t number_of_bins = 36;
double min_dist_to_source;

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

namespace xrbf = spurt::RBF;

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
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, linear_type, fast_solver_type>     fast_linear_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, linear_type, medium_solver_type>   medium_linear_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, linear_type, slow_solver_type>     slow_linear_rbf_type;
// cubic
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, cubic_type, fast_solver_type>      fast_cubic_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, cubic_type, medium_solver_type>    medium_cubic_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, cubic_type, slow_solver_type>      slow_cubic_rbf_type;
// quintic
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, quintic_type, fast_solver_type>    fast_quintic_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, quintic_type, medium_solver_type>  medium_quintic_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, quintic_type, slow_solver_type>    slow_quintic_rbf_type;
// gaussian
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, gaussian_type, fast_solver_type>   fast_gaussian_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, gaussian_type, medium_solver_type> medium_gaussian_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, gaussian_type, slow_solver_type>   slow_gaussian_rbf_type;

// RBF interpolators with compact support
// Wendland C2
typedef xrbf::CompactSupportRBFInterpolator<nvis::vec1, double, 2, wendland_type>            wendland_rbf_type;
// truncated gaussian
typedef xrbf::CompactSupportRBFInterpolator<nvis::vec1, double, 2, truncated_gaussian_type>  trunc_gaussian_rbf_type;

class raster_helper {
public:
    raster_helper(const nvis::bbox2& bounds, const nvis::ivec3& size)
        : _bounds(bounds), _size(size) 
    {
        _data = (float *)calloc(size[0]*size[1]*size[2], sizeof(float));
        _step = bounds.size()/nvis::vec2(size[1]-1, size[2]-1);
        _dtheta = 2*M_PI/(double)size[0];
    }
    
    ~raster_helper() {
        delete[] _data;
    }
    
    const nvis::bbox2& bounds() const { return _bounds; }
    const nvis::vec2& spacing() const { return _step; }
    const nvis::ivec3& resolution() const { return _size; }
    
    int nb_positions() const {
        return _size[1]*_size[2];
    }
    
    nvis::vec2 pos(int n) const {
        return pos(n%_size[1], n/_size[1]);
    }
    
    float value(int n, int bin) const {
        int i = n%_size[1];
        int j = n/_size[1];
        return _data[_size[0]*(i+_size[1]*j)+bin];
    }
    
    float& value(int n, int bin) {
        int i = n%_size[1];
        int j = n/_size[1];
        return _data[_size[0]*(i+_size[1]*j)+bin];
    }
    
    float& bin(int n, const nvis::vec2& v) {
        double theta = atan2(v[1], v[0]);
        if (theta < 0) theta += 2*M_PI;
        int b = floor(theta / _dtheta);
        return value(n, b);
    }

    nvis::vec2 pos(int i, int j) const {
        return _bounds.min() + nvis::vec2(i,j)*_step;
    }
    
    float* get_data() { return _data; }

private:    
    float* _data;
    nvis::bbox2 _bounds;
    nvis::ivec3 _size;
    nvis::vec2 _step;
    double _dtheta;
};

std::vector<std::pair<nvis::vec2, nvis::vec2> > measures_at_102_29;

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
    
    if (nvis::norm(source - nvis::vec2(102, 29)) > min_dist_to_source) {
        derivative_type g = interpolator.derivative(nvis::vec2(102, 29));
        nvis::vec2 v(g[0][0], g[1][0]);
        v /= nvis::inner(v, v);
        v *= spurt::maarten::km_per_angle_ratio(nvis::vec2(102, 29), v);
        measures_at_102_29.push_back(std::pair<nvis::vec2, nvis::vec2>(source, v));
    }
    
    progress_message msg(nb_jobs, "interpolations");
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (int n=0 ; n<number_of_samples ; ++n) {
            
            nvis::vec2 x = raster.pos(n);
            if (nvis::norm(x-source) < min_dist_to_source) continue;
            
            derivative_type nabla = interpolator.derivative(x);
                
            nvis::vec2 g(nabla[0][0], nabla[1][0]); // dt/drx, dt/dry
            raster.bin(n, g) = 1/nvis::norm(g);
            
            if (verbose) {
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
    spurt::to_vector(data, nin);
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
              << " -m | --min <float>         Min distance to source (default: 10\% of domain diameter)\n"
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
        std::vector<nvis::vec1> weights;
        nvis::vec2 source = read_nrrd(points, times, weights, filenames[i]);
        Int_ int_(points, times, weights, fun, verbose);
        std::cout << "Processing source #" << i+1 << "/" 
                  << filenames.size() << " at " << source << '\n';
        anisotropy(int_, raster, source);
    }
    
    std::fstream out("test_aniso.csv", std::ios::out);
    for (int i=0 ; i<measures_at_102_29.size() ; ++i) {
        const nvis::vec2& at = measures_at_102_29[i].first;
        const nvis::vec2& v = measures_at_102_29[i].second;
	if (nvis::norm(at-nvis::vec2(102,29)) < min_dist_to_source) continue;
        out << at[0] << ',' << at[1] << ',' 
            << atan2(v[1], v[0]) << ',' << nvis::norm(v) << '\n';
    }
    out.close();
}

int main(int argc, char* argv[])
{
    std::string input_name="", output_name="", path="";
    me = argv[0];
    bounds.reset();
    min_dist_to_source = 0;
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
            resolution[0] = atoi(argv[++i]);
            resolution[1] = atoi(argv[++i]);
        } else if (arg == "-m" || arg == "--min") {
            if (i == argc-1) {
               usage("missing minimum distance to source");
            }
            std::string tmp_str = argv[++i];
            if (*tmp_str.rbegin() == '%') {
                min_dist_to_source = atof(tmp_str.substr(0, tmp_str.size()-1).c_str());
                relative_distance = true;
            } else {
                min_dist_to_source = atof(tmp_str.c_str());
            }
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
            verbose = true;
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
    if (verbose) {
        std::cout << "input:  " << input_name << '\n';
        std::cout << "output: " << output_name << '\n';
        std::cout << "bounds: " << bounds << '\n';
    }
    
    nvis::ivec3 size(number_of_bins, resolution[0], resolution[1]);
    raster_helper raster(bounds, size);
    
    nvis::vec2 domain_size = bounds.size();
    if (min_dist_to_source == 0) {
        min_dist_to_source = 0.1*nvis::norm(domain_size);
    }
    else if (relative_distance) {
        min_dist_to_source *= nvis::norm(domain_size);
    }
    if (verbose) {
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
        if (verbose) {
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
        if (verbose) {
            std::cout << "Selected kernel: truncated gaussian: sigma^2="
                      << sigma
                      << ", radius = " << radius << '\n';
        }
        truncated_gaussian_type tgauss(1./(2.*sigma*sigma), radius);
        reconstruct<trunc_gaussian_rbf_type>(filenames, raster, tgauss);
    } else {
        usage("Unrecognized kernel type: " + kernel_name);
    }
    
    spurt::nrrd_params<float, 3> params;
    params.sizes()[0] = number_of_bins;
    params.sizes()[1] = resolution[0];
    params.sizes()[2] = resolution[1];
    params.mins()[0] = 0;
    params.mins()[1] = bounds.min()[0];
    params.mins()[2] = bounds.min()[1];
    nvis::vec2 spacing = bounds.size() / nvis::vec2(resolution - nvis::ivec2(1,1));
    params.spacings()[0] = TWO_PI/(double)number_of_bins;
    params.spacings()[1] = spacing[0];
    params.spacings()[2] = spacing[1];
    spurt::writeNrrd(raster.get_data(), output_name, params);
    std::cout << output_name << " has been exported\n";
    
    return 0;
}
