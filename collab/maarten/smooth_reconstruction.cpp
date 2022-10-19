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
#define BPO_WRAPPER_IS_BROKEN
#ifndef BPO_WRAPPER_IS_BROKEN
#   include <misc/option_parse.hpp>
#endif
#include <format/format.hpp>
#include <boost/shared_ptr.hpp>

// global variables
nvis::ivec2 resolution(800, 800);
nvis::bbox2 bounds;
bool verbose = false;
std::string kernel_name = "r3";

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
        if (_what.size()) _os << " " << _what;
        if (elapsed > 0)
            _os << " in " << elapsed << " seconds (" << (float)n/elapsed << " Hz)";
        _os << "                       \r";
        return _os.str();
    }

    void reset(const std::string& what = "") {
        if (what.size()) _what = what;
    }
};
namespace xrbf = xavier::RBF;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>  matrix_type;
typedef Eigen::HouseholderQR<matrix_type>                      solver_type;

typedef xrbf::linear_function<double>   linear_type;
typedef xrbf::cubic_function<double>    cubic_type;
typedef xrbf::quintic_function<double>  quintic_type;
typedef xrbf::wendland_function<double> wendland_type;

typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, linear_type, solver_type >  linear_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, cubic_type, solver_type>    cubic_rbf_type;
typedef xrbf::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, quintic_type, solver_type>  quintic_rbf_type;
typedef xrbf::CompactSupportRBFInterpolator<nvis::vec1, double, 2, wendland_type>               wendland_rbf_type;

template<typename _Interpolator>
void check_solution(const _Interpolator& f) {
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

void save_rbf(const std::vector<nvis::vec2>& points,
              const std::vector<nvis::vec1>& values,
              const std::vector<nvis::vec1>& weights,
              const std::string& kernel_name,
              const std::string& filename);

template<typename Int_>
float* reconstruct(const std::vector<nvis::vec2>& points,
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

// custom reader for KCD* files
// format:
// source_lon; source_lat; receiver_lon; receiver_lat; travel_time; \
//             <optional: source-receiver distance>; receiver index.
nvis::bbox2 read_text(std::vector<nvis::vec2>& points, std::vector<nvis::vec1>& times,
                      const std::string& file_name) {
    static const double invalid = std::numeric_limits<double>::max();
    points.clear();
    times.clear();
    nvis::bbox2 _bounds;
    double min=invalid, max=-invalid;
    std::fstream file(file_name.c_str(), std::ios::in);
    
    if (!file) {
        throw std::runtime_error("unable to open " + file_name);
    }
    
    // check if file contains information about source-receiver distance
    bool has_7_terms = false;
    {
        std::string first_line;
        std::getline(file, first_line);
        std::istringstream iss(first_line);
        size_t nwords;
        for (nwords=0 ; true ; ++nwords) {
            std::string word;
            iss >> word;
            if (word.empty()) break;
        }
        has_7_terms = ( nwords == 7 );
        if (verbose) {
            std::cout << "data file has " << nwords << " terms per row\n";
        }
        file.seekg(0);
    }
    
    while (!file.eof()) {
        double blah, x, y, t, sentinel=invalid;
        file >> blah >> blah >> x >> y >> t >> sentinel;
        if (has_7_terms) file >> sentinel;
        if (sentinel == invalid) break;
        points.push_back(nvis::vec2(x,y));
        times.push_back(nvis::vec1(t));
        min = std::min(t, min);
        max = std::max(t, max);
        _bounds.add(points.back());
    }
    file.close();
    
    return _bounds;
}

nvis::bbox2
read_nrrd(std::vector<nvis::vec2>& points, 
          std::vector<nvis::vec1>& times,
          std::vector<nvis::vec1>& weights,
          const std::string& filename) {
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, filename.c_str(), NULL)) {
        std::cerr << "read_data_from_nrrd: " << biffGetDone(NRRD) << std::endl;
        throw;
    }
    std::vector<float> data;
    xavier::to_vector(data, nin);
    size_t N = nin->axis[0].size;
    bool has_weights = (N == 4);
    size_t nb_pts = data.size()/N;
    points.resize(nb_pts);
    times.resize(nb_pts);
    if (has_weights) weights.resize(nb_pts);
    if (verbose && has_weights) 
        std::cout << "input data contains precomputed weights\n";
    nvis::bbox2 _bounds;
    for (size_t i=0 ; i<nb_pts ; ++i) {
        points[i][0] = data[N*i  ];
        points[i][1] = data[N*i+1];
        times[i][0]  = data[N*i+2];
        if (has_weights) weights[i][0] = data[N*i+3];
        _bounds.add(points[i]);
    }
    nrrdNuke(nin);
    return _bounds;
}

// save full information about RBF reconstruction
void save_rbf(const std::vector<nvis::vec2>& points,
              const std::vector<nvis::vec1>& values,
              const std::vector<nvis::vec1>& weights,
              const std::string& kernel_name,
              const std::string& filename) {
    size_t npts = points.size();
    double *data = (double*)calloc(npts*4, sizeof(double));
    for (size_t i=0 ; i<npts ; ++i) {
        data[4*i  ] = points[i][0];
        data[4*i+1] = points[i][1];
        data[4*i+2] = values[i][0];
        data[4*i+3] = weights[i][0];
    }
    xavier::nrrd_params<double, 2> params;
    params.sizes()[0] = 4;
    params.sizes()[1] = npts;
    params.description() = "RBF reconstruction data for kernel " + kernel_name;
    params.labels()[0] = "x_rec;y_rec;time;weight";
    params.labels()[1] = "RBF centers";
    xavier::writeNrrd(data, filename, params);
    std::cout << filename << " has been exported\n";
}

#ifdef BPO_WRAPPER_IS_BROKEN
std::string me;
void usage(const std::string& message="") {
    if (!message.empty()) {
        std::cerr << "ERROR: " << message << '\n';
    }
    std::cout << '\n'
              << "DESCRIPTION: Compute a smooth reconstruction of the\n"
              << "tomographic travel time associated with a single source\n"
              << "based on the data available at a discrete set of\n"
              << "stations. Radial basis functions are used for the\n"
              << "reconstruction, whose weight coefficients can be stored\n"
              << "for reuse. In addition to the travel time, the gradient\n"
              << "can be smoothly reconstructed using the same basis.\n"
              << '\n'
              << "USAGE: " << me << " [parameters] [options]\n"
              << '\n'
              << "PARAMETERS:\n"
              << " -i | --input <string>       Input file name\n"
              << " -o | --output <string>      Output base name\n"
              << '\n'                          
              << "OPTIONS:\n"                  
              << " -h | --help                 Print this information\n"
              << " -p | --path <string>        Path to prepend to all file names\n"
              << " -w | --weights <string>     Filename to export computed weights\n"
              << " -k | --kernel <string>      Reconstruction kernel\n"
              << "                             (\"r\", \"r3\", \"r5\", or \"wendland:<radius>\")\n"
              << " -r | --resolution <int> x2  Resampling resolution\n"
              << " -g | --gradient             Compute gradient\n"
              << " -v | --verbose              Turn on verbose mode\n"
              << std::endl;
    exit(!message.empty());
}
#endif

int main(int argc, char* argv[]) {
    std::string input_name="", output_name="", weights_name="", path="";
    bool export_weights = false;
    bool import_weights = false;

#ifndef BPO_WRAPPER_IS_BROKEN
    namespace bpo = boost::program_options;
    namespace clt = xavier::cmdline_tools;

    bpo::options_description desc(
        "Compute anisotropy at vertices of a raster grid\n"
        "based on the smooth reconstruction of the mapping\n"
        "from (x_source, y_source) x (x_receiver, y_receiver)\n"
        "to the associated travel time");

    desc.add_options()
    ("help,h", "Display this message")
    ("input,i",      clt::value(&input_name, false)->required(),   "File containing list of input file names")
    ("output,o",     clt::value(&output_name, false)->required(),  "Output file name")
    ("weights,w",    clt::value(&weights_name, false)->required(), "RBF weights base name")
    ("export,e",     clt::value(&export_name),                     "Single data filename")
    ("path,p",       clt::value(&path),                            "Path to be prepended to all file names")
    ("kernel,k",     clt::value(&kernel_name),                     "RBF kernel")
    ("resolution,r", clt::array(&resolution),                      "Resampling resolution")
    ("bins,b",       clt::value(&number_of_bins),                  "Number of bins to discretize orientation spectrum")
    ("verbose,v",    bpo::bool_switch(&verbose),                   "Turn ON verbose mode");

    clt::parse_command_line(argc, argv, desc);
#else
    me = argv[0];
    bool do_gradient = false;
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
        } else if (arg == "-w" || arg == "--weights") {
            if (i == argc-1) {
                usage("missing weights filename");
            }
            weights_name = argv[++i];
            export_weights = true;
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
        } else if (arg == "-g" || arg == "--gradient") {
            do_gradient = true;
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else {
            usage("unrecognized argument: " + arg);
        }
    }
    if (input_name.empty()) usage("Missing input filename");
    if (output_name.empty()) usage("Missing output filename");
#endif

    if (path != "") {
        if (*path.rbegin() != '/') path.push_back('/');
        input_name = path + input_name;
        output_name = path + output_name;
        if (export_weights) weights_name = path + weights_name;
    }
    if (verbose) {
        std::cout << "input:  " << input_name << '\n';
        std::cout << "output: " << output_name << '\n';
    }
    weights_name = xavier::get_basename(weights_name);

    size_t number_of_samples = resolution[0]*resolution[1];

    std::vector<nvis::vec2> points;
    std::vector<nvis::vec1> values, weights;
    
    // determine type of file to import
    std::string ext = xavier::get_extension(input_name);
    if (ext == "nrrd" || ext == ".nhdr") {
        bounds = read_nrrd(points, values, weights, input_name);
        import_weights = !weights.empty();
    }
    else {
        bounds = read_text(points, values, input_name);
    }
     
    if (import_weights && export_weights) {
        std::cerr << "WARNING: odd request to import and export RBF weights\n";
        std::cerr << "         export request will be ignored\n";
        export_weights = false;
    }
    
    if (verbose) {
        std::cout << "imported " << points.size() << " data points\n";
        std::cout << "bounding box: " << bounds << '\n';
    }

    nvis::vec2 domain_size = bounds.size();

    float* result;
    if (kernel_name == "r") 
        result = reconstruct<linear_rbf_type>(points, values, weights,
                                              weights_name, do_gradient);
    else if (kernel_name == "r3") 
        result = reconstruct<cubic_rbf_type>(points, values, weights, 
                                             weights_name, do_gradient);
    else if (kernel_name == "r5" )
        result = reconstruct<quintic_rbf_type>(points, values, weights, 
                                               weights_name, do_gradient);
    else if (kernel_name.substr(0, 8) == "wendland") {
        double radius = 1;
        if (kernel_name[8] == ':') {
            radius = atof(kernel_name.substr(9).c_str());
        }
        else usage("Syntax error in kernel definition: " + kernel_name);
        if (verbose) std::cout << "kernel radius = " << radius << '\n';
        result = reconstruct<wendland_rbf_type>(points, values, weights,
                                                weights_name, do_gradient, 
                                                wendland_type(radius));
    }
    else {
#ifndef BPO_WRAPPER_IS_BROKEN
        std::cerr << "ERROR: Unrecognized kernel type";
        std::cerr << desc << '\n';
#else
        usage("Unrecognized kernel type: " + kernel_name);
#endif
    }
    
    xavier::nrrd_params<float, 3> params;
    params.sizes()[0] = (do_gradient ? 3 : 1);
    params.sizes()[1] = resolution[0];
    params.sizes()[2] = resolution[1];
    params.mins()[0] = 0;
    params.mins()[1] = bounds.min()[0];
    params.mins()[2] = bounds.min()[1];
    nvis::vec2 spacing = bounds.size() / nvis::vec2(resolution - nvis::ivec2(1,1));
    params.spacings()[0] = 1;
    params.spacings()[1] = spacing[0];
    params.spacings()[2] = spacing[1];
    if (do_gradient) params.labels()[0] = "time;dtime/dx;dtime;dy";
    else params.labels()[0] = "time";
    params.labels()[1] = "X";
    params.labels()[2] = "Y";
    xavier::writeNrrd(result, output_name, params);
    std::cout << output_name << " has been exported\n";

    return 0;
}
