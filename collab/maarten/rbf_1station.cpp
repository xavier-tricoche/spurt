#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
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
#include <image/nrrd_wrapper.hpp>

bool verbose;

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
        size_t length = _os.str().size();
        _os << std::string(' ', 20);
        _os << '\r';
        return _os.str();
    }

    void reset(const std::string& what = "") {
        if (what.size()) _what = what;
    }
};

std::string me;
void printUsage(const std::string& msg) {
    if (!msg.empty()) std::cerr << "ERROR: " << msg << '\n';
    std::cout
            << "USAGE: " << me << " [parameters] [options]\n"
            << "DESCRIPTION: Resample scattered scalar 2D dataset on\n"
            << "             regular lattice grid using polyharmonic RBF\n"
            << "             reconstruction\n"
            << "PARAMETERS:\n"
            << " -i | --input <string>        Input file name\n"
            << " -o | --output <string>       Output file name\n"
            << "OPTIONS:\n"
            << " -h | --help                  Print this information\n"
            << " -p | --path <string>         Path to be prepended to all file names (default: \"\")\n"
            << " -k | --kernel <string>       RBF (\"r\", \"r3\", \"r5\")                  (default: \"r3\")\n"
            << " -r | --resolution <int> (x2) Resampling resolution                  (default: 200^2)\n"
            << " -v | --verbose               Activate verbose mode                  (default: false)\n";
    exit(1);
}

struct Linear {
    double operator()(double r) const {
        return r;
    }
};

struct Cubic {
    double operator()(double r) const {
        return r*r*r;
    }
};

struct Quintic {
    double operator()(double r) const {
        double r2 = r*r;
        return r2*r2*r;
    }
};

typedef spurt::RBF::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, Linear>    linear_rbf_type;
typedef spurt::RBF::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, Cubic>     cubic_rbf_type;
typedef spurt::RBF::InfiniteSupportRBFInterpolator<nvis::vec1, double, 2, Quintic>   quintic_rbf_type;

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
float* reconstruct(const nvis::ivec2& resolution,
                   const nvis::bbox2& bounds,
                   const _Interpolator& interpolator) {
    nvis::vec2 spacing = bounds.size() / nvis::vec2(resolution - nvis::ivec2(1, 1));
    size_t number_of_samples = resolution[0]*resolution[1];
    float* result = (float*)calloc(number_of_samples, sizeof(float));

    size_t number_of_threads = 1;
#ifdef _OPENMP
    number_of_threads = omp_get_max_threads();
#endif
    size_t counter = 0;
    size_t nb_jobs = number_of_samples;
    nvis::timer _timer;
    
    typedef typename boost::shared_ptr<progress_message> msg_ptr;
    std::vector<msg_ptr> messages(number_of_threads);
    for (size_t i=0 ; i<number_of_threads ; ++i) {
        messages[i].reset(new progress_message(nb_jobs, "interpolations"));
    }

    progress_message msg(number_of_samples, "interpolations");

    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic, 1)
        for (size_t n=0 ; n<number_of_samples ; ++n) {
            int i = n%resolution[0];
            int j = n/resolution[0];
            nvis::vec2 x = bounds.min() + nvis::vec2(i, j)*spacing;
            result[n] = interpolator(x)[0];

            if (verbose) {
                #pragma omp atomic
                ++counter;
                size_t my_counter = counter;
                double elapsed = _timer.elapsed();
                std::cout << msg(my_counter, elapsed) << std::flush;
            }
        }
    }
    std::cout << "\nRBF reconstruction completed in " << _timer.elapsed() << " seconds ("
              << (float)number_of_samples/_timer.elapsed() << " Hz)\n";

    if (verbose) {
        std::cout << "reconstruction min: " << *std::min_element(result, &result[number_of_samples]) << '\n';
        std::cout << "reconstruction max: " << *std::max_element(result, &result[number_of_samples]) << '\n';
        check_solution(interpolator);
    }

    return result;
}

// custom reader for KCD* files
// format:
// source_lon; source_lat; receiver_lon; receiver_lat; travel_time; \
//             <optional: source-receiver distance>; receiver index.
nvis::bbox2 read_file(std::vector<nvis::vec2>& points, std::vector<nvis::vec1>& times,
                      const std::string& file_name) {
    static const double invalid = std::numeric_limits<double>::max();
    points.clear();
    times.clear();
    nvis::bbox2 bounds;
    double min=invalid, max=-invalid;
    std::fstream file(file_name.c_str(), std::ios::in);
    
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
            if (word.empty() || !iss.gcount()) break;
        }
        has_7_terms = ( nwords == 7 );
        if (verbose) {
            std::cout << "data file has " << nwords << " terms per row\n";
        }
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
        bounds.add(points.back());
    }
    file.close();

    return bounds;
}

int main(int argc, char* argv[]) {
    std::string input_name, output_name, path="";
    input_name = output_name = "undefined!";
    me = argv[0];
    nvis::ivec2 resolution(200, 200);
    verbose = false;
    std::string kernel_name = "r3";

    for (int i=1 ; i<argc ; ++i) {
        std::string argument = argv[i];
        if (argument == "-h") printUsage("");
        else if (argument == "-i" || argument == "--input") {
            if (i == argc-1) printUsage("missing input name");
            input_name = argv[++i];
        }
        else if (argument == "-o" || argument == "--output") {
            if (i == argc-1) printUsage("missing output name");
            output_name = argv[++i];
        }
        else if (argument == "-p" || argument == "--path") {
            if (i == argc-1) printUsage("missing path name");
            path = argv[++i];
        }
        else if (argument == "-k" || argument == "--kernel") {
            if (i == argc-1) printUsage("missing kernel name");
            kernel_name = argv[++i];
        }
        else if (argument == "-r" || argument == "--res") {
            if (i >= argc-2) printUsage("missing resolution parameters");
            resolution[0] = atoi(argv[++i]);
            resolution[1] = atoi(argv[++i]);
        }
        else if (argument == "-v" || argument == "--verbose") {
            verbose = true;
        }
        else printUsage("unrecognized input argument");
    }

    if (input_name == "undefined!" || output_name == "undefined!") {
        printUsage("Missing input or output file name");
    }
    if (path != "") {
        if (*path.rbegin() != '/') path.push_back('/');
        input_name = path + input_name;
        output_name = path + output_name;
    }
    if (verbose) {
        std::cout << "input:  " << input_name << '\n';
        std::cout << "output: " << output_name << '\n';
    }

    size_t number_of_samples = resolution[0]*resolution[1];
    std::vector<nvis::vec2> points;
    std::vector<nvis::vec1> values;
    nvis::bbox2 bounds = read_file(points, values, input_name);
    if (verbose) {
        std::cout << "bounding box: " << bounds << '\n';
    }

    float* result;
    nvis::timer _timer;
    if (kernel_name == "r") {
        linear_rbf_type interpolator(points, values, Linear());
        result = reconstruct<linear_rbf_type>(resolution, bounds, interpolator);
    }
    else if (kernel_name == "r3") {
        cubic_rbf_type interpolator(points, values, Cubic());
        result = reconstruct<cubic_rbf_type>(resolution, bounds, interpolator);
    }
    else if (kernel_name == "r5") {
        quintic_rbf_type interpolator(points, values, Quintic());
        result = reconstruct<quintic_rbf_type>(resolution, bounds, interpolator);
    }
    else {
        printUsage("Unrecognized kernel type");
    }

    spurt::nrrd_params<float, 2> params;
    params.sizes()[0] = resolution[0];
    params.sizes()[1] = resolution[1];
    params.mins()[0] = bounds.min()[0];
    params.mins()[1] = bounds.min()[1];
    nvis::vec2 spacing = bounds.size() / nvis::vec2(resolution - nvis::ivec2(1,1));
    params.spacings()[0] = spacing[0];
    params.spacings()[1] = spacing[1];
    spurt::writeNrrd(result, output_name, params);
    std::cout << output_name << " has been exported\n";

    return 0;
}
