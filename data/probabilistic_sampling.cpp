#include <sstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <cctype>
#include <regex>
#include <stdexcept>
#include <locale>
#include <iomanip>
#include <algorithm>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/filesystem.hpp>

#include <data/field_wrapper.hpp>
#include <data/raster.hpp>
#include <format/filename.hpp>
#include <image/nrrd_wrapper.hpp>
#include <image/nrrd_field.hpp>
#include <image/probe.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>

#include <VTK/vtk_utils.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>

#ifdef _OPENMP
#include <omp.h>
#endif

std::string name_image, name_feature, name_out;
size_t npts;
bool verbose;

void initialize(int argc, const char* argv[])
{
    namespace xcl = xavier::command_line;
        
    xcl::option_traits 
            required(true, false, "Required Options"), 
            optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Probabilistic sampling of an image based on feature magnitude");
    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("image", name_image, "Image filename", required);
        parser.add_value("feature", name_feature, "Feature filename", required);
        parser.add_value("output", name_out, "Output file basename ", required);
        parser.add_value("number", npts, 1000, "Number of samples", optional);
        parser.add_flag("verbose", verbose, "Toggle verbose mode", optional);
        
        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

typedef std::pair<double, size_t> value_type;
typedef std::array<double, 2> pos_t;

std::ostream& operator<<(std::ostream& os, const value_type& value) {
    os << "(" << value.first << ", " << value.second << ")";
    return os;
}

size_t select(const std::vector<long double>& cdf, double value) {
    auto iter = std::lower_bound(cdf.begin(), cdf.end(), value);
    // std::cerr << "idx(" << value << ")=" << std::distance(cdf.begin(), iter) << " ("
        // << (double)std::distance(cdf.begin(), iter) / (double)cdf.size() << ")\n";
    return std::distance(cdf.begin(), iter);
}

int main(int argc, const char* argv[]) {
    initialize(argc, argv);
    
    Nrrd* image_nrrd = xavier::readNrrd(name_image);
    
    if (image_nrrd && verbose) std::cout << name_image << " successfully imported\n";
    else if (!verbose) {
        std::cerr << "ERROR: unable to open " << name_image << '\n';
        exit(1);
    } 
    
    Nrrd* feature_nrrd = xavier::readNrrd(name_feature);
    if (feature_nrrd && verbose) std::cout << name_feature << " successfully imported\n";
    else if (!verbose) {
        std::cerr << "ERROR: unable to open " << name_feature << '\n';
        exit(1);
    } 
    
    xavier::nrrd_data_wrapper<int> image(image_nrrd);
    xavier::nrrd_data_wrapper<double> feature(feature_nrrd);
    
    size_t size[2] = { image_nrrd->axis[0].size, image_nrrd->axis[1].size };
    
    std::vector<value_type> pixels(size[0]*size[1]);
    xavier::ProgressDisplay progress(verbose);
    progress.start(size[0]*size[1], "Reading feature values");
    for (size_t i=0; i<size[0]*size[1]; ++i) {
        pixels[i] = value_type(std::abs(feature[i]), i);
        progress.update(i);
    }
    progress.end();
    if (verbose) std::cout << "feature values read\n";
    
    struct less {
        bool operator()(const value_type& c0, const value_type& c1) {
            return c0.first < c1.first;
        }
    };
    
    typedef std::numeric_limits<long double> ldbl_limits;
    std::cout.precision(ldbl_limits::max_digits10);
    
    std::sort(pixels.begin(), pixels.end(), less());
    
    std::vector<long double> cdf(pixels.size());
    cdf[0] = 0;
    long double sum = 0;
    
    progress.start(pixels.size(), "Computing CDF");
    for (size_t i=0; i<pixels.size(); ++i) {
        cdf[i] = pixels[i].first + sum;
        sum += pixels[i].first;
        progress.update(i);
    }
    progress.end();
    
    if (verbose) std::cout << "Cumulative probability density function computed, sum=" << sum << "\n";
    
    size_t ii=0;
    std::for_each(cdf.begin(), cdf.end(), [&](long double& v) { 
        v = v/sum;
    });
    
    std::set<size_t> selected;
    srand48(time(0));
    progress.start(npts, "Random sampling");
    while (selected.size() < npts) {
        double d = drand48();
        selected.insert(select(cdf, d));
        progress.update(selected.size());
    }
    progress.end();
    
    std::vector<pos_t> data_points;
    std::vector<int> values;
    std::vector<pos_t> tcoords;
    for (auto it=selected.begin(); it!=selected.end(); ++it) {
        size_t idx = pixels[*it].second;
        double x = idx % size[0];
        double y = size[1] - 1 - idx / size[0];
        data_points.push_back({{x,y}});
        tcoords.push_back({{x/(double)size[0], y/(double)size[1]}});
        values.push_back(image[idx]);
    }
    
    vtkPolyData* data = vtk_utils::make_points(data_points);
    vtk_utils::add_scalars(data, values, true, "altitude");
    vtk_utils::add_tcoords(data, tcoords);
    vtk_utils::saveVTK(data, name_out);
    if (verbose) std::cout << "discrete samples successfully exported\n";
    
    return 0;
}

