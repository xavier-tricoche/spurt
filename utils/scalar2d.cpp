#include <array>
#include <cctype>
#include <fstream>
#include <string>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <format/filename.hpp>
#include <misc/option_parse.hpp>
#include <math/types.hpp>

#include <utils/functions.hpp>

using namespace spurt;

constexpr size_t res = 512;

inline std::string lower_case(const std::string& str) {
    std::string r(str);
    std::transform(r.begin(), r.end(), r.begin(),
                   [](unsigned char c){ return std::tolower(c); }
                  );
    return r;
}

inline bool valid_bounds(const bounding_box<vec2>& bounds) {
    return all(bounds.min() <= bounds.max());
}

inline vec2 random_pos(const bounding_box<vec2>& bounds) {
    return bounds.min() + bounds.size() * vec2(drand48(), drand48());
}

size_t npts = 100;
std::string fname = "circle", name_out;
size_t seedval = time(0);
std::array<double, 4> bnds({-1, 1, -1, 1});
bool verbose = false;
bool ground_truth = false;

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;
        
    xcl::option_traits 
            required_group(true, false, "Required Options"), 
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Create scattered samples of 2D scalar field");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("output", name_out, "Output filename", required_group);
        parser.add_value("n", npts, npts, "Number of sample points", optional_group);
        parser.add_value("f", fname, fname, "function name", optional_group);
        parser.add_value("s", seedval, seedval, "Random seed value", optional_group);
        parser.add_flag("g", ground_truth, "Export ground truth", optional_group);
        parser.add_tuple<4>("bounds", bnds, bnds, "Sampling bounds", optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional_group);
        
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

int main(int argc, char* argv[]) {
    
    initialize(argc, const_cast<const char**>(argv));
    
    srand48(seedval);
    
    std::string filename = spurt::filename::remove_extension(name_out);
    
    fname = lower_case(fname);
    
    bbox2 bbox(vec2(bnds[0], bnds[2]), vec2(bnds[1], bnds[3]));
    
    std::fstream output(filename + ".xyz", std::ios::out);
    output << "#" << npts << " rows - \"x y f(x,y)\"\n";
    for (size_t i=0; i<npts; ++i) {
        vec2 p = random_pos(bbox);
        double f;
        
        if (fname == "sinc" || fname == "s") 
            f = functions2d::sinc(p);
        else if (fname == "real complex sinc" || fname == "real complex" || fname == "real" || fname == "r")
            f = std::real(functions2d::complex_sinc(p));
        else if (fname == "imaginary complex sinc" || fname == "imaginary complex" || fname == "imaginary" || fname == "imag complex" || fname == "imag" || fname == "i")
            f = std::imag(functions2d::complex_sinc(p));
        else if (fname == "marschner-lobb" || fname == "marschner" || fname == "marschner lobb" || fname == "ml") 
            f = functions2d::z_marschner_lobb(p);
        else if (fname == "circle" || fname == "c") 
            f = functions2d::circle(p);
        else if (fname == "disk" || fname == "d")
            f = functions2d::disk(p);
        else if (fname == "radial") 
            f = functions2d::radial(p);
        else if (fname == "sincos") 
            f = functions2d::sincos(p);
        else if (fname == "mathworks" || fname == "mathworks1")
            f = functions2d::mathworks1(p);
        else {
            std::cerr << "Unrecognized function name: " << fname << '\n';
            exit(1);
        }
        
        output << p << " " << f << '\n';
    }
    output.close();
    
    if (ground_truth) {
        std::vector<double> values(res*res);
        svec2 dims(res, res);
        vec2 spc = bbox.size() / (dims-1);
        
        for (size_t i=0; i<res*res; ++i) {
            int v = i/res;
            int u = i%res;
            vec2 p = bbox.min() + vec2(u,v) * spc;
            double& f = values[i];
            if (fname == "sinc" || fname == "s") 
                f = functions2d::sinc(p);
            else if (fname == "real complex sinc" || fname == "real complex" || fname == "real" || fname == "r")
                f = std::real(functions2d::complex_sinc(p));
            else if (fname == "imaginary complex sinc" || fname == "imaginary complex" || fname == "imaginary" || fname == "imag complex" || fname == "imag" || fname == "i")
                f = std::imag(functions2d::complex_sinc(p));
            else if (fname == "marschner-lobb" || fname == "marschner" || fname == "marschner lobb" || fname == "ml") 
                f = functions2d::z_marschner_lobb(p);
            else if (fname == "circle" || fname == "c") 
                f = functions2d::circle(p);
            else if (fname == "disk" || fname == "d")
                f = functions2d::disk(p);
            else if (fname == "radial") 
                f = functions2d::radial(p);
            else if (fname == "sincos") 
                f = functions2d::sincos(p);
            else if (fname == "mathworks" || fname == "mathworks1")
                f = functions2d::mathworks1(p);
            else {
                std::cerr << "Unrecognized function name: " << fname << '\n';
                exit(1);
            }
        }
        
        std::cout << "computed values: min=" << *std::min_element(values.begin(), values.end()) << ", max=" << *std::max_element(values.begin(), values.end()) << '\n';
        
        spurt::nrrd_utils::writeNrrdFromContainers(&values[0], filename + ".nrrd", dims, spc, bbox.min());
    }   
    
    return 0;
}