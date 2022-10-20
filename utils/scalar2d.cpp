#include <array>
#include <cctype>
#include <fstream>
#include <string>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <format/filename.hpp>
#include <misc/option_parse.hpp>

#include <boost/math/special_functions/sinc.hpp>

typedef double value_t;
typedef std::complex<value_t> complex_t;

constexpr value_t PI = 3.14159265358979323846264338;
constexpr value_t TWO_PI=6.28318530717958647688;
constexpr size_t res = 512;

inline std::string lower_case(const std::string& str) {
    std::string r(str);
    std::transform(r.begin(), r.end(), r.begin(),
                   [](unsigned char c){ return std::tolower(c); }
                  );
    return r;
}

inline value_t sqr(const value_t& x) {
    return x*x;
}
inline value_t cube(const value_t& x) {
    return x*x*x;
}

inline value_t mathworks1(value_t x, value_t y) {
    return 3.*exp(-sqr(y+1)-sqr(x))*sqr(x-1)-exp(-sqr(x-1)-sqr(y))/3 + exp(-sqr(x)-sqr(y))*(10*cube(x)-2*x+10*pow(y,5));
}

inline value_t l2norm(value_t x, value_t y) {
    return sqrt(x*x + y*y);
}

inline value_t gauss(value_t r) {
    return exp(-r*r/2)/sqrt(TWO_PI);
}

inline value_t sinc(value_t x, value_t y) {
    value_t r = sqrt(x*x + y*y);
    if (r==0) return 1;
    else return sin(r)/r;
}

inline value_t sincos(value_t x, value_t y) {
    return sin(TWO_PI*x)*cos(TWO_PI*y);
}

inline complex_t complex_sinc(value_t x, value_t y) {
    complex_t z(x, y);
    return boost::math::sinc_pi(z);
}

inline value_t z_marschner_lobb(value_t x, value_t y) {
    value_t r = sqrt(x*x + y*y);
    return 2/PI*asin(cos(12*PI*cos(PI*r/2))/4);
}

inline value_t circle(value_t x, value_t y) {
    return fabs(l2norm(x, y) - 1);
}

inline value_t radial(value_t x, value_t y) {
    return l2norm(x,y);
}

inline value_t disk(value_t x, value_t y) {
    value_t r = l2norm(x,y);
    if (r < 1) return 0;
    else return gauss(r - 1);
}

typedef std::array<value_t, 2> pos_t;
typedef std::array<value_t, 4> bnds_t;

inline bool valid_bounds(const bnds_t& bounds) {
    for (unsigned int i=0; i<2; ++i) {
        if (bounds[2*i+1] <= bounds[2*i]) return false;
    }
    return true;
}

inline pos_t random_pos(const bnds_t& bounds) {
    pos_t r;
    r[0] = bounds[0] + drand48()*(bounds[1]-bounds[0]);
    r[1] = bounds[2] + drand48()*(bounds[3]-bounds[2]);
    return r;
}

size_t npts = 100;
std::string fname = "circle", name_out;
size_t seedval = time(0);
bnds_t bnds({-1, 1, -1, 1});
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
    
    std::fstream output(filename + ".xyz", std::ios::out);
    output << "#" << npts << " rows - \"x y f(x,y)\"\n";
    for (size_t i=0; i<npts; ++i) {
        pos_t p = random_pos(bnds);
        const value_t& x = p[0];
        const value_t& y = p[1];
        value_t f;
        
        if (fname == "sinc" || fname == "s") 
            f = sinc(x,y);
        else if (fname == "real complex sinc" || fname == "real complex" || fname == "real" || fname == "r")
            f = std::real(complex_sinc(x,y));
        else if (fname == "imaginary complex sinc" || fname == "imaginary complex" || fname == "imaginary" || fname == "imag complex" || fname == "imag" || fname == "i")
            f = std::imag(complex_sinc(x,y));
        else if (fname == "marschner-lobb" || fname == "marschner" || fname == "marschner lobb" || fname == "ml") 
            f = z_marschner_lobb(x,y);
        else if (fname == "circle" || fname == "c") 
            f = circle(x,y);
        else if (fname == "disk" || fname == "d")
            f = disk(x,y);
        else if (fname == "radial") 
            f = radial(x,y);
        else if (fname == "sincos") 
            f = sincos(x,y);
        else if (fname == "mathworks" || fname == "mathworks1")
            f = mathworks1(x,y);
        else {
            std::cerr << "Unrecognized function name: " << fname << '\n';
            exit(1);
        }
        
        output << x << " " << y << " " << f << '\n';
    }
    output.close();
    
    if (ground_truth) {
        std::vector<value_t> values(res*res);
        std::array<size_t, 2> dims({res, res});
        std::array<value_t, 2> spc({(bnds[1]-bnds[0])/(value_t)(res-1), (bnds[3]-bnds[2])/(value_t)(res-1)});
        std::array<value_t, 2> mins({bnds[0], bnds[2]});
        
        for (size_t i=0; i<res*res; ++i) {
            int v = i/res;
            int u = i%res;
            value_t x = mins[0] + u*spc[0];
            value_t y = mins[1] + v*spc[1];
            value_t& f = values[i];
            if (fname == "sinc" || fname == "s") 
                f = sinc(x,y);
            else if (fname == "real complex sinc" || fname == "real complex" || fname == "real" || fname == "r")
                f = std::real(complex_sinc(x,y));
            else if (fname == "imaginary complex sinc" || fname == "imaginary complex" || fname == "imaginary" || fname == "imag complex" || fname == "imag" || fname == "i")
                f = std::imag(complex_sinc(x,y));
            else if (fname == "marschner-lobb" || fname == "marschner" || fname == "marschner lobb" || fname == "ml") 
                f = z_marschner_lobb(x,y);
            else if (fname == "circle" || fname == "c") 
                f = circle(x,y);
            else if (fname == "disk" || fname == "d")
                f = disk(x,y);
            else if (fname == "radial") 
                f = radial(x,y);
            else if (fname == "sincos") 
                f = sincos(x,y);
            else if (fname == "mathworks" || fname == "mathworks1")
                f = mathworks1(x,y);
            else {
                std::cerr << "Unrecognized function name: " << fname << '\n';
                exit(1);
            }
        }
        
        std::cout << "computed values: min=" << *std::min_element(values.begin(), values.end()) << ", max=" << *std::max_element(values.begin(), values.end()) << '\n';
        
        
        spurt::nrrd_utils::writeNrrdFromContainers(&values[0], filename + ".nrrd", dims, spc, mins);
    }   
    
    return 0;
}