#include <array>
#include <cctype>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <format/filename.hpp>
#include <misc/option_parse.hpp>
#include <math/types.hpp>

#include <utils/functions.hpp>

using namespace spurt;

size_t npts = 100;
std::string fname = "sphere", name_out;
size_t seedval = time(0);
vec6 bnds({-1, 1, -1, 1, -1, 1});
bool verbose = false;
int order = 0;
double h=1.0e-6;

inline std::string lower_case(const std::string& str) {
    std::string r(str);
    std::transform(r.begin(), r.end(), r.begin(),
                   [](unsigned char c){ return std::tolower(c); }
                  );
    return r;
}

inline bool valid_bounds(const vec6& bounds) {
    for (unsigned int i=0; i<3; ++i) {
        if (bounds[2*i+1] <= bounds[2*i]) return false;
    }
    return true;
}

inline vec3 random_pos(const vec6& bounds) {
    vec3 r;
    r[0] = bounds[0] + drand48()*(bounds[1]-bounds[0]);
    r[1] = bounds[2] + drand48()*(bounds[3]-bounds[2]);
    r[2] = bounds[4] + drand48()*(bounds[5]-bounds[4]);
    return r;
}

inline vec3 first_derivative(double (*f)(const vec3&), const vec3& p) {
    vec3 r;
    vec3 dx(h, 0, 0);
    vec3 dy(0, h, 0);
    vec3 dz(0, 0, h);
    r[0] = (f(p+dx)-f(p-dx))/(2*h);
    r[1] = (f(p+dy)-f(p-dy))/(2*h);
    r[2] = (f(p+dz)-f(p-dz))/(2*h);
    return r;
}

inline mat3 second_derivative(double (*f)(const vec3&),
    const vec3& p) {
    mat3 r;
    vec3 dx(h, 0, 0);
    vec3 dy(0, h, 0);
    vec3 dz(0, 0, h);
    r(0,0) = (f(p-dx) - 2*f(p) + f(p+dx))/(h*h);
    r(1,0) = r(0,1) = (f(p+dx+dy) - f(p-dx+dy) - f(p+dx-dy) + f(p-dx-dy))/(4*h*h);
    r(2,0) = r(0,2) = (f(p+dx+dz) - f(p-dx+dz) - f(p+dx-dz) + f(p-dx-dz))/(4*h*h);
    r(1,1) = (f(p-dy) - 2*f(p) + f(p+dy))/(h*h);
    r(2,1) = r(1,2) = (f(p+dy+dz) - f(p-dy+dz) - f(p+dy-dz) + f(p-dy-dz))/(4*h*h);
    r(2,2) = (f(p-dz) - 2*f(p) + f(p+dz))/(h*h);
    return r;
}

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;

    xcl::option_traits
            required_group(true, false, "Required Options"),
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Create scattered samples of 3D scalar field");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("output", name_out, "Output filename", required_group);
        parser.add_value("number", npts, npts, "Number of sample points", optional_group);
        parser.add_value("func", fname, fname, "function name", optional_group);
        parser.add_value("seed", seedval, seedval, "Random seed value", optional_group);
        parser.add_value("deriv", order, order, "Derivative order", optional_group);
        parser.add_value("step", h, h, "Discretization step", optional_group);
        parser.add_tuple<6>("bounds", bnds, bnds, "Sampling bounds", optional_group);
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

void evaluate(double& v, vec3& df, mat3& d2f, int order, double (*f)(const vec3&), const vec3& p) {
    v = f(p);
    if (order >= 1) df = first_derivative(f, p);
    if (order >= 2) d2f = second_derivative(f, p);
}


int main(int argc, char* argv[]) {

    initialize(argc, const_cast<const char**>(argv));

    srand48(seedval);

    std::string filename = spurt::filename::remove_extension(name_out);

    fname = lower_case(fname);

    std::fstream output(filename + ".xyz", std::ios::out);
    output << "#" << npts << " rows - \"x y z f(x,y,z)";
    if (order>=1) output << " df/dx df/dy df/dz";
    if (order==2) output << " d2f/dx2 d2f/dxdy d2f/dxdz d2f/dy2 d2f/dydz d2f/dz2";
    output << "\"\n";
    for (size_t i=0; i<npts; ++i) {
        vec3 p = random_pos(bnds);

        double f;
        vec3 df;
        mat3 d2f;

        if (fname == "sphere" || fname == "s")
            evaluate(f, df, d2f, order, functions3d::sphere, p);
        else if (fname == "ball" || fname == "b")
            evaluate(f, df, d2f, order, functions3d::ball, p);
        else if (fname == "ellipsoid" || fname == "e")
            evaluate(f, df, d2f, order, functions3d::ellipsoid, p);
        else if (fname == "torus" || fname == "t")
            evaluate(f, df, d2f, order, functions3d::torus, p);
        else if (fname == "helix" || fname == "t")
            evaluate(f, df, d2f, order, functions3d::helix, p);
        else {
            std::cerr << "Unrecognized function name: " << fname << '\n';
            exit(1);
        }

        output << p << " " << f << '\n';
        if (order >=1) output << df << '\n';
        if (order >=2) output << '\n' << d2f << '\n';
    }
    output.close();

    return 0;
}
