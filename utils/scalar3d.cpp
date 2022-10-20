#include <array>
#include <cctype>
#include <fstream>
#include <limits>
#include <string>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <format/filename.hpp>
#include <misc/option_parse.hpp>

typedef double value_t;

constexpr value_t PI = 3.14159265358979323846264338;
constexpr value_t TWO_PI=6.28318530717958647688;

typedef std::array<value_t, 3> pos_t;
typedef std::array<value_t, 3> vec_t;
typedef std::array<value_t, 6> bnds_t;
typedef std::array<value_t, 6> smat_t;

size_t npts = 100;
std::string fname = "sphere", name_out;
size_t seedval = time(0);
bnds_t bnds({-1, 1, -1, 1, -1, 1});
bool verbose = false;
int order = 0;
value_t h=1.0e-6;

template<unsigned long N>
inline std::array<value_t, N> operator+(const std::array<value_t, N>& a, const std::array<value_t, N>& b) {
    std::array<value_t, N> r;
    for (unsigned int i=0; i<N; ++i) r[i] = a[i] + b[i];
    return r;
}
template<unsigned long N>
inline std::array<value_t, N> operator-(const std::array<value_t, N>& a, const std::array<value_t, N>& b) {
    std::array<value_t, N> r;
    for (unsigned int i=0; i<N; ++i) r[i] = a[i] - b[i];
    return r;
}

inline std::string lower_case(const std::string& str) {
    std::string r(str);
    std::transform(r.begin(), r.end(), r.begin(),
                   [](unsigned char c){ return std::tolower(c); }
                  );
    return r;
}

value_t gauss(value_t r) {
    return exp(-r*r/2)/sqrt(TWO_PI);
}

inline bool valid_bounds(const bnds_t& bounds) {
    for (unsigned int i=0; i<3; ++i) {
        if (bounds[2*i+1] <= bounds[2*i]) return false;
    }
    return true;
}

inline pos_t random_pos(const bnds_t& bounds) {
    pos_t r;
    r[0] = bounds[0] + drand48()*(bounds[1]-bounds[0]);
    r[1] = bounds[2] + drand48()*(bounds[3]-bounds[2]);
    r[2] = bounds[4] + drand48()*(bounds[5]-bounds[4]);
    return r;
}

inline value_t l2norm(value_t x, value_t y, value_t z) {
    return sqrt(x*x + y*y + z*z);
}
inline value_t l2norm(const pos_t& p) {
    return l2norm(p[0], p[1], p[2]);
}

inline value_t sphere(value_t x, value_t y, value_t z) {
    return l2norm(x, y, z);
}
inline value_t sphere(const pos_t& p) {
    return sphere(p[0], p[1], p[2]);
}

inline value_t ball(value_t x, value_t y, value_t z) {
    value_t r = l2norm(x,y,z);
    if (r > 1) return 0;
    else return gauss(r-1);
}

inline value_t ball(const pos_t& p) {
    return ball(p[0], p[1], p[2]);
}

inline value_t ellipsoid(value_t x, value_t y, value_t z) {
    return l2norm(sqrt(3)*x, sqrt(2)*y, sqrt(5)*z);
}

inline value_t ellipsoid(const pos_t& p) {
    return ellipsoid(p[0], p[1], p[2]);
}

inline value_t torus(value_t x, value_t y, value_t z) {
    double theta = atan2(y, x);
    double px = cos(theta);
    double py = sin(theta);
    return l2norm(x-px, y-py, z*z);
}

inline value_t torus(const pos_t& p) {
    return torus(p[0], p[1], p[2]);
}

inline value_t helix(value_t x, value_t y, value_t z) {
    std::vector<double> values;
    double theta = atan2(y, x);
    double px = cos(theta);
    double py = sin(theta);
    double pz = 0.5*theta;
    return l2norm(x-px, y-py, z-pz);
}

inline value_t helix(const pos_t& p) {
    return helix(p[0], p[1], p[2]);
}

inline vec_t first_derivative(value_t (*f)(const pos_t&), const pos_t& p) {
    vec_t r;
    pos_t dx = {h, 0, 0};
    pos_t dy = {0, h, 0};
    pos_t dz = {0, 0, h};
    r[0] = (f(p+dx)-f(p-dx))/(2*h);
    r[1] = (f(p+dy)-f(p-dy))/(2*h);
    r[2] = (f(p+dz)-f(p-dz))/(2*h);
    return r;
}

inline std::array<value_t, 6> second_derivative(value_t (*f)(const pos_t&),
    const pos_t& p) {
    smat_t r;
    pos_t dx = {h, 0, 0};
    pos_t dy = {0, h, 0};
    pos_t dz = {0, 0, h};
    r[0] = (f(p-dx) - 2*f(p) + f(p+dx))/(h*h);
    r[1] = (f(p+dx+dy) - f(p-dx+dy) - f(p+dx-dy) + f(p-dx-dy))/(4*h*h);
    r[2] = (f(p+dx+dz) - f(p-dx+dz) - f(p+dx-dz) + f(p-dx-dz))/(4*h*h);
    r[3] = (f(p-dy) - 2*f(p) + f(p+dy))/(h*h);
    r[4] = (f(p+dy+dz) - f(p-dy+dz) - f(p+dy-dz) + f(p-dy-dz))/(4*h*h);
    r[5] = (f(p-dz) - 2*f(p) + f(p+dz))/(h*h);
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

void evaluate(value_t& v, vec_t& df, smat_t& d2f, int order, value_t (*f)(const pos_t&), const pos_t& p) {
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
        pos_t p = random_pos(bnds);
        const value_t& x = p[0];
        const value_t& y = p[1];
        const value_t& z = p[2];

        double f;
        vec_t df;
        smat_t d2f;

        if (fname == "sphere" || fname == "s")
            evaluate(f, df, d2f, order, sphere, p);
        else if (fname == "ball" || fname == "b")
            evaluate(f, df, d2f, order, ball, p);
        else if (fname == "ellipsoid" || fname == "e")
            evaluate(f, df, d2f, order, ellipsoid, p);
        else if (fname == "torus" || fname == "t")
            evaluate(f, df, d2f, order, torus, p);
        else if (fname == "helix" || fname == "t")
            evaluate(f, df, d2f, order, helix, p);
        else {
            std::cerr << "Unrecognized function name: " << fname << '\n';
            exit(1);
        }

        output << x << " " << y << " " << z << " " << f;
        if (order >=1) output << " " << df[0] << " " << df[1] << " " << df[2];
        if (order >=2) output << " " << d2f[0] << " "<< d2f[1] << " "<< d2f[2] << " " << d2f[3] << " "<< d2f[4] << " "<< d2f[5];
        output << '\n';
    }
    output.close();

    return 0;
}
