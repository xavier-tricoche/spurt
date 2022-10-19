#include <array>
#include <iostream>
#include <fstream>
#include <image/nrrd_wrapper.hpp>
#include <string>
#include <misc/option_parse.hpp>
#include <vector>
#include <math/fixed_vector.hpp>
#include <format/filename.hpp>
#include <sstream>
#include <limits>

std::array<size_t, 3> res;
std::array<double, 6> bounds;
std::array<double, 3> range;
std::string function_name, name_out;
bool verbose;

constexpr double PI=3.14159265358979323844;
constexpr double TWO_PI=6.28318530717958647688;


struct ABC_field {
    ABC_field(double a=sqrt(3), double b=sqrt(2), double c=1) : A(a), B(b), C(c) {}
        
    nvis::vec3 operator()(const nvis::vec3& p, double t=0) const {
		nvis::vec3 dxdt;
		const double& x = p[0];
		const double& y = p[1];
		const double& z = p[2];
        dxdt[0] = (A + 0.5*t*sin(PI*t))*sin(z) + C*cos(y);
        dxdt[1] = B*sin(x) + (A + 0.5*t*sin(PI*t))*cos(z);
        dxdt[2] = C*sin(y) + B*cos(x);
		
		return dxdt;
    }
    
    double A, B, C;
};

bool valid_bounds() {
	return (bounds[0]<bounds[1] &&
			bounds[2]<bounds[3] && 
			bounds[4]<bounds[5]);
}

bool valid_range() {
	return (range[2]>0 && range[0]<range[1] && range[2]<=(range[1]-range[0]));
}

void initialize(int argc, const char* argv[])
{
    namespace xcl = xavier::command_line;
        
    xcl::option_traits 
            required(true, false, "Required Options"), 
            optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Create a (possibly time-dependent) flow dataset");

    res = { 64, 64, 64 };
    bounds = { 1, -1, 1, -1, 1, -1 }; // invalid bounds
    verbose = false;
	range = {0, 0, -1};

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("name", function_name, "Function name", required);
        parser.add_value("output", name_out, "Output base name", required);
        parser.add_tuple<3>("range", range, "Time interval (min, max, dt)", optional);
        parser.add_tuple<3>("res", res, res, "Sampling resolution", optional);
        parser.add_tuple<6>("bounds", bounds, "Sampling bounds", optional);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional);
        
        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

void make_abc(const std::string& name, const nvis::bbox3& domain, const nvis::vec3& step, double t) {
	ABC_field ABC;

	double* data = (double *)calloc(3*res[0]*res[1]*res[2], sizeof(double));
	size_t npts = res[0]*res[1]*res[2];
	for (size_t n=0; n<npts; ++n) {
		// n = i + res[0]*(j + res[1]*k)
		size_t k(n);
		size_t i=k%res[0];
		k /= res[0];
		size_t j = k%res[1];
		k /= res[1];
		nvis::vec3 x = domain.min() + nvis::vec3(i,j,k)*step;
		nvis::vec3 f = ABC(x, t);
		data[3*n  ] = f[0];
		data[3*n+1] = f[1];
		data[3*n+2] = f[2];
	}
	
	std::vector<size_t> size(4);
	std::vector<double> spacing(4);
	std::vector<double> mins(4);
	size[0] = 3;
	size[1] = res[0];
	size[2] = res[1];
	size[3] = res[2];
	spacing[0] = std::numeric_limits<double>::quiet_NaN();
	spacing[1] = step[0];
	spacing[2] = step[1];
	spacing[3] = step[2];
	mins[0] = std::numeric_limits<double>::quiet_NaN();
	mins[1] = 0;
	mins[2] = 0;
	mins[3] = 0;
	xavier::nrrd_utils::writeNrrdFromContainers(data, name, size, spacing, mins);
}

int main(int argc, const char* argv[]) {
	
	initialize(argc, argv);
	
	if (function_name == "ABC") {
		nvis::bbox3 domain;
		domain.min() = nvis::vec3(0,0,0);
		domain.max() = nvis::vec3(TWO_PI, TWO_PI, TWO_PI);
		nvis::vec3 step = (domain.max() - domain.min())/nvis::vec3(res[0]-1, res[1]-1, res[2]-1);
				
		double mint, maxt, dt;
		if (!valid_range()) {
			make_abc(xavier::filename::replace_extension(name_out, "nrrd"), domain, step, 0);
		}
		else {
			for (double t=range[0]; t<=range[1]; t+=range[2]) {
				std::ostringstream os;
				os << xavier::filename::remove_extension(name_out) << "_t=" << std::setprecision(3) << t << ".nrrd";
				make_abc(os.str(), domain, step, t);
			}
		}
	}
	
	return 0;
}
