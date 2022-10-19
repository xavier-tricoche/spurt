#include <iostream>
#include <sstream>

#include <boost/numeric/odeint.hpp>

#include <image/nrrd_wrapper.hpp>
#include <data/field_wrapper.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <data/raster.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <format/filename.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace odeint = boost::numeric::odeint;

typedef double value_t;

constexpr value_t PI=3.14159265358979323844;
constexpr value_t TWO_PI=6.28318530717958647688;

typedef nvis::fixed_vector< value_t, 3> state_t;
typedef nvis::fixed_vector< value_t, 3 > pos_t;
typedef nvis::bounding_box< pos_t > bbox_t;

std::string name_out;
value_t t_max=10., eps=1.0e-8;
std::array< size_t, 3 > res={ 256, 256, 256 };
std::array< value_t, 6 > bnds={ -PI, PI, -PI, PI, -0.001, 0.001};
std::array< value_t, 3 > abc={ sqrt(3), sqrt(2), 1};
bool verbose=false;


bbox_t to_bbox(const std::array<value_t, 6>& array) {    
    bbox_t b;
    b.min()=pos_t(array[0], array[2], array[4]);
    b.max()=pos_t(array[1], array[3], array[5]);
    return b;
}

template<typename T, size_t N>
nvis::fixed_vector<T, N> to_vec(const std::array<T, N>& array) {
    nvis::fixed_vector<T, N> v;
    for (size_t i=0; i<N; ++i) v[i]=array[i];
    return v;
}

void initialize(int argc, const char* argv[])
{
    namespace xcl = xavier::command_line;
        
    xcl::option_traits 
            required_group(true, false, "Required Options"), 
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compute flow map of ABC flow");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("output", name_out, "Output filename", required_group);
        parser.add_value("T", t_max, t_max, "Integration length", optional_group);
        parser.add_value("eps", eps, eps, "Integration precision", optional_group);
        parser.add_tuple<3>("cst", abc, abc, "ABC constants", optional_group);
        parser.add_tuple<3>("res", res, res, "Sampling resolution", optional_group);
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

struct ABC_field {
    ABC_field(value_t a, value_t b, value_t c) : A(a), B(b), C(c) {}
        
    void operator()(const state_t& x, state_t& dxdt, value_t) const {
        dxdt[0] = A*sin(x[2]) + C*cos(x[1]);
        dxdt[1] = B*sin(x[0]) + A*cos(x[2]);
        dxdt[2] = C*sin(x[1]) + B*cos(x[0]);
    }
    
    value_t A, B, C;
};

int main(int argc, const char* argv[])
{
    using namespace xavier;
    using namespace odeint;
    
    initialize(argc, argv);
    
    name_out=xavier::filename::remove_extension(name_out);
    
    ABC_field rhs(abc[0], abc[1], abc[2]);
    
    if (verbose) std::cout << "Resolution = " << res[0] << "x" << res[1] << "x" << res[2] << std::endl;
    xavier::raster_grid<3> sampling_grid(to_vec(res), to_bbox(bnds));
          
    size_t npoints = sampling_grid.size();
        
    if (verbose) std::cout << "nb points = " << npoints << '\n';
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    value_t* fmap = (value_t*)calloc(3 * npoints, sizeof(value_t));
    
    int counter = 0;
    
    // create a stepper
    runge_kutta_dopri5<state_t> stepper;
    value_t dt = 1.0e-2;
    
    xavier::ProgressDisplay progress(true);
    
    progress.start(npoints, "Computing flow map");
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (size_t n = 0 ; n < npoints ; ++n) {
        
            #pragma omp atomic
            ++counter;
            
#if _OPENMP
            const int thread=omp_get_thread_num();
#else
            const int thread=0;
#endif
            if (!thread) progress.update(counter);
            
            nvis::ivec3 c = sampling_grid.coordinates(n);
            state_t x = sampling_grid(c);
            
            integrate_adaptive(make_controlled(eps, eps, stepper), rhs, x,
                               static_cast<value_t>(0), t_max, dt);
                            
            fmap[3*n  ]=x[0];
            fmap[3*n+1]=x[1];
            fmap[3*n+2]=x[2];
        }
    }
    progress.end();
    
    std::vector<size_t> size(4);
    std::vector<double> step(4);
    std::vector<double> mins(4);
    const nvis::vec3& s = sampling_grid.spacing();
    step[0] = AIR_NAN;
    mins[0] = AIR_NAN;
    for (int i = 0 ; i < 3 ; ++i) {
        size[i+1] = res[i];
        step[i+1] = s[i];
        mins[i+1] = bnds[2*i];
    }
    
    std::ostringstream os;
    
    size[0] = 3;
    os << name_out << "-flowmap-e=" << eps << "-T=" << std::setw(4) << std::setfill('0') << t_max << ".nrrd";
    xavier::nrrd_utils::writeNrrdFromContainers(fmap, os.str(), /*nrrd_value_traits_from_type<value_t>::index,*/ size, step, mins);
    
    return 0;
}
