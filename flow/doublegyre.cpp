#include <iostream>
#include <sstream>

#include <boost/numeric/odeint.hpp>

#include <image/nrrd_wrapper.hpp>
#include <data/field_wrapper.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <data/raster.hpp>
#include <misc/option_parse.hpp>
#include <misc/time_helper.hpp>

#include <math/vector_manip.hpp>

#include <flow/ftle_rhs.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace odeint = boost::numeric::odeint;

typedef double value_t;

constexpr value_t PI=3.14159265358979323844;
constexpr value_t TWO_PI=6.28318530717958647688;

typedef spurt::fixed_vector< value_t, 2> state_t;
typedef spurt::fixed_vector< value_t, 2> pos_t;
typedef nvis::bounding_box< pos_t > bbox_t;

std::string name_out;
value_t t_max=25., eps=1.0e-8;
std::array< size_t, 2 > res={ 512, 256 };
std::array< value_t, 4 > bnds={ 0, 2, 0, 1 };
bool verbose=false;

bbox_t to_bbox(const std::array<value_t, 4>& array) {    
    bbox_t b;
    b.min()=pos_t(array[0], array[2]);
    b.max()=pos_t(array[1], array[3]);
    return b;
}

template<typename T, size_t N>
spurt::fixed_vector<T, N> to_vec(const std::array<T, N>& array) {
    spurt::fixed_vector<T, N> v;
    for (size_t i=0; i<N; ++i) v[i]=array[i];
    return v;
}

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;
        
    xcl::option_traits 
            required_group(true, false, "Required Options"), 
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compute flow map of double gyre flow");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("output", name_out, "Output filename", required_group);
        parser.add_value("T", t_max, t_max, "Integration length", optional_group);
        parser.add_value("eps", eps, eps, "Integration precision", optional_group);
        parser.add_tuple<2>("res", res, res, "Sampling resolution", optional_group);
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

int main(int argc, const char* argv[])
{
    using namespace spurt;
    using namespace odeint;
    
    typedef spurt::raster_grid<2> grid_t;
    typedef grid_t::coord_type coord_t;
    
    initialize(argc, argv);
    
    spurt::DoubleGyre<value_t, state_t> rhs;
    
    if (verbose) std::cout << "Resolution = " << res[0] << "x" << res[1] << std::endl;
     grid_t sampling_grid(to_vec(res), to_bbox(bnds));
          
    size_t npoints = sampling_grid.size();
        
    if (verbose) std::cout << "nb points = " << npoints << '\n';
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    value_t* fmap = (value_t*)calloc(2 * npoints, sizeof(value_t));
    
    int counter = 0;
    
    
    // runge_kutta_dopri5<state_t> stepper;
    value_t dt = 1.0e-2;
    
    spurt::progress_display progress(true);
    
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
            
            coord_t c = sampling_grid.coordinates(n);
            state_t x = sampling_grid(c);
            
            // create a stepper
            auto stepper = make_controlled(eps, eps, runge_kutta_dopri5<spurt::vec2>());
            
            integrate_adaptive(stepper, rhs, x,
                               static_cast<value_t>(0), t_max, dt);
                            
            fmap[2*n  ]=x[0];
            fmap[2*n+1]=x[1];
        }
    }
    progress.stop();
            
    spurt::nrrd_params<value_t, 3> nrrd_params;
    nrrd_params.mins()[0] = AIR_NAN;
    spurt::vector::copy(sampling_grid.bounds().min(), nrrd_params.mins(), 0, 1);
    nrrd_params.spacings()[0] = AIR_NAN;
    spurt::vector::copy(sampling_grid.spacing(), nrrd_params.spacings(), 0, 1);
    nrrd_params.sizes()[0]=2;
    spurt::vector::copy(sampling_grid.resolution(), nrrd_params.sizes(), 0, 1);
    nrrd_params.centers().fill(nrrdCenterNode);
    nrrd_params.centers()[0]=nrrdCenterUnknown;
    nrrd_params.labels() = {{"flow map", "x", "y"}};
    std::ostringstream os;
    os << "Double gyre flow map: RK45 DOPRI, eps=" << eps << ", Tmax=" << t_max;
    nrrd_params.comments().push_back(os.str());
    
    os.clear();
    os.str("");
    
    os << name_out << "-flowmap-e=" << eps << "-T=" << std::setw(4) << std::setfill('0') << t_max << ".nrrd";
    spurt::writeNrrdFromParams(fmap, os.str(), nrrd_params);
    
    return 0;
}
