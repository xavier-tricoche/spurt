#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <fstream>
#include <new>

#include <flow/time_dependent_field.hpp>
#include <flow/vector_field.hpp>
#include <misc/progress.hpp>
#include <misc/strings.hpp>
#include <misc/option_parse.hpp>

#include <math/types.hpp>
#include <flow/vector_field.hpp>
#include <data/vtk_field.hpp>
#include <data/raster.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/filesystem.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif


std::string name_in, name_out, path;
double length, eps, T;
size_t dim, mem;
std::array<size_t, 3> res;
std::array<bool, 3> periodic;
std::array<double, 6> bounds;
bool verbose;
size_t max_rhs_evals = 1000000;

namespace odeint = boost::numeric::odeint;

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;

    xcl::option_traits
            required(true, false, "Required Options"),
            optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compute flow map in a vector field given as VTK file");

    mem = 1024;
    eps = 1.0e-6;
    res = { 64, 64, 64 };
    periodic = { false, false, false };
    bounds = { 1, -1, 1, -1, 1, -1 }; // invalid bounds
    verbose = false;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input filename", required);
        parser.add_value("output", name_out, "Output base name", required);
        parser.add_value("eps", eps, eps, "Integration precision", optional);
        parser.add_value("T", T, "Integration time", required);
        parser.add_value("maxeval", max_rhs_evals, max_rhs_evals, "Maximum number of RHS evaluations", optional);
        parser.add_tuple<3>("res", res, res, "Sampling resolution", optional);
        parser.add_tuple<3>("periodic", periodic, periodic, "Periodic boundary conditions", optional);
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

using namespace spurt;

typedef vtk_field                 field_type;  // steady 3D vector field
typedef field_type::point_type    point_type;  // position in 3D space
typedef field_type::vector_type   vector_type; // 3D vector
typedef field_type::scalar_type   scalar_type; // scalar (e.g., time)

struct RHS {
    RHS(const field_type& field, size_t max_evals=1000)
        : m_field(field), m_counter(0), m_max_evals(max_evals) {
    }

    void operator()(const point_type& x, vector_type& dxdt, scalar_type t) const {
        dxdt = m_field(x /*, t*/); // exception will be passed along if position is invalid
        ++m_counter;
    }

    const field_type& m_field;
    mutable size_t m_counter;
    size_t m_max_evals;
};
typedef RHS rhs_type;

struct Observer {
    Observer(point_type& seed, double& t, double& d) : last_p(seed), last_t(t), distance(d) {}
    void operator()(const point_type& p, double t) {
        distance += norm(last_p-p);
        last_p = p;
        last_t = t;
        // std::cout << "distance = " << distance << '\n';
    }
    
    point_type& last_p;
    double& last_t;
    double& distance;
};
typedef Observer observer_t;

int main(int argc, const char* argv[])
{
    using namespace spurt;
    using namespace odeint;

    initialize(argc, argv);

    size_t nb_threads = 1;
#if _OPENMP
    nb_threads = omp_get_max_threads();
#endif
    std::cout << nb_threads << " threads available\n";

    field_type field(name_in);
    RHS rhs(field);
    
    std::cerr << "Resolution = " << res[0] << " x " << res[1] << " x " << res[2] << std::endl;
    spurt::raster_grid<3> sampling_grid(res, field.bounds());

    std::cout << "sampling grid bounds are: " << sampling_grid.bounds().min()
    << " -> " << sampling_grid.bounds().max() << '\n';
    int npoints = sampling_grid.size();

    float* flowmap = (float*)calloc(3*npoints, sizeof(float));

    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';

    // initialize coordinates
#pragma openmp parallel
    for (int n=0 ; n<npoints ; ++n) {
        ivec3 c = sampling_grid.coordinates(n);
        vec3 x = sampling_grid(c);
        flowmap[3*n  ] = x[0];
        flowmap[3*n+1] = x[1];
        flowmap[3*n+2] = x[2];
    }
    std::vector<bool> stopped(npoints, false);

    int incr = (T > 0) ? +1 : -1;

   size_t nbcomputed=0;
   spurt::ProgressDisplay progress(true);

   size_t nb_lost = 0;
   progress.start(npoints);

#pragma omp parallel
   {
   #pragma omp for schedule(dynamic,1)
   for (int n = 0 ; n < npoints ; ++n) {
       #if _OPENMP
       const int thread=omp_get_thread_num();
       #else
       const int thread=0;
       #endif

       if (!thread) progress.update(n);
       
       progress.update(n);

       if (stopped[n]) {
           continue;
       }

       // create a stepper
       auto stepper = make_controlled(eps, eps, runge_kutta_dopri5<point_type>());

       point_type x(flowmap[3*n], flowmap[3*n+1], flowmap[3*n+2]);
       double t=0, d=0;
       point_type p(x);
       observer_t obs(p, t, d);
       try {
           point_type y(x);
           integrate_adaptive(stepper, rhs, y, static_cast<double>(0), T, 1.0e-4, obs);
           flowmap[3*n  ] = y[0];
           flowmap[3*n+1] = y[1];
           flowmap[3*n+2] = y[2];
           std::cout << "integration successful: final position=" << y << ", (" << obs.last_p << ", " << obs.distance << ")\n";
       }
       catch (std::exception& e) {
           if (verbose) {
               std::cerr << "caught exception while integrating from "
                     << x << ":" << e.what() << '\n';
               std::cout << "last position reached was " << obs.last_p << " at time " << obs.last_t << " for a total distance of " << obs.distance << '\n';
           }
           ++nb_lost;
           stopped[n] = true;
           flowmap[3*n  ] = obs.last_p[0];
           flowmap[3*n+1] = obs.last_p[1];
           flowmap[3*n+2] = obs.last_p[2];
       }
   }
   }

   progress.end();
   
   std::cout << "number of successful cell searches: " << field.n_found << "\n"
       << "number of failed searches: " << field.n_failed << '\n';

   std::vector<size_t> size(4);
   std::vector<double> step(4), mins(4);
   step[0] = std::numeric_limits<double>::quiet_NaN();
   mins[0] = std::numeric_limits<double>::quiet_NaN();
   size[0] = 3;
   for (int i = 0 ; i < 3 ; ++i) {
       size[i+1] = res[i];
       step[i+1] = sampling_grid.spacing()[i];
       mins[i+1] = sampling_grid.bounds().min()[i];
       mins[i+1] = sampling_grid.bounds().min()[i];
   }

   std::ostringstream os;
   os << name_out << "-flowmap-deltaT=" << T << ".nrrd";
   spurt::nrrd_utils::writeNrrdFromContainers(flowmap, os.str(), size, step, mins);

   std::cout << "done (5)\n";

   delete[] flowmap;

   std::cout << "done (6)\n";

   return 0;
}
