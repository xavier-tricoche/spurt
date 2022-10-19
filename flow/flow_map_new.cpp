#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <fstream>
#include <new>

#include <image/nrrd_wrapper.hpp>
#include <data/field_wrapper.hpp>
#include <flow/time_dependent_field.hpp>
#include <flow/vector_field.hpp>
#include <misc/time_helper.hpp>
#include <misc/strings.hpp>
#include <misc/option_parse.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <flow/vector_field.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/filesystem.hpp>

#include <data/raster.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

std::string name_in, name_out, path;
double length, eps, t0, T;
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
            "Compute flow map in time-dependent vector field given as stack of NRRD files");

    mem = 1024;
    eps = 1.0e-6;
    res = { 64, 64, 64 };
    periodic = { false, false, false };
    bounds = { 1, -1, 1, -1, 1, -1 }; // invalid bounds
    verbose = false;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "List of input files", required);
        parser.add_value("path", path, "Path to input files (to be prepended to filenames)", optional);
        parser.add_value("output", name_out, "Output base name", required);
        parser.add_value("t0", t0, "Initial time", required);
        parser.add_value("T", T, "Integration length", required);
        parser.add_value("eps", eps, eps, "Integration precision", optional);
        parser.add_value("mem", mem, mem, "Available memory (in MB)", optional);
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

typedef gage_vector_field                         interpolator_type;   // single thread steady 3D vector field
typedef time_dependent_field< interpolator_type > unsteady_field_type; // single thread unsteady (3+1D) vector field
typedef unsteady_field_type::point_type           point_type;          // position in 3D space
typedef unsteady_field_type::vector_type          vector_type;         // 3D vector
typedef unsteady_field_type::scalar_type          scalar_type;         // scalar (e.g., time)

struct RHS {
    RHS(std::shared_ptr<unsteady_field_type> field, size_t max_evals=1000)
        : m_field(field), m_counter(0), m_max_evals(max_evals) {
        std::ostringstream os;
    }

    void operator()(const point_type& x, vector_type& dxdt, scalar_type t) const {
        dxdt = m_field->operator()(x, t); // exception will be passed along if time or position are invalid
        ++m_counter;
    }

    std::shared_ptr<unsteady_field_type> m_field;
    mutable size_t m_counter;
    size_t m_max_evals;
};
typedef RHS rhs_type;

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

    std::vector<std::pair< double, std::string> > data_files;

    if (!path.empty() && path.back() != '/') path = path + '/';

    // read file names
    std::fstream in(name_in, std::ios::in);
    size_t nb_files;
    in >> nb_files;
    for (int i=0 ; i<nb_files ; ++i) {
        std::string name;
        double t;
        in >> t >> name;
        if (!path.empty()) name = path + name;
        data_files.push_back(std::make_pair(t, name));
        if (verbose) std::cerr << name << " at t=" << t << std::endl;
    }
    in.close();

    Nrrd* nin = spurt::readNrrd(data_files.begin()->second);
    std::cerr << "Resolution = " << res[0] << " x " << res[1] << " x " << res[2] << std::endl;
    spurt::raster_grid<3> sampling_grid(res, spurt::bounds<3>(nin));
    nrrdNuke(nin);

    std::cout << "sampling grid bounds are: " << sampling_grid.bounds().min()
    << " -> " << sampling_grid.bounds().max() << '\n';
    int npoints = sampling_grid.size();

    const double min_time = T > 0 ? t0 : t0+T;
    const double max_time = T > 0 ? t0+T : t0;

    std::vector< std::shared_ptr<spurt::nrrd_wrapper> > input_nrrds;
    for (auto name : data_files) {
        input_nrrds.push_back(std::make_shared<spurt::nrrd_wrapper>(name.second));
    }
    const size_t nb_steps = input_nrrds.size();
    
    std::vector< std::vector< std::shared_ptr<gage_vector_field> > > steps(nb_threads);
    std::vector< double > times(input_nrrds.size());

    for (size_t n=0; n<nb_threads; ++n) {
        steps[n].resize(nb_steps);
        for (size_t i=0; i<nb_steps; ++i) {
            times[i] = data_files[i].first;
            steps[n][i] = std::make_shared<gage_vector_field>(input_nrrds[i]->pointer(), data_files[i].second, false, periodic);
        }
    }
    
    std::vector< std::shared_ptr<rhs_type> > per_thread_rhs(nb_threads);
    for (size_t n=0; n<nb_threads; ++n) {
        std::shared_ptr<unsteady_field_type> uf = std::make_shared<unsteady_field_type>(steps[n], times);
        per_thread_rhs[n] = std::make_shared<rhs_type>(uf, max_rhs_evals);
    }

    double max_available_time = per_thread_rhs.front()->m_field->time_range().second;
    double min_available_time = per_thread_rhs.front()->m_field->time_range().first;

    float* flowmap = (float*)calloc(3*npoints, sizeof(float));

    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';

    double target_time = T > 0 ? max_time : min_time;
    double initial_time = T > 0 ? min_time : max_time;

    // initialize coordinates
#pragma openmp parallel
    for (int n=0 ; n<npoints ; ++n) {
        spurt::ivec3 c = sampling_grid.coordinates(n);
        spurt::vec3 x = sampling_grid(c);
        flowmap[3*n  ] = x[0];
        flowmap[3*n+1] = x[1];
        flowmap[3*n+2] = x[2];
    }
    std::vector<bool> stopped(npoints, false);

    int incr = (T > 0) ? +1 : -1;

   size_t nbcomputed=0;
   spurt::progress_display progress(true);

   size_t nb_lost = 0;
   progress.start(npoints);

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
       try {
           rhs_type& rhs = *(per_thread_rhs[thread]);
           
           point_type y(x);
           integrate_adaptive(stepper, rhs, y, initial_time, target_time, 1.0e-4);
           flowmap[3*n  ] = y[0];
           flowmap[3*n+1] = y[1];
           flowmap[3*n+2] = y[2];

           size_t nsamples = rhs.m_field->nb_samples();
           std::ostringstream os;
           // os << "integration from " << y << " required " << nsamples << " gage evaluations corresponding to " << rhs.m_counter << " samples" << std::endl;
           // std::cout << os.str();
           rhs.m_field->reset_counter();
       }
       catch (std::exception& e) {
           std::cerr << "caught exception while integrating from "
                     << x << ":" << e.what() << '\n';
           ++nb_lost;
           stopped[n] = true;
       }
   }
   }

   progress.stop();

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
   os << name_out << "-flowmap-t0=" << t0 << "-T=" << T << ".nrrd";
   spurt::writeNrrdFromContainers(flowmap, os.str(), size, step, mins);
   delete[] flowmap;

   return 0;
}
