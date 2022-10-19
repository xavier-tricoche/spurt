#include <iostream>
#include <map>
#include <queue>
#include <sstream>
#include <fstream>
#include <new>
#include <iterator>

// #include <flow/time_dependent_field.hpp>
#include <flow/vector_field.hpp>
#include <misc/progress.hpp>
#include <misc/strings.hpp>
#include <misc/option_parse.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <flow/vector_field.hpp>
#include <data/raster.hpp>
#include <VTK/vtk_interpolator.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/filesystem.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
tbb::atomic<size_t> tbb_progress_counter;


std::string name_in, name_out, path, cmdline;
double length, eps, T;
size_t dim, mem;
std::array<size_t, 3> res;
std::array<bool, 3> periodic;
std::array<double, 6> bounds;
bool verbose;
size_t max_rhs_evals = 1000000;
size_t nb_threads = 1;
bool save_lines = false;
bool monitor = false;
bool confine = false;
nvis::bbox3 the_bounds;
bool in_parallel=false;

namespace odeint = boost::numeric::odeint;

template< typename T, typename T1, typename T2=void, typename T3=void, typename T4=void >
struct is_one_of_those_types : public std::false_type {};

template<typename T, typename T1>
struct is_one_of_those_types< T, T1, typename std::enable_if< std::is_same<T, T1>::value >::type > : public std::true_type {};

template<typename T, typename T1, typename T2>
struct is_one_of_those_types< T, T1, T2, typename std::enable_if< std::is_same<T, T1>::value || std::is_same<T, T2>::value >::type > : public std::true_type {};

template<typename T, typename T1, typename T2, typename T3>
struct is_one_of_those_types< T, T1, T2, T3,
            typename std::enable_if<
                std::is_same<T, T1>::value ||
                std::is_same<T, T2>::value ||
                std::is_same<T, T3>::value >::type > : public std::true_type {};


void initialize(int argc, const char* argv[])
{
    namespace xcl = xavier::command_line;

    cmdline = "Command line: " + std::string(argv[0]);
    for (int i=1; i<argc; i++) {
        cmdline += " " + std::string(argv[i]);
    }

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
        parser.add_value("parallel", in_parallel, in_parallel, "Compute flow map in parallel", optional);
        parser.add_tuple<3>("periodic", periodic, periodic, "Periodic boundary conditions", optional);
        parser.add_tuple<6>("bounds", bounds, "Sampling bounds", optional);
        parser.add_value("confine", confine, confine, "Confine integration to prescribed bounds", optional);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional);
        parser.add_value("geometry", save_lines, save_lines, "Save streamlines geometry", optional);
        parser.add_value("monitor", monitor, monitor, "Save geometry of problematic streamlines", optional);

        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n"
            << e.what() << "\n"
                << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

using namespace xavier;
using namespace vtk_utils;

struct vec3 : public Eigen::Vector3d {
    typedef Eigen::Vector3d base_type;
    typedef base_type::value_type value_type;

    vec3() : base_type(0,0,0) {}
    vec3(const base_type& v) : base_type(v) {}
    vec3(double x, double y, double z) : base_type(x,y,z) {}

#if EIGEN_VERSION_AT_LEAST(3,4,0)
    // starting with version 3.4, Eigen offers STL-type iterators
    typedef base_type::iterator iterator;
    typedef base_type::const_iterator const_iterator;
#else
    typedef const double* const_iterator;
    typedef double* iterator;

    const_iterator begin() const {
        return &this->operator()(0);
    }
    iterator begin() {
        return &this->operator()(0);
    }
    const_iterator end() const {
        return begin()+3;
    }
    iterator end() {
        return begin()+3;
    }
#endif

    nvis::vec3& as_nvis_vec3() { return *((nvis::vec3*)(&base_type::operator()(0))); }
    const nvis::vec3& as_nvis_vec3() const { return *((const nvis::vec3*)(&base_type::operator()(0))); }
    // nvis::vec3 as_nvis_vec3() const { return nvis::vec3(this->operator()(0), this->operator()(1), this->operator()(2)); }
};


inline std::string to_str(const vec3& p) {
    std::ostringstream os;
    os << "[" << p(0) << ", " << p(1) << ", " << p(2) << "]";
    return os.str();
}

typedef vec3 point_type;
typedef vec3 vector_type;
typedef double scalar_type;

typedef interpolator<vtkImageData, double, 3, vector_type>  img_intp_t;
typedef interpolator<vtkRectilinearGrid, double, 3, vector_type> rect_intp_t;
typedef interpolator<vtkUnstructuredGrid, double, 3, vector_type> unst_intp_t;
typedef interpolator<vtkStructuredGrid, double, 3, vector_type> curv_intp_t;

struct point_out_of_bounds_exception : public std::runtime_error {
    typedef std::runtime_error base_type;
    explicit point_out_of_bounds_exception( const std::string& what_arg) : base_type(what_arg) {}
    explicit point_out_of_bounds_exception( const char* what_arg ) : base_type(what_arg) {}
};

template<typename Interpolator, typename Enable=void>
struct rhs_type {};

template<typename Interpolator>
struct rhs_type<
    Interpolator,
    typename std::enable_if<
        std::is_same<typename Interpolator::dataset_type, vtkRectilinearGrid>::value ||
        std::is_same<typename Interpolator::dataset_type, vtkImageData>::value
    >::type > {
    typedef Interpolator intp_t;

    rhs_type(const intp_t& intp, size_t max_evals=1000, bool forward=true, bool _verbose=false)
    : m_intp(intp), m_counter(0), m_max_evals(max_evals),
      m_orient(forward ? 1 : -1), m_verbose(_verbose) {
        if (m_verbose) {
            std::ostringstream os;
            os << "RHS initialized" << std::endl;
            std::cout << os.str();
        }
    }

    void operator()(const vec3& x, vec3& dxdt, scalar_type t) const {
        std::ostringstream os;
        if (m_verbose) {
            os << "rhs called at " << to_str(x) << "at t=" << t << '\n';
        }
        if ((confine && !the_bounds.inside(x.as_nvis_vec3())) ||
            !m_intp.interpolate(dxdt, x, m_verbose)) {
            throw point_out_of_bounds_exception("invalid position: " + to_str(x) + " at t=" + to_string(t));
        }
        if (m_verbose) {
            os << "interpolated value: " << to_str(dxdt) << std::endl;
            std::cout << os.str();
        }
        dxdt *= m_orient;
        ++m_counter;
    }

    const intp_t& m_intp;
    mutable size_t m_counter;
    const size_t m_max_evals;
    int m_orient;
    bool m_verbose;
};

template<typename Interpolator>
struct rhs_type<
    Interpolator,
    typename std::enable_if<
        std::is_same<typename Interpolator::dataset_type, vtkStructuredGrid >::value ||
        std::is_same<typename Interpolator::dataset_type, vtkUnstructuredGrid>::value
    >::type > {
    typedef Interpolator intp_t;

    rhs_type(const intp_t& intp, size_t max_evals=1000, bool forward=true,
             bool _verbose=false)
    : m_intp(intp), m_counter(0), m_max_evals(max_evals), m_verbose(_verbose),
      m_orient(forward ? 1 : -1) {
        m_cell = vtkGenericCell::New();
    }

    rhs_type(const intp_t& intp, vtkGenericCell* acell, size_t max_evals=1000, bool forward=true, bool _verbose=false)
        : m_intp(intp), m_cell(acell), m_counter(0), m_max_evals(max_evals),
          m_orient(forward ? 1 : -1), m_verbose(_verbose) {}

    rhs_type(const rhs_type& other)
        : m_intp(other.m_intp), m_counter(other.m_counter),
          m_cell(other.m_cell), m_max_evals(other.m_max_evals),
          m_orient(other.m_orient), m_verbose(other.m_verbose) {}

    void operator()(const vec3& x, vec3& dxdt, scalar_type t) const {
        if ((confine && !the_bounds.inside(x.as_nvis_vec3())) ||
            !m_intp.interpolate(dxdt, x)) {
            throw point_out_of_bounds_exception("invalid position: " + to_str(x) + " at t=" + to_string(t));
        }
        dxdt *= m_orient;
        ++m_counter;
    }

    const intp_t& m_intp;
    mutable size_t m_counter;
    const size_t m_max_evals;
    bool m_verbose;
    int m_orient;
    VTK_SMART(vtkGenericCell) m_cell;
};



namespace boost { namespace numeric { namespace odeint {

template<
class State ,
class Value = double ,
class Deriv = State ,
class Time = Value ,
class Algebra = typename algebra_dispatcher< State >::algebra_type ,
class Operations = typename operations_dispatcher< State >::operations_type ,
class Resizer = initially_resizer
>
class boundary_aware_stepper : public runge_kutta_dopri5<State, Value, Deriv, Time, Algebra, Operations, Resizer> {
public:

    typedef runge_kutta_dopri5<State, Value, Deriv, Time, Algebra, Operations, Resizer> base_type;
    typedef typename base_type::stepper_base_type stepper_base_type;
    typedef State state_type;
    typedef Deriv deriv_type;
    typedef Value value_type;
    typedef Time time_type;
    typedef Algebra algebra_type;
    typedef Operations operations_type;
    typedef Resizer resizer_type;
    typedef unsigned short order_type;
    typedef typename base_type::stepper_type stepper_type;
    typedef typename base_type::wrapped_state_type wrapped_state_type;
    typedef typename base_type::wrapped_deriv_type wrapped_deriv_type;


    boundary_aware_stepper() : base_type() {}

    template< class System , class StateIn , class DerivIn , class StateOut , class DerivOut , class Err >
    void do_step_impl( System system , const StateIn &in , const DerivIn &dxdt_in , time_type t ,
            StateOut &out , DerivOut &dxdt_out , time_type dt , Err &xerr )
    {
        try {
            base_type::template do_step_impl<System, StateIn, DerivIn, StateOut, DerivOut, Err>(out, dxdt_out, dt, xerr);
        }
        catch (point_out_of_bounds_exception& e) {
            // integrator stepped out of domain: use error to force step size reduction
            xerr *= 1000;
        }
    }
};

template< class State , class Value , class Deriv , class Time , class Algebra , class Operations , class Resize >
struct get_controller< boundary_aware_stepper< State , Value , Deriv , Time , Algebra , Operations , Resize > >
{
    typedef boundary_aware_stepper< State , Value , Deriv , Time , Algebra , Operations , Resize > stepper_type;
    typedef controlled_runge_kutta< stepper_type > type;
};

}}}


struct Observer {
    Observer(vec3& seed, double& t, double& d, std::vector<vec3>& c,
    std::vector<double>& ts, bool _verbose=false)
    : last_p(seed), last_t(t), distance(d), curve(c), times(ts), m_verbose(_verbose) {
        if (monitor || save_lines) {
            curve.push_back(seed);
            ts.push_back(t);
        }
        if (m_verbose) {
            std::ostringstream os;
            os << "Observer: seeding at " << to_str(seed) << " at time " << t << std::endl;
            std::cout << os.str();
        }
    }
    void operator()(const vec3& p, double t) {
        distance += (last_p-p).norm();
        last_p = p;
        last_t = t;
        if (monitor || save_lines) {
            curve.push_back(p);
            times.push_back(t);
        }
        if (m_verbose) {
            std::ostringstream os;
            os << "\nObserver: p=" << to_str(p) << ", t=" << t << std::endl;
            std::cout << os.str();
        }
    }

    vec3& last_p;
    double& last_t;
    double& distance;
    std::vector<vec3>& curve;
    std::vector<double>& times;
    bool m_verbose;
};
typedef Observer observer_t;

template<typename DataSet>
int run(VTK_SMART(DataSet) dataset) {
    typedef interpolator<DataSet> intp_t;
    typedef odeint::boundary_aware_stepper<vec3> stepper_type;

    if (dataset->GetPointData()->GetVectors() == nullptr) {
        std::cerr << "This dataset does not contain a vector field\n";
        std::cerr << "Exiting Now.\n";
        exit(1);
    }

    intp_t intp(dataset);

    double global_bounds[6];
    dataset->GetBounds(global_bounds);
    nvis::bbox3 bnds;

    if (::bounds[0]<::bounds[1] && ::bounds[2]<::bounds[3] && ::bounds[4]<::bounds[5] &&
        ::bounds[0] >= global_bounds[0] && ::bounds[1] <= global_bounds[1] &&
        ::bounds[2] >= global_bounds[2] && ::bounds[3] <= global_bounds[3] &&
        ::bounds[4] >= global_bounds[4] && ::bounds[5] <= global_bounds[5]) {
        // valid bounds supplied by user
        bnds.min() = nvis::vec3(::bounds[0], ::bounds[2], ::bounds[4]);
        bnds.max() = nvis::vec3(::bounds[1], ::bounds[3], ::bounds[5]);
    }
    else {
        bnds.min() = nvis::vec3(global_bounds[0], global_bounds[2], global_bounds[4]);
        bnds.max() = nvis::vec3(global_bounds[1], global_bounds[3], global_bounds[5]);
    }

    the_bounds = bnds;

    std::cout << "Resolution = " << res[0] << " x " << res[1] << " x " << res[2] << std::endl;
    xavier::raster_grid<3> sampling_grid(res, bnds);

    std::cout << "sampling grid bounds are: " << sampling_grid.bounds().min()
        << " -> " << sampling_grid.bounds().max() << '\n';
    int npoints = sampling_grid.size();

    float* flowmap = (float*)calloc(3*npoints, sizeof(float));
    float* flowtimes = (float*)calloc(npoints, sizeof(float));

    std::vector<VTK_SMART(vtkGenericCell)> thread_cell(nb_threads);
    for (int i=0; i<nb_threads; ++i) thread_cell[i] = vtkGenericCell::New();

    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';

    // initialize coordinates
#pragma openmp parallel
    for (int n=0 ; n<npoints ; ++n) {
        nvis::ivec3 c = sampling_grid.coordinates(n);
        nvis::vec3 x = sampling_grid(c);
        flowmap[3*n  ] = x[0];
        flowmap[3*n+1] = x[1];
        flowmap[3*n+2] = x[2];
    }

    size_t nbcomputed=0;
    xavier::ProgressDisplay progress(true);

    size_t nb_lost = 0;
    progress.start(npoints);

    std::vector< std::vector< std::vector< vec3 > > > thread_lines(nb_threads);
    std::vector< std::vector< std::vector< double > > > thread_times(nb_threads);
    std::vector<bool> stopped(npoints, false);

    std::vector< std::vector<vec3> > lines(npoints);
    std::vector< std::vector<double> > times(npoints);

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
            // auto stepper = make_controlled(eps, eps, odeint::runge_kutta_dopri5<point_type>());
            auto stepper = make_controlled(eps, eps, stepper_type());

            point_type x = {flowmap[3*n], flowmap[3*n+1], flowmap[3*n+2]};
            double t=0, d=0;
            point_type p(x);
            std::vector<vec3>& sline = lines[n];
            std::vector<double>& slinetimes = times[n];

            rhs_type<intp_t> rhs(intp, 1000, (T>0));

            observer_t obs(p, t, d, sline, slinetimes);
            try {
                point_type y(x);
                integrate_adaptive(stepper, rhs, y, static_cast<double>(0),
                                   fabs(T), 1.0e-4, obs);
                flowmap[3*n  ] = y[0];
                flowmap[3*n+1] = y[1];
                flowmap[3*n+2] = y[2];
                flowtimes[n] = obs.last_t;
                if (verbose)
                    std::cout << "integration successful: final position="
                        << to_str(y) << ", (last sampled location: " << to_str(obs.last_p) << ", distance: "
                            << obs.distance << ", final time: " << slinetimes.back() << ")\n";
            }
            catch (std::exception& e) {
                if (true || verbose) {
                    std::cerr << "caught exception while integrating from "
                        << to_str(x) << ":" << e.what() << '\n';
                    std::cout << "last position reached was "
                        << to_str(obs.last_p) << " at time "
                            << obs.last_t << " for a total distance of "
                                << obs.distance << '\n';
                }
                ++nb_lost;
                flowtimes[n] = obs.last_t;
                flowmap[3*n  ] = obs.last_p[0];
                flowmap[3*n+1] = obs.last_p[1];
                flowmap[3*n+2] = obs.last_p[2];
                stopped[n] = true;

            }
        }
    }

    progress.end();

    if (save_lines) {
        std::vector<nvis::ivec2> dummy;
        VTK_SMART(vtkPolyData) pdata = vtk_utils::make_polylines(lines,dummy, 0);
        std::vector<double> all_times;
        std::for_each(times.begin(), times.end(),
            [&](const std::vector<double>& ts) {
                all_times.insert(all_times.end(), ts.begin(), ts.end());
            });
        vtk_utils::add_scalars(pdata, all_times);
        {
            std::ostringstream os;
            os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "streamlines-deltaT=" << fabs(T) << ".vtp";
            vtk_utils::saveVTK(pdata, os.str());
        }
    }

    std::vector<size_t> size(4);
    std::vector<double> step(4), mins(4);
    std::vector<int> ctrs(4);
    std::fill(ctrs.begin(), ctrs.end(), nrrdCenterNode);
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
    xavier::nrrd_utils::writeNrrdFromContainers(flowmap, os.str(), size, step, mins, ctrs, cmdline);

    size[0] = 1;
    os.clear();
    os.str("");
    os << name_out << "-flowmap-deltaT=" << T << "-time.nrrd";
    xavier::nrrd_utils::writeNrrdFromContainers(flowtimes, os.str(), size, step, mins, ctrs, cmdline);

    delete[] flowmap;
    delete[] flowtimes;

    return 0;
}

template<typename DataSet>
int runTBB(VTK_SMART(DataSet) dataset) {
    typedef interpolator<DataSet> intp_t;
    typedef odeint::boundary_aware_stepper<vec3> stepper_type;

    if (dataset->GetPointData()->GetVectors() == nullptr) {
        std::cerr << "This dataset does not contain a vector field\n";
        std::cerr << "Exiting Now.\n";
        exit(1);
    }

    double global_bounds[6];
    dataset->GetBounds(global_bounds);
    nvis::bbox3 bnds;

    vec3 target(0.548387, 0.0266667, 0.114714);

    if (::bounds[0]<::bounds[1] && ::bounds[2]<::bounds[3] && ::bounds[4]<::bounds[5] &&
        ::bounds[0] >= global_bounds[0] && ::bounds[1] <= global_bounds[1] &&
        ::bounds[2] >= global_bounds[2] && ::bounds[3] <= global_bounds[3] &&
        ::bounds[4] >= global_bounds[4] && ::bounds[5] <= global_bounds[5]) {
        // valid bounds supplied by user
        bnds.min() = nvis::vec3(::bounds[0], ::bounds[2], ::bounds[4]);
        bnds.max() = nvis::vec3(::bounds[1], ::bounds[3], ::bounds[5]);
    }
    else {
        bnds.min() = nvis::vec3(global_bounds[0], global_bounds[2], global_bounds[4]);
        bnds.max() = nvis::vec3(global_bounds[1], global_bounds[3], global_bounds[5]);
    }

    the_bounds = bnds;

    std::cout << "Resolution = " << res[0] << " x " << res[1] << " x " << res[2] << std::endl;
    xavier::raster_grid<3> sampling_grid(res, bnds);

    std::cout << "sampling grid bounds are: " << sampling_grid.bounds().min()
        << " -> " << sampling_grid.bounds().max() << '\n';
    int npoints = sampling_grid.size();

    float* flowmap = (float*)calloc(3*npoints, sizeof(float));
    float* flowtimes = (float*)calloc(npoints, sizeof(float));

    std::vector<VTK_SMART(vtkGenericCell)> thread_cell(nb_threads);
    for (int i=0; i<nb_threads; ++i) thread_cell[i] = vtkGenericCell::New();

    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';

    // initialize coordinates

    tbb_progress_counter = 0;

    tbb::parallel_for(tbb::blocked_range<int>(0,npoints),
                       [&](tbb::blocked_range<int> r) {

        for (int n=r.begin(); n!=r.end(); ++n) {
        // for (int n=0 ; n<npoints ; ++n) {
            nvis::ivec3 c = sampling_grid.coordinates(n);
            nvis::vec3 x = sampling_grid(c);
            flowmap[3*n  ] = x[0];
            flowmap[3*n+1] = x[1];
            flowmap[3*n+2] = x[2];

            ++tbb_progress_counter;
        }
    });

    size_t nbcomputed=0;
    xavier::ProgressDisplay progress(true);
    progress.fraction_on();

    size_t nb_lost = 0;
    progress.start(npoints);

    std::vector<bool> stopped(npoints, false);
    std::vector< std::vector<vec3> > lines(npoints);
    std::vector< std::vector<double> > times(npoints);

    tbb_progress_counter = 0;
    tbb::parallel_for(tbb::blocked_range<int>(0,npoints),
                       [&](tbb::blocked_range<int> r) {
        intp_t intp(dataset);

        for (int n=r.begin(); n!=r.end(); ++n) {
            ++tbb_progress_counter;

            if (!(tbb_progress_counter%10)) progress.update(tbb_progress_counter);

            if (stopped[n]) {
                continue;
            }

            // create a stepper
            auto stepper = make_controlled(eps, eps, stepper_type());

            point_type x = {flowmap[3*n], flowmap[3*n+1], flowmap[3*n+2]};
            double t=0, d=0;
            point_type p(x);
            std::vector<vec3>& sline = lines[n];
            std::vector<double>& slinetimes = times[n];

            bool do_verbose = (x-target).norm() < 1.0e-6;
            rhs_type<intp_t> rhs(intp, 1000, (T>0), do_verbose);

            observer_t obs(p, t, d, sline, slinetimes, do_verbose);
            try {
                point_type y(x);
                integrate_adaptive(stepper, rhs, y, static_cast<double>(0),
                                   fabs(T), 1.0e-4, obs);
                flowmap[3*n  ] = y[0];
                flowmap[3*n+1] = y[1];
                flowmap[3*n+2] = y[2];
                flowtimes[n] = obs.last_t;
                if (verbose) {
                    std::ostringstream os;
                    os << n << ": successful: final position: "
                        << to_str(y) << ", (last sampled: " << to_str(obs.last_p) << ", distance: "
                            << obs.distance << ", final time: " << slinetimes.back() << ")\n";
                    std::cout << os.str();
                }
                if (flowtimes[n] >= 0.99*fabs(T))
                {
                    sline.clear();
                    slinetimes.clear();
                }
            }
            catch (std::exception& e) {
                if (verbose) {
                    std::ostringstream os;
                    os << "\n\ncaught exception while integrating from "
                        << to_str(x) << ":" << e.what() << '\n'
                        << "last position reached was "
                        << to_str(obs.last_p) << " at time "
                        << obs.last_t << " for a total distance of "
                        << obs.distance << std::endl;
                    std::cerr << os.str();
                }
                ++nb_lost;
                flowtimes[n] = obs.last_t;
                flowmap[3*n  ] = obs.last_p[0];
                flowmap[3*n+1] = obs.last_p[1];
                flowmap[3*n+2] = obs.last_p[2];
                stopped[n] = true;
            }
        }
    }); // end of tbb::parallel_for

    progress.end();

    if (save_lines || monitor) {
        if (monitor) {
            std::vector< std::vector<vec3> > newlines;
            std::vector< std::vector<double> > newtimes;
            for (size_t i=0; i<lines.size(); i++) {
                if (lines[i].size() > 2) {
                    newlines.push_back(lines[i]);
                    newtimes.push_back(times[i]);
                }
            }
            lines.swap(newlines);
            times.swap(newtimes);
        }
        std::vector<nvis::ivec2> dummy;
        VTK_SMART(vtkPolyData) pdata = vtk_utils::make_polylines(lines,dummy, 0);
        std::vector<double> all_times;
        std::for_each(times.begin(), times.end(),
            [&](const std::vector<double>& ts) {
                all_times.insert(all_times.end(), ts.begin(), ts.end());
            });
        vtk_utils::add_scalars(pdata, all_times);
        if (monitor) {
            VTK_CREATE(vtkCellArray, vertices);
            vertices->InitTraversal();
            size_t offset=0;
            for (int i=0; i<lines.size(); ++i) {
                vertices->InsertNextCell(1);
                vertices->InsertCellPoint(offset + lines[i].size() - 1);
                offset += lines[i].size();
            }
            pdata->SetVerts(vertices);
        }
        {
            std::ostringstream os;
            os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "streamlines-deltaT=" << fabs(T) << ".vtp";
            vtk_utils::saveVTK(pdata, os.str());
        }
    }

    std::vector<size_t> size(4);
    std::vector<double> step(4), mins(4);
    step[0] = std::numeric_limits<double>::quiet_NaN();
    mins[0] = std::numeric_limits<double>::quiet_NaN();
    std::vector<int> ctrs(4);
    std::fill(ctrs.begin(), ctrs.end(), nrrdCenterNode);
    size[0] = 3;
    for (int i = 0 ; i < 3 ; ++i) {
        size[i+1] = res[i];
        step[i+1] = sampling_grid.spacing()[i];
        mins[i+1] = sampling_grid.bounds().min()[i];
        mins[i+1] = sampling_grid.bounds().min()[i];
    }

    std::ostringstream os;
    os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "flowmap-deltaT=" << fabs(T) << ".nrrd";
    xavier::nrrd_utils::writeNrrdFromContainers(flowmap, os.str(), size, step, mins, ctrs, cmdline);

    size[0] = 1;
    os.clear();
    os.str("");
    os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "flowtime-deltaT=" << fabs(T) << ".nrrd";
    xavier::nrrd_utils::writeNrrdFromContainers(flowtimes, os.str(), size, step, mins, ctrs, cmdline);

    delete[] flowmap;
    delete[] flowtimes;

    return 0;
}

int main(int argc, const char* argv[])
{
    using namespace xavier;
    using namespace odeint;

    initialize(argc, argv);

#if _OPENMP
    nb_threads = omp_get_max_threads();
#endif
    std::cout << nb_threads << " threads available\n";

    VTK_SMART(vtkDataSet) dataset;
    dataset = vtk_utils::readVTK(name_in);

    vtkUnstructuredGrid* ugrid = vtkUnstructuredGrid::SafeDownCast(dataset);
    vtkStructuredGrid* sgrid = vtkStructuredGrid::SafeDownCast(dataset);
    vtkImageData* img = vtkImageData::SafeDownCast(dataset);
    vtkRectilinearGrid* rgrid = vtkRectilinearGrid::SafeDownCast(dataset);

    if (!in_parallel) {
        if (ugrid != nullptr) run<vtkUnstructuredGrid>(ugrid);
        else if (sgrid != nullptr) run<vtkStructuredGrid>(sgrid);
        else if (img != nullptr) run<vtkImageData>(img);
        else if (rgrid != nullptr) run<vtkRectilinearGrid>(rgrid);
        else {
            std::cerr << "Unsupported dataset type\n";
            return 1;
        }
    }
    else {
        if (ugrid != nullptr) runTBB<vtkUnstructuredGrid>(ugrid);
        else if (sgrid != nullptr) runTBB<vtkStructuredGrid>(sgrid);
        else if (img != nullptr) runTBB<vtkImageData>(img);
        else if (rgrid != nullptr) runTBB<vtkRectilinearGrid>(rgrid);
        else {
            std::cerr << "Unsupported dataset type\n";
            return 1;
        }
    }

    return 0;
}
