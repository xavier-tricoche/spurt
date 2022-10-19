#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <new>
#include <queue>
#include <sstream>
#include <thread>

#include <boost/filesystem.hpp>
#include <boost/numeric/odeint.hpp>

#include <data/raster.hpp>
#include <flow/time_dependent_field.hpp>
#include <flow/vector_field.hpp>
#include <flow/vector_field.hpp>
#include <format/DLRreader.hpp>
#include <math/bounding_box.hpp>
#include <math/dopri5.hpp>
#include <math/fixed_vector.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <misc/strings.hpp>
#include <VTK/vtk_interpolator.hpp>

// #define NO_TBB

#ifndef NO_TBB
#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
tbb::atomic<size_t> progress_counter;
#else
size_t progress_counter;
#endif


std::string name_in, name_out, seed_name, path, cmdline;
double length, eps, T, t0, deltaT;
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

void initialize(int argc, const char* argv[]) {
    namespace xcl = xavier::command_line;

    cmdline = "Command line: " + std::string(argv[0]);
    for (int i=1; i<argc; i++) {
        cmdline += " " + std::string(argv[i]);
    }

    xcl::option_traits
        required(true, false, "Required Options"),
    optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
    "Compute flow map across a time-dependent vector field");

    mem = 1024;
    eps = 1.0e-6;
    res = { 64, 64, 64 };
    periodic = { false, false, false };
    bounds = { 1, -1, 1, -1, 1, -1 }; // invalid bounds
    verbose = false;
    deltaT = -1;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input filename for list of file / time", required);
        parser.add_value("seeds", seed_name, "File containing seeds", optional);
        parser.add_value("output", name_out, "Output base name", required);
        parser.add_value("eps", eps, eps, "Integration precision", optional);
        parser.add_value("T", T, "Integration time", required);
        parser.add_value("t0", t0, "Integration start time", optional);
        parser.add_value("dT", deltaT, deltaT, "Time between intermediate exported results", optional);
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

struct vec3 : public Eigen::Vector3d
{
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
};

vec3& as_vec3(nvis::vec3& v) {
    return *((vec3*)(&v[0]));
}

const vec3& as_vec3(const nvis::vec3& v) {
    return *((const vec3*)(&v[0]));
}

inline std::string to_str(const vec3& p) {
    std::ostringstream os;
    os << "[" << p(0) << ", " << p(1) << ", " << p(2) << "]";
    return os.str();
}

namespace xavier {
template<>
struct data_traits< vec3 > {
    typedef vec3::value_type value_type;
    typedef vec3 data_type;

    constexpr static int size() { return 3; }
    static const value_type& value(const data_type& v, size_t i) {
        return v[i];
    }
    static value_type& value(data_type& v, size_t i) { return v[i]; }
    static data_type& assign(data_type& inout, value_type val) {
        inout.setConstant(val);
        return inout;
    }
    static value_type norm(const data_type& v) {
        return v.norm();
    }
};
} // xavier

typedef vec3 point_type;
typedef vec3 vector_type;
typedef double scalar_type;

typedef interpolator<vtkImageData, double, 3, vector_type>  img_intp_t;
typedef interpolator<vtkRectilinearGrid, double, 3, vector_type> rect_intp_t;
typedef interpolator<vtkUnstructuredGrid, double, 3, vector_type> unst_intp_t;
typedef interpolator<vtkStructuredGrid, double, 3, vector_type> curv_intp_t;

typedef vtk_utils::point_locator<vtkUnstructuredGrid, double, 3, vec3> locator_type;
typedef xavier::fixed_mesh_time_dependent_field<locator_type, std::vector<vec3> > field_type;

struct RKDOPRI5 {
    typedef nvis::vec3 vec_t; // vector type required by nvis::dopri5
    typedef nvis::dopri5<vec_t> ode_solver_type;
    typedef ode_solver_type::step step_type;

    template<typename RHS>
    struct dopri_rhs {
        dopri_rhs(const RHS& rhs) : m_rhs(rhs) {}

        vec_t operator()(double t, const vec_t& y) const {
            vector_type dydt;
            try {
                m_rhs(as_vec3(y), dydt, t);
            }
            catch(...) {
                throw nvis::invalid_position_exception("Unable to interpolate");
            }
            return dydt.as_nvis_vec3();
        }
        const RHS& m_rhs;
    };

    template<typename RHS, typename Obs>
    static void integrate(const RHS& rhs, point_type& y, double t0,
                          double T, double hinit, double hmax, double eps,
                          Obs& obs) {
        ode_solver_type solver;
        dopri_rhs<RHS> my_rhs(rhs);
        if (eps > 0) {
            solver.reltol = solver.abstol = eps;
        }
        if (hmax > 0) {
            solver.h_max = hmax;
        }
        solver.t = t0;
        solver.t_max = t0+T;
        solver.y = y.as_nvis_vec3();
        step_type step;
        try {
            while (true) {
                auto result = solver.do_step(my_rhs, step);
                if (result == ode_solver_type::OK || result ==
                    ode_solver_type::T_MAX_REACHED) {
                    obs(as_vec3(step.y1()), step.t1());
                    if (result == ode_solver_type::T_MAX_REACHED) break;
                }
                else {
                    std::ostringstream os;
                    os << "unable to complete integration: ";
                    if (result == ode_solver_type::STIFFNESS_DETECTED) {
                        os << "stiffness detected";
                    }
                    else if (result == ode_solver_type::STEPSIZE_UNDERFLOW) {
                        os << "step size underflow";
                    }
                    os << " at (" << to_str(as_vec3(solver.y)) << ", "
                       << solver.t << ")\n";
                    break;
                }
            }
        }
        catch(std::runtime_error& e) {
            std::ostringstream os;
            os << "RKDOPRI5::integrate: exception caught: " << e.what() << '\n';
            if (verbose) std::cerr << os.str() << std::flush;
            throw nvis::invalid_position_exception(os.str());
        }
    }
};

template<typename Field, typename Enable=void>
struct rhs_type {};

template<typename Field>
struct rhs_type<
    Field,
    typename std::enable_if<
        std::is_same<typename Field::dataset_type, vtkStructuredGrid >::value ||
        std::is_same<typename Field::dataset_type, vtkUnstructuredGrid>::value
    >::type > {
    typedef Field field_type;

    rhs_type(std::shared_ptr<const field_type> field, size_t max_evals=1000,
             bool _verbose=false)
    : m_field(field), m_counter(0), m_max_evals(max_evals), m_verbose(_verbose),
      m_cell(VTK_SMART(vtkGenericCell)::New()) {}

    rhs_type(std::shared_ptr<const field_type> field, vtkGenericCell* acell, size_t max_evals=1000, bool forward=true, bool _verbose=false)
        : m_field(field), m_cell(acell), m_counter(0), m_max_evals(max_evals),
          m_verbose(_verbose) {}

    rhs_type(const rhs_type& other)
        : m_field(other.m_field), m_counter(other.m_counter),
          m_cell(other.m_cell), m_max_evals(other.m_max_evals),
          m_verbose(other.m_verbose) {}

    void operator()(const vec3& x, vec3& dxdt, scalar_type t) const {
        if (confine && !the_bounds.inside(x.as_nvis_vec3())) {
            throw nvis::invalid_position_exception("invalid position: " + to_str(x) + " at t=" + to_string(t));
        }
        std::ostringstream os;
        if (false) {
            std::cout << "About to run interpolator\n";
            os << "rhs(" << to_str(x) << ", " << t << ")=";
        }
        dxdt = (*m_field)(x, t); // may throw
        if (false) {
            os << to_str(dxdt) << '\n';
            std::cout << os.str() << std::flush;
        }
        ++m_counter;
    }

    std::shared_ptr<const field_type> m_field;
    mutable size_t m_counter;
    const size_t m_max_evals;
    VTK_SMART(vtkGenericCell) m_cell;
    bool m_verbose;
};

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


int runTBB(shared_ptr<field_type> field) {
    double global_bounds[6];
    field->get_locator()->get_dataset()->GetBounds(global_bounds);
    nvis::bbox3 bnds;

    std::vector<vec3> seeds;

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
    int npoints;
    xavier::raster_grid<3> sampling_grid(res, bnds);

    if (seed_name.empty()) {
        std::cout << "Resolution = " << res[0] << " x " << res[1] << " x " << res[2] << std::endl;

        std::cout << "sampling grid bounds are: " << sampling_grid.bounds().min()
            << " -> " << sampling_grid.bounds().max() << '\n';
        npoints = sampling_grid.size();
        seeds.resize(npoints);
        for (int i=0; i<npoints; i++) {
            seeds[i] = as_vec3(sampling_grid(sampling_grid.coordinates(i)));
        }
    }
    else {
        VTK_SMART(vtkDataSet) dataset = vtk_utils::readVTK(seed_name);
        npoints = dataset->GetNumberOfPoints();
        seeds.resize(npoints);
        for (int i=0; i<npoints; i++) {
            dataset->GetPoint(i, (double*)(&(seeds[i][0])));
        }
    }

    float* flowmap = (float*)calloc(3*npoints, sizeof(float));
    float* flowtimes = (float*)calloc(npoints, sizeof(float));
    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';

    // initialize coordinates

    progress_counter = 0;

#ifndef NO_TBB
    tbb::parallel_for(tbb::blocked_range<int>(0,npoints),
                       [&](tbb::blocked_range<int> r) {
        for (int n=r.begin(); n!=r.end(); ++n) {
#else
        for (int n=0 ; n<npoints ; ++n) {
#endif
            // nvis::ivec3 c = sampling_grid.coordinates(n);
            // nvis::vec3 x = sampling_grid(c);
            flowmap[3*n  ] = seeds[n][0];
            flowmap[3*n+1] = seeds[n][1];
            flowmap[3*n+2] = seeds[n][2];

            ++progress_counter;
        }
#ifndef NO_TBB
    });
#endif

    size_t nbcomputed=0;
    xavier::ProgressDisplay progress(true);
    progress.fraction_on();

    size_t nb_lost = 0;
    progress.start(npoints);

    std::vector<bool> stopped(npoints, false);
    std::vector< std::vector<vec3> > lines(npoints);
    std::vector< std::vector<double> > times(npoints);

    progress_counter = 0;
#ifndef NO_TBB
    tbb::parallel_for(tbb::blocked_range<int>(0,npoints),
                       [&](tbb::blocked_range<int> r) {
        for (int n=r.begin(); n!=r.end(); ++n) {
#else
        for (int n=0; n<npoints; ++n) {
#endif
            ++progress_counter;

            if (!(progress_counter%10)) progress.update(progress_counter);

            if (stopped[n]) {
                continue;
            }
            point_type x = {flowmap[3*n], flowmap[3*n+1], flowmap[3*n+2]};
            double t=t0, d=0;
            point_type p(x);
            std::vector<vec3>& sline = lines[n];
            std::vector<double>& slinetimes = times[n];

            bool do_verbose = false;

            rhs_type<field_type> rhs(field, 1000, do_verbose);
            observer_t obs(p, t, d, sline, slinetimes, do_verbose);
            try {
                point_type y(x);
                RKDOPRI5::integrate(rhs, y, t0, T, 0, 0, eps, obs);
                flowmap[3*n  ] = obs.last_p[0];
                flowmap[3*n+1] = obs.last_p[1];
                flowmap[3*n+2] = obs.last_p[2];
                flowtimes[n] = obs.last_t;
                if (verbose && flowtimes[n] < t0 + 0.99*T) {
                    std::ostringstream os;
                    os << n << ": successful: final position: "
                        << to_str(y) << ", (last sampled: " << to_str(obs.last_p) << ", distance: "
                            << obs.distance << ", final time: " << slinetimes.back() << ")\n";
                    std::cout << os.str();
                }
                if (!save_lines && (!monitor || (fabs(flowtimes[n]-(t0+T)) < 1.0e-3))) {
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
                if (!monitor) {
                    sline.clear();
                    slinetimes.clear();
                }
            }
        }
#ifndef NO_TBB
    }); // end of tbb::parallel_for
#endif

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
            if (save_lines)
              os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "pathlines-deltaT=" << fabs(T) << "-t0=" << t0 << ".vtp";
            else {
              os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "interrupted_pathlines-deltaT=" << fabs(T) << "-t0=" << t0 << ".vtp";
            }
            vtk_utils::saveVTK(pdata, os.str());
        }
    }

    if (seed_name.empty()) {
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
        os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "flowmap-deltaT=" << fabs(T) << "_t0=" << t0 << ".nrrd";
        xavier::nrrd_utils::writeNrrdFromContainers(flowmap, os.str(), size, step, mins, ctrs, cmdline);

        size[0] = 1;
        os.clear();
        os.str("");
        os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "flowtime-deltaT=" << fabs(T) << "_t0=" << t0 << ".nrrd";
        xavier::nrrd_utils::writeNrrdFromContainers(flowtimes, os.str(), size, step, mins, ctrs, cmdline);
    }

    delete[] flowmap;
    delete[] flowtimes;

    return 0;
}

shared_ptr<field_type> load_DLR_time_steps() {
    std::ifstream info_file(name_in.c_str());
    if (!info_file) {
        std::cerr << "Unable to open input file " << name_in << '\n';
        exit(1);
    }

    std::string mesh_name;
    std::vector<std::string> steps;
    std::vector<double> times;
    std::string buffer;
    while (!info_file.eof() && info_file.good()) {
        std::getline(info_file, buffer);
        if (buffer[0] == '#') continue;
        else if (buffer.empty()) break;
        std::istringstream iss(buffer);
        if (mesh_name.empty())
            iss >> mesh_name;
        else {
            std::string name;
            double t;
            iss >> name >> t;
            steps.push_back(name);
            times.push_back(t);
        }
    }
    info_file.close();

    xavier::DLRreader reader(mesh_name, "");
    std::vector<nvis::fvec3> vertices;
    std::vector<long> cell_indices;
    std::vector<std::pair<DLRreader::cell_type, long> > cell_types;
    reader.read_mesh(false, vertices, cell_indices, cell_types);
    size_t ncells = cell_types.size()-1; // last entry is not an actual cell
    VTK_CREATE(vtkUnstructuredGrid, grid);
    VTK_SMART(vtkPoints) points = vtk_utils::make_vtkpoints(vertices);
    grid->SetPoints(points);
    VTK_CREATE(vtkCellArray, cells);
    cells->SetNumberOfCells(ncells);
    for (long cell_id=0 ; cell_id<ncells ; ++cell_id) {
        size_t start = cell_types[cell_id].second;
        size_t end = cell_types[cell_id+1].second;
        size_t size = end-start;
        cells->InsertNextCell(size);
        for (long i=0 ; i<size ; ++i) {
            cells->InsertCellPoint(cell_indices[start+i]);
        }
    }
    VTK_CREATE(vtkUnsignedCharArray, types);
    VTK_CREATE(vtkIdTypeArray, locations);
    types->SetNumberOfComponents(1);
    types->SetNumberOfTuples(cell_types.size());
    locations->SetNumberOfComponents(1);
    locations->SetNumberOfTuples(ncells);
    for (long cell_id=0 ; cell_id<ncells ; ++cell_id) {
        unsigned char type_name;
        switch(cell_types[cell_id].first) {
            case DLRreader::TRIANGLE:
                type_name = VTK_TRIANGLE;
                break;
            case DLRreader::QUADRILATERAL:
                type_name = VTK_QUAD;
                break;
            case DLRreader::TETRAHEDRON:
                type_name = VTK_TETRA;
                break;
            case DLRreader::HEXAHEDRON:
                type_name = VTK_HEXAHEDRON;
                break;
            case DLRreader::PRISM:
                type_name = VTK_WEDGE;
                break;
            case DLRreader::PYRAMID:
                type_name = VTK_PYRAMID;
                break;
            default:
                std::cerr << "cell #" << cell_id << " has type " << cell_types[cell_id].first << "\n";
                throw std::runtime_error("invalid cell type");
        }
        types->SetValue(cell_id, type_name);
        locations->SetValue(cell_id, cell_types[cell_id].second);
    }
    grid->SetCells(types, locations, cells);
    std::shared_ptr<locator_type> locator(new locator_type(grid, false, true));

    std::vector< std::shared_ptr<std::vector<vector_type> > > data(steps.size());
    for (int i=0; i<steps.size(); i++) {
        const std::string ext = xavier::filename::extension(steps[i]);
        data[i] = std::shared_ptr< std::vector< vector_type > >(new std::vector< vector_type >());
        if (ext == "nrrd" || ext == "nhdr") {
            Nrrd* nin = xavier::nrrd_utils::readNrrd(steps[i]);
            data[i]->resize(nin->axis[1].size);
            // std::cout << "nin->axis[1].size=" << nin->axis[1].size << '\n';
            // std::cout << "data[i]->size()=" << data[i]->size() << '\n';
            xavier::nrrd_utils::to_vector<vector_type,float>(*data[i], nin->data, data[i]->size());
        }
        else {
            std::vector<double> vx, vy, vz;
            reader.read_data_from_file(steps[i], "x_velocity", vx);
            reader.read_data_from_file(steps[i], "y_velocity", vy);
            reader.read_data_from_file(steps[i], "z_velocity", vz);
            data[i]->resize(vx.size());
            for (size_t j=0; j<data[i]->size(); j++) {
                (*data[i])[j] = vec3(vx[j], vy[j], vz[j]);
            }
        }
        if (verbose) std::cout << steps[i] << " imported\n";
    }

    return shared_ptr<field_type>(new field_type(locator, data, times, verbose));
}

int main(int argc, const char* argv[])
{
    using namespace xavier;
    using namespace odeint;

    initialize(argc, argv);

#if _OPENMP
    nb_threads = omp_get_max_threads();
#else
    nb_threads = std::thread::hardware_concurrency();
#endif
    std::cout << nb_threads << " threads available\n";

    std::shared_ptr<field_type> field = load_DLR_time_steps();
    runTBB(field);

    return 0;
}
