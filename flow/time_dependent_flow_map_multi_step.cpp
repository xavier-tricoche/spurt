#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <new>
#include <queue>
#include <sstream>
#include <thread>
#include <mutex>

#include <boost/filesystem.hpp>
#include <boost/numeric/odeint.hpp>

#include <math/types.hpp>
#include <data/image.hpp>
#include <flow/ode_observer.hpp>

#include <flow/time_dependent_field.hpp>
#include <flow/vector_field.hpp>
#include <format/dlr_reader.hpp>
#include <format/filename.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <misc/strings.hpp>
#include <misc/meta_utils.hpp>
#include <vtk/vtk_interpolator.hpp>

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>

using namespace spurt;

std::atomic<size_t> progress_counter;
std::string name_in, name_out, seed_name, path, cmdline;
double length, eps, eps_refined, T, t0, t1, deltaT;
size_t dim, mem;
svec3 res;
std::array<bool, 3> periodic;
std::array<double, 6> bounds;
bool verbose=false, do_log=true;
size_t max_rhs_evals = 10000;
size_t nb_threads;
bool save_lines = false;
bool monitor = false;
bool confine = false;
bool split_prisms = false;
bbox3 the_bounds;
std::ofstream ofs;
std::mutex output_mutex;
std::mutex progress_mutex;
double threshold=1.01;
bool use_bdry_aware = true;

int prism_rotation[6][6] = {
    {0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3},
    {2, 0, 1, 5, 3, 4}, {3, 5, 4, 0, 2, 1},
    {4, 3, 5, 1, 0, 2}, {5, 4, 3, 2, 1, 0}
};
int prism_cases[2][4] = {
    {1, 5, 2, 4}, {2, 4, 1, 5}
};
int prism_tets[2][3][4] = {
        { {0, 1, 2, 5}, {0, 1, 5, 4}, {0, 4, 5, 3} },
        { {0, 1, 2, 4}, {0, 4, 2, 5}, {0, 4, 5, 3} }
};

void write_to_ostream(std::ostream& os, const std::string& str) {
    {
        std::scoped_lock lock(output_mutex);
        os << str << '\n';
    }
}

void update_progress(spurt::ProgressDisplay& progress) {
    {
        std::scoped_lock lock(progress_mutex);
        progress.update(progress_counter);
    }
}

namespace odeint = boost::numeric::odeint;

void initialize(int argc, const char* argv[]) {
    namespace xcl = spurt::command_line;

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
    eps_refined = 1.0e-10;
    res = { 64, 64, 64 };
    periodic = { false, false, false };
    bounds = { 1, -1, 1, -1, 1, -1 }; // invalid bounds
    deltaT = -1; // valid deltaT (time discretization) is always >= 0
    t1=0;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input filename for list of file / time", required);
        parser.add_value("seeds", seed_name, "File containing seeds", optional);
        parser.add_value("output", name_out, "Output base name", required);
        parser.add_value("eps", eps, eps, "Integration precision", optional);
        parser.add_value("eps2", eps_refined, eps_refined, "Finest integration precision", optional);
        parser.add_value("T", T, "Integration time", required);
        parser.add_value("t0", t0, "Integration start time", optional);
        parser.add_value("t1", t1, t1, "Shortest integration time", optional);
        parser.add_value("dT", deltaT, deltaT, "Time between intermediate exported results", optional);
        parser.add_value("maxeval", max_rhs_evals, max_rhs_evals, "Maximum number of RHS evaluations", optional);
        parser.add_tuple<3>("res", res, res, "Sampling resolution", optional);
        parser.add_tuple<3>("periodic", periodic, periodic, "Periodic boundary conditions", optional);
        parser.add_tuple<6>("bounds", bounds, "Sampling bounds", optional);
        parser.add_value("confine", confine, confine, "Confine integration to prescribed bounds", optional);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional);
        parser.add_value("geometry", save_lines, save_lines, "Save streamlines geometry", optional);
        parser.add_value("monitor", monitor, monitor, "Save geometry of problematic streamlines", optional);
        parser.add_value("split", split_prisms, split_prisms, "Split prisms into three tetrahedra", optional);
        parser.add_value("log", do_log, do_log, "Log timings in a file", optional);
        parser.add_value("threshold", threshold, threshold, "Verbose activation threshold", optional);
        parser.add_value("bdry", use_bdry_aware, use_bdry_aware, "Use boundary-aware for Dana's rectilinear reconstruction", optional);

        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n"
            << e.what() << "\n"
                << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

using namespace vtk_utils;

typedef vec3 point_type;
typedef vec3 vector_type;
typedef double scalar_type;

typedef vtk_utils::point_locator<vtkUnstructuredGrid, double, 3> 
    unstructured_locator_type;
typedef spurt::fixed_mesh_time_dependent_field<unstructured_locator_type, std::vector<vec3> > unstructured_field_type;
typedef vtk_utils::point_locator<vtkStructuredGrid, double, 3> curvilinear_locator_type;
typedef spurt::fixed_mesh_time_dependent_field<curvilinear_locator_type, std::vector<vec3> > curvilinear_field_type;
typedef spurt::structured_mesh_time_dependent_field< vtkRectilinearGrid, double, 3, std::vector<vec3> > rectilinear_field_type;
typedef spurt::structured_mesh_time_dependent_field< vtkImageData, double, 3, std::vector<vec3> > image_field_type;
typedef spurt::tp_bspline_time_dependent_field< boundaryAwareRectGrid, double, 3, std::vector<vec3> > BARG_field_type;  //BARG edit

template<typename Field>
struct rhs_base {
    typedef Field field_type;

    rhs_base(std::shared_ptr<const field_type> field, size_t max_evals=1000,
             bool _verbose=false, int seedid=-1, size_t counter=0)
    : m_field(field), m_counter(counter), m_max_evals(max_evals),
      m_verbose(_verbose), m_seedid(seedid) {}

    rhs_base(const rhs_base& other)
        : m_field(other.m_field), m_counter(other.m_counter),
          m_max_evals(other.m_max_evals), m_verbose(other.m_verbose),
          m_seedid(other.m_seedid) {}

    std::shared_ptr<const field_type> m_field;
    mutable size_t m_counter;
    const size_t m_max_evals;
    bool m_verbose;
    int m_seedid;
};

template<typename Field, typename Enable=void>
struct rhs_type {};

template<typename Field>
struct rhs_type<
    Field,
    typename std::enable_if<
        std::is_same<typename Field::dataset_type, vtkStructuredGrid >::value ||
        std::is_same<typename Field::dataset_type, vtkRectilinearGrid>::value ||
        std::is_same<typename Field::dataset_type, vtkImageData>::value ||
		std::is_same<typename Field::dataset_type, boundaryAwareRectGrid>::value   //BARG edit #3
    >::type > : public rhs_base<Field> {
    typedef Field field_type;
    typedef rhs_base<Field> base_type;

    rhs_type(std::shared_ptr<const field_type> field, size_t max_evals=1000,
             bool _verbose=false, int seedid=-1)
        : base_type(field, max_evals, _verbose, seedid) {}

    rhs_type(const rhs_type& other)
        : base_type(other.m_field, other.m_max_evals,
                    other.m_verbose, other.m_seedid, other.m_counter) {}

    void operator()(const vec3& x, vec3& dxdt, scalar_type t) const {
        if (confine && !the_bounds.inside(x)) {
            if (verbose && do_log) {
                std::ostringstream os;
                os << "invalid position: " << x << " at t=" << t << '\n';
                write_to_ostream(ofs, os.str());
            }
            throw invalid_position_exception("invalid position: " + to_str(x) + " at t=" + to_string(t));
        }
        dxdt = (*base_type::m_field)(x, t); // may throw
        ++base_type::m_counter;
        if (base_type::m_counter >= base_type::m_max_evals) {
            std::ostringstream os;
            os << "Max RHS evaluations reached at " << x
               << " at t=" << t << '\n';
            write_to_ostream(ofs, os.str());
            // will trigger a invalid_position_exception:
            throw std::runtime_error(os.str());
        }
    }
};

template<typename Field>
struct rhs_type<
    Field,
    typename std::enable_if<
        std::is_same<typename Field::dataset_type, vtkUnstructuredGrid >::value
    >::type > : public rhs_base<Field> {
    typedef Field field_type;
    typedef rhs_base<Field> base_type;

    rhs_type(std::shared_ptr<const field_type> field, size_t max_evals=1000,
             bool _verbose=false, int seedid=-1, size_t counter=0)
        : base_type(field, max_evals, _verbose, seedid, counter) {
        VTK_INIT(vtkGenericCell, m_cell);
    }

    rhs_type(const rhs_type& other)
        : base_type(other.m_field, other.m_max_evals, other.m_verbose, other.m_seedid, other.m_counter) {
        VTK_INIT(vtkGenericCell, m_cell);
    }

    void operator()(const vec3& x, vec3& dxdt, scalar_type t) const {
        if (confine && !the_bounds.inside(x)) {
            if (verbose && do_log) {
                std::ostringstream os;
                os << "invalid position: " << x << " at t=" << t << '\n';
                write_to_ostream(ofs, os.str());
            }
            throw invalid_position_exception("invalid position: " + to_str(x) + " at t=" + to_string(t));
        }
        dxdt = (*base_type::m_field)(x, t, m_cell);
        ++base_type::m_counter;
        if (base_type::m_counter >= base_type::m_max_evals) {
            std::ostringstream os;
            os << "Max RHS evaluations reached at "
               << to_str(x) << " at t=" << t << '\n';
            write_to_ostream(ofs, os.str());
            // will trigger a invalid_position_exception:
            throw std::runtime_error(os.str());
        }
    }

    VTK_SMART(vtkGenericCell) m_cell;
};
typedef spurt::Observer<vec3> observer_t;

const vtkDataSet* get_dataset(const std::shared_ptr<unstructured_field_type> field) {
    return field->get_locator()->get_dataset();
}

const vtkDataSet* get_dataset(const std::shared_ptr<curvilinear_field_type> field) {
    return field->get_locator()->get_dataset();
}

const vtkDataSet* get_dataset(const std::shared_ptr<rectilinear_field_type> field) {
    return field->get_interpolator()->get_dataset();
}

const vtkDataSet* get_dataset(const std::shared_ptr<image_field_type> field) {
    return field->get_interpolator()->get_dataset();
}

const vtkDataSet* get_dataset(const std::shared_ptr<BARG_field_type> field) {
    return field->get_interpolator()->get_dataset();
}

template<typename Field>
int runTBB(shared_ptr<Field> field) {
    double global_bounds[6];
    const_cast<vtkDataSet*>(get_dataset(field))->GetBounds(global_bounds);
    bbox3 bnds;

    std::cout << "dataset imported.\n";

    // serial portion of the algorithm
    std::vector<vec3> seeds;

    if (::bounds[0]<::bounds[1] &&
        ::bounds[2]<::bounds[3] &&
        ::bounds[4]<::bounds[5] &&
        ::bounds[0] >= global_bounds[0] &&
        ::bounds[1] <= global_bounds[1] &&
        ::bounds[2] >= global_bounds[2] &&
        ::bounds[3] <= global_bounds[3] &&
        ::bounds[4] >= global_bounds[4] &&
        ::bounds[5] <= global_bounds[5]) {
        // valid bounds supplied by user
        bnds.min() = vec3(::bounds[0], ::bounds[2], ::bounds[4]);
        bnds.max() = vec3(::bounds[1], ::bounds[3], ::bounds[5]);
    }
    else {
        bnds.min() = vec3(global_bounds[0], global_bounds[2], global_bounds[4]);
        bnds.max() = vec3(global_bounds[1], global_bounds[3], global_bounds[5]);
    }

    the_bounds = bnds;
    int npoints;
    typedef spurt::raster_grid<long, double, 3> grid_t;
    typedef grid_t::coord_type coord_t;
    grid_t sampling_grid(res, bnds);

    if (seed_name.empty()) {
        std::cout << "initializing seed locations..." << std::flush;
        std::cout << "Resolution = " << res[0] << " x " << res[1] << " x " << res[2] << std::endl;

        std::cout << "sampling grid bounds are: " 
            << (sampling_grid.bounds().min())
            << " -> " 
            << (sampling_grid.bounds().max()) << '\n';

        if (do_log) {
            ofs << "Resolution = " << res[0] << " x " << res[1] << " x " << res[2] << std::endl;
            ofs << "sampling grid bounds are: " 
                << (sampling_grid.bounds().min())
                << " -> " 
                << (sampling_grid.bounds().max()) << '\n';
        }
        npoints = sampling_grid.size();
        seeds.resize(npoints);
        for (int i=0; i<npoints; i++) {
            seeds[i] = sampling_grid(sampling_grid.coordinates(i));
        }
        std::cout << " done\n";
    }
    else {
        std::cout << "importing seed locations from file..." << std::flush;
        std::string ext = spurt::filename::extension(seed_name);
        if (ext == "vtk" || ext == "VTK" || ext == "vtp" || ext == "VTP") {
            VTK_SMART(vtkDataSet) dataset = vtk_utils::readVTK(seed_name);
            npoints = dataset->GetNumberOfPoints();
            seeds.resize(npoints);
            for (int i=0; i<npoints; i++) {
                dataset->GetPoint(i, (double*)(&(seeds[i][0])));
            }
        }
        else if (ext == "txt" || ext == "TXT") {
            std::ifstream ifs(seed_name);
            ifs >> npoints;
            seeds.resize(npoints);
            for (int i=0; i<npoints; ++i) {
                ifs >> seeds[i][0] >> seeds[i][1] >> seeds[i][2];
            }
            ifs.close();
        }
        std::cout << " done\n";
    }

    float* flowmap = (float*)calloc(3*npoints, sizeof(float));
    float* flowtimes = (float*)calloc(npoints, sizeof(float));
    int lastpct = -1;
    std::cout << "nb points = " << npoints << '\n';
    if (do_log) {
        ofs << "nb points = " << npoints << '\n';
    }
    // initialize coordinates

    progress_counter = 0;

    // serial portion of the algorithm ends

    // parallel portion of the algorithm:
    // initialization of initial conditions
    std::cout << "initializing flow map with seeds..." << std::flush;
#ifndef NO_TBB
    tbb::parallel_for(tbb::blocked_range<int>(0,npoints),
                       [&](tbb::blocked_range<int> r) {
        for (int n=r.begin(); n!=r.end(); ++n) {
#else
        for (int n=0 ; n<npoints ; ++n) {
#endif
            flowmap[3*n  ] = seeds[n][0];
            flowmap[3*n+1] = seeds[n][1];
            flowmap[3*n+2] = seeds[n][2];

            ++progress_counter;
        }
#ifndef NO_TBB
    });
#endif
    std::cout << " done\n";
    using namespace boost::numeric::odeint;

    // parallel portion of the algorithm
    std::fill(flowtimes, flowtimes + npoints, t0);

    size_t nbcomputed=0;
    spurt::ProgressDisplay progress(true);
    progress.fraction_on();
    size_t nb_lost = 0;
    if (deltaT < 0) {
        deltaT = T;
    }
    else {
        deltaT *= sign(T);
    }

    double current_time = t0;
    if (t1==0) t1=deltaT;
    else {
        t1 *= sign(T);
    }

    int n_rounds=0;
    if (fabs(t1) >= fabs(T)) n_rounds = 1;
    else if (fabs(T)-fabs(t1) <= deltaT) n_rounds = 2;
    else n_rounds = 1 + std::ceil((fabs(T)-fabs(t1))/fabs(deltaT));

    std::cout << "n_rounds = " << n_rounds << '\n';
    if (do_log) {
        ofs << "n_rounds = " << n_rounds << '\n';
    }

    std::cout << "verbose = " << (verbose ? "true" : "false") << '\n';

    std::vector<bool> failed(npoints, false);
    std::vector<bool> stopped(npoints, false);
    std::vector< std::vector<vec3> > lines(npoints);
    std::vector< std::vector<double> > times(npoints);

    for (int round=0; round<n_rounds; ++round) {
        double next_time = current_time + (!round ? t1 : deltaT);
        if (T>0 && next_time>t0+T) {
            next_time = t0+T;
        }
        else if (T<0 && next_time<t0+T) {
            next_time = t0+T;
        }
        progress.start(npoints);
        std::cout << "Round #" << round << ": Flow map computation from " << current_time << " to " << next_time << '\n';
        if (do_log) {
            ofs << "Starting round #" << round << ": Flow map computation from " << current_time << " to " << next_time << '\n';
        }

        lines.resize(npoints);
        times.resize(npoints);

        std::cout << "lines initially contains " << lines.size() << " lines\n";

        progress_counter = 0;
        tbb::parallel_for(tbb::blocked_range<int>(0,npoints),
                           [&](tbb::blocked_range<int> r) {
            for (int n=r.begin(); n!=r.end(); ++n) {
                ++progress_counter;

                update_progress(progress);

                bool super_verbose = !seed_name.empty();
                point_type x = {flowmap[3*n], flowmap[3*n+1], flowmap[3*n+2]};
                if (stopped[n]) {
                    if (do_log && verbose) {
                        std::ostringstream oss;
                        oss << "pathline " << n << ": already stopped at "
                        << x;
                        write_to_ostream(ofs, oss.str());
                    }
                    continue;
                }

                std::vector<vec3>& sline = lines[n];
                std::vector<double>& slinetimes = times[n];

                bool do_verbose = verbose;

                rhs_type<Field> rhs(field, max_rhs_evals, do_verbose, n);
                double t;
                double d;
                point_type p;
                bool too_short = true;
                try {
                    for (int k=0; k<2; ++k) {
                        double local_eps = (!k) ? eps : eps_refined;
                        point_type y(x);
                        t=current_time; // obs.last_t
                        d=0; // obs.dist
                        p=x; // obs.last_p
                        sline.clear();
                        slinetimes.clear();
                        observer_t obs(p, t, d, sline, slinetimes /*, do_verbose*/);
                        auto stepper = make_controlled(local_eps, local_eps, runge_kutta_dopri5<point_type>());
                        integrate_adaptive(stepper, rhs, y, current_time, 
                                           next_time-current_time, 1.0e-4, obs);
                        if (fabs(t-next_time)/fabs(next_time-current_time) < 0.001) {
                            too_short = false;
                            break;
                        }
                        too_short = true;
                        if (do_log) {
                            std::ostringstream os;
                            os << n << ": unable to complete requested integration time. Advected from " << to_str(x) << " to " << to_str(p) << ", distance: "
                                    << d << ", final time: " << t << ")\n";
                            os << "refining step size: " << local_eps << " -> " << local_eps/2 << '\n';
                            write_to_ostream(ofs, os.str());
                        }
                        // will refine step size if this is first attempt
                    }
                    flowmap[3*n  ] = p[0];
                    flowmap[3*n+1] = p[1];
                    flowmap[3*n+2] = p[2];
                    flowtimes[n]   = t;
                    if (do_log && verbose && !too_short) {
                        std::ostringstream os;
                        os << "Seed #" << n << " advected from " << to_str(x) << " to "
                            << to_str(p) << ", distance: "
                                << d << ", final time: " << t << ")\n";
                        os << "line contains " << sline.size() << " (" << lines[n].size() << ") points" << '\n';
                        if (do_log) write_to_ostream(ofs, os.str());
		            }
                    if (too_short) {
                        if (do_log) {
                            std::ostringstream os;
                            os << "This pathline (id: " << n << ", seed: "
                            << seeds[n]
                            << ") stopped prematurely. Last recorded time: "
                            << t << " (" << fabs(next_time-flowtimes[n])
                            << ", "
                            << 100.*fabs((flowtimes[n]-next_time)/(next_time-current_time))
                            << "% missing) for a distance of " << d
                            << ". There are " << lines[n].size() << "points in this curve\n";
                            write_to_ostream(ofs, os.str());
                        }
                        failed[n] = true;
                        stopped[n] = true;
                    }
                }
                catch (std::exception& e) {
                    if (do_log) {
                        std::ostringstream os;
                        os << "\n\ncaught exception while integrating from "
                           << x << ". reason: " << e.what() << '\n'
                           << "last position reached was "
                           << p << " at time "
                           << t << " for a total distance of "
                           << d;
                        write_to_ostream(ofs, os.str());
                    }
                    flowtimes[n]   = t;
                    flowmap[3*n  ] = p[0];
                    flowmap[3*n+1] = p[1];
                    flowmap[3*n+2] = p[2];
                    stopped[n] = true;
                }
            }
        }); // end of tbb::parallel_for

        progress.end();
        if (do_log) {
            ofs << progress << '\n';
        }

        if ((save_lines || monitor) && verbose) {
            std::cout << "before filtering, lines contains " << lines.size()
                      << " and their lengths are ";

            for (int i=0; i<lines.size(); ++i)
                std::cout << lines[i].size() << " ";
            std::cout << "\n";
        }

        if (save_lines || monitor) {
            if (true || monitor) {
                std::vector< std::vector<vec3> > newlines;
                std::vector< std::vector<double> > newtimes;
                for (size_t i=0; i<lines.size(); i++) {
                    if (!save_lines && !stopped[i]) continue;
                    if (lines[i].size() > 2) {
                        newlines.push_back(lines[i]);
                        newtimes.push_back(times[i]);
                    }
                    else if (verbose) {
                        std::cout << "rejected line of length " << lines[i].size() << '\n';
                    }
                }
                lines.swap(newlines);
                times.swap(newtimes);
            }
            if (verbose) {
                std::cout << "after filtering lines contains "
                          << lines.size() << " lines.\n";
                if (do_log) {
                    std::ostringstream os;
                    os << "after filtering, lines contains " << lines.size() << " lines\n";
                    write_to_ostream(ofs, os.str());
                }
            }
            std::vector<ivec2> dummy;
            VTK_SMART(vtkPolyData) pdata = vtk_utils::make_polylines(lines, dummy, 0);
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
                  os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "pathlines_from_" << current_time << "_to_" << next_time << "_round_" << round << "_of_" << n_rounds << ".vtp";
                else {
                  os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "interrupted_pathlines_from_" << current_time << "_to_" << next_time << "_round_" << round << "_of_" << n_rounds << ".vtp";
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
            os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "flowmap-deltaT=" << fabs(next_time-t0) << "_t0=" << t0 << ".nrrd";
            spurt::nrrd_utils::writeNrrdFromContainers(flowmap, os.str(), size, step, mins, ctrs, cmdline);
            if (do_log) {
                ofs << "exported: " << os.str() << '\n';
            }

            size[0] = 1;
            os.clear();
            os.str("");
            os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "flowtime-deltaT=" << fabs(next_time-t0) << "_t0=" << t0 << ".nrrd";
            spurt::nrrd_utils::writeNrrdFromContainers(flowtimes, os.str(), size, step, mins, ctrs, cmdline);
            if (do_log) {
                ofs << "exported: " << os.str() << '\n';
            }
        }

        current_time = next_time;
    } // n_rounds

    if (!failed.empty()) {
        if (do_log) {
            std::ostringstream os;
            os << "\nFailed integrations started at following locations at time " << t0 << '\n';
            for (int i=0; i<failed.size(); ++i) {
                if (failed[i]) {
                    os << "id: " << i << ", pos: " << seeds[i] << '\n';
                }
            }
            write_to_ostream(ofs, os.str());
        }
        if (monitor) {
            std::vector<point_type> new_seeds;
            std::vector<int> seed_ids;
            for (int i=0; i<failed.size(); ++i) {
                if (failed[i]) {
                    new_seeds.push_back(seeds[i]);
                    seed_ids.push_back(i);
                }
            }
            if (do_log) {
                std::ostringstream os;
                os << "There were " << new_seeds.size() << " failed integrations overall\n";
                write_to_ostream(ofs, os.str());
            }
            lines.resize(new_seeds.size());
            times.resize(new_seeds.size());
            progress_counter = 0;
            tbb::parallel_for(tbb::blocked_range<int>(0,new_seeds.size()),
            [&](tbb::blocked_range<int> r) {
                for (int n=r.begin(); n!=r.end(); ++n) {
                    ++progress_counter;

                    double t=t0, d=0;
                    const point_type& x = new_seeds[n];
                    point_type p(x);

                    std::vector<vec3>& sline = lines[n];
                    std::vector<double>& slinetimes = times[n];
                    sline.clear();
                    slinetimes.clear();

                    bool do_verbose = verbose;

                    rhs_type<Field> rhs(field, max_rhs_evals, do_verbose, n);
                    observer_t obs(p, t, d, sline, slinetimes, do_verbose);

                    try {
                        point_type y(x);
                        // attempt integration over the full time interval, knowing that it will
                        // stop somewhere along the way

                        auto stepper = make_controlled(eps_refined, eps_refined, runge_kutta_dopri5<point_type>());
                        integrate_adaptive(stepper, rhs, y, t0,
                                           T, 1.0e-4, obs);
                        if (do_log) {
                            std::ostringstream os;
                            os << "failed streamline #" << seed_ids[n]
                               << " advected from " << to_str(x) << " to "
                               << to_str(obs.last_p) << ", distance: "
                               << obs.distance << ", final time: " << obs.last_t << ")\n";
                            write_to_ostream(ofs, os.str());
                        }
                    }
                    catch (std::exception& e) {
                        if (do_log) {
                            std::ostringstream os;
                            os << "\n\ncaught exception while integrating from "
                            << to_str(x) << " (failed streamline #"
                            << seed_ids[n]
                            << "). reason: " << e.what() << '\n'
                            << "last position reached was "
                            << to_str(obs.last_p) << " at time "
                            << obs.last_t << " for a total distance of "
                            << obs.distance;
                            write_to_ostream(ofs, os.str());
                        }
                    }
                } // for loop
#ifndef NO_TBB
            });
#endif

            VTK_SMART(vtkPolyData) failed_lines = vtk_utils::make_polylines(lines, 0, false);
            VTK_CREATE(vtkFloatArray, line_times);
            line_times->SetNumberOfComponents(1);
            line_times->SetNumberOfTuples(failed_lines->GetNumberOfPoints());
            vtkIdType _count=0;
            for (int n=0; n<times.size(); n++) {
                for (int i=0; i<times[n].size(); ++i, ++_count) {
                    line_times->InsertTuple(_count, &(times[n][i]));
                }
            }
            failed_lines->GetPointData()->SetScalars(line_times);
            std::ostringstream os;
            os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "failed_pathlines_from_" << t0 << "_to_" << t0+T << ".vtp";
            vtk_utils::saveVTK(failed_lines, os.str());
            // end of monitor section
        }
    }

    delete[] flowmap;
    delete[] flowtimes;

    return 0;
}

void get_filenames(std::vector<std::string>& steps, std::vector<double>& times) {
    std::ifstream info_file(name_in.c_str());
    if (!info_file) {
        std::cerr << "Unable to open input file " << name_in << '\n';
        if (do_log) {
            ofs << "Unable to open input file " << name_in << '\n';
            ofs.close();
        }
        exit(1);
    }
    std::string buffer;
    while (!info_file.eof() && info_file.good()) {
        std::getline(info_file, buffer);
        if (buffer[0] == '#') continue;
        else if (buffer.empty()) break;
        std::istringstream iss(buffer);
        std::string filename;
        double t;
        iss >> filename >> t;
        steps.push_back(filename);
        times.push_back(t);
    }
}

template<typename Dataset>
VTK_SMART(Dataset)
import_data(std::vector< std::shared_ptr<std::vector<vector_type> > >& data,
            std::vector<bool>& out_of_bounds,
            const std::vector<std::string>& steps) {
    data.resize(steps.size());
    VTK_SMART(Dataset) the_grid;
    for (int i=0; i<steps.size(); i++) {
        VTK_SMART(Dataset) a_grid = Dataset::SafeDownCast(vtk_utils::readVTK(steps[i]));
        if (!a_grid) throw std::runtime_error("Invalid grid type");
        if (!i) the_grid = a_grid;
        data[i] = std::shared_ptr< std::vector< vector_type > >(new std::vector< vector_type >());
        data[i]->resize(the_grid->GetNumberOfPoints());
        // std::cout << "grid contains " << the_grid->GetNumberOfPoints() << " points\n";
        // the_grid->PrintSelf(std::cout, vtkIndent(0));
        for (size_t j=0; j<data[i]->size(); j++) {
            double vec[3];
            the_grid->GetPointData()->GetVectors()->GetTuple(j, vec);
            (*data[i])[j] = vec3(vec[0], vec[1], vec[2]);
        }
        if (!i) {
            out_of_bounds.resize(the_grid->GetNumberOfCells());
            VTK_CREATE(vtkIdList, ptsids);
            for (int j=0; j<out_of_bounds.size(); ++j) {
                the_grid->GetCellPoints(j, ptsids);
                bool is_zero=true;
                for (int n=0; n<ptsids->GetNumberOfIds() && is_zero; ++n) {
                    const vec3& v = (*data[i])[ptsids->GetId(n)];
                    if (v[0]!=0 || v[1]!=0 || v[2]!=0) {
                        is_zero=false;
                    }
                }
                out_of_bounds[j] = is_zero;
            }
        }
        if (verbose) std::cout << steps[i] << " imported\n";
        size_t ninvalid;
        if (!i) ninvalid = std::count(out_of_bounds.begin(), out_of_bounds.end(), true);
        if (verbose && !i) std::cout << ninvalid << " cells out of bounds\n";
        if (do_log) {
            ofs << steps[i] << " imported\n";
            if (!i) ofs << ninvalid << " cells out of bounds\n";
        }
    }
    return the_grid;
}

shared_ptr<rectilinear_field_type>
load_rectilinear_VTK_time_steps(const std::vector<std::string>& steps,
                                const std::vector<double>& times) {
    assert(!steps.empty());
    std::vector< std::shared_ptr<std::vector<vector_type> > > data;
    std::vector<bool> out_of_bounds;
    VTK_SMART(vtkRectilinearGrid) grid =
        import_data<vtkRectilinearGrid>(data, out_of_bounds, steps);

    return shared_ptr<rectilinear_field_type>(
        new rectilinear_field_type(grid, data, times, out_of_bounds, verbose));
}

shared_ptr<image_field_type>
load_image_VTK_time_steps(const std::vector<std::string>& steps,
                          const std::vector<double>& times) {
    assert(!steps.empty());
    std::vector< std::shared_ptr<std::vector<vector_type> > > data;
    std::vector<bool> out_of_bounds;
    VTK_SMART(vtkImageData) grid =
        import_data<vtkImageData>(data, out_of_bounds, steps);

    return shared_ptr<image_field_type>(
        new image_field_type(grid, data, times, out_of_bounds/*, verbose*/));
}

shared_ptr<curvilinear_field_type>
load_curvilinear_VTK_time_steps(const std::vector<std::string>& steps,
                         const std::vector<double>& times) {
    assert(!steps.empty());
    std::vector< std::shared_ptr<std::vector<vector_type> > > data;
    std::vector<bool> out_of_bounds;
    VTK_SMART(vtkStructuredGrid) grid =
        import_data<vtkStructuredGrid>(data, out_of_bounds, steps);

    shared_ptr<curvilinear_locator_type>
        locator(new curvilinear_locator_type(grid, verbose));
    return shared_ptr<curvilinear_field_type>(
        new curvilinear_field_type(locator, data, times, verbose));
}

shared_ptr<unstructured_field_type>
load_unstructured_VTK_time_steps(const std::vector<std::string>& steps,
                                 const std::vector<double>& times) {
    assert(!steps.empty());
    std::vector< std::shared_ptr<std::vector<vector_type> > > data;
    std::vector<bool> out_of_bounds;
    VTK_SMART(vtkUnstructuredGrid) grid =
        import_data<vtkUnstructuredGrid>(data, out_of_bounds, steps);

    shared_ptr<unstructured_locator_type>
        locator(new unstructured_locator_type(grid, verbose));
    return shared_ptr<unstructured_field_type>(
        new unstructured_field_type(locator, data, times, verbose));
}

shared_ptr<BARG_field_type>
load_BARG_time_steps(const std::vector<std::string>& steps,
                                const std::vector<double>& times) {
    assert(!steps.empty());
    std::vector< std::shared_ptr<std::vector<vector_type> > > data;
    std::vector<bool> out_of_bounds;

	VTK_SMART(vtkRectilinearGrid) rgrid =
		import_data<vtkRectilinearGrid>(data, out_of_bounds, steps);

	VTK_CREATE(boundaryAwareRectGrid, grid) ;
	grid->ShallowCopy(rgrid);

    return shared_ptr<BARG_field_type>(
        new BARG_field_type(grid, data, times, out_of_bounds, false/*, verbose*/, use_bdry_aware));
}

shared_ptr<unstructured_field_type> load_DLR_time_steps() {
    std::ifstream info_file(name_in.c_str());
    if (!info_file) {
        std::cerr << "Unable to open input file " << name_in << '\n';
        if (do_log) {
            ofs << "Unable to open input file " << name_in << '\n';
            ofs.close();
        }
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

    // select subset of files that is needed for requested computation
    {
        std::vector<double> sorted_times(times.begin(), times.end());
        std::sort(sorted_times.begin(), sorted_times.end());
        double min = t0;
        double max = t0+T;
        if (T<0) std::swap(min, max);
        // first element that is not less than min
        auto lower = std::lower_bound(sorted_times.begin(), sorted_times.end(), min);
        // first element that is greater than max
        auto upper = std::upper_bound(sorted_times.begin(), sorted_times.end(), max);
        if (lower == sorted_times.end()) {
            // all time steps lie before min time
            std::cerr << "ERROR (1): interval [" << min << ", " << max << "] out of bounds\n";
            exit(1);
        }
        else if (lower == sorted_times.begin()) {
            if (*lower > min) {
                std::cerr << "ERROR (2): min time " << min << " out of bounds\n";
                exit(1);
            }
        }
        else {
            --lower;
        }
        min = *lower;
        if (upper == sorted_times.end()) {
            // all time steps lie before max time
            std::cout << "ERROR (3): max time " << max << " out of bounds\n";
            exit(1);
        }
        max = *upper;
        std::cout << "min=" << min << ", max=" << max << '\n';
        std::vector<double> new_times;
        std::vector<std::string> new_steps;
        for (int i=0; i<times.size(); ++i) {
            double t = times[i];
            if (t >= min && t <= max) {
                std::cout << "added time t=" << t << '\n';
                new_times.push_back(times[i]);
                new_steps.push_back(steps[i]);
            }
            else {
                std::cout << "rejecting time t=" << t << '\n';
            }
        }
        times.swap(new_times);
        steps.swap(new_steps);
    }

    spurt::dlr_reader reader(mesh_name, "");
    std::vector<fvec3> vertices;
    std::vector<long> cell_indices;
    std::vector<std::pair<dlr_reader::cell_type, long> > cell_types;
    reader.read_mesh(false, vertices, cell_indices, cell_types);

    if (split_prisms) {
        std::vector<long> new_cell_indices;
        std::vector<std::pair<dlr_reader::cell_type, long> > new_cell_types;
        for (long n=0; n<cell_types.size()-1; ++n) {
            auto cell_type = cell_types[n];
            long first = cell_type.second;
            long next = cell_types[n+1].second;
            if (cell_type.first != dlr_reader::PRISM) {
                new_cell_types.push_back(cell_type);
                new_cell_types.back().second = new_cell_indices.size();
                for (long k=first; k<next; ++k) {
                    new_cell_indices.push_back(cell_indices[k]);
                }
            }
            else {
                /*
                int prism_rotation[6][6] = {
                    {0, 1, 2, 3, 4, 5}, {1, 2, 0, 4, 5, 3},
                    {2, 0, 1, 5, 3, 4}, {3, 5, 4, 0, 2, 1},
                    {4, 3, 5, 1, 0, 2}, {5, 4, 3, 2, 1, 0} };
                int prism_cases[2][4] = {
                    {1, 5, 2, 4}, {2, 4, 1, 5}
                };
                int prism_tets[2][3][4] = {
                        { {0, 1, 2, 5}, {0, 1, 5, 4}, {0, 4, 5, 3} },
                        { {0, 1, 2, 4}, {0, 4, 2, 5}, {0, 4, 5, 3} }
                };
                */
                std::vector<long> prism_ids(6);
                for (long k=first; k<next; ++k) {
                    prism_ids[k-first] = cell_indices[k];
                }
                long minid=std::distance(prism_ids.begin(), std::min_element(prism_ids.begin(), prism_ids.end()));
                std::vector<int> rotated_ids(6);
                for (int i=0; i<6; ++i) {
                    rotated_ids[i] = prism_ids[prism_rotation[minid][i]];
                }
                if (std::min(rotated_ids[prism_cases[0][0]], rotated_ids[prism_cases[0][1]]) <
                    std::min(rotated_ids[prism_cases[0][2]], rotated_ids[prism_cases[0][3]])) {
                    for (int t=0; t<3; ++t) {
                        new_cell_types.push_back(std::make_pair(dlr_reader::TETRAHEDRON, new_cell_indices.size()));
                        for (int id=0; id<4; ++id) {
                            new_cell_indices.push_back(rotated_ids[prism_tets[0][t][id]]);
                        }
                    }
                }
                else {
                    // tets case 1
                    for (int t=0; t<3; ++t) {
                        new_cell_types.push_back(std::make_pair(dlr_reader::TETRAHEDRON, new_cell_indices.size()));
                        for (int id=0; id<4; ++id) {
                            new_cell_indices.push_back(rotated_ids[prism_tets[1][t][id]]);
                        }
                    }
                }
            }
        }
        std::cout << "after conversion, there are "<< cell_types.size()-1 << " tetrahedra in the mesh\n";
        cell_types.swap(new_cell_types);
        cell_indices.swap(new_cell_indices);
    }

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
        // reorder vertices if cell is a prism
        if (cell_types[cell_id].first == dlr_reader::PRISM) {
            cells->InsertCellPoint(cell_indices[start]);
            cells->InsertCellPoint(cell_indices[start+2]);
            cells->InsertCellPoint(cell_indices[start+1]);
            cells->InsertCellPoint(cell_indices[start+3]);
            cells->InsertCellPoint(cell_indices[start+5]);
            cells->InsertCellPoint(cell_indices[start+4]);
        }
        else {
            for (long i=0 ; i<size ; ++i) {
                cells->InsertCellPoint(cell_indices[start+i]);
            }
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
            case dlr_reader::TRIANGLE:
                type_name = VTK_TRIANGLE;
                break;
            case dlr_reader::QUADRILATERAL:
                type_name = VTK_QUAD;
                break;
            case dlr_reader::TETRAHEDRON:
                type_name = VTK_TETRA;
                break;
            case dlr_reader::HEXAHEDRON:
                type_name = VTK_HEXAHEDRON;
                break;
            case dlr_reader::PRISM:
                type_name = VTK_WEDGE;
                break;
            case dlr_reader::PYRAMID:
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
    std::shared_ptr<unstructured_locator_type>
        locator(new unstructured_locator_type(grid, false, true));

    std::vector< std::shared_ptr<std::vector<vector_type> > > data(steps.size());
    for (int i=0; i<steps.size(); i++) {
        const std::string ext = spurt::filename::extension(steps[i]);
        data[i] = std::shared_ptr< std::vector< vector_type > >(new std::vector< vector_type >());
        if (ext == "nrrd" || ext == "nhdr") {
            Nrrd* nin = spurt::nrrd_utils::readNrrd(steps[i]);
            data[i]->resize(nin->axis[1].size);
            // std::cout << "nin->axis[1].size=" << nin->axis[1].size << '\n';
            // std::cout << "data[i]->size()=" << data[i]->size() << '\n';
            spurt::nrrd_utils::to_vector<vector_type,float>(*data[i], nin->data, data[i]->size());
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
        if (do_log) {
            ofs << steps[i] << " imported\n";
        }
    }

    return shared_ptr<unstructured_field_type>(new unstructured_field_type(locator, data, times, verbose));
}

int main(int argc, const char* argv[])
{
    using namespace spurt;
    using namespace odeint;

    initialize(argc, argv);

    if (do_log) {
        ofs.open(name_out + ".log", std::ios::app);
    }

#if _OPENMP
    nb_threads = omp_get_max_threads();
#else
    nb_threads = std::thread::hardware_concurrency();
#endif
    std::cout << nb_threads << " threads available\n";
    const std::string ext = spurt::filename::extension(name_in);

    if (ext == "dlr") {
        std::shared_ptr<unstructured_field_type> field = load_DLR_time_steps();
        runTBB(field);
    }
    else if (ext == "tdvtk") {
        std::vector<std::string> steps;
        std::vector<double> times;
        get_filenames(steps, times);
        assert(!steps.empty());
        const std::string& name = steps[0];
        std::string ext = spurt::filename::extension(name);
        spurt::lower_case(ext);
        if (ext == "vtk") {
            // we need to load the dataset to figure out the type...
            VTK_SMART(vtkDataSet) mystery = vtk_utils::readVTK(name);
            if (vtkImageData::SafeDownCast(mystery)) {
                std::shared_ptr<image_field_type> field = load_image_VTK_time_steps(steps, times);
                runTBB(field);
            }
            else if (vtkRectilinearGrid::SafeDownCast(mystery)) {
                std::shared_ptr<rectilinear_field_type> field = load_rectilinear_VTK_time_steps(steps, times);
                runTBB(field);
            }
            else if (vtkStructuredGrid::SafeDownCast(mystery)) {
                std::shared_ptr<curvilinear_field_type> field = load_curvilinear_VTK_time_steps(steps, times);
                runTBB(field);
            }
            else if (vtkUnstructuredGrid::SafeDownCast(mystery)) {
                std::shared_ptr<unstructured_field_type> field = load_unstructured_VTK_time_steps(steps, times);
                runTBB(field);
            }
            else {
                std::cerr << "Unsupported VTK dataset type\n";
                exit(1);
            }
        }
        else if (ext == "vti") {
            std::shared_ptr<image_field_type> field = load_image_VTK_time_steps(steps, times);
            runTBB(field);
        }
        else if (ext == "vtu") {
            std::shared_ptr<unstructured_field_type> field = load_unstructured_VTK_time_steps(steps, times);
            runTBB(field);
        }
        else if (ext == "vtr") {
            std::shared_ptr<rectilinear_field_type> field = load_rectilinear_VTK_time_steps(steps, times);
            runTBB(field);
        }
        else {
            std::cerr << "Unsupported VTK dataset type\n";
            exit(1);
        }
    }
	else if (ext == "atdvtk")
	{
		//std::shared_ptr<BARG_field_type> field = load_VTK_time_steps_BARG();
        std::vector<std::string> steps;
        std::vector<double> times;
        get_filenames(steps, times);
        assert(!steps.empty());
        const std::string& name = steps[0];
        std::string ext = spurt::filename::extension(name);
        spurt::lower_case(ext);
		assert(ext == "vtr");
		std::shared_ptr<BARG_field_type> field = load_BARG_time_steps(steps, times);
        runTBB(field);
	}

    else {
        std::cerr << "Unrecognized input file extension: " << ext << '\n';
        exit(1);
    }

    if (do_log) ofs.close();

    return 0;
}
