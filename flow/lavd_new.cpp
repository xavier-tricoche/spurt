#include <flow/lavd_new.hpp>

#include <sstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <cctype>
#include <stdexcept>
#include <locale>
#include <iomanip>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <boost/numeric/odeint.hpp>
#include <boost/filesystem.hpp>

#include <data/field_wrapper.hpp>
#include <data/raster.hpp>
#include <format/filename.hpp>
#include <image/nrrd_wrapper.hpp>
#include <image/probe.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <misc/log_helper.hpp>
#include <VTK/vtk_utils.hpp>

#include <Eigen/Core>
#include <Eigen/SVD>

// #include <vtkPoints.h>
// #include <vtkCellArray.h>
// #include <vtkPolyLine.h>
// #include <vtkPolyData.h>
// #include <vtkDataSetWriter.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace odeint = boost::numeric::odeint;

using namespace xavier::lavd;


std::string name_in, name_out;
std::vector<std::string> t_init_str(1, "0");
std::string t_max_str = "2d";
std::string t_expstp_str = "3h";
std::string dt_str="30";
std::string t_skip_str="0";
std::string me;
std::string log_name = "log.txt";
std::string start_time_str;
std::string border_mask_name;
std::ofstream log_file;

xavier::log::dual_ostream xavier::lavd::_log_(log_file, std::cout, 1, 0, 0, true);

value_t t_max=2*xavier::lavd::DAY;
value_t t_init=0;
value_t t_skip=0;
value_t t_export_step=3*xavier::lavd::HOUR;
value_t t_between_files=3*xavier::lavd::HOUR;
value_t eps=1.0e-8;
value_t dt=30;
std::array<size_t, 2> res( { 512, 512 } );
std::array< value_t, 4 > bnds( { 0, 0, 0, 0 } );
bbox_t region; // invalid bounds by default
bbox_t domain;
value_t step_x, step_y;
size_t nslices = 8*10; // 10 days of data processed simultaneously
size_t last_time_step;
size_t nb_samples;
size_t nb_lost=0;
size_t support_radius;
size_t nb_threads;
std::vector< int > verbose_vec;
int verbose_cout=0;
int verbose_log=1;
bool export_region = false;
bool long_name = false;
bool export_trajectories = false;
value_t next_export_time, current_start_time;
size_t max_rhs_evals = 1000000;
vec2 input_spc;

typedef const void* address_t;
const Nrrd* watched_nrrd;
double* stored_address;

inline vec3 append(const vec2& v, value_t t) {
    return vec3(v[0], v[1], t);
}

void set_verbose_values() {
    if (verbose_vec.size()<1 || verbose_vec.size()>2) {
        throw std::runtime_error("Invalid verbose values");
    }
    verbose_cout = verbose_vec[0];
    if (verbose_vec.size()==1) {
        verbose_log = verbose_cout;
    }
    else verbose_log = verbose_vec[1];
    _log_(0) << "verbose_cout=" << verbose_cout << ", verbose_log="
        << verbose_log << std::endl;
}

void initialize(int argc, const char* argv[])
{
    namespace xcl = xavier::command_line;

    xcl::option_traits
            required_group(true, false, "Required Options"),
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compute Lagrangian averaged vorticity deviation");

    verbose_vec.resize(2);
    verbose_vec[0] = verbose_cout;
    verbose_vec[1] = verbose_log;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input info filename", required_group);
        parser.add_value("output", name_out, "Output file basename ", required_group);
        parser.add_sequence("tinit", t_init_str, t_init_str, "Integration starting time, expressed either in time since initial time step (in seconds; for minutes append \"m\", for hours \"h\", for days \"d\", for weeks \"w\") or as a date MM/DD", required_group, "<float>[<char>]");
        parser.add_value("tmax", t_max_str, "Integration length (in seconds; for minutes append \"m\", for hours \"h\", for days \"d\", for weeks \"w\")", required_group, "<float>[<char>]");
        parser.add_value("step", t_expstp_str, t_expstp_str, "Time interval between intermediate result output", optional_group);
        parser.add_value("skip", t_skip_str, t_skip_str, "Time skipped before first output", optional_group);
        parser.add_value("dt", dt_str, dt_str, "Time sampling of integral curve for LAVD computation (in seconds, for minutes add \"m\", for hours \"h\", for times \"x\")", optional_group, "<float>[<char>]");
        parser.add_tuple<2>("res", res, res, "Sampling resolution", optional_group);
        parser.add_tuple<4>("bounds", bnds, bnds, "Sampling bounds (min longitude and latitude followed by max's)", optional_group);
        parser.add_value("nslices", nslices, nslices, "Number of time steps to process simultaneously", optional_group);
        parser.add_value("eps", eps, eps, "Integration precision", optional_group);
        parser.add_value("mask", border_mask_name, border_mask_name, "Mask filename", optional_group);
        parser.add_value("max", max_rhs_evals, max_rhs_evals, "Maximum authorized number of rhs evaluations during intermediate steps of the integration");
        parser.add_value("export", export_region, export_region, "Export index range of selected region", optional_group);
        parser.add_value("traj", export_trajectories, export_trajectories, "Export trajectories", optional_group);
        parser.add_value("log", log_name, log_name, "Log export filename", optional_group);
        parser.add_sequence("verbose", verbose_vec, verbose_vec, "Verbosity level(s) on stdout and log file", optional_group, "<uint> [<uint>]");
        parser.add_value("lname", long_name, long_name, "Export result with long name indicating computation parameters", optional_group);

        parser.parse(argc, argv);
        validate_value_with_time_unit(t_max, t_max_str, false);
        validate_value_with_time_unit(dt, dt_str, true);
        validate_value_with_time_unit(t_export_step, t_expstp_str, false);
        std::string t_init_single_str;
        for (auto it=t_init_str.begin(); it!=t_init_str.end(); ++it) {
            t_init_single_str.append(*it);
        }
        try {
            validate_value_with_time_unit(t_init, t_init_single_str);
        }
        catch(...) {
            validate_date_string(t_init, t_init_single_str);
        }
        validate_value_with_time_unit(t_skip, t_skip_str, false);
        set_verbose_values();
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

// global static data storage for parallel computation
Nrrd* current_velocity_volume=0;
Nrrd* current_vorticity_volume=0;
Nrrd* border_mask=0;
value_t current_t_min; // starting time
value_t current_t_max; // stopping time

inline void nuke_current_volumes() {
    if (current_velocity_volume != NULL) {
        nrrdNuke(current_velocity_volume);
    }
    if (current_vorticity_volume != NULL) {
        nrrdNuke(current_vorticity_volume);
    }
}

void clean_exit(int n) {
    if (log_file.is_open()) log_file.close();
    nuke_current_volumes();
    exit(n);
}

struct observer
{
    orbit_type& m_orbit;
    NrrdScalarField<3>& m_vort;

    value_t vorticity(const vec2& p, value_t t)
    {
        return m_vort( vec3(p[0], p[1], t) );
    }

    integration_state state(const vec2& p, value_t t)
    {
        return integration_state(p, t, vorticity(p, t));
    }

    observer(orbit_type& orbit, NrrdScalarField<3>& vorticity)
        : m_orbit(orbit), m_vort(vorticity) {}

    void initialize(const vec2& p0, value_t t0)
    {
        m_orbit.clear();
        m_orbit.push_back(state(p0, t0));
    }

    void operator()(const vec2 &x, value_t t)
    {
        m_orbit.push_back(state(x, t));
    }
};

std::vector< orbit_type > all_orbits;
std::vector< lavd_state > all_states;

void import_data(const std::vector< std::string >& vel_filenames,
                 const std::vector< std::string >& vor_filenames,
                 size_t start_id)
{
    size_t n_available = std::min(nslices, vel_filenames.size()-start_id+support_radius);
    size_t n_used = std::min(n_available, static_cast<size_t>((t_max-t_init)/t_between_files+2*support_radius));
    _log_(1)
        << "import_data: start_id=" << start_id
        << ", #filenames=" << vel_filenames.size()
        << ", n_available=" << n_available
        << ", (tmax-tinit)/t_between_files=" << (t_max-t_init)/t_between_files
        << ", number of time steps needed to complete entire integration length="
        << static_cast<size_t>((t_max-t_init)/t_between_files+2*support_radius)
        << ", n_used=" << n_used << std::endl;

    nuke_current_volumes();

    // initialize velocity and vorticity volumes
    std::vector< Nrrd* > vel_tsteps(n_used);
    std::vector< Nrrd* > vor_tsteps(n_used);
    size_t first = ( start_id >= support_radius ? start_id-support_radius : static_cast<size_t>(0) );
    for (size_t i=first; i<first+n_used; ++i) {
        _log_(1) << "Reading " << vel_filenames[i] << "... " << std::flush;
        vel_tsteps[i-first]=xavier::nrrd_utils::readNrrd(vel_filenames[i]);
        _log_(1) << "done" << std::endl;
        print(vel_tsteps[i-first]);

        _log_(1) << "Reading " << vor_filenames[i] << "... " << std::flush;
        vor_tsteps[i-first]=xavier::nrrd_utils::readNrrd(vor_filenames[i]);
        _log_(1) << "done" << std::endl;
    }
    last_time_step = first + n_used - 1;
    value_t outer_tmin = /*t_init +*/ first*t_between_files;
    _log_(2) << "outer_tmin=" << outer_tmin << std::endl;
    current_velocity_volume = xavier::lavd::create_nrrd_volume(vel_tsteps, outer_tmin, t_between_files);
    _log_(2) << "after velocity volume creation:" << std::endl;
    print(current_velocity_volume);

    // clean up
    for (size_t i=0; i<vel_tsteps.size(); ++i) {
        nrrdNuke(vel_tsteps[i]);
    }

    current_vorticity_volume = xavier::lavd::create_nrrd_volume(vor_tsteps, outer_tmin, t_between_files);
    current_t_min = start_id*t_between_files;
    current_t_max = last_time_step*t_between_files;
    if (current_t_max < t_max) current_t_max -= support_radius*t_between_files;
    else current_t_max = t_max;
    _log_(1)
        << "n_used=" << n_used << '\n'
        << "support_radius=" << support_radius << '\n'
        << "current_t_max=" << current_t_max << '\n';

    _log_(1)
        << "selected region = " << region << '\n'
        << "x spacing: " << input_spc[0] << '\n'
        << "y spacing: " << input_spc[1] << std::endl;


    // clean up
    for (size_t i=0; i<vor_tsteps.size(); ++i) {
        nrrdNuke(vor_tsteps[i]);
    }
}

void export_results(double current_time, double wall_time, double cpu_time, bool final) {

    std::ostringstream os;
    ios::fmtflags default_settings = os.flags();
    if (long_name) {
        os << name_out << "_started_at_" << start_time_str << "_flowmap+lavd_" << std::setw(5) << std::setfill('0') << current_time/xavier::lavd::HOUR << "h_";
        os << res[0] << "x" << res[1] << "_" << std::setprecision(2)
           << std::scientific << eps;
        os.flags(default_settings);
    }
    else {
        os << name_out << "_fmap+lavd_" << std::setw(5) << std::setfill('0') << current_time/xavier::lavd::HOUR << "h";
    }
    std::string basename = os.str();

    os.clear();
    os << "This file was produced by " << me
       << ". 1st axis corresponds to (flowmap_x, flowmap_y, time, lavd), where "
       << "lavd stands for Lagrangian averaged vorticity deviation. ";
    os << "Computation parameters: ";
    os << "*input file=" << name_in << "; ";
    os << "*tmax=" << current_time << "; ";
    os << "*epsilon=" << eps << "; ";
    os << "*resolution=" << res[0] << "x" << res[1] << "; ";
    os << "*bounds=" << region.min() << "->" << region.max() << "; ";
    os << "*dt=" << dt << "; ";
    os << "*kernel size=" << support_radius << '\n';
    std::string comment = os.str();

    xavier::lavd::export_results< float >(
        current_time,
        wall_time, cpu_time,
        final,
        nb_lost,
        res[0], res[1],
        step_x, step_y,
        export_trajectories,
        region,
        basename,
        std::vector<std::string>(1, comment),
        all_states,
        all_orbits
    );
}

int main(int argc, const char* argv[])
{
    using namespace xavier;
    using namespace odeint;

    me=argv[0];

    nb_threads=1;

#if _OPENMP
    nb_threads = omp_get_max_threads();
    _log_.set_nthreads(nb_threads);
    _log_(1) << nb_threads << " threads available" << std::endl;
#else
    _log_.set_nthreads(1);
#endif

    initialize(argc, argv);

    log_file.open(log_name.c_str());
    _log_.set_thresholds(verbose_cout, verbose_log);

    if (t_max < 0) {
        _log_(0) << "Invalid time duration symbol in definition of tmax" << std::endl;
        clean_exit(1);
    }

    if (dt < 0) {
        dt = (t_max - t_init)/(fabs(dt)-1);
    }

    nb_samples = res[0]*res[1];

    name_out=xavier::filename::remove_extension(name_out);

    std::vector<std::string> velocity_filenames;
    std::vector<std::string> vorticity_filenames;

    std::fstream info_file(name_in, std::ios::in);
    if (!info_file) {
        _log_(0) << "Unable to open info file named " << name_in << std::endl;
        exit(1);
    }

    boost::filesystem::path p(name_in);
    std::string parent_dir=p.parent_path().string();
    while (info_file.good()) {
        std::string velname, vortname;
        value_t avg;
        info_file >> velname >> vortname >> avg;
        velocity_filenames.push_back(parent_dir + '/' + velname);
        vorticity_filenames.push_back(parent_dir + '/' + vortname);
    }
    info_file.close();
    assert(!velocity_filenames.empty());

    _log_(1) << velocity_filenames.size() << " time steps in input" << std::endl;

    // compute bounds of entire domain
    get_spatial_info(domain, input_spc, velocity_filenames[0], 1);
    border_mask = get_border_mask(border_mask_name, velocity_filenames[0]);

    if ( !( bnds[0]==bnds[1]==bnds[2]==bnds[3] ) ) {
        region.min()=vec2( bnds[0], bnds[1] );
        region.max()=vec2( bnds[2], bnds[3] );
    }
    else {
        region.min() = domain.min();
        region.max() = domain.max();
    }

    _log_(1) << "region=\n" << region << std::endl;

    // storage for LAVD states and trajectories
    all_orbits.resize(nb_samples);
    all_states.resize(nb_samples);

    support_radius = xavier::lavd::compute_support_radius();
    _log_(1) << "support radius = " << support_radius << std::endl;

    // initialize velocity and vorticity volumes
    import_data(velocity_filenames, vorticity_filenames, t_init/HOUR/3);

    size_t counter = 0;
    size_t nb_early = 0;

    step_x = region.size()[0] / static_cast<value_t>(res[0]-1);
    step_y = region.size()[1] / static_cast<value_t>(res[1]-1);

    _log_(1) << "bounds=" << region << std::endl;

    if (export_region)
        xavier::lavd::export_mask(vorticity_filenames[0], region, true);

    std::vector<size_t> sample_counter(nb_threads);
    std::vector< shared_ptr< NrrdODERHS > > rhs_copies(nb_threads);
    std::vector< shared_ptr< NrrdVectorField > > vf_copies(nb_threads);
    for (int i=0; i<nb_threads; ++i) {
        vf_copies[i] = shared_ptr< NrrdVectorField >(new NrrdVectorField(current_velocity_volume, std::string("velocity volume for thread #") + std::to_string(i), false));
        rhs_copies[i] = shared_ptr< NrrdODERHS >(new NrrdODERHS(*vf_copies[i], sample_counter[i], domain, max_rhs_evals));
    }

    start_time_str=sec2time(t_init);
    std::for_each(start_time_str.begin(), start_time_str.end(),
                 [](char& c) { if (c==' ') c='_'; });
    _log_(1) << "start_time_string=" << start_time_str << std::endl;

    xavier::ProgressDisplay progress(false), total_progress(false);

    total_progress.start(1); // we only care for the timer function
    total_progress.set_active(false);

    std::clock_t clock_begin=std::clock();
    auto timer_start = std::chrono::high_resolution_clock::now();

    bool initial_loop = true;
    size_t niter = 0;
    current_start_time = t_init;
    next_export_time = t_init + t_export_step + t_skip;
    while ( true ) {
        _log_(0, "\n\n") << "next target time for this loop: "
            << sec2time(next_export_time)
            << " (" << next_export_time << " s.)" << std::endl;

        progress.start(nb_samples, "flow map");
        progress.set_active(true);
        counter = 0;
        nb_lost = 0;
        size_t nb_early = 0;
        value_t current_target_time = std::min(current_t_max, next_export_time);

        std::clock_t loop_clock_begin=std::clock();
        auto loop_timer_start = std::chrono::high_resolution_clock::now();
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for (size_t n = 0 ; n < nb_samples ; ++n) {

                #if _OPENMP
                const int thread=omp_get_thread_num();
                #else
                const int thread=0;
                #endif

                if (!thread) progress.update(n);

                if (all_states[n].stopped) {
                    #pragma omp atomic
                    ++nb_lost;
                    continue;
                }

                vec2 x0;
                if (initial_loop) {
                    size_t i = n % res[0];
                    size_t j = n / res[0];
                    x0 = region.min() + vec2(i*step_x, j*step_y);
                    _log_(2) << "\n\n\n\nstarting integration #" << n << " at " << x0 << std::endl;
                }
                else {
                    x0 = all_states[n].position;
                }

                #pragma omp atomic
                ++counter;

                // create a stepper
                auto stepper = make_controlled(eps, eps, runge_kutta_dopri5<vec2>());

                NrrdODERHS& rhs=*rhs_copies[thread];
                rhs.m_counter=0;

                NrrdScalarField<3> my_vorticityf(current_vorticity_volume);
                observer an_observer(all_orbits[n], my_vorticityf);

                if (initial_loop) {
                    _log_(2) << "initializing observer at " << x0 << std::endl;
                    try {
                        an_observer.initialize(x0, current_t_min);
                    }
                    catch( std::exception& e ) {
                        _log_(1) << "caught exception while attempting "
                            << "to initialize observer at " << x0 << std::endl;
                        all_states[n].stopped = true;
                        all_states[n].position = x0;
                        #pragma omp atommic
                        ++nb_lost;
                        continue;
                    }
                }

                try {
                    integrate_const(stepper, rhs, x0, current_start_time,
                        current_target_time, dt, an_observer);
                }
                catch(std::exception& e) {
                    _log_(1) << "caught exception while integrating from "
                            << x0 << ":" << e.what() << '\n'
                            << "integration started at " << x0
                            << " and ended at " << all_orbits[n].back().position()
                            << '\n' << std::flush;
                    all_states[n].stopped = true;
                    #pragma omp atomic
                    ++nb_lost;
                }
                catch(...) {
                    _log_(0) << "unknown exception thrown!" << std::endl;
                }
                _log_(2) << rhs.m_counter << " RHS evaluations for this IC"
                    << std::endl;
            }

            initial_loop = false;
        }
	    progress.end();

        if (current_target_time == next_export_time) {
            ++niter;
            _log_(0)
                << "Compute time this iteration: cpu: " << human_readable_duration(progress.cpu_time())
                    << " (" << progress.cpu_time() << " ms.) | wall: "
                << human_readable_duration(progress.wall_time())
                << " (" << progress.wall_time() << " ms.)\n"
                << "Overall compute time so far (" << niter << " iterations): cpu: ";
            double tcpu = total_progress.instant_cpu_time();
            double twal = total_progress.instant_wall_time();
             _log_(0)
                << human_readable_duration(tcpu) << " (" << tcpu << " ms.) | wall: "
                << human_readable_duration(twal) << " (" << twal << " ms.)"
                << std::endl;
        }

        if (current_target_time == t_max) break;
        else if (current_target_time == current_t_max) {
            // we have reached the end of the data currently available in
            // RAM. Import next time section.
            // [last_time_step-2*R, last_time_step-2*R+nslices-1]
            import_data(velocity_filenames, vorticity_filenames,
                        last_time_step-2*support_radius);
            for (int i=0; i<nb_threads; ++i) {
                vf_copies[i] = shared_ptr< NrrdVectorField >(new NrrdVectorField(current_velocity_volume, std::string("velocity volume for thread #") + std::to_string(i), false));
                rhs_copies[i] = shared_ptr< NrrdODERHS >(new NrrdODERHS(*vf_copies[i], sample_counter[i], domain, max_rhs_evals));
            }
        }
        if (current_target_time == next_export_time) {
            std::clock_t loop_clock_end=std::clock();
            auto loop_timer_end = std::chrono::high_resolution_clock::now();
            export_results(next_export_time,
                           std::chrono::duration<double, std::milli>(loop_timer_end-loop_timer_start).count(),
                           1000.*(loop_clock_end-loop_clock_begin)/CLOCKS_PER_SEC, false);
            next_export_time += t_export_step;
            if (next_export_time > t_max) next_export_time = t_max;
            // reset orbits
            for (size_t n=0; n<all_orbits.size(); ++n) {
                integration_state save = all_orbits[n].back();
                all_orbits[n].clear();
                all_orbits[n].push_back(save);
            }
        }
        // advance starting time
        current_start_time = current_target_time;
        _log_(1) << "current start time=" << current_start_time << std::endl;
    }
    std::clock_t clock_end=std::clock();
    auto timer_end = std::chrono::high_resolution_clock::now();
    export_results(t_max, std::chrono::duration<double, std::milli>(timer_end-timer_start).count(),
                   1000.*(clock_end-clock_begin)/CLOCKS_PER_SEC, true);

    return 0;
}
