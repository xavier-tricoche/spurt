#include <flow/lavd.hpp>

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

#ifdef _OPENMP
#include <omp.h>
#endif

namespace odeint = boost::numeric::odeint;

using namespace xavier::lavd;

typedef nvis::fixed_vector<value_t, 6> vec6;

std::string name_in, name_mask, name_out;
std::vector<std::string> t_init_str(1, "0");
std::string t_max_str = "2d";
std::string t_expstp_str = "3h";
std::string dt_str="30";
std::string t_skip_str="0";
std::string me;
std::string log_name;
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
value_t seed_ratio=0.01;
value_t min_dist=0.001;
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
bool resume = false;
value_t next_export_time, current_start_time;
size_t max_rhs_evals = 1000000;
vec2 input_spc;

typedef const void* address_t;

std::vector< std::vector< vec4 > > failed_paths;

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
        parser.add_value("output", name_out, "Output file basename", required_group);
        parser.add_value("mask", border_mask_name, border_mask_name, "Mask filename", optional_group);
        parser.add_value("ratio", seed_ratio, seed_ratio, "Seeding ratio", optional_group);
        parser.add_value("mind", min_dist, min_dist, "Min step size along trajectories", optional_group);
        parser.add_flag("resume", resume, "Use input file as starting point", optional_group);
        parser.add_sequence("tinit", t_init_str, t_init_str, "Integration starting time, expressed either in time since initial time step (in seconds; for minutes append \"m\", for hours \"h\", for days \"d\", for weeks \"w\") or as a date MM-DD", optional_group, "<float>[<char>]");
        parser.add_value("tmax", t_max_str, "Integration length (in seconds; for minutes append \"m\", for hours \"h\", for days \"d\", for weeks \"w\")", optional_group, "<float>[<char>]");
        parser.add_value("step", t_expstp_str, t_expstp_str, "Time interval between intermediate result output", optional_group);
        parser.add_value("dt", dt_str, dt_str, "Time sampling of integral curve for LAVD computation (in seconds, for minutes add \"m\", for hours \"h\", for times \"x\")", optional_group, "<float>[<char>]");
        parser.add_tuple<2>("res", res, res, "Sampling resolution", optional_group);
        parser.add_tuple<4>("bounds", bnds, bnds, "Sampling bounds (min longitude and latitude followed by max's)", optional_group);
        parser.add_value("nslices", nslices, nslices, "Number of time steps to process simultaneously", optional_group);
        parser.add_value("eps", eps, eps, "Integration precision", optional_group);
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
        // validate_value_with_time_unit(t_skip, t_skip_str, false);
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
Nrrd* current_average_vorticity_vector=0;
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
    if (current_average_vorticity_vector != NULL) {
        nrrdNuke(current_average_vorticity_vector);
    }
}

void clean_exit(int n) {
    if (log_file.is_open()) log_file.close();
    nuke_current_volumes();
    exit(n);
}

size_t seed_id;
struct observer
{
    LAVD_state&           m_state;
    std::vector< vec4 >&  m_steps;
    NrrdScalarField<1>&   m_avg_vort;
    NrrdScalarField<3>&   m_vort;

    value_t vorticity_deviation(const vec2& x, value_t t) {
        try {
            value_t w = m_vort( vec3(x[0], x[1], t ) );
            value_t aw = m_avg_vort( vec1(t) );
            return fabs(w-aw);
        }
        catch (std::exception& e) {
            _log_(1) << "exception caught in vorticity_deviation:\n"
                    << e.what() << std::endl;
            throw;
        }
    }
    
    observer(LAVD_state& state, std::vector<vec4>& steps, 
             NrrdScalarField<3>& vorticity, NrrdScalarField<1>& average_vorticity) 
        : m_state(state), m_steps(steps), m_vort(vorticity), 
		  m_avg_vort(average_vorticity) {}
    
    void initialize(const vec2& p0, value_t t0) {
        value_t vd0 = vorticity_deviation(p0, t0);
        m_state = LAVD_state(p0, vd0, t0);
        m_steps.clear();
        m_steps.push_back(vec4(p0[0], p0[1], t0, 0.));
    }

    void operator()(const vec2 &x, value_t t)
    {
        value_t vd = vorticity_deviation(x, t);
        m_state.update(x, vd, t);
		
        // append latest position
        m_steps.push_back(vec4(x[0], x[1], t, m_state.m_acc_vd));
        
        /*_log_(1) << "LAVD(" << t << ")=" << m_state.evaluate() 
            << " (vd=" << vd << ")" << std::endl;*/
    }
};

typedef std::vector<vec4> trajectory_t;
std::vector< trajectory_t > all_trajectories; // x, y, t, lavd
std::vector< LAVD_state > all_states;

void filter_trajectory(trajectory_t& out, const trajectory_t& in) {
    out.clear();
    if (in.empty()) {
        _log_(0) << "WARNING: empty trajectory in input" << std::endl;
        return;
    }
    double mean_dist=0;
    out.push_back(in[0]);
    vec2 last = nvis::subv<0, 2, value_t, 4>(in[0]);
    for (size_t i=1; i<in.size()-1; ++i) {
        const vec4& p = in[i];
        vec2 cur = nvis::subv<0, 2, value_t, 4>(p);
        // mean_dist += nvis::norm(cur - nvis::subv<0, 2, value_t, 4>(in[i-1]));
        if (nvis::norm(cur-last) > min_dist) {
            out.push_back(p);
            last = cur;
        }
    }
    if (in.size()>1) {
        out.push_back(in.back());
        // mean_dist /= (value_t)(in.size()-1);
        // std::cout << "mean step size on trajectory: " << mean_dist << '\n';
    }
}

void extract_parameters() {
    Nrrd* nin = xavier::nrrd_utils::readNrrd(name_in);
    int nbcomments = nin->cmtArr->len;
    value_t start_time;
    for (int i=0; i<nbcomments; ++i) {
        std::vector<std::string> strs;
        xavier::tokenize(strs, nin->cmt[i], " =[],");
        if (strs[0] == "input") {
            name_in = strs[2];
        }
        else if (strs[0] == "initial") {
            start_time = std::stod(strs[2]);
        }
        else if (strs[0] == "current") {
            t_init = std::stod(strs[2]);
        }
        else if (strs[0] == "epsilon") {
            eps = std::stod(strs[1]);
        }
        else if (strs[0] == "resolution") {
            std::vector<std::string> _strs;
            xavier::tokenize(_strs, strs[1], "x");
            res[0] = std::stoi(_strs[0]);
            res[1] = std::stoi(_strs[1]);
        }
        else if (strs[0] == "bounds") {
            std::vector<std::string> _strs;
            bnds[0] = std::stod(strs[1]);
            bnds[1] = std::stod(strs[2]);
            bnds[2] = std::stod(strs[4]);
            bnds[3] = std::stod(strs[5]);
        }
        else if (strs[0] == "dt") {
            dt = std::stod(strs[1]);
        }
    }

    size_t nb = res[0]*res[1];
    all_trajectories.resize(nb);
    all_states.resize(nb);

    xavier::ProgressDisplay progress(false);
    progress.start(nb, "importing seeds and states", 100);
    progress.set_active(true);

    #pragma omp parallel
    {
        #pragma omp for schedule(static,1)
        for (size_t n=0; n<nb ; ++n) {
            #if _OPENMP
            const int thread=omp_get_thread_num();
            #else
            const int thread=0;
            #endif

            if (!thread) progress.update(n);

            all_trajectories[n].resize(1);
            vec4& p = all_trajectories[n].back();
            p[0] = nrrd_value<value_t>(nin, 4*n  );
            p[1] = nrrd_value<value_t>(nin, 4*n+1);
            p[2] = nrrd_value<value_t>(nin, 4*n+2);
            p[3] = nrrd_value<value_t>(nin, 4*n+3);

            LAVD_state& state = all_states[n];
            value_t lavd = nrrd_value<value_t>(nin, 4*n+3);
            state.m_acc_vd = lavd;
            state.m_acc_time = t_init - start_time;
            state.m_time = t_init;
            state.m_stopped = (state.m_acc_vd == -1);
        }
    }
    nrrdNuke(nin);
}

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
        << static_cast<size_t>((t_max-t_init)/t_between_files-start_id+2*support_radius)
        << ", n_used=" << n_used << std::endl;

    nuke_current_volumes();

    // initialize velocity and vorticity volumes
    std::vector< Nrrd* > vel_tsteps(n_used);
    std::vector< Nrrd* > vor_tsteps(n_used);
    size_t first = ( start_id >= support_radius ?
                     start_id-support_radius :
                     static_cast<size_t>(0) );
    for (size_t i=first; i<first+n_used; ++i) {
        vel_tsteps[i-first]=xavier::nrrd_utils::readNrrd(vel_filenames[i]);
        _log_(1) << "Read " << vel_filenames[i] << std::endl;
        print(vel_tsteps[i-first]);

        vor_tsteps[i-first]=xavier::nrrd_utils::readNrrd(vor_filenames[i]);
        _log_(1) << "Read " << vor_filenames[i] << std::endl;
    }
    last_time_step = first + n_used - 1;
    value_t outer_tmin = /*t_init +*/ first*t_between_files;
    _log_(2) << "outer_tmin=" << outer_tmin << std::endl;
    current_velocity_volume = xavier::lavd::create_nrrd_volume(vel_tsteps, outer_tmin, t_between_files);
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

    // compute average vorticity
    value_t* av_array = (value_t *)calloc(n_used, sizeof(value_t));
    value_t region_size = domain.size()[0]*domain.size()[1];
    _log_(1) << "sampling vorticity slices" << std::endl;

    size_t nx = vor_tsteps[0]->axis[0].size;
    size_t ny = vor_tsteps[0]->axis[1].size;
    value_t minx = domain.min()[0];
    value_t miny = domain.min()[1];
    value_t dx = input_spc[0];
    value_t dy = input_spc[1];

    _log_(1) << "nx=" << nx << ", ny=" << ny << '\n';
    _log_(1) << "dx=" << dx << ", dy=" << dy << '\n';

    size_t mini = 0;
    size_t minj = 0;
    size_t maxi = nx-1;
    size_t maxj = ny-1;
    size_t nsamples = (maxi-mini+1)*(maxj-minj+1);

    _log_(1)
        << "index range for selected region: ["
        << mini << ", " << maxi << "] x ["
        << minj << ", " << maxj << "]" << std::endl;

    xavier::ProgressDisplay progress(false);
    progress.start(n_used*nsamples, "avg. vort.");
    progress.set_active(true);
    size_t base_n=0;
    size_t ncols=maxi-mini+1;
    std::vector< value_t > sum_copies(nb_threads, 0);
    std::vector< size_t > nvalid_copies(nb_threads, 0);
    for (size_t n=0; n<n_used; ++n, base_n+=nsamples) {
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for (size_t k=0; k<nsamples; ++k) {
                #if _OPENMP
                const int thread=omp_get_thread_num();
                #else
                const int thread=0;
                #endif

                if (!thread) progress.update(base_n+k);
                int j=k/ncols;
                int i=k%ncols;

                size_t id=(minj+j)*nx+(mini+i);

                if (nrrd_value(border_mask, id) != 0) {
                    sum_copies[thread] += nrrd_value(vor_tsteps[n], id);
                    ++nvalid_copies[thread];
                }
            }
        }
        value_t sum = std::accumulate(sum_copies.begin(), sum_copies.end(), static_cast<value_t>(0));
        size_t nvalid = std::accumulate(nvalid_copies.begin(), nvalid_copies.end(), static_cast<size_t>(0));
        av_array[n] = sum/static_cast<value_t>(nvalid);
        _log_(1) << std::setprecision(12) << "average vorticity[" << n << "]=" << av_array[n] << std::endl;
        _log_(1) << "nsamples=" << nsamples << ", nvalid_samples=" << nvalid << std::endl;
    }
	progress.end();

    Nrrd* av=nrrdNew();
    if (nrrdWrap_nva(av, av_array, xavier::nrrd_utils::nrrd_value_traits_from_type<value_t>::index,
                     1, &n_used)) {
        throw std::runtime_error(xavier::nrrd_utils::error_msg("unable to create 1D average vorticity nrrd"));
    }
    // Per-axis info
    int center = nrrdCenterNode;
    nrrdAxisInfoSet_nva(av, nrrdAxisInfoSpacing, &t_between_files);
    nrrdAxisInfoSet_nva(av, nrrdAxisInfoMin, &outer_tmin);
    nrrdAxisInfoSet_nva(av, nrrdAxisInfoCenter, &center);

    // upsampling average vorticity
    // (much cheaper than to average upsampled vorticity!)
    current_average_vorticity_vector = NrrdScalarField<1>::make3d(av);

    // clean up
    for (size_t i=0; i<vor_tsteps.size(); ++i) {
        nrrdNuke(vor_tsteps[i]);
    }
    nrrdNuke(av);
}

void export_results(double current_time, double wall_time, double cpu_time, bool final) {

    std::ostringstream os;
    ios::fmtflags default_settings = os.flags();

    boost::filesystem::path out_path(name_out);
    std::string parent = out_path.parent_path().string();
    std::string filename = out_path.filename().string();
    std::string stem = out_path.stem().string();

    std::string prefix = parent + "/trajectories";

    if (long_name) {
        os << '_' << stem << "_started_at_" << start_time_str << "_"
           << std::setw(5) << std::setfill('0')
           << current_time/xavier::lavd::HOUR << "h_";
        os << std::setprecision(2) << std::scientific << eps;
    }
    else {
        os << '_' << stem << "_" << std::setw(5) << std::setfill('0')
           << current_time/xavier::lavd::HOUR << "h";
    }
    os.flags(default_settings);
    std::string suffix = os.str() + ".nrrd";

    std::string name = prefix + suffix;

    std::vector<std::string> comments;
    os.clear(); os.str("");
    os << "This file was produced by " << me << '.';
    comments.push_back(os.str()); os.clear(); os.str("");
    comments.push_back("1st axis: (flowmap_x, flowmap_y, t, lavd, line_id),");
    comments.push_back("lavd: \"Lagrangian averaged vorticity deviation\"");
    comments.push_back("2nd axis: consecutive positions along time axis of particles.");
    comments.push_back("Computation parameters: ");
    comments.push_back("input file: " + name_in);
    comments.push_back("mask file: " + name_mask);
    os << "Considered region: min=" << region.min() << " -- max=" << region.max();
    comments.push_back(os.str()); os.clear(); os.str("");
    os << std::setprecision(10);
    os << "initial time=" << t_init;
    comments.push_back(os.str()); os.clear(); os.str("");
    os << "current time=" << current_time;
    comments.push_back(os.str()); os.clear(); os.str("");
    os << "epsilon=" << eps;
    comments.push_back(os.str()); os.clear(); os.str("");
    os << "dt=" << dt << " s";
    comments.push_back(os.str()); os.clear(); os.str("");
    os << "kernel size=" << support_radius << '\n';
    comments.push_back(os.str());

    typedef nvis::fixed_vector<float, 5> fvec5;
    std::vector< fvec5 > out_array;

    for (size_t i=0; i<all_trajectories.size(); ++i) {
        auto& orig = all_trajectories[i];
        trajectory_t t;
        filter_trajectory(t, orig);
        for (size_t n=0; n<t.size(); ++n) {
            const auto& pt = t[n];
            out_array.push_back(fvec5(pt[0], pt[1], pt[2], pt[3], i));
        }
        auto pt = orig.back();
        orig.clear();
        orig.push_back(pt);
    }

    std::array<size_t, 2> dims( {5, out_array.size()} );
    std::array<double, 2> spc({ AIR_NAN, AIR_NAN }), min({ AIR_NAN, AIR_NAN });
    std::array<int, 2> ctr({ nrrdCenterUnknown, nrrdCenterNode });
    std::cout << "exporting " << name << "...\n";
    try {
        xavier::nrrd_utils::writeNrrdFromContainers((float*)&out_array[0], name, dims, spc, min, ctr, comments);
    }
    catch(std::exception& e) {
        std::cout << "WARNING: exception caught: " << e.what() << '\n';
    }
    std::cout << "done\n";

    prefix = parent + "/failed_trajectories";
    name = prefix + suffix;

    out_array.clear();

    for (size_t i=0; i<failed_paths.size(); ++i) {
        const auto& t = failed_paths[i];
        for (size_t n=0; n<t.size(); ++n) {
            const auto& pt = t[n];
            out_array.push_back(fvec5(pt[0], pt[1], pt[2], pt[3], i));
        }
    }
    dims[1] = out_array.size();
    if (dims[1] > 0) {
        std::cout << "exporting " << name << "...\n";
        xavier::nrrd_utils::writeNrrdFromContainers((float*)&out_array[0], name, dims, spc, min, ctr, comments);
        std::cout << "done\n";
    }
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

    if (resume) {
        extract_parameters();
    }

    if (t_max < 0) {
        _log_(0) << "Invalid time duration symbol in definition of tmax" << std::endl;
        clean_exit(1);
    }

    if (dt < 0) {
        dt = (t_max - t_init)/(fabs(dt)-1);
    }

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
	
    // if user provided bounds for computation, use them, otherwise use the
    // domain's bounds
    if ( !( bnds[0]==bnds[1]==bnds[2]==bnds[3] ) ) {
        region.min()=vec2( bnds[0], bnds[1] );
        region.max()=vec2( bnds[2], bnds[3] );
    }
    else {
        region = domain;
    }
	std::vector<nvis::vec2> seeds;
    nvis::vec2 span = region.max() - region.min();
    nvis::vec2 spc = span / nvis::vec2(res[0]-1, res[1]-1);
    seeds.resize(res[0]*res[1]);
    for (int j=0; j<res[1]; ++j) {
        for (int i=0; i<res[0]; ++i) {
            seeds[j*res[0]+i] = region.min() + nvis::vec2(i,j)*spc;
        }
    }
    nb_samples = seeds.size();

    // storage for LAVD states and trajectories
    if (!resume) {
        all_states.resize(nb_samples);
        all_trajectories.resize(nb_samples);
    }

    _log_(0) << nb_samples << " seed point in input\n";

    support_radius = xavier::lavd::compute_support_radius();
    _log_(1) << "support radius = " << support_radius << std::endl;

    // initialize velocity and vorticity volumes
    import_data(velocity_filenames, vorticity_filenames, t_init/HOUR/3);

    // current_velocity_field->save_as_nrrd("upsampled_data.nrrd");

    size_t counter = 0;
    size_t nb_early = 0;

    _log_(1) << "bounds=" << domain << std::endl;

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

    std::vector<size_t> nb_lost_copies(nb_threads);
    std::vector<size_t> counter_copies(nb_threads);

    bool initial_loop = !resume;
    size_t niter = 0;
    current_start_time = t_init;
    next_export_time = t_init + t_export_step;
    while ( true ) {

        _log_(0, "\n\n") << "next target time for this loop: "
            << sec2time(next_export_time)
            << " (" << next_export_time << " s.)" << std::endl;

        // counter = 0;
        std::fill(counter_copies.begin(), counter_copies.end(), 0);
        // nb_lost = 0;
        std::fill(nb_lost_copies.begin(), nb_lost_copies.end(), 0);
        size_t nb_early = 0;
        value_t current_target_time = std::min(current_t_max, next_export_time);

        progress.start(nb_samples, "flow map");
        progress.set_active(true);
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for (size_t n = 0 ; n < nb_samples ; ++n) {

                LAVD_state& _state = all_states[n];
                trajectory_t& _traj = all_trajectories[n];

                #if _OPENMP
                const int thread=omp_get_thread_num();
                #else
                const int thread=0;
                #endif

                seed_id = n;

                if (!thread) progress.update(n);

                if (_state.m_stopped) {
                    ++nb_lost_copies[thread];
                    continue;
                }

                vec2 x0;
                if (initial_loop) {
                    x0 = seeds[n];
#ifndef __SKIP_LOGS__
                    _log_(2) << "\n\n\n\nstarting integration #" << n << " at " << x0 << std::endl;
#endif
                }
                else {
                    x0 = _state.m_pos;
                }

                ++counter_copies[thread];

                // create a stepper
                auto stepper = make_controlled(eps, eps, runge_kutta_dopri5<vec2>());

                NrrdODERHS& rhs=*rhs_copies[thread];
                rhs.m_counter=0;
                
                NrrdScalarField<1> my_avg_vorticityf(current_average_vorticity_vector, std::string("average vorticity field for thread #") + std::to_string(thread));
                NrrdScalarField<3> my_vorticityf(current_vorticity_volume, std::string("vorticity field for thread #" + std::to_string(thread))); 

                // empty contents of trajectory and retain only last element
                // since this is all we need to proceed
                if (!initial_loop && !_traj.empty()) {
                    vec4 last_p = _traj.back();
                    _traj.empty();
                    _traj.push_back(last_p);
                }
                observer an_observer(_state, _traj, my_vorticityf, my_avg_vorticityf);

                if (initial_loop) {
#ifndef __SKIP_LOGS__
                    _log_(2) << "initializing observer at " << x0 << std::endl;
#endif
                    try {
                        an_observer.initialize(x0, current_t_min);
                    }
                    catch( std::exception& e ) {
                        _log_(1) << "caught exception while attempting "
                            << "to initialize observer at " << x0 << std::endl;
                        _state.stop();
                        _traj.push_back(vec4(x0[0], x0[1], current_t_min, all_states[n].m_acc_vd));
                        ++nb_lost_copies[thread];
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
                            << " and ended at " << _traj.back()
                            << " (" << _traj.size()-1 << " steps)"
                            << '\n' << std::flush;
                    all_states[n].stop();
                    ++nb_lost_copies[thread];

                    // DEBUG
                    failed_paths.push_back(std::vector<vec4>());
                    std::vector<vec4>& path = failed_paths.back();
                    LAVD_state another_state;
                    trajectory_t another_traj;
                    observer another_observer(another_state, another_traj,
                                              my_vorticityf,
                                              my_avg_vorticityf);
                    if (initial_loop) {
                        try {
                            another_observer.initialize(x0, current_t_min);
                        }
                        catch( std::exception& e ) {

                        }
                    }
                    rhs.m_counter=0;
                }
                catch(...) {
                    _log_(0) << "unknown exception thrown!" << std::endl;
                }
#ifndef __SKIP_LOGS__
                _log_(2) << rhs.m_counter << " RHS evaluations for this IC"
                    << std::endl;
#endif
            }

            initial_loop = false;
        }
	    progress.end();

        nb_lost = std::accumulate(nb_lost_copies.begin(), nb_lost_copies.end(), 0);

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
            export_results(next_export_time,
                           progress.wall_time(),
                           progress.cpu_time(), false);
            next_export_time += t_export_step;
            if (next_export_time > t_max) next_export_time = t_max;
        }
        // advance starting time
        current_start_time = current_target_time;
        _log_(1) << "current start time=" << current_start_time << std::endl;
    }
    total_progress.end();
    export_results(t_max, total_progress.wall_time(),
                   total_progress.cpu_time(), true);

    return 0;
}
