#include <flow/lavd.hpp>

#include <sstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <cctype>
#include <stdexcept>
#include <locale>
#include <iomanip>

#include <math/types.hpp>
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
#include <vtk/vtk_utils.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace odeint = boost::numeric::odeint;

using namespace spurt::lavd;

typedef spurt::image< vec3, 3, value_t, size_t >    vector_field_t;
typedef spurt::image< value_t, 3, value_t, size_t > scalar_field_t;
typedef spurt::image< value_t, 1, value_t, size_t > scalar_field_1d_t;

inline vec3 append(const vec2& v, value_t t) {
    return vec3(v[0], v[1], t);
}

inline vec2 crop(const vec3& v) {
    return vec2(v[0], v[1]);
}


struct lavd_window {
    bbox_t   m_bounds;
    lvec2    m_res;
    lvec2    m_offset;
    value_t* m_lavd;
    value_t* m_times;
    
    static const value_t invalid_time;
    
    lavd_window(const bbox_t& bounds=bbox_t(), 
                const lvec2& res=lvec2(0), const lvec2& offset=lvec2(0)) 
            : m_bounds(bounds), m_res(res), 
              m_offset(offset), m_lavd(0), m_times(0) {}
    
    void initialize(const bbox_t& bounds, const lvec2& res, 
                    const lvec2& offset) {
        m_bounds = bounds;
        m_res = res; 
        m_offset = offset; 
        allocate();
    }
    
    void allocate() {
        m_lavd = (value_t *)calloc(m_res[0]*m_res[1], sizeof(value_t));
        m_times = (value_t *)calloc(m_res[0]*m_res[1], sizeof(value_t));
        std::fill(m_times, m_times+(m_res[0]*m_res[1]), lavd_window::invalid_time); // invalid default time
    }
    
    ~lavd_window() {
        if (m_lavd != NULL) free(m_lavd);
        if (m_times != NULL) free(m_times);
    }
    
    value_t step_x() const {
        return m_bounds.size()[0]/static_cast<value_t>(m_res[0]-1);
    }
    
    value_t step_y() const {
        return m_bounds.size()[1]/static_cast<value_t>(m_res[1]-1);
    }
    
    std::pair<bool, size_t> index(const lvec2& coord) const {
        lvec2 c(coord);
        c -= m_offset;
        if ((c[0] >= m_res[0]) || (c[1] >= m_res[1])) {
            return std::pair<bool, size_t>(false, 0);
        }
        else return std::pair<bool, size_t>(true, c[0] + c[1]*m_res[0]);
    }
    
    value_t& lavd(const lvec2& coord, bool is_global=true) {
        lvec2 c(coord);
        if (is_global) c -= m_offset;
        return m_lavd[c[0] + m_res[0]*c[1]];
    }
    
    const value_t& lavd(const lvec2& coord, bool is_global=true) const {
        lvec2 c(coord);
        if (is_global) c -= m_offset;
        return m_lavd[c[0] + m_res[0]*c[1]];
    }
    
    value_t& time(const lvec2& coord, bool is_global=true) {
        lvec2 c(coord);
        if (is_global) c -= m_offset;
        return m_times[c[0] + m_res[0]*c[1]];
    }
    
    const value_t& time(const lvec2& coord, bool is_global=true) const {
        lvec2 c(coord);
        if (is_global) c -= m_offset;
        return m_times[c[0] + m_res[0]*c[1]];
    }
    
    void update(size_t n, const std::vector<vec3>& pos, 
                const std::vector<value_t>& vort,
                scalar_field_1d_t& av_vort) {
        if (m_times[n] == invalid_time) {
            m_times[n] = pos[0][2];
        }
        for (size_t i=0; i<pos.size(); ++i) {
            m_lavd[n] += (pos[i][2] - m_times[n])*fabs(vort[i] - av_vort.value(vec1(pos[i][2])));
            m_times[n] = pos[i][2];
        }
    }
};

std::ostream& operator<<(std::ostream& os, const lavd_window& w) {
    os << "[res=" << w.m_res 
        << ", bounds=" << w.m_bounds.min() << "-" << w.m_bounds.max()
        << ", offset=" << w.m_offset << "]";
    return os;
}

const value_t lavd_window::invalid_time = -1;

struct ode_state {
    ode_state(std::vector<vec3>& pos)
        : m_pos(pos), m_stopped(false) {
        m_pos.clear();
    }
    
    ode_state(std::vector<vec3>& pos, const vec2& p0, value_t t0) 
        : m_pos(pos), m_stopped(false) {
        m_pos.clear();
        m_pos.push_back(append(p0, t0));
    }
    
    ode_state(std::vector<vec3>& pos, const vec3& x0) 
        : m_pos(pos), m_stopped(false) {
        m_pos.clear();
        m_pos.push_back(x0);
    }
    
    void initialize(const vec2& p, value_t t) {
        m_pos.clear();
        this->update(p, t);
        m_stopped = false;
    }
    
    void update(const vec2& p, value_t t) {
        m_pos.push_back(append(p, t));
    }
    
    void stop() {
        m_stopped = true;
    }
    
    bool m_stopped;
    std::vector< vec3 >&  m_pos;
};

namespace params {
    
std::vector< lavd_window > lavd_windows;
size_t nb_windows=0; // invariant: nb_windows == lavd_windows.size()

std::string name_in, name_out;
std::string me;
std::string log_name = "log.txt";
std::string border_mask_name;
std::ofstream log_file;

} // params

spurt::log::dual_ostream spurt::lavd::_log_(params::log_file, std::cout, 1, 0, 0, true);

size_t nb_threads;

namespace params {

// time-related
value_t t_max=2*spurt::lavd::DAY;
value_t t_init=0;
value_t t_skip=0;
value_t t_export_step=3*spurt::lavd::HOUR;
value_t t_between_files=3*spurt::lavd::HOUR;
value_t eps=1.0e-8;
value_t dt=30;
size_t nslices = 8*10; // 10 days of data processed simultaneously
size_t last_time_step;
value_t next_export_time, current_start_time;
value_t current_t_min; // starting time
value_t current_t_max; // stopping time
std::string start_time_str;

// space-related
std::array<size_t, 2> res( { 512, 512 } );
bbox_t region; // invariant: region == bbox(\sigma lavd_windows)
bbox_t domain; // invariant: region C domain
value_t step_x, step_y;
size_t nb_samples; // invariant: nb_samples = res[0]*res[1]
size_t nb_lost=0;
size_t support_radius;
svec3 up(1);
vec2 input_spc;

// i/o-related
std::vector< int > verbose_vec;
int verbose_cout=0;
int verbose_log=1;
bool export_region = false;
bool long_name = false;
bool export_trajectories = false;
size_t max_rhs_evals = 1000000;

// storage-related
Nrrd* border_mask=0;
std::shared_ptr< vector_field_t >    velocity_field;
std::shared_ptr< scalar_field_t >    vorticity_field;
std::vector< std::shared_ptr< scalar_field_1d_t > > average_vorticity_fields;
std::vector< std::vector< vec2 > > trajectories;
std::vector< std::shared_ptr< ode_state > > states;
std::vector< std::vector< vec3 > > last_steps;
std::vector< std::vector< value_t > > last_vorticities;

} // params

void set_verbose_values() {
    if (params::verbose_vec.size()<1 || params::verbose_vec.size()>2) {
        throw std::runtime_error("Invalid verbose values");
    }
    params::verbose_cout = params::verbose_vec[0];
    if (params::verbose_vec.size()==1) {
        params::verbose_log = params::verbose_cout;
    }
    else params::verbose_log = params::verbose_vec[1];
    _log_(0) << "verbose_cout=" << params::verbose_cout << ", verbose_log="
        << params::verbose_log << std::endl;
}

void validate_window_bounds(const std::vector<value_t>& bnds) {
    if (bnds.size()%4) {
        throw std::runtime_error("Invalid windows definitions");
    }
    params::region.reset();
    if (!bnds.empty()) {
        params::lavd_windows.resize(bnds.size()/4);
        params::nb_windows = params::lavd_windows.size();
        for (size_t i=0; i<params::nb_windows; ++i) {
            bbox_t& b = params::lavd_windows[i].m_bounds;
            b.min()[0] = bnds[4*i  ];
            b.min()[1] = bnds[4*i+1];
            b.max()[0] = bnds[4*i+2];
            b.max()[1] = bnds[4*i+3];
            params::region.add(b.min());
            params::region.add(b.max());
        }
    }
    else {
        params::lavd_windows.clear();
        params::nb_windows = 0;
        // TBD
        // single window will be set to global domain bounds (== region)
    }
}

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;
        
    xcl::option_traits 
            required_group(true, false, "Required Options"), 
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compute Lagrangian averaged vorticity deviation");
    
    using namespace params;
    
    verbose_vec.resize(2);    
    verbose_vec[0] = verbose_cout;
    verbose_vec[1] = verbose_log;
    
    std::vector<std::string> t_init_str(1, "0");
    std::string t_max_str = "2d";
    std::string t_expstp_str = "3h";
    std::string dt_str="30";
    std::string t_skip_str="0";
    std::vector< value_t > windows_bnds;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input info filename", required_group);
        parser.add_value("output", name_out, "Output file basename ", required_group);
        parser.add_sequence("tinit", t_init_str, t_init_str, "Integration starting time, expressed either in time since initial time step (in seconds; for minutes append \"m\", for hours \"h\", for days \"d\", for weeks \"w\") or as a date MM/DD", required_group, "<float>[<char>]");
        parser.add_value("tmax", t_max_str, "Integration length (in seconds; for minutes append \"m\", for hours \"h\", for days \"d\", for weeks \"w\")", required_group, "<float>[<char>]");
        parser.add_value("step", t_expstp_str, t_expstp_str, "Time interval between intermediate result output", optional_group);
        parser.add_tuple<3>("up", up, up, "Upsampling factor before interpolation", optional_group);
        parser.add_value("skip", t_skip_str, t_skip_str, "Time skipped before first output", optional_group);
        parser.add_value("dt", dt_str, dt_str, "Time sampling of integral curve for LAVD computation (in seconds, for minutes add \"m\", for hours \"h\", for times \"x\")", optional_group, "<float>[<char>]");
        parser.add_tuple<2>("res", res, res, "Sampling resolution", optional_group);
        parser.add_sequence("windows", windows_bnds, windows_bnds, "LAVD windows (series of bounds definitions)", optional_group, "(4x<float>)+");
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
        validate_window_bounds(windows_bnds);
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

inline void nuke_current_volumes() {
    params::velocity_field.reset();
    params::vorticity_field.reset();
    std::for_each(params::average_vorticity_fields.begin(),
                 params::average_vorticity_fields.end(),
                 [](std::shared_ptr<scalar_field_1d_t> sptr){ sptr.reset(); });
}

void clean_exit(int n) {
    if (params::log_file.is_open()) params::log_file.close();
    nuke_current_volumes();
    exit(n);
}

// lightweight observer that only records space time coordinates along the way
struct observer
{   
    ode_state& m_state;
    
    observer(ode_state& state) 
        :  m_state(state) {}
    
    void initialize(const vec2& p0, value_t t0) {
        m_state.initialize(p0, t0);
    }

    void operator()(const vec2 &x, value_t t)
    {  
        m_state.m_pos.push_back(append(x, t));
    }
};

void import_data(const std::vector< std::string >& vel_filenames,
                 const std::vector< std::string >& vor_filenames,
                 size_t start_id)
{
    size_t n_available = std::min(params::nslices, vel_filenames.size() - start_id + params::support_radius);
    size_t n_used = std::min(n_available, static_cast<size_t>((params::t_max-params::t_init)/params::t_between_files+2*params::support_radius));
    _log_(1) 
        << "import_data: start_id=" << start_id
        << ", #filenames=" << vel_filenames.size()
        << ", n_available=" << n_available
        << ", (tmax-tinit)/t_between_files=" << (params::t_max-params::t_init)/params::t_between_files
        << ", number of time steps needed to complete entire integration length="
        << static_cast<size_t>((params::t_max-params::t_init)/params::t_between_files - start_id + 2*params::support_radius)
        << ", n_used=" << n_used << std::endl;
    
    nuke_current_volumes();
    
    // initialize velocity and vorticity volumes
    std::vector< Nrrd* > vel_tsteps(n_used);
    std::vector< Nrrd* > vor_tsteps(n_used);
    size_t first = ( start_id >= params::support_radius ? 
                     start_id - params::support_radius : 
                     static_cast<size_t>(0) );
    for (size_t i=first; i<first+n_used; ++i) {
        vel_tsteps[i-first]=spurt::nrrd_utils::readNrrd(vel_filenames[i]);
        _log_(1) << "Read " << vel_filenames[i] << std::endl;
        print(vel_tsteps[i-first]);
        
        vor_tsteps[i-first]=spurt::nrrd_utils::readNrrd(vor_filenames[i]);
        _log_(1) << "Read " << vor_filenames[i] << std::endl;
    }
    params::last_time_step = first + n_used - 1;
    value_t outer_tmin = /*t_init +*/ first*params::t_between_files;
    _log_(2) << "outer_tmin=" << outer_tmin << std::endl;
    Nrrd* vel3d = spurt::lavd::create_nrrd_volume(vel_tsteps, outer_tmin, params::t_between_files);
    params::velocity_field.reset(spurt::lavd::upsample_vector(vel3d, params::up, "velocity"));
    
    // clean up
    for (size_t i=0; i<vel_tsteps.size(); ++i) {
        nrrdNuke(vel_tsteps[i]);
    }
    nrrdNuke(vel3d);
    
    Nrrd* vort3d = spurt::lavd::create_nrrd_volume(vor_tsteps, outer_tmin, params::t_between_files);
    params::vorticity_field.reset(spurt::lavd::upsample_scalar(vort3d, params::up, "vorticity"));
    
    // clean up
    nrrdNuke(vort3d);
    
    params::current_t_min = start_id*params::t_between_files;
    params::current_t_max = params::last_time_step*params::t_between_files;
    if (params::current_t_max < params::t_max) params::current_t_max -= params::support_radius*params::t_between_files;
    else params::current_t_max = params::t_max;
    _log_(1) 
        << "n_used=" << n_used << '\n'
        << "support_radius=" << params::support_radius << '\n'
        << "current_t_max=" << params::current_t_max << '\n';
    
    size_t nx = vor_tsteps[0]->axis[0].size;
    size_t ny = vor_tsteps[0]->axis[1].size;
    value_t minx = params::domain.min()[0];
    value_t miny = params::domain.min()[1];
    value_t dx = params::input_spc[0];
    value_t dy = params::input_spc[1];
    
    _log_(1) << "nx=" << nx << ", ny=" << ny << '\n';
    _log_(1) << "minx=" << minx << ", miny=" << miny << '\n';
    _log_(1) << "dx=" << dx << ", dy=" << dy << '\n';
    
    params::average_vorticity_fields.resize(params::nb_windows);
    for (size_t i=0; i<params::nb_windows; ++i) {
        lavd_window& window = params::lavd_windows[i];
        bbox_t& subregion = window.m_bounds;
        
        _log_(1)
            << "selected subregion #" << i << " = " << subregion << '\n'
            << "x spacing: " << params::input_spc[0] << '\n'
            << "y spacing: " << params::input_spc[1] << std::endl;
        
        // compute average vorticity for selected region
        value_t* av_array = (value_t *)calloc(n_used, sizeof(value_t));
        value_t subregion_size = subregion.size()[0]*subregion.size()[1];
        _log_(1) << "sampling vorticity slices" << std::endl;
        
        double min_offset_x = subregion.min()[0]-minx;
        double min_offset_y = subregion.min()[1]-miny;
        double max_offset_x = subregion.max()[0]-minx;
        double max_offset_y = subregion.max()[1]-miny;
        _log_(1) << "min offset x=" << min_offset_x << '\n';
        _log_(1) << "min offset y=" << min_offset_y << '\n';
        _log_(1) << "max offset x=" << max_offset_x << '\n';
        _log_(1) << "max offset y=" << max_offset_y << std::endl;
        
        size_t mini = min_offset_x <= 0 ? 0 : static_cast<size_t>(min_offset_x/dx);
        size_t minj = min_offset_y <= 0 ? 0 : static_cast<size_t>(min_offset_y/dy);
        size_t maxi = max_offset_x <= 0 ? 0 : std::min(static_cast<size_t>(max_offset_x/dx), nx-1);
        size_t maxj = max_offset_y <= 0 ? 0 : std::min(static_cast<size_t>(max_offset_y/dy), ny-1);
        
        if (mini>=nx || mini>maxi) _log_(1) << "invalid mini value: " << mini << std::endl;
        if (minj>=ny || minj>maxj) _log_(1) << "invalid minj value: " << minj << std::endl;
        if (maxi>=nx) _log_(1) << "invalid maxi value: " << maxi << std::endl;
        if (maxj>=ny) _log_(1) << "invalid maxj value: " << maxj << std::endl;
        size_t nsamples = (maxi-mini+1)*(maxj-minj+1);

        _log_(1) 
            << "index range for selected region: [" 
            << mini << ", " << maxi << "] x [" 
            << minj << ", " << maxj << "]" << std::endl;
        
        spurt::ProgressDisplay progress(false);
        progress.start(n_used*nsamples, "avg. vort.");
        progress.set_active(true);
        size_t base_n=0;
        size_t ncols=maxi-mini+1;
        std::vector< value_t > sum_copies(nb_threads, 0);
        std::vector< size_t > nvalid_copies(nb_threads, 0);
        for (size_t n=0; n<n_used; ++n, base_n+=nsamples) {
            #pragma omp parallel
            {
                #pragma omp for schedule(static,1)      
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
                    
                    if (nrrd_value(params::border_mask, id) != 0) {
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
        if (nrrdWrap_nva(av, av_array, spurt::nrrd_utils::nrrd_value_traits_from_type<value_t>::index, 
                         1, &n_used)) {
            throw std::runtime_error(spurt::nrrd_utils::error_msg("unable to create 1D average vorticity nrrd"));
        }
        // Per-axis info
        int center = nrrdCenterNode;
        nrrdAxisInfoSet_nva(av, nrrdAxisInfoSpacing, &params::t_between_files);
        nrrdAxisInfoSet_nva(av, nrrdAxisInfoMin, &outer_tmin);
        nrrdAxisInfoSet_nva(av, nrrdAxisInfoCenter, &center);
        
        // upsampling average vorticity
        // (much cheaper than to average upsampled vorticity!)
        params::average_vorticity_fields[i].reset(spurt::lavd::upsample<value_t, NrrdScalarField<1>, 1>(av, lvec1(params::up[2]), "average vorticity #" + std::to_string(i)));
        nrrdNuke(av);
    }
    
    // clean up
    for (size_t i=0; i<vor_tsteps.size(); ++i) {
        nrrdNuke(vor_tsteps[i]);
    }
}

template< typename T = float > 
void export_results(value_t current_time, value_t wall_time, value_t cpu_time, 
                    bool final) 
{
    typedef T export_value_t;
    
    std::string qualifier;
    if (final) qualifier = "Overall";
    else qualifier = "Intermediate";
    
    _log_(1) << "\nEntering export_results\n";
    _log_(0) << "\n" << qualifier << " wall time was " << wall_time << " ms.\n";
    _log_(0) << qualifier << " cpu time was " << cpu_time << " ms.\n";
    _log_(0) << params::nb_lost << " interrupted trajectories ("
        << static_cast<value_t>(params::nb_lost)/static_cast<value_t>(params::nb_samples)*100.
        << "\%)\n";

    for (size_t i=0; i<params::nb_windows; ++i) { 
        std::ostringstream os;
        ios::fmtflags default_settings = os.flags();

        const lavd_window& window = params::lavd_windows[i];
        std::string basename;
        std::string comment;
        if (params::long_name) {
            os << params::name_out 
               << "_window_#" << std::setw(3) << std::setfill('0') << i 
               << "_started_at_" << params::start_time_str 
               << "_lavd_" << std::setw(5) << std::setfill('0') 
               << current_time/spurt::lavd::HOUR << "h_";
            os << window.m_res[0] << "x" << window.m_res[1]  // reset default width / fill values
               << "_" << std::setprecision(2) << std::scientific << params::eps;
        }
        else {
            os << params::name_out 
               << "_window_#" << std::setw(3) << std::setfill('0') << i
               << "_lavd_" << std::setw(5) << std::setfill('0') 
               << current_time/spurt::lavd::HOUR << "h";
        }
        os.flags(default_settings);
        basename = os.str();
        std::string filename = basename + ".nrrd";
        
        os.clear();
        os.str("");
        os << "This file was produced by " << params::me
           << "; data = lavd. ";
        os << "Computation parameters: ";
        os << "input file=" << params::name_in << "; ";
        os << "tmax=" << current_time << "; ";
        os << "epsilon=" << params::eps << "; ";
        os << "resolution=" << window.m_res[0] << "x" << window.m_res[1] << "; ";
        os << "bounds=" << window.m_bounds.min() << "->" << window.m_bounds.max() << "; ";
        os << "dt=" << params::dt << "; ";
        os << "kernel size=" << params::support_radius << '\n';
        comment = os.str();
    
    
        _log_(1) << "Setting Nrrd header values... " << std::flush;
        std::vector<size_t> __res(2);
        __res[0] = window.m_res[0];
        __res[1] = window.m_res[1];
    
        std::vector<double> __mins(2);
        __mins[0] = window.m_bounds.min()[0];
        __mins[1] = window.m_bounds.min()[1];
    
        std::vector<double> __spc(2);
        __spc[0] = window.step_x();
        __spc[1] = window.step_y();
    
        std::vector<int> __ctr(2);
        __ctr[0] = nrrdCenterNode;
        __ctr[1] = nrrdCenterNode;
        
        _log_(1) << "done\n";
    
        _log_(1) << "Writing NRRD file under " << filename << "... " << std::flush;
        spurt::nrrd_utils::writeNrrdFromContainers(window.m_lavd, filename, 
                                        __res, __spc, __mins, __ctr, comment);
    }
    
    // export flow map separately
    {
        std::ostringstream os;
        ios::fmtflags default_settings = os.flags();
        
        std::string basename;
        std::string comment;
        if (params::long_name) {
            os << params::name_out 
               << "_started_at_" << params::start_time_str 
               << "_flowmap_" << std::setw(5) << std::setfill('0') 
               << current_time/spurt::lavd::HOUR << "h_";
            os << params::res[0] << "x" << params::res[1]  // reset default width / fill values
               << "_" << std::setprecision(2) << std::scientific << params::eps;
        }
        else {
            os << params::name_out 
               << "_flowmap_" << std::setw(5) << std::setfill('0') 
               << current_time/spurt::lavd::HOUR << "h";
        }
        os.flags(default_settings);
        basename = os.str();
        
        os.clear();
        os << "This file was produced by " << params::me
           << "; data = flow map; ";
        os << "Computation parameters: ";
        os << "input file=" << params::name_in << "; ";
        os << "tmax=" << current_time << "; ";
        os << "epsilon=" << params::eps << "; ";
        os << "resolution=" << params::res[0] << "x" << params::res[1] << "; ";
        os << "bounds=" << params::region.min() << "->" << params::region.max() << "; ";
        os << "dt=" << params::dt << "; ";
        os << "kernel size=" << params::support_radius << '\n';
        comment = os.str();
        
        std::string filename = basename + ".nrrd";
        
        _log_(1) << "Setting Nrrd header values... " << std::flush;
        std::vector<size_t> __res(3);
        __res[0] = 2;
        __res[1] = params::res[0];
        __res[2] = params::res[1];
        
        std::vector<double> __mins(3);
        __mins[0] = AIR_NAN;
        __mins[1] = params::region.min()[0];
        __mins[2] = params::region.min()[1];
        
        std::vector<double> __spc(3);
        __spc[0] = AIR_NAN;
        __spc[1] = params::step_x;
        __spc[2] = params::step_y;
        
        std::vector<int> __ctr(3);
        __ctr[0] = nrrdCenterUnknown;
        __ctr[1] = nrrdCenterNode;
        __ctr[2] = nrrdCenterNode;
        
        export_value_t* fmap = (export_value_t*)calloc(params::nb_samples*2, sizeof(export_value_t));
        for (size_t i=0; i<params::trajectories.size(); ++i) {
            fmap[2*i  ] = params::trajectories[i].back()[0];
            fmap[2*i+1] = params::trajectories[i].back()[1];
        }
        
        _log_(1) << "done\n";
        
        _log_(1) << "Writing NRRD file under " << filename << "... " << std::flush;
        spurt::nrrd_utils::writeNrrdFromContainers(fmap, filename, 
                                        __res, __spc, __mins, __ctr, comment);
        free(fmap);
    }
    _log_(1) << "Leaving export_results\n\n";
}

int main(int argc, const char* argv[])
{
    using namespace spurt;
    using namespace odeint;
    
    params::me=argv[0];
    
    nb_threads=1;
    
#if _OPENMP
    nb_threads = omp_get_max_threads();
    _log_.set_nthreads(nb_threads);
    _log_(1) << nb_threads << " threads available" << std::endl;
#else
    _log_.set_nthreads(1);
#endif
    
    initialize(argc, argv);
    
    params::log_file.open(params::log_name.c_str());
    _log_.set_thresholds(params::verbose_cout, params::verbose_log);
    
    params::region.reset();
    for (size_t i=0; i<params::lavd_windows.size(); ++i) {
        params::region.add(params::lavd_windows[i].m_bounds.min());
        params::region.add(params::lavd_windows[i].m_bounds.max());
    }
    if (!params::region.empty()) {
        _log_(1) << "overall sampling region: " << params::region << std::endl;
    }
    
    if (params::t_max < 0) {
        _log_(0) << "Invalid time duration symbol in definition of tmax" << std::endl;
        clean_exit(1);
    }
    
    if (params::dt < 0) {
        params::dt = (params::t_max - params::t_init)/(fabs(params::dt)-1);
    }
    
    params::nb_samples = params::res[0]*params::res[1];
    
    std::cout << "params::res = [" << params::res[0] << "," << params::res[1] << "]" 
              << ", params::nb_samples=" << params::nb_samples << '\n';
    
    params::name_out = spurt::filename::remove_extension(params::name_out);
    
    std::vector<std::string> velocity_filenames;
    std::vector<std::string> vorticity_filenames;
    
    std::fstream info_file(params::name_in, std::ios::in);
    if (!info_file) {
        _log_(0) << "Unable to open info file named " << params::name_in << std::endl;
        exit(1);
    }

    boost::filesystem::path p(params::name_in);
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
    get_spatial_info(params::domain, params::input_spc, velocity_filenames[0], 1);
    params::border_mask = get_border_mask(params::border_mask_name, velocity_filenames[0]);
    
    if ( params::region.empty() ) {
        params::region.min() = params::domain.min();
        params::region.max() = params::domain.max();
        params::lavd_windows.push_back(lavd_window());
        params::lavd_windows[0].m_bounds = params::region;
        params::lavd_windows[0].m_res[0] = params::res[0];
        params::lavd_windows[0].m_res[1] = params::res[1];
        params::lavd_windows[0].m_offset = lvec2(0);
        params::lavd_windows[0].allocate();
        params::nb_windows = 1;
        
        std::cout << "after initialization: 1st window=" << params::lavd_windows[0] << '\n';
    }
    else {
        // assign each window to its actual position on the overall region
        // raster since we have a prescribed resolution for that raster.
        value_t _dx = (params::region.max()[0] - params::region.min()[0])/(params::res[0]-1);
        value_t _dy = (params::region.max()[1] - params::region.min()[1])/(params::res[1]-1);
        for (size_t i=0; i<params::lavd_windows.size(); ++i) {
            lavd_window& window = params::lavd_windows[i];
            bbox_t& bounds = window.m_bounds;
            window.m_offset[0] = std::max(static_cast<long>((bounds.min()[0]-params::region.min()[0])/_dx), 
                                          static_cast<long>(0));
            window.m_offset[1] = std::max(static_cast<long>((bounds.min()[1]-params::region.min()[1])/_dy), 
                                          static_cast<long>(0));
            window.m_res[0] = static_cast<size_t>(bounds.size()[0]/_dx) + 1;
            window.m_res[1] = static_cast<size_t>(bounds.size()[1]/_dy) + 1;
            
            // adjust bounds as necessary
            bounds.min()[0] = params::region.min()[0] + window.m_offset[0]*_dx;
            bounds.min()[1] = params::region.min()[1] + window.m_offset[1]*_dx;
            bounds.max()[0] = bounds.min()[0] + (window.m_res[0]-1)*_dx;
            bounds.max()[1] = bounds.min()[1] + (window.m_res[1]-1)*_dy;
            
            // sanity check TBD
        }
    }
    
    _log_(1) << "region=\n" << params::region << std::endl;
    
    // storage for ODE states and trajectories
    params::trajectories.resize(params::nb_samples);
    params::states.resize(params::nb_samples);
    params::last_steps.resize(params::nb_samples);
    for (size_t i=0; i<params::states.size(); ++i) {
        params::states[i].reset(new ode_state(params::last_steps[i]));
    }
    params::last_vorticities.resize(params::nb_samples);
    
    params::support_radius = spurt::lavd::compute_support_radius();
    _log_(1) << "support radius = " << params::support_radius << std::endl;
    
    // initialize velocity and vorticity volumes
    import_data(velocity_filenames, vorticity_filenames, params::t_init/HOUR/3);
    
    size_t counter = 0;
    size_t nb_early = 0;
    
    params::step_x = params::region.size()[0] / static_cast<value_t>(params::res[0]-1);
    params::step_y = params::region.size()[1] / static_cast<value_t>(params::res[1]-1);
    
    _log_(1) << "bounds=" << params::region << std::endl;
    
    if (params::export_region) 
        spurt::lavd::export_mask(vorticity_filenames[0], params::region, true);
    
    std::vector<size_t> sample_counter(nb_threads);
    std::vector< shared_ptr< RasterODERHS > > rhs_copies(nb_threads);
    for (int i=0; i<nb_threads; ++i) {
        rhs_copies[i] = shared_ptr< RasterODERHS >(
            new RasterODERHS(*params::velocity_field, 
                             sample_counter[i],
                             params::domain, 
                             params::max_rhs_evals) );
    }
    
    params::start_time_str=sec2time(params::t_init);
    std::for_each(params::start_time_str.begin(), params::start_time_str.end(), 
                 [](char& c) { if (c==' ') c='_'; });
    _log_(1) << "start_time_string=" << params::start_time_str << std::endl;
    
    spurt::ProgressDisplay progress(false), step_progress(false), total_progress(false);
    
    total_progress.start(1); // we only care for the timer function
    total_progress.set_active(false);
    
    std::vector<size_t> nb_lost_copies(nb_threads);
    std::vector<size_t> counter_copies(nb_threads);
    
    bool initial_loop = true; 
    size_t niter = 0;
    params::current_start_time = params::t_init;
    params::next_export_time = params::t_init + params::t_export_step + params::t_skip;
    while ( true ) {
        step_progress.start(1);
        step_progress.set_active(false);
        
        _log_(0, "\n\n") << "next target time for this loop: " 
            << sec2time(params::next_export_time) 
            << " (" << params::next_export_time << " s.)" << std::endl;
        
        // counter = 0;
        std::fill(counter_copies.begin(), counter_copies.end(), 0);
        // nb_lost = 0;
        std::fill(nb_lost_copies.begin(), nb_lost_copies.end(), 0);
        size_t nb_early = 0;
        value_t current_target_time = std::min(params::current_t_max, 
                                               params::next_export_time);

        progress.start(params::nb_samples, "integrate flow map");
        progress.set_active(true);
        #pragma omp parallel
        {
            #pragma omp for schedule(dynamic,1)
            for (size_t n = 0 ; n < params::nb_samples ; ++n) {
                
                #if _OPENMP
                const int thread=omp_get_thread_num();
                #else
                const int thread=0;
                #endif
                
                if (!thread) progress.update(n);
                
                if (params::states[n]->m_stopped) {
                    ++nb_lost_copies[thread];
                    continue;
                }
                
                vec2 x0;
                if (initial_loop) {
                    size_t i = n % params::res[0];
                    size_t j = n / params::res[0];
                    x0 = params::region.min() + vec2(i*params::step_x, 
                                                     j*params::step_y);
#ifndef __SKIP_LOGS__
                    _log_(2) << "\n\n\n\nstarting integration #" << n << " at " << x0 << std::endl;
#endif
                }
                else {
                    x0 = params::trajectories[n].back();
                }
                params::last_steps[n].clear();
            
                ++counter_copies[thread];
        
                // create a stepper
                auto stepper = make_controlled(params::eps, params::eps, runge_kutta_dopri5<vec2>());
                                
                RasterODERHS& rhs=*rhs_copies[thread];
                rhs.m_counter=0;
                
                observer an_observer(*params::states[n]);
                
                if (initial_loop) {
#ifndef __SKIP_LOGS__
                    _log_(2) << "initializing observer at " << x0 << std::endl;
#endif
                    try {
                        params::trajectories[n].push_back(x0);
                        an_observer.initialize(x0, params::current_t_min);
                    }
                    catch( std::exception& e ) {
                        _log_(1) << "caught exception while attempting "
                            << "to initialize observer at " << x0 << std::endl;
                        params::states[n]->stop();
                        ++nb_lost_copies[thread];
                        continue;
                    }
                }
                            
                try { 
                    integrate_const(stepper, rhs, x0, params::current_start_time, 
                        current_target_time, params::dt, an_observer);
                    params::trajectories[n].push_back(crop(params::last_steps[n].back()));
                }
                catch(std::exception& e) {
                    _log_(1) << "caught exception while integrating from " 
                            << x0 << ":" << e.what() << '\n'
                            << "integration started at " << x0  
                            << " and ended at " << params::states[n]->m_pos.back() 
                            << '\n' << std::flush;
                    params::states[n]->stop();
                    ++nb_lost_copies[thread];
                    if (!params::last_steps[n].empty()) {
                        params::trajectories[n].push_back(crop(params::last_steps[n].back()));
                    }
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
        params::nb_lost = std::accumulate(nb_lost_copies.begin(), nb_lost_copies.end(), 0);

        progress.start(params::nb_samples, "sample vorticity");
        progress.set_active(true);
        #pragma omp parallel
        {
            #pragma omp for schedule(static,1)
            for (size_t n = 0; n < params::nb_samples; ++n) {
                #if _OPENMP
                const int thread=omp_get_thread_num();
                #else
                const int thread=0;
                #endif
                
                if (!thread) progress.update(n);
                
                const std::vector<vec3>& steps = params::last_steps[n];
                std::vector< value_t >& vorts = params::last_vorticities[n];
                vorts.clear();
                
                if (steps.empty()) continue;
                
                vorts.resize(steps.size());
                for (size_t i=0; i<steps.size(); ++i) {
                    vorts[i] = params::vorticity_field->value(steps[i]);
                }
            }
        }
        progress.end();
        
        progress.start(params::nb_samples, "update lavd");
        progress.set_active(true);
        #pragma omp parallel
        {
            #pragma omp for schedule(static,1)
            for (size_t n = 0; n < params::nb_samples; ++n) {
                #if _OPENMP
                const int thread=omp_get_thread_num();
                #else
                const int thread=0;
                #endif
                
                if (!thread) progress.update(n);
                
                const std::vector< vec3 >& steps = params::last_steps[n];
                const std::vector< value_t >& vorts = params::last_vorticities[n];
                if (vorts.empty()) continue;
                
                size_t i = n % params::res[0];
                size_t j = n / params::res[0];
                
                std::pair<bool, size_t> idx;
                for (size_t k=0; k<params::lavd_windows.size(); ++k) {
                    idx = params::lavd_windows[k].index(lvec2(i,j));
                    if (!idx.first) continue;
                    params::lavd_windows[k].update(idx.second, steps, vorts, *params::average_vorticity_fields[k]);
                }
            }
        }
        progress.end();
        
        step_progress.end();
        if (current_target_time == params::next_export_time) {
            ++niter;
            _log_(0) 
                << "Compute time this iteration: cpu: " << human_readable_duration(step_progress.cpu_time())
                    << " (" << step_progress.cpu_time() << " ms.) | wall: "
                << human_readable_duration(step_progress.wall_time()) 
                << " (" << step_progress.wall_time() << " ms.)\n"
                << "Overall compute time so far (" << niter << " iterations): cpu: ";
            double tcpu = total_progress.instant_cpu_time();
            double twal = total_progress.instant_wall_time();
             _log_(0)
                << human_readable_duration(tcpu) << " (" << tcpu << " ms.) | wall: "
                << human_readable_duration(twal) << " (" << twal << " ms.)"
                << std::endl;
        }
        
        if (current_target_time == params::t_max) break;
        else if (current_target_time == params::current_t_max) {
            // we have reached the end of the data currently available in
            // RAM. Import next time section.
            // [last_time_step-2*R, last_time_step-2*R+nslices-1]
            import_data(velocity_filenames, vorticity_filenames, 
                        params::last_time_step-2*params::support_radius);
            for (int i=0; i<nb_threads; ++i) {
                rhs_copies[i] = shared_ptr< RasterODERHS >(new RasterODERHS(*params::velocity_field, sample_counter[i], params::domain, params::max_rhs_evals));            
            }
        }
        if (current_target_time == params::next_export_time) {
            export_results(params::next_export_time,
                           progress.wall_time(),
                           progress.cpu_time(), false);
            params::next_export_time += params::t_export_step;
            if (params::next_export_time > params::t_max) params::next_export_time = params::t_max;
        }
        // advance starting time
        params::current_start_time = current_target_time;
        _log_(1) << "current start time=" << params::current_start_time << std::endl;
    }
    total_progress.end();
    export_results(params::t_max, total_progress.wall_time(),
                   total_progress.cpu_time(), true);

    return 0;
}
