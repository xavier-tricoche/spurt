#include <flow/lavd.hpp>

#include <sstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <cctype>
#include <regex>
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

template<typename T, size_t N>
std::ostream& operator<<(std::ostream& os, const std::array<T, N>& a) {
	os << "[" << a[0];
	for (int i=1; i<N; ++i) {
		os << ", " << a[i];
	}
	os << "]";
	return os;
}

namespace odeint = boost::numeric::odeint;

using namespace spurt::lavd;


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
std::string restart_name;
std::ofstream log_file;

spurt::log::dual_ostream spurt::lavd::_log_(log_file, std::cout, 1, 0, 0, true);

value_t t_max=2*spurt::lavd::DAY;
value_t t_init=0;
value_t t_skip=0;
value_t t_export_step=3*spurt::lavd::HOUR;
value_t t_between_files=3*spurt::lavd::HOUR;
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

void parse_restart_file() {
	_log_(2) << "parsing restart file" << std::endl;
	Nrrd* nin = spurt::nrrd_utils::readNrrd(restart_name);
	spurt::nrrd_utils::nrrd_traits traits(nin);
	res[0] = traits.sizes()[1];
	res[1] = traits.sizes()[2];

	std::cout << "resolution: " << res[0] << " x " << res[1] << std::endl;

	const std::vector<std::string>& comments = traits.comments();
	std::cout << "There are " << comments.size() << " lines of comments" << std::endl;

	std::string param, val_as_str;

	std::regex line_regex("[ ]*\\*[ ]*(.*)[ ]*=[ ]*(.*)$"); // "* <what>=<value>"
	std::regex res_regex("(.*)[ ]*x[ ]*(.*)$");
	std::regex bnd_regex("(.*)[ ]*->[ ]*(.*)$");
	std::smatch line_match, res_match, bnd_match;
	for (int i=0; i<comments.size(); ++i) { // skip first 2 lines of comments
		const std::string& c = comments[i];
		if (std::regex_search(c, line_match, line_regex)) {
			param = line_match[1].str();
			val_as_str = line_match[2].str();
		}
		else {
			std::cout << "skipping comment line: " << c << std::endl;
			continue;
		}
		std::cout << "parameter: " << param << ", value: " << val_as_str << std::endl;
		if (param == "input file") name_in = val_as_str;
		else if (param == "current time" || param == "tmax") t_init = std::stod(val_as_str);
		else if (param == "epsilon") eps = std::stod(val_as_str);
		else if (param == "resolution") {
			if (std::regex_search(val_as_str, res_match, res_regex)) {
				res[0] = std::stoi(res_match[1].str());
				res[1] = std::stoi(res_match[2].str());
			}
			else throw std::runtime_error("invalid resolution syntax in restart file");
		}
		else if (param == "bounds") {
			if (std::regex_search(val_as_str, bnd_match, bnd_regex)) {
				std::istringstream iss1(bnd_match[1].str());
				spurt::vec2 x;
				iss1 >> x;
				bnds[0] = x[0];
				bnds[1] = x[1];
				std::istringstream iss2(bnd_match[2].str());
				iss2 >> x;
				bnds[2] = x[0];
				bnds[3] = x[1];
			}
			else throw std::runtime_error("invalid bounds syntax in restart file");
		}
		else if (param == "dt") dt = std::stod(val_as_str);
		else {
			std::cout << "skipping value assignment for " << param << std::endl;
			continue;
		}
	}

	std::cout << "imported parameters:\n";
	std::cout << "name_in=" << name_in << '\n';
	std::cout << "t_init=" << t_init << '\n';
	std::cout << "epsilon=" << eps << '\n';
	std::cout << "resolution=" << res << '\n';
	std::cout << "bounds=" << bnds << '\n';
	std::cout << "dt=" << dt << '\n';

	nrrdNuke(nin); // we will reopen the file later to access particle coordinates
}

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;

    xcl::option_traits
            required(true, false, "Required Options"),
            optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Compute Lagrangian averaged vorticity deviation");

    verbose_vec.resize(2);
    verbose_vec[0] = verbose_cout;
    verbose_vec[1] = verbose_log;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input info filename", optional);
        parser.add_value("output", name_out, "Output file basename ", required);
		parser.add_value("restart", restart_name, "Restart filename", optional);
        parser.add_sequence("tinit", t_init_str, t_init_str, "Integration starting time, expressed either in time since initial time step (in seconds; for minutes append \"m\", for hours \"h\", for days \"d\", for weeks \"w\") or as a date MM/DD", optional, "<float>[<char>]");
        parser.add_value("tmax", t_max_str, "Integration ending time (in seconds; for minutes append \"m\", for hours \"h\", for days \"d\", for weeks \"w\")", optional, "<float>[<char>]");
        parser.add_value("step", t_expstp_str, t_expstp_str, "Time interval between intermediate result outputs", optional);
        parser.add_value("skip", t_skip_str, t_skip_str, "Time skipped before first output", optional);
        parser.add_value("dt", dt_str, dt_str, "Time sampling of integral curve for LAVD computation (in seconds, for minutes add \"m\", for hours \"h\", for times \"x\")", optional, "<float>[<char>]");
        parser.add_tuple<2>("res", res, res, "Sampling resolution", optional);
        parser.add_tuple<4>("bounds", bnds, bnds, "Sampling bounds (min longitude and latitude followed by max's)", optional);
        parser.add_value("nslices", nslices, nslices, "Number of time steps to process simultaneously", optional);
        parser.add_value("eps", eps, eps, "Integration precision", optional);
        parser.add_value("mask", border_mask_name, border_mask_name, "Mask filename", optional);
        parser.add_value("max", max_rhs_evals, max_rhs_evals, "Maximum authorized number of rhs evaluations during intermediate steps of the integration");
        parser.add_value("export", export_region, export_region, "Export index range of selected region", optional);
        parser.add_value("traj", export_trajectories, export_trajectories, "Export trajectories", optional);
        parser.add_value("log", log_name, log_name, "Log export filename", optional);
        parser.add_sequence("verbose", verbose_vec, verbose_vec, "Verbosity level(s) on stdout and log file", optional, "<uint> [<uint>]");
        parser.add_value("lname", long_name, long_name, "Export result with long name indicating computation parameters", optional);

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
		if (!restart_name.empty()) { // input file will override existing values
			_log_(1) << "Restart file provided: " << restart_name << std::endl;
			parse_restart_file();
		}
		if (name_in.empty()) {
			std::cerr << "ERROR: " << argv[0] << ": Neither input file nor restart file have been provided.\n";
			exit(1);
		}
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n"
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

struct observer
{
    LAVD_state& m_state;
    std::vector< vec2 >& m_pos;
    NrrdScalarField<1>& m_avg_vort;
    NrrdScalarField<3>& m_vort;

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

    observer(LAVD_state& state, std::vector<vec2>& pos,
             NrrdScalarField<3>& vorticity, NrrdScalarField<1>& average_vorticity)
        : m_state(state), m_pos(pos), m_vort(vorticity),
          m_avg_vort(average_vorticity) {}

    void initialize(const vec2& p0, value_t t0) {
        value_t vd0 = vorticity_deviation(p0, t0);
        m_state = LAVD_state(p0, vd0, t0);
        m_pos.clear();
        m_pos.push_back(p0);
    }

    void operator()(const vec2 &x, value_t t)
    {
        // save latest position: if trajectories are not exported
        // we only retain the latest position for flow map
        // computation purposes
        if (export_trajectories) m_pos.push_back(x);
        else m_pos.back() = x;

        value_t vd = vorticity_deviation(x, t);
        m_state.update(x, vd, t);
        /*_log_(1) << "LAVD(" << t << ")=" << m_state.evaluate()
            << " (vd=" << vd << ")" << std::endl;*/
    }
};

std::vector< std::vector< vec2 > > all_trajectories;
std::vector< LAVD_state > all_states;


void import_data(const std::vector< std::string >& vel_filenames,
                 const std::vector< std::string >& vor_filenames,
                 size_t start_id)
{
    // for (size_t i=0; i<vel_filenames.size(); ++i) {
    //     std::cout << "filename[" << i << "]=" << vel_filenames[i] << '\n';
    // }

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
        vel_tsteps[i-first]=spurt::nrrd_utils::readNrrd(vel_filenames[i]);
        _log_(1) << "done" << std::endl;
        print(vel_tsteps[i-first]);

        _log_(1) << "Reading " << vor_filenames[i] << "... " << std::flush;
        vor_tsteps[i-first]=spurt::nrrd_utils::readNrrd(vor_filenames[i]);
		_log_(1) << "printing contents of vorticity input Nrrd #" << i << "\n";
		print(vor_tsteps[i-first]);
        _log_(1) << "done" << std::endl;
    }
    last_time_step = first + n_used - 1;
    value_t outer_tmin = /*t_init +*/ first*t_between_files;
    _log_(2) << "outer_tmin=" << outer_tmin << std::endl;
    current_velocity_volume = spurt::lavd::create_nrrd_volume(vel_tsteps, outer_tmin, t_between_files);
    _log_(2) << "after velocity volume creation:" << std::endl;
    print(current_velocity_volume);

    // clean up
    for (size_t i=0; i<vel_tsteps.size(); ++i) {
        nrrdNuke(vel_tsteps[i]);
    }

    current_vorticity_volume = spurt::lavd::create_nrrd_volume(vor_tsteps, outer_tmin, t_between_files);
    current_t_min = start_id*t_between_files;
    // if (outer_tmin > t_init) {
    //     current_t_min += support_radius*t_between_files;
    // }
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

    // compute average vorticity for selected region
    value_t* av_array = (value_t *)calloc(n_used, sizeof(value_t));
    value_t region_size = region.size()[0]*region.size()[1];
    _log_(1) << "sampling vorticity slices" << std::endl;

    size_t nx = vor_tsteps[0]->axis[0].size;
    size_t ny = vor_tsteps[0]->axis[1].size;
    value_t minx = domain.min()[0];
    value_t miny = domain.min()[1];
    value_t dx = input_spc[0];
    value_t dy = input_spc[1];

    _log_(1) << "nx=" << nx << ", ny=" << ny << '\n';
    _log_(1) << "minx=" << minx << ", mimy=" << miny << '\n';
    _log_(1) << "dx=" << dx << ", dy=" << dy << '\n';

    double min_offset_x = region.min()[0]-minx;
    double min_offset_y = region.min()[1]-miny;
    double max_offset_x = region.max()[0]-minx;
    double max_offset_y = region.max()[1]-miny;
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
    for (size_t n=0; n<n_used; ++n, base_n+=nsamples) {
        size_t nvalid_samples=0;
        value_t sum=0;
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
                    value_t vort=nrrd_value(vor_tsteps[n], id);
                    #pragma omp atomic
                    sum += vort;
                    #pragma omp atomic
                    ++nvalid_samples;
                }
            }
        }
        av_array[n] = sum/static_cast<value_t>(nvalid_samples);
        _log_(1) << std::setprecision(12) << "average vorticity[" << n << "]=" << av_array[n] << std::endl;
        _log_(1) << "nsamples=" << nsamples << ", nvalid_samples=" << nvalid_samples << std::endl;
    }
	progress.end();

    Nrrd* av=nrrdNew();
    if (nrrdWrap_nva(av, av_array, spurt::nrrd_utils::nrrd_value_traits_from_type<value_t>::index,
                     1, &n_used)) {
        throw std::runtime_error(spurt::nrrd_utils::error_msg("unable to create 1D average vorticity nrrd"));
    }
    // Per-axis info
    int center = nrrdCenterNode;
    nrrdAxisInfoSet_nva(av, nrrdAxisInfoSpacing, &t_between_files);
    nrrdAxisInfoSet_nva(av, nrrdAxisInfoMin, &outer_tmin);
    nrrdAxisInfoSet_nva(av, nrrdAxisInfoCenter, &center);
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
    if (long_name) {
        os << name_out << "_started_at_" << start_time_str << "_flowmap+lavd_" << std::setw(5) << std::setfill('0') << current_time/spurt::lavd::HOUR << "h_";
        os << res[0] << "x" << res[1] << "_" << std::setprecision(2)
           << std::scientific << eps;
        os.flags(default_settings);
    }
    else {
        os << name_out << "_fmap+lavd_" << std::setw(5) << std::setfill('0') << current_time/spurt::lavd::HOUR << "h";
    }
    std::string basename = os.str();

	std::vector<std::string> comments;
    os.clear(); os.str("");
    os << "This file was produced by " << me
       << ". 1st axis corresponds to (flowmap_x, flowmap_y, lavd), where "
       << "lavd stands for Lagrangian averaged vorticity deviation.\n";
	comments.push_back(os.str()); os.clear(); os.str("");
    os << "Computation parameters:\n";
	comments.push_back(os.str()); os.clear(); os.str("");
    os << " * input file=" << name_in << '\n';
	comments.push_back(os.str()); os.clear(); os.str("");
	os << " * current time=" << current_time << "\n";
	comments.push_back(os.str()); os.clear(); os.str("");
    os << " * tmax=" << current_time << '\n';
	comments.push_back(os.str()); os.clear(); os.str("");
    os << " * epsilon=" << eps << '\n';
	comments.push_back(os.str()); os.clear(); os.str("");
    os << " * resolution=" << res[0] << "x" << res[1] << '\n';
	comments.push_back(os.str()); os.clear(); os.str("");
    os << " * bounds=" << region.min() << "->" << region.max() << '\n';
	comments.push_back(os.str()); os.clear(); os.str("");
    os << " * dt=" << dt << '\n';
	comments.push_back(os.str()); os.clear(); os.str("");
    os << " * kernel size=" << support_radius << '\n';
	comments.push_back(os.str());

    spurt::lavd::export_results< float >(
        current_time,
        wall_time, cpu_time,
        final,
        nb_lost,
        res[0], res[1],
        step_x, step_y,
        export_trajectories,
        region,
        basename,
        comments,
        all_trajectories,
        all_states
    );
}

int main(int argc, const char* argv[])
{
    using namespace spurt;
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

    name_out=spurt::filename::remove_extension(name_out);

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
    get_spatial_info(domain, input_spc, velocity_filenames[t_init/HOUR/3], 1);
    border_mask = get_border_mask(border_mask_name, velocity_filenames[t_init/HOUR/3]);

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
    all_trajectories.resize(nb_samples);
    all_states.resize(nb_samples);

	// if restarting from last output, initialize these arrays with available information

	std::vector<double> restart_values;
	if (!restart_name.empty()) {
		Nrrd* nin = spurt::nrrd_utils::readNrrd(restart_name);
		spurt::nrrd_utils::to_vector(restart_values, nin);
		for (size_t i=0; i<all_trajectories.size(); ++i) {
			all_trajectories[i].push_back(spurt::vec2(restart_values[3*i], restart_values[3*i+1]));
		}
	}

    support_radius = spurt::lavd::compute_support_radius();
    _log_(1) << "support radius = " << support_radius << std::endl;

    // initialize velocity and vorticity volumes
    import_data(velocity_filenames, vorticity_filenames, t_init/HOUR/3);

    size_t counter = 0;
    size_t nb_early = 0;

    step_x = region.size()[0] / static_cast<value_t>(res[0]-1);
    step_y = region.size()[1] / static_cast<value_t>(res[1]-1);

    _log_(1) << "bounds=" << region << std::endl;

    if (export_region)
        spurt::lavd::export_mask(vorticity_filenames[t_init/HOUR/3], region, true);

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

    spurt::ProgressDisplay progress(false), total_progress(false);

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

                if (all_states[n].m_stopped) {
                    #pragma omp atomic
                    ++nb_lost;
                    continue;
                }

                vec2 x0;
                if (initial_loop && restart_name.empty()) {
                    size_t i = n % res[0];
                    size_t j = n / res[0];
                    x0 = region.min() + vec2(i*step_x, j*step_y);
                    _log_(2) << "\n\n\n\nstarting integration #" << n << " at " << x0 << std::endl;
                }
                else {
                    x0 = all_trajectories[n].back();
					if (initial_loop && !restart_name.empty()) {
						_log_(2) << "\n\n\n\nresuming integration #" << n << " at " << x0 << std::endl;
					}
                }

                #pragma omp atomic
                ++counter;

                // create a stepper
                auto stepper = make_controlled(eps, eps, runge_kutta_dopri5<vec2>());

                NrrdODERHS& rhs=*rhs_copies[thread];
                rhs.m_counter=0;

                NrrdScalarField<1> my_avg_vorticityf(current_average_vorticity_vector, std::string("average vorticity field for thread #") + std::to_string(thread));
                NrrdScalarField<3> my_vorticityf(current_vorticity_volume, std::string("vorticity field for thread #" + std::to_string(thread)));
                observer an_observer(all_states[n], all_trajectories[n],
                                     my_vorticityf, my_avg_vorticityf);

                if (initial_loop) {
                    _log_(2) << "initializing observer at " << x0 << std::endl;
                    try {
                        an_observer.initialize(x0, current_t_min);
						if (!restart_name.empty()) {
							all_states[n].m_acc_vd = restart_values[3*n+2];
						}
                    }
                    catch( std::exception& e ) {
                        _log_(1) << "caught exception while attempting "
                            << "to initialize observer at " << x0 << std::endl;
                        all_states[n].stop();
                        if (restart_name.empty()) all_trajectories[n].push_back(x0); // otherwise trajectory already contains seed point
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
                            << " and ended at " << all_trajectories[n].back()
                            << '\n' << std::flush;
                    all_states[n].stop();
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
                << " (" << progress.wall_time() << " ms.)\n";
			_log_(0) << "Interrupted trajectories: " << nb_lost << " (" << nb_lost*100/nb_samples << "%)" << std::endl;
            _log_(0) << "Overall compute time so far (" << niter << " iterations): cpu: ";

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
