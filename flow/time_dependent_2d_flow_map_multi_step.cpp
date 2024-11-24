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

#include <math/types.hpp>
#include <flow/time_dependent_field.hpp>
#include <flow/vector_field.hpp>  //BARG
#include <flow/ode_observer.hpp>
#include <data/image.hpp>

#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <misc/strings.hpp>
#include <vtk/vtk_interpolator.hpp>
#include <tbb/parallel_for.h>

#include <mutex>
#include <atomic>

using namespace spurt;

std::atomic<size_t> progress_counter;
std::string name_in, name_out, seed_name, path, cmdline;
double length, eps, eps_refined, T, t0, t1, deltaT;
size_t dim, mem;
std::array<size_t, 2> res;
std::array<bool, 2> periodic;
std::array<double, 4> bounds;
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
    res = { 64, 64 };
    periodic = { false, false };
    bounds = { 1, -1, 1, -1}; // invalid bounds
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
        parser.add_tuple<2>("res", res, res, "Sampling resolution", optional);
        parser.add_tuple<2>("periodic", periodic, periodic, "Periodic boundary conditions", optional);
        parser.add_tuple<4>("bounds", bounds, "Sampling bounds", optional);
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

using namespace spurt;
using namespace vtk_utils;

typedef vec3 point_t;
typedef vec3 vec_t;
typedef double scalar_type;

typedef interpolator<vtkImageData, double, 2, vec_t>  img_intp_t;
typedef interpolator<vtkRectilinearGrid, double, 2, vec_t> rect_intp_t;
typedef interpolator<vtkUnstructuredGrid, double, 2, vec_t> unst_intp_t;
typedef interpolator<vtkStructuredGrid, double, 2, vec_t> curv_intp_t;
typedef interpolator<boundaryAwareRectGrid, double, 2, vec_t> BARG_intp_t; //BARG edit #1


typedef vtk_utils::point_locator<vtkUnstructuredGrid, double, 2> unstr_locator_type;
typedef spurt::fixed_mesh_time_dependent_field<unstr_locator_type, std::vector<vec_t> > unstr_field_type;
typedef vtk_utils::point_locator<vtkStructuredGrid, double, 2> curv_locator_type;
typedef spurt::fixed_mesh_time_dependent_field<curv_locator_type, std::vector<vec_t> > curv_field_type;
typedef spurt::structured_mesh_time_dependent_field< vtkRectilinearGrid, double, 2, std::vector<vec_t> > rect_field_type;
typedef spurt::structured_mesh_time_dependent_field< vtkImageData, double, 2, std::vector<vec_t> > img_field_type;
typedef spurt::tp_bspline_time_dependent_field< boundaryAwareRectGrid, double, 2, std::vector<vec_t> > BARG_field_type;  //BARG edit


template<typename Field, typename Enable = void>
struct rhs_type {};

template<typename Field>
struct rhs_type<
	Field,
	typename std::enable_if<
	std::is_same<typename Field::dataset_type, vtkStructuredGrid >::value ||
	std::is_same<typename Field::dataset_type, vtkUnstructuredGrid>::value ||
	std::is_same<typename Field::dataset_type, vtkRectilinearGrid>::value ||
	std::is_same<typename Field::dataset_type, vtkImageData>::value ||
	std::is_same<typename Field::dataset_type, boundaryAwareRectGrid>::value
	>::type > {
	typedef Field field_type;

	rhs_type(std::shared_ptr<const field_type> field, size_t max_evals = 1000,
		bool _verbose = false, int seedid = -1)
		: m_field(field), m_counter(0), m_max_evals(max_evals), m_verbose(_verbose),
		m_seedid(seedid), m_cell(VTK_SMART(vtkGenericCell)::New()) {}

	rhs_type(std::shared_ptr<const field_type> field, vtkGenericCell* acell, size_t max_evals = 1000, bool forward = true, bool _verbose = false, int seedid = -1)
		: m_field(field), m_cell(acell), m_counter(0), m_max_evals(max_evals),
		m_seedid(seedid), m_verbose(_verbose) {}

	rhs_type(const rhs_type& other)
		: m_field(other.m_field), m_counter(other.m_counter),
		m_cell(other.m_cell), m_max_evals(other.m_max_evals),
		m_verbose(other.m_verbose), m_seedid(other.m_seedid) {}

	void operator()(const vec3& x, vec3& dxdt, scalar_type t) const {
		if (confine && !the_bounds.inside(x)) {
			if (verbose && do_log) {
				std::ostringstream os;
				os << "invalid position: " << (x) << " at t=" << t << '\n';
				write_to_ostream(ofs, os.str());
			}
			throw invalid_position_exception("invalid position: " + to_str(x) + " at t=" + to_string(t));
		}
		std::ostringstream os;
		if (verbose && do_log) {
			os << "rhs(" << (x) << ", " << t << ")=";
		}
		dxdt = (*m_field)(x, t); // may throw
		if (verbose && do_log) {
			os << to_str(dxdt);
		}
		if (verbose && do_log) {
			write_to_ostream(ofs, os.str());
		}
		++m_counter;
		if (m_counter >= m_max_evals) {
			std::ostringstream os;
			os << "Max RHS evaluations reached at "
				<< (x) << " at t=" << t << '\n';
			write_to_ostream(ofs, os.str());
			// will trigger a invalid_position_exception:
			throw std::runtime_error(os.str());
		}
	}

	std::shared_ptr<const field_type> m_field;
	mutable size_t m_counter;
	const size_t m_max_evals;
	VTK_SMART(vtkGenericCell) m_cell;
	bool m_verbose;
	int m_seedid;
};

typedef Observer<vec_t> observer_t;

const vtkDataSet* get_dataset(const std::shared_ptr<unstr_field_type> field) {
	return field->get_locator()->get_dataset();
}

const vtkDataSet* get_dataset(const std::shared_ptr<curv_field_type> field) {
	return field->get_locator()->get_dataset();
}

const vtkDataSet* get_dataset(const std::shared_ptr<rect_field_type> field) {
	return field->get_interpolator()->get_dataset();
}

const vtkDataSet* get_dataset(const std::shared_ptr<img_field_type> field) {
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

	// serial portion of the algorithm
	std::vector<vec3> seeds;
	typedef spurt::raster_grid<long, double, 3, lvec3, vec3> grid_t;

	if (::bounds[0] < ::bounds[1] &&
		::bounds[2] < ::bounds[3] &&
		::bounds[0] >= global_bounds[0] &&
		::bounds[1] <= global_bounds[1] &&
		::bounds[2] >= global_bounds[2] &&
		::bounds[3] <= global_bounds[3]) {
		// valid bounds supplied by user
		bnds.min() = vec3(::bounds[0], ::bounds[2], 0);
		bnds.max() = vec3(::bounds[1], ::bounds[3], 0);
	}
	else {
		bnds.min() = vec3(global_bounds[0], global_bounds[2], 0);
		bnds.max() = vec3(global_bounds[1], global_bounds[3], 0);
	}

	the_bounds = bnds;
	int npoints;
	grid_t sampling_grid(as<grid_t::coord_type>(res), bnds);

	if (seed_name.empty()) {
		std::cout << "Resolution = " << res[0] << " x " << res[1] << std::endl;

		std::cout << "sampling grid bounds are: " 
				  << sampling_grid.bounds().min()
			      << " -> " 
				  << sampling_grid.bounds().max() << '\n';

		if (do_log) {
			ofs << "Resolution = " << res[0] << " x " << res[1] << std::endl;
			ofs << "sampling grid bounds are: " << sampling_grid.bounds().min()
				<< " -> " << sampling_grid.bounds().max() << '\n';
		}
		npoints = sampling_grid.size();
		seeds.resize(npoints);
		for (int i = 0; i < npoints; i++) {
			auto p = sampling_grid(sampling_grid.coordinates(i));
			seeds[i] = vec3(p[0], p[1], 0);
		}
	}
	else {
		std::string ext = spurt::filename::extension(seed_name);
		if (ext == "vtk" || ext == "VTK" || ext == "vtp" || ext == "VTP") {
			VTK_SMART(vtkDataSet) dataset = vtk_utils::readVTK(seed_name);
			npoints = dataset->GetNumberOfPoints();
			seeds.resize(npoints);
			for (int i = 0; i < npoints; i++) {
				double pt[3];
				dataset->GetPoint(i, (double*)(pt));
				seeds[i] = vec3(pt[0], pt[1], 0);
			}
		}
		else if (ext == "txt" || ext == "TXT") {
			std::ifstream ifs(seed_name);
			ifs >> npoints;
			seeds.resize(npoints);
			for (int i = 0; i < npoints; ++i) {
				ifs >> seeds[i][0] >> seeds[i][1];
				seeds[i][2] = 0;
			}
			ifs.close();
		}
	}

	float* flowmap = (float*)calloc(2 * npoints, sizeof(float));
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
	// initialization of I.C.
	tbb::parallel_for(tbb::blocked_range<int>(0, npoints),
		[&](tbb::blocked_range<int> r) {
		for (int n = r.begin(); n != r.end(); ++n) {
		flowmap[2 * n] = seeds[n][0];
		flowmap[2 * n + 1] = seeds[n][1];

		++progress_counter;
		}
	});
	// parallel portion of the algorithm
	std::fill(flowtimes, flowtimes + npoints, t0);

	size_t nbcomputed = 0;
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
	if (t1 == 0) t1 = deltaT;
	else {
		t1 *= sign(T);
	}

	int n_rounds = 0;
	if (fabs(t1) >= fabs(T)) n_rounds = 1;
	else if (fabs(T) - fabs(t1) <= deltaT) n_rounds = 2;
	else n_rounds = 1 + std::ceil((fabs(T) - fabs(t1)) / fabs(deltaT));

	std::cout << "n_rounds = " << n_rounds << '\n';
	if (do_log) {
		ofs << "n_rounds = " << n_rounds << '\n';
	}

	std::cout << "verbose = " << (verbose ? "true" : "false") << '\n';

	std::vector<bool> failed(npoints, false);
	std::vector<bool> stopped(npoints, false);
	std::vector< std::vector<vec3> > lines(npoints);
	std::vector< std::vector<double> > times(npoints);

	for (int round = 0; round < n_rounds; ++round) {
		double next_time = current_time + (!round ? t1 : deltaT);
		if (T > 0 && next_time > t0 + T) {
			next_time = t0 + T;
		}
		else if (T < 0 && next_time < t0 + T) {
			next_time = t0 + T;
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
		tbb::parallel_for(tbb::blocked_range<int>(0, npoints),
			[&](tbb::blocked_range<int> r) {
			for (int n = r.begin(); n != r.end(); ++n) {
			++progress_counter;

			update_progress(progress);

			bool super_verbose = !seed_name.empty();
			point_t x = { flowmap[2 * n], flowmap[2 * n + 1], 0 };
			if (stopped[n]) {
				if (do_log && verbose) {
					std::ostringstream oss;
					oss << "pathline " << n << ": already stopped at " << (x);
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
			double local_eps = eps;
			point_t p;
			try {
				while (local_eps >= eps_refined) {
					point_t y(x);
					t = flowtimes[n]; // obs.last_t
					d = 0; // obs.dist
					p = x; // obs.last_p
					sline.clear();
					slinetimes.clear();
					using namespace boost::numeric::odeint;
					observer_t obs(p, t, d, sline, slinetimes, do_verbose);
					auto stepper = make_controlled(local_eps, local_eps, runge_kutta_dopri5<point_t>());
					integrate_adaptive(stepper, rhs, y, t,
									   next_time - flowtimes[n], 1.0e-4, obs);
					if (true || fabs(t - next_time) / fabs(next_time - current_time) < 0.001) {
						break;
					}
					local_eps /= 2.;
				}
				flowmap[2 * n] = p[0];
				flowmap[2 * n + 1] = p[1];
				flowtimes[n] = t;
				if (do_log && verbose) {
					std::ostringstream os;
					os << n << ": successful: advected from " << x << " to "
						<< p << ", distance: "
						<< d << ", final time: " << t << ")\n";
					os << "line contains " << sline.size() << " (" << lines[n].size() << ") points" << '\n';
					if (do_log) write_to_ostream(ofs, os.str());
				}
				if (fabs((flowtimes[n] - next_time) / (next_time - current_time)) > 0.001) {
					if (do_log) {
						std::ostringstream os;
						os << "This pathline (id: " << n << ", seed: "
							<< sampling_grid(sampling_grid.coordinates(n))
							<< ") stopped prematurily. Last recorded time: "
							<< t << " (" << fabs(next_time - flowtimes[n])
							<< ", "
							<< 100.*fabs((flowtimes[n] - next_time) / (next_time - current_time))
							<< "%) for a distance of " << d << '\n';
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
						<< to_str(x) << ". reason: " << e.what() << '\n'
						<< "last position reached was "
						<< to_str(p) << " at time "
						<< t << " for a total distance of "
						<< d;
					write_to_ostream(ofs, os.str());
				}
				flowtimes[n] = t;
				flowmap[2 * n] = p[0];
				flowmap[2 * n + 1] = p[1];
				stopped[n] = true;
			}
		}
		}); // end of tbb::parallel_for

		progress.end();
		if (do_log) {
			ofs << progress << '\n';
		}

		if (save_lines || monitor) {
			if (verbose) {
				std::cout << "before filtering, lines contains " 
							<< lines.size()
							<< " and their lengths are ";
				for (int i = 0; i < lines.size(); ++i)
					std::cout << lines[i].size() << " ";
				std::cout << "\n";
			}
			if (true || monitor) {
				std::vector< std::vector<vec3> > newlines;
				std::vector< std::vector<double> > newtimes;
				for (size_t i = 0; i < lines.size(); i++) {
					if (lines[i].size() > 2) {
						newlines.push_back(lines[i]);
						newtimes.push_back(times[i]);
					}
					else {
						std::cout << "rejected line of length " << lines[i].size() << '\n';
					}
				}
				lines.swap(newlines);
				times.swap(newtimes);
			}
			std::cout << "after filtering lines contains " << lines.size() << " lines.\n";
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
				size_t offset = 0;
				for (int i = 0; i < lines.size(); ++i) {
					vertices->InsertNextCell(1);
					vertices->InsertCellPoint(offset + lines[i].size() - 1);
					offset += lines[i].size();
				}
				pdata->SetVerts(vertices);
			}
			{
				std::ostringstream os;
				if (save_lines)
					os << name_out << (T > 0 ? "-fwd_" : "-bwd_") << "pathlines_from_" << current_time << "_to_" << next_time << "_round_" << round << "_of_" << n_rounds << ".vtp";
				else {
					os << name_out << (T > 0 ? "-fwd_" : "-bwd_") << "interrupted_pathlines_from_" << current_time << "_to_" << next_time << "_round_" << round << "_of_" << n_rounds << ".vtp";
				}
				vtk_utils::saveVTK(pdata, os.str());
			}
		}

		if (seed_name.empty()) {
			std::vector<size_t> size(3);
			std::vector<double> step(3), mins(3);
			step[0] = std::numeric_limits<double>::quiet_NaN();
			mins[0] = std::numeric_limits<double>::quiet_NaN();
			std::vector<int> ctrs(3);
			std::fill(ctrs.begin(), ctrs.end(), nrrdCenterNode);
			size[0] = 2;
			for (int i = 0; i < 2; ++i) {
				size[i + 1] = res[i];
				step[i + 1] = sampling_grid.spacing()[i];
				mins[i + 1] = sampling_grid.bounds().min()[i];
				mins[i + 1] = sampling_grid.bounds().min()[i];
			}

			std::ostringstream os;
			os << name_out << (T > 0 ? "-fwd_" : "-bwd_") << "flowmap-deltaT=" << fabs(next_time - t0) << "_t0=" << t0 << ".nrrd";
			spurt::nrrd_utils::writeNrrdFromContainers(flowmap, os.str(), size, step, mins, ctrs, cmdline);
			if (do_log) {
				ofs << "exported: " << os.str() << '\n';
			}

			size[0] = 1;
			os.clear();
			os.str("");
			os << name_out << (T > 0 ? "-fwd_" : "-bwd_") << "flowtime-deltaT=" << fabs(next_time - t0) << "_t0=" << t0 << ".nrrd";
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
			for (int i = 0; i < failed.size(); ++i) {
				if (failed[i]) {
					os << "id: " << i << ", pos: " << sampling_grid(sampling_grid.coordinates(i)) << '\n';
				}
			}
			write_to_ostream(ofs, os.str());
		}
		if (monitor) {
			std::vector<vec3> seeds;
			std::vector<int> seed_ids;
			for (int i = 0; i < failed.size(); ++i) {
				if (failed[i]) {
					auto p = sampling_grid(sampling_grid.coordinates(i));
					seeds.push_back(vec3(p[0], p[1], 0));
					seed_ids.push_back(i);
				}
			}
			if (do_log) {
				std::ostringstream os;
				os << "There were " << seeds.size() << " failed integrations overall\n";
				write_to_ostream(ofs, os.str());
			}
			lines.resize(seeds.size());
			times.resize(seeds.size());
			progress_counter = 0;

			tbb::parallel_for(tbb::blocked_range<int>(0, seeds.size()),
				[&](tbb::blocked_range<int> r) {
			for (int n = r.begin(); n != r.end(); ++n) 
			{
				++progress_counter;

				double t = t0, d = 0;
				const point_t& x = seeds[n];
				point_t p(x);

				std::vector<vec3>& sline = lines[n];
				std::vector<double>& slinetimes = times[n];

				bool do_verbose = verbose;

				rhs_type<Field> rhs(field, 1000, do_verbose, n);
				observer_t obs(p, t, d, sline, slinetimes, do_verbose);

				try {
					point_t y(x);
					// attempt integration over the full time interval, knowing that it will
					// stop somewhere along the way
					using namespace boost::numeric::odeint;
					observer_t obs(p, t, d, sline, slinetimes, do_verbose);
					auto stepper = make_controlled(eps_refined, eps_refined, runge_kutta_dopri5<vec3>());
					integrate_adaptive(stepper, rhs, y, t0, T, 1.0e-4, obs);
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
							<< to_str(seeds[n])
							<< "). reason: " << e.what() << '\n'
							<< "last position reached was "
							<< to_str(obs.last_p) << " at time "
							<< obs.last_t << " for a total distance of "
							<< obs.distance;
						write_to_ostream(ofs, os.str());
					}
				}
			} // for loop
			}); // tbb parallel_for

			VTK_SMART(vtkPolyData) failed_lines =
				vtk_utils::make_polylines(lines, 0, false);
			VTK_CREATE(vtkFloatArray, line_times);
			line_times->SetNumberOfComponents(1);
			line_times->SetNumberOfTuples(failed_lines->GetNumberOfPoints());
			vtkIdType _count = 0;
			for (int n = 0; n < times.size(); n++) {
				for (int i = 0; i < times[n].size(); ++i, ++_count) {
					line_times->InsertTuple(_count, &(times[n][i]));
				}
			}
			failed_lines->GetPointData()->SetScalars(line_times);
			std::ostringstream os;
			os << name_out << (T > 0 ? "-fwd_" : "-bwd_") << "failed_pathlines_from_" << t0 << "_to_" << t0 + T << ".vtp";
			vtk_utils::saveVTK(failed_lines, os.str());
		} // end of monitor section
	} // end of failed processing

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
	import_data(std::vector< std::shared_ptr<std::vector<vec_t> > >& data,
		std::vector<bool>& out_of_bounds,
		const std::vector<std::string>& steps) {
	data.resize(steps.size());
	VTK_SMART(Dataset) the_grid;
	for (int i = 0; i < steps.size(); i++) {
		VTK_SMART(Dataset) a_grid = Dataset::SafeDownCast(vtk_utils::readVTK(steps[i]));
		if (!a_grid) throw std::runtime_error("Invalid grid type");
		if (!i) the_grid = a_grid;
		data[i] = std::shared_ptr< std::vector< vec_t > >(new std::vector< vec_t >());
		data[i]->resize(the_grid->GetNumberOfPoints());
		std::cout << "grid contains " << the_grid->GetNumberOfPoints() << " points\n";
		the_grid->PrintSelf(std::cout, vtkIndent(0));
		for (size_t j = 0; j < data[i]->size(); j++) {
			double vec[3];
			the_grid->GetPointData()->GetVectors()->GetTuple(j, vec);
			(*data[i])[j] = vec3(vec[0], vec[1], 0);
		}
		if (!i) {
			out_of_bounds.resize(the_grid->GetNumberOfCells());
			VTK_CREATE(vtkIdList, ptsids);
			for (int j = 0; j < out_of_bounds.size(); ++j) {
				the_grid->GetCellPoints(j, ptsids);
				bool is_zero = true;
				for (int n = 0; n < ptsids->GetNumberOfIds() && is_zero; ++n) {
					const vec3& v = (*data[i])[ptsids->GetId(n)];
					if (v[0] != 0 || v[1] != 0) {
						is_zero = false;
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

shared_ptr<rect_field_type>
	load_rectilinear_VTK_time_steps(const std::vector<std::string>& steps,
		const std::vector<double>& times) {
	assert(!steps.empty());
	std::vector< std::shared_ptr<std::vector<vec_t> > > data;
	std::vector<bool> out_of_bounds;
	VTK_SMART(vtkRectilinearGrid) grid =
		import_data<vtkRectilinearGrid>(data, out_of_bounds, steps);

	return shared_ptr<rect_field_type>(
		new rect_field_type(grid, data, times, out_of_bounds, verbose));
}

shared_ptr<img_field_type>
	load_image_VTK_time_steps(const std::vector<std::string>& steps,
		const std::vector<double>& times) {
	assert(!steps.empty());
	std::vector< std::shared_ptr<std::vector<vec_t> > > data;
	std::vector<bool> out_of_bounds;
	VTK_SMART(vtkImageData) grid =
		import_data<vtkImageData>(data, out_of_bounds, steps);

	return shared_ptr<img_field_type>(
		new img_field_type(grid, data, times, out_of_bounds, verbose));
}

shared_ptr<curv_field_type>
	load_curvilinear_VTK_time_steps(const std::vector<std::string>& steps,
		const std::vector<double>& times) {
	assert(!steps.empty());
	std::vector< std::shared_ptr<std::vector<vec_t> > > data;
	std::vector<bool> out_of_bounds;
	VTK_SMART(vtkStructuredGrid) grid =
		import_data<vtkStructuredGrid>(data, out_of_bounds, steps);

	shared_ptr<curv_locator_type>
		locator(new curv_locator_type(grid, verbose));
	return shared_ptr<curv_field_type>(
		new curv_field_type(locator, data, times, verbose));
}

shared_ptr<unstr_field_type>
	load_unstructured_VTK_time_steps(const std::vector<std::string>& steps,
		const std::vector<double>& times) {
	assert(!steps.empty());
	std::vector< std::shared_ptr<std::vector<vec_t> > > data;
	std::vector<bool> out_of_bounds;
	VTK_SMART(vtkUnstructuredGrid) grid =
		import_data<vtkUnstructuredGrid>(data, out_of_bounds, steps);

	shared_ptr<unstr_locator_type>
		locator(new unstr_locator_type(grid, verbose));
	return shared_ptr<unstr_field_type>(
		new unstr_field_type(locator, data, times, verbose));
}

// added for BARG
shared_ptr<BARG_field_type>
	load_BARG_time_steps(const std::vector<std::string>& steps,
		const std::vector<double>& times) {
	assert(!steps.empty());
	std::vector< std::shared_ptr<std::vector<vec_t> > > data;
	std::vector<bool> out_of_bounds;

	VTK_SMART(vtkRectilinearGrid) rgrid =
		import_data<vtkRectilinearGrid>(data, out_of_bounds, steps);

	VTK_CREATE(boundaryAwareRectGrid, grid);
	grid->ShallowCopy(rgrid);

	return shared_ptr<BARG_field_type>(
		new BARG_field_type(grid, data, times, out_of_bounds, verbose, use_bdry_aware));
}

int main(int argc, const char* argv[])
{
	using namespace spurt;
	using namespace odeint;

	initialize(argc, argv);

	if (do_log) {
		ofs.open(name_out + ".log", std::ios::app);
	}
	nb_threads = std::thread::hardware_concurrency();
	std::cout << nb_threads << " threads available\n";
	const std::string ext = spurt::filename::extension(name_in);

	if (ext == "tdvtk") {
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
				std::shared_ptr<img_field_type> field = load_image_VTK_time_steps(steps, times);
				runTBB(field);
			}
			else if (vtkRectilinearGrid::SafeDownCast(mystery)) {
				std::shared_ptr<rect_field_type> field = load_rectilinear_VTK_time_steps(steps, times);
				runTBB(field);
			}
			else if (vtkStructuredGrid::SafeDownCast(mystery)) {
				std::shared_ptr<curv_field_type> field = load_curvilinear_VTK_time_steps(steps, times);
				runTBB(field);
			}
			else if (vtkUnstructuredGrid::SafeDownCast(mystery)) {
				std::shared_ptr<unstr_field_type> field = load_unstructured_VTK_time_steps(steps, times);
				runTBB(field);
			}
			else {
				std::cerr << "Unsupported VTK dataset type\n";
				exit(1);
			}
		}
		else if (ext == "vti") {
			std::shared_ptr<img_field_type> field = load_image_VTK_time_steps(steps, times);
			runTBB(field);
		}
		else if (ext == "vtu") {
			std::shared_ptr<unstr_field_type> field = load_unstructured_VTK_time_steps(steps, times);
			runTBB(field);
		}
		else if (ext == "vtr") {
			std::shared_ptr<rect_field_type> field = load_rectilinear_VTK_time_steps(steps, times);
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
	else
	{
		std::cerr << "Unrecognized input file extension: " << ext << '\n';
		exit(1);
	}

	if (do_log) ofs.close();

	return 0;
}
