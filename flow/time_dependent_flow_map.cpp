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
#include <math/small_vector.hpp>
#include <math/bounding_box.hpp>

#include <flow/ode_observer.hpp>
#include <flow/time_dependent_field.hpp>
#include <format/dlr_reader.hpp>

#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <misc/strings.hpp>
#include <misc/meta_utils.hpp>

#include <vtk/vtk_interpolator.hpp>

#include <tbb/parallel_for.h>
#include <tbb/tbb.h>
std::atomic<size_t> progress_counter;

#include <boost/numeric/odeint.hpp>
#include <boost/filesystem.hpp>
namespace odeint = boost::numeric::odeint;


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
spurt::bbox3 the_bounds;
bool in_parallel=false;

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

using namespace spurt;
using namespace vtk_utils;

typedef double scalar_type; 
typedef size_t size_type;
typedef spurt::small_vector<scalar_type, 3> point_type;
typedef spurt::small_vector<scalar_type, 3> vector_type;
typedef spurt::small_vector<size_type, 3> coord_type;
typedef spurt::bounding_box<point_type> bounds_type;

typedef interpolator<vtkImageData, scalar_type, 3, vector_type>  img_intp_t;
typedef interpolator<vtkRectilinearGrid, scalar_type, 3, vector_type> rect_intp_t;
typedef interpolator<vtkUnstructuredGrid, scalar_type, 3, vector_type> unst_intp_t;
typedef interpolator<vtkStructuredGrid, scalar_type, 3, vector_type> curv_intp_t;

typedef vtk_utils::point_locator<vtkUnstructuredGrid, scalar_type, 3, point_type> locator_type;
typedef spurt::fixed_mesh_time_dependent_field<locator_type, std::vector<vector_type> > field_type;

typedef Observer<point_type> observer_t;


int runTBB(shared_ptr<field_type> field) {
    scalar_type global_bounds[6];
    field->get_locator()->get_dataset()->GetBounds(global_bounds);
    bbox3 bnds;

    std::vector<point_type> seeds;

    if (::bounds[0]<::bounds[1] && ::bounds[2]<::bounds[3] && ::bounds[4]<::bounds[5] &&
        ::bounds[0] >= global_bounds[0] && ::bounds[1] <= global_bounds[1] &&
        ::bounds[2] >= global_bounds[2] && ::bounds[3] <= global_bounds[3] &&
        ::bounds[4] >= global_bounds[4] && ::bounds[5] <= global_bounds[5]) {
        // valid bounds supplied by user
        bnds.min() = point_type(::bounds[0], ::bounds[2], ::bounds[4]);
        bnds.max() = point_type(::bounds[1], ::bounds[3], ::bounds[5]);
    }
    else {
        bnds.min() = point_type(global_bounds[0], global_bounds[2], global_bounds[4]);
        bnds.max() = point_type(global_bounds[1], global_bounds[3], global_bounds[5]);
    }

    the_bounds = bnds;
    int npoints;
    spurt::raster_grid<size_t, scalar_type, 3> sampling_grid(res, bnds);

    if (seed_name.empty()) {
        std::cout << "Resolution = " << res[0] << " x " << res[1] << " x " << res[2] << std::endl;

        std::cout << "sampling grid bounds are: " << sampling_grid.bounds().min()
            << " -> " << sampling_grid.bounds().max() << '\n';
        npoints = sampling_grid.size();
        seeds.resize(npoints);
        for (int i=0; i<npoints; i++) {
            seeds[i] = sampling_grid(sampling_grid.coordinates(i));
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

    tbb::parallel_for(tbb::blocked_range<int>(0,npoints),
                       [&](tbb::blocked_range<int> r) 
    {
        for (int n=r.begin(); n!=r.end(); ++n) {
            flowmap[3*n  ] = seeds[n][0];
            flowmap[3*n+1] = seeds[n][1];
            flowmap[3*n+2] = seeds[n][2];

            ++progress_counter;
        }
    });


    size_t nbcomputed=0;
    spurt::ProgressDisplay progress(true);
    progress.fraction_on();

    size_t nb_lost = 0;
    progress.start(npoints);

    std::vector<bool> stopped(npoints, false);
    std::vector< std::vector<point_type> > lines(npoints);
    std::vector< std::vector<scalar_type> > times(npoints);

    progress_counter = 0;
    tbb::parallel_for(tbb::blocked_range<int>(0,npoints),
                       [&](tbb::blocked_range<int> r) 
    {
        for (int n=r.begin(); n!=r.end(); ++n) {
            ++progress_counter;

            if (!(progress_counter%10)) progress.update(progress_counter);

            if (stopped[n]) {
                continue;
            }
            point_type x = {flowmap[3*n], flowmap[3*n+1], flowmap[3*n+2]};
            scalar_type t=t0, d=0;
            point_type p(x);
            std::vector<point_type>& sline = lines[n];
            std::vector<scalar_type>& slinetimes = times[n];

            bool do_verbose = false;

            // create a stepper
            auto stepper = odeint::make_controlled(eps, eps, odeint::runge_kutta_dopri5<point_type>());
            observer_t obs(p, t, d, sline, slinetimes, do_verbose);

            try {
                point_type y(x);
                integrate_adaptive(stepper, *field, y, t0, T, 1.0e-4, obs);
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
    }); // end of tbb::parallel_for

    progress.end();

    if (save_lines || monitor) {
        if (monitor) {
            std::vector< std::vector<point_type> > newlines;
            std::vector< std::vector<scalar_type> > newtimes;
            for (size_t i=0; i<lines.size(); i++) {
                if (lines[i].size() > 2) {
                    newlines.push_back(lines[i]);
                    newtimes.push_back(times[i]);
                }
            }
            lines.swap(newlines);
            times.swap(newtimes);
        }
        std::vector<spurt::ivec2> dummy;
        VTK_SMART(vtkPolyData) pdata = vtk_utils::make_polylines(lines, dummy, 0);
        std::vector<scalar_type> all_times;
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
        std::vector<size_type> size(4);
        std::vector<scalar_type> step(4), mins(4);
        step[0] = std::numeric_limits<scalar_type>::quiet_NaN();
        mins[0] = std::numeric_limits<scalar_type>::quiet_NaN();
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
        spurt::nrrd_utils::writeNrrdFromContainers(flowmap, os.str(), size, step, mins, ctrs, cmdline);

        size[0] = 1;
        os.clear();
        os.str("");
        os << name_out << (T>0 ? "-fwd_" : "-bwd_") << "flowtime-deltaT=" << fabs(T) << "_t0=" << t0 << ".nrrd";
        spurt::nrrd_utils::writeNrrdFromContainers(flowtimes, os.str(), size, step, mins, ctrs, cmdline);
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
    std::vector<scalar_type> times;
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
            scalar_type t;
            iss >> name >> t;
            steps.push_back(name);
            times.push_back(t);
        }
    }
    info_file.close();

    spurt::dlr_reader reader(mesh_name, "");
    std::vector<fvec3> vertices;
    std::vector<long> cell_indices;
    std::vector<std::pair<spurt::dlr_reader::cell_type, long> > cell_types;
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
            case spurt::dlr_reader::TRIANGLE:
                type_name = VTK_TRIANGLE;
                break;
            case spurt::dlr_reader::QUADRILATERAL:
                type_name = VTK_QUAD;
                break;
            case spurt::dlr_reader::TETRAHEDRON:
                type_name = VTK_TETRA;
                break;
            case spurt::dlr_reader::HEXAHEDRON:
                type_name = VTK_HEXAHEDRON;
                break;
            case spurt::dlr_reader::PRISM:
                type_name = VTK_WEDGE;
                break;
            case spurt::dlr_reader::PYRAMID:
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
    }

    return shared_ptr<field_type>(new field_type(locator, data, times, verbose));
}

int main(int argc, const char* argv[])
{
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
