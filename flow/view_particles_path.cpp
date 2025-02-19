#include <algorithm>
#include <random>
#include <regex>
#include <string>

#include <vtkTextActor.h>
#include <vtkTextProperty.h>

#include <vtk/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <math/stat.hpp>
#include <graphics/colors.hpp>
#include <format/filename.hpp>

#include <math/small_vector.hpp>
#include <math/bounding_box.hpp>

#include <boost/filesystem.hpp>

typedef spurt::small_vector< double, 2 > vec2d;
typedef spurt::small_vector< double, 3 > vec3d;
typedef spurt::small_vector< int, 2 >    vec2i;
typedef spurt::small_vector< int, 3 >    vec3i;
typedef spurt::bounding_box< vec2d >    bbox2d;
typedef spurt::vec4                       vec4;

typedef vec4                            pos_t;
typedef std::vector<pos_t>       trajectory_t;



std::string name_in, name_out, name_mask, name_cmap;
bool verbose;
bool follow_camera=false;
bool do_spheres=false;
std::string name_camera;
vec2i res(800, 800), data_res;
vec2d spc;
size_t npts, natt;
float down_factor=1;
bbox2d bounds, gbounds;
std::array< double, 4 > bounds_as_array;
Nrrd *mask=0, *fmap;
int nsteps=-1;
double tinit=-1;
double radius=0.1;
double scale=0.001;
int color_method=0; // constant
size_t nframes=10;
vec2d range(0, -1);
vec2i start(0, 0);
float point_size = 2;
bool do_colorbar = false;
bool normalize = false;
int quality = 75;

constexpr double invalid=-30000;



void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;

    xcl::option_traits
            required_group(true, false, "Required Options"),
            optional_group(false, false, "Optional Group");

    xcl::option_parser parser(argv[0],
            "Visualize flow of uniform seeding distribution");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input filename", required_group);
        parser.add_value("output", name_out, "Output base name", optional_group);
        parser.add_value("tinit", tinit, tinit, "Initial integration time", optional_group);
        parser.add_value("scalar", color_method, color_method, "Scalar to use for color mapping (0: none, 1: lavd, 2: time, 3: id)", optional_group);
        parser.add_value("mask", name_mask, "Mask filename", optional_group);
        parser.add_tuple< 2 >("range", range, "Value range to consider for color mapping", optional_group);
        parser.add_tuple< 2 >("res", res, res, "Output image resolution", optional_group);
        parser.add_value("nframes", nframes, nframes, "Number of frames to create", optional_group);
        parser.add_tuple< 2 >("start", start, start, "Number of files to skip + first frame index", optional_group);
        parser.add_flag("normalize", normalize, "Normalize lavd values", optional_group);
        parser.add_tuple< 4 >("bounds", bounds_as_array, bounds_as_array, "Sampling bounds (min longitude and latitude followed by max's)", optional_group);
        parser.add_value("down", down_factor, down_factor, "Downsampling factor", optional_group);
        parser.add_value("length", nsteps, nsteps, "Number of files to consider", optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbosity level", optional_group);
        parser.add_value("radius", radius, radius, "Tube radius", optional_group);
        parser.add_value("scale", scale, scale, "Time axis scaling factor", optional_group);
        parser.add_value("sphere", do_spheres, do_spheres, "Use sphere representation for particles", optional_group);
        parser.add_flag("track", follow_camera, "Track camera settings", optional_group);
        parser.add_flag("bar", do_colorbar, "Display color bar", optional_group);
        parser.add_value("camera", name_camera, "Import camera settings", optional_group);
        parser.add_value("cmap", name_cmap, "Import color map definition", optional_group);
        parser.add_value("ptsize", point_size, point_size, "Point size", optional_group);

        parser.parse(argc, argv);
        bounds.min() = vec2d(bounds_as_array[0], bounds_as_array[1]);
        bounds.max() = vec2d(bounds_as_array[2], bounds_as_array[3]);
        if (res[0] == -1 && res[1] == -1) {
            std::cout << "invalid dimensions: " << res << '\n';
            exit(0);
        }
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

template<typename T = double>
inline double nrrd_value(const Nrrd* nin, size_t n) {
    return spurt::nrrd_utils::nrrd_data_wrapper<T>(nin)[n];
}

inline vec2d pos(const vec2i& c, const vec2i& res) {
    vec2d r;
    r[0] = static_cast<double>(c[0])/static_cast<double>(res[0]);
    r[1] = static_cast<double>(c[1])/static_cast<double>(res[1]);
    return r;
}

inline vec3i color(const vec2d& x, const bbox2d& bounds) {
    double u = x[0];
    double v = x[1];
    vec3d col =
        (1.-u)*(1.-v)*vec3d(0,0,0) + // black
        u*(1.-v)*vec3d(255,0,0) +    // red
        u*v*vec3d(255,255,255) +     // white
        (1.-u)*v*vec3d(0,0,255);     // blue
    return vec3i(static_cast<int>(col[0]),
                 static_cast<int>(col[1]),
                 static_cast<int>(col[2]));
}

inline vec2d where(size_t id) {
    size_t i = id % data_res[0];
    size_t j = id / data_res[0];
    return gbounds.min() + vec2d(i,j)*spc;
}

vtkSmartPointer<vtkColorTransferFunction> create_ctf(const std::vector<double>& values) {
    VTK_CREATE(vtkColorTransferFunction, ctf);

    if (name_cmap.empty()) {
        std::vector<spurt::fvec3> scale;

        // heat map:
        std::vector<spurt::fvec3> heat_scale;
        heat_scale.push_back(spurt::black);
        heat_scale.push_back(spurt::red);
        heat_scale.push_back(spurt::yellow);
        heat_scale.push_back(spurt::white);

        if (verbose) std::cout << "heat_scale created\n";

        std::vector<double> copy;
        bool has_range = (range[0] < range[1]);
        for (double d : values) {
            if (!has_range || (has_range && d>=range[0] && d<=range[1])) {
                copy.push_back(d);
            }
        }

        spurt::adaptive_color_map<double> cmap(copy, heat_scale, true, 20);

        if (verbose) std::cout << "heat color map created\n";
        for (int i=0; i<cmap.t.size() ; ++i) {
            ctf->AddRGBPoint(cmap.t[i], cmap.colors[i][0], cmap.colors[i][1], cmap.colors[i][2]);
        }
    }
    else {
        vtk_utils::import_colormap(name_cmap, ctf);
    }

    return ctf;
}

bool where(vec4& p, const trajectory_t& traj, double t) {
    if (traj.empty() || t < traj.front()[2] || t > traj.back()[2]) {
        return false;
    }

    // binary search
    size_t lo=0, hi=traj.size()-1;
    while (lo <= hi) {
        size_t mid = lo + (hi-lo)/2;
        if (traj[mid][2] < t) lo=mid+1;
        else if (traj[mid][2] > t) hi=mid-1;
        else {
            p = traj[mid];
            return true;
        }
    }
    double u = (t - traj[hi][2])/(traj[lo][2]-traj[hi][2]);
    p = (1.-u)*traj[hi] + u*traj[lo];
    return true;
}

VTK_SMART(vtkActor)
create_frame_actor(const std::vector<trajectory_t>& trajs, double t) {
    std::vector<double> values;
    std::vector<vec2d> points;
    size_t id=0;
    for (const trajectory_t& tr : trajs) {
        vec4 x;
        ++id;
        if (!where(x, tr, t)) continue;
        vec2d p;
        p[0] = x[0];
        p[1] = x[1];
        points.push_back(p);
        switch (color_method) {
            case 1: values.push_back(x[3]); break;
            case 2: values.push_back(x[2]); break;
            case 3: values.push_back(id); break;
            default: break;
        }
    }
    vtkSmartPointer<vtkPolyData> poly = vtk_utils::make_points(points);
    if (color_method > 0) vtk_utils::add_scalars(poly, values);

    if (verbose) {
        std::cout << "there are " << points.size() << " points in this frame\n";
    }
    VTK_CREATE(vtkPolyDataMapper, mapper);
    if (do_spheres) {
        VTK_SMART(vtkPolyData) spheres = vtk_utils::make_spheres(poly, radius, 6, 6);
        mapper->SetInputData(spheres);
    }
    else {
        vtk_utils::add_vertices(poly);
        mapper->SetInputData(poly);
    }
    if (color_method > 0) {
        mapper->ScalarVisibilityOn();
        mapper->SetLookupTable(create_ctf(values));
    }
    else mapper->ScalarVisibilityOff();

    VTK_CREATE(vtkActor, actor);
    actor->SetMapper(mapper);
    if (!color_method) actor->GetProperty()->SetColor(1,0,0);
    if (!do_spheres) actor->GetProperty()->SetPointSize(point_size);

    return actor;
}

int main(int argc, char* argv[]) {
    initialize(argc, (const char**)argv);

    int skip = start[0];
    if (verbose) std::cout << "skip = " << skip << '\n';


    if (!name_mask.empty()) {
        mask = spurt::nrrd_utils::readNrrd(name_mask);
        if (verbose) std::cout << "mask has been imported\n";
    }

    std::regex time_regex("([0-9]+)h");
    std::smatch time_match;

    boost::filesystem::path p(name_in);
    std::string parent_dir=p.parent_path().string();
    std::vector<Nrrd*> datasets;
    std::vector<int> times;
    std::fstream input(name_in, std::ios::in);
    for (size_t n=0, nused=0; (nsteps <=0 || nused<nsteps) && !input.eof() && input.good(); ++n) {
        std::string name;
        input >> name;
        if (input.fail()) break;
        if (n>=skip) {
            datasets.push_back(spurt::nrrd_utils::readNrrd(parent_dir + '/' + name));
            if (verbose) std::cout << "just imported: " << name << '\n';
            if (std::regex_search(name, time_match, time_regex)) {
                int t = std::stoi(time_match[1].str());
                times.push_back(t);
                if (verbose) std::cout << "corresponding time stamp: " << t << '\n';
            }
            ++nused;
        }
    }
    input.close();

    std::string extension;
    if (!name_out.empty()) {
        boost::filesystem::path p(name_out);
        if (!p.parent_path().empty())
            name_out = p.parent_path().string() + '/';
        name_out += p.stem().string();
        extension = p.extension().string();
        if (verbose) std::cout << "name_out = " << name_out << '\n';
    }

    // VTK_CREATE(vtkActor, mask_actor);
    // if (!name_mask.empty()) {
    //     auto mask = vtk_utils::load_nrrd(name_mask);
    //     VTK_CREATE(vtkDataSetMapper, mask_mapper);
    //     mask_mapper->SetInputData(mask);
    //     mask_mapper->ScalarVisibilityOn();
    //     VTK_CREATE(vtkColorTransferFunction, mask_ctf);
    //     mask_ctf->AddRGBPoint(0, 1, 1, 1);
    //     mask_ctf->AddRGBPoint(1, 0, 0, 0);
    //     mask_mapper->SetLookupTable(mask_ctf);
    //     mask_actor->SetMapper(mask_mapper);
    // }

    std::vector<trajectory_t> trajectories;
    std::map<int, size_t> id2pos;
    int natt = datasets[0]->axis[0].size;
    spurt::bbox2 bounds;
    size_t total_npts=0;
    std::vector<double> timesteps;

    // input file may be of two types:
    // * trajectory file: (x,y,t,lavd,id) x npts
    // * lavd file: (x,y,lavd) x X x Y
    bool is_lavd_file = (natt == 3 && datasets[0]->dim == 3);

    if (verbose) std::cout << "this is " << (is_lavd_file ? "" : "not ") << "a lavd file\n";

    size_t npts = datasets[0]->axis[1].size;
    if (is_lavd_file) npts *= datasets[0]->axis[2].size;
    std::vector<size_t> indices(npts);
    std::iota(indices.begin(), indices.end(), 0);

    if (down_factor<1.) {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle(indices.begin(), indices.end(), g);
        npts = (size_t)((float)npts * down_factor);
        indices.resize(npts);
        std::sort(indices.begin(), indices.end());
    }

    std::cout << "there are " << datasets.size() << " datasets in input\n";

    for (size_t i=0 ; i<datasets.size(); ++i) {
        Nrrd* nin = datasets[i];
        size_t npts = nin->axis[1].size;
        if (is_lavd_file) npts *= nin->axis[2].size;
        for (size_t n=0; n<npts; ++n, ++total_npts) {
            pos_t pt;
            int id;
            pt[0] = nrrd_value<float>(nin, n*natt);    // x
            pt[1] = nrrd_value<float>(nin, n*natt+1);  // y

            if (!is_lavd_file) {
                pt[2] = nrrd_value<float>(nin, n*natt+2);  // t
                pt[3] = nrrd_value<float>(nin, n*natt+3);  // lavd
                id =    nrrd_value<float>(nin, n*natt+4);  // curve id
            }
            else {
                pt[2] = times[i];                          // t
                pt[3] = nrrd_value<float>(nin, n*natt+2);  // lavd
                id = n;                                    // curve id
                if (tinit > 0 && pt[3] > 0) {
                    pt[3] /= (double)(times[i]-tinit);
                }

                if (mask != NULL && nrrd_value(mask, n) < 0.5) {
                    continue;
                }
            }

            bounds.add(vec2d(pt[0], pt[1]));

            timesteps.push_back(pt[2]);

            int entry_id;
            if (id2pos.find(id) == id2pos.end()) {
                trajectories.push_back(trajectory_t());
                id2pos[id] = trajectories.size()-1;
            }
            trajectories[id2pos[id]].push_back(pt);
        }
        nrrdNuke(datasets[i]);
    }
    if (verbose) {
        std::cout << "we have imported " << trajectories.size()
                  << " trajectories and " << total_npts << " points\n";
    }

    double tmin = *std::min_element(timesteps.begin(), timesteps.end());
    double tmax = *std::max_element(timesteps.begin(), timesteps.end());
    if (verbose) {
        std::cout << "time range: " << std::setw(12) << tmin << " - " << tmax << '\n';
        std::cout << "bounds are\n" << bounds << '\n';
    }

    VTK_CREATE(vtkRenderer, renderer);
    renderer->SetBackground(0, 0, 0);
    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
    window->SetSize(res[0], res[1]);

    if (color_method > 0 && !name_cmap.empty() && do_colorbar) {
        vtk_utils::colorbar_param param;
        if (color_method == 1) param.title = "Normalized\nLAVD";
        else if (color_method == 2) param.title = "time";
        else if (color_method == 3) param.title = "ID";
        param.width = res[0]/15;
        param.height = res[1]/3;
        param.pos[0] = 0.9;
        param.pos[1] = 0.65;
        param.font_size = 22;
        param.nlabels = 8;

        VTK_SMART(vtkScalarBarActor) baractor;
        std::vector<double> dummy;
        VTK_SMART(vtkColorTransferFunction) _ctf = create_ctf(dummy);
        baractor = vtk_utils::colorbar(_ctf, param);
        baractor->SetBarRatio(0.5);
        renderer->AddActor(baractor);
    }

    // if (!name_mask.empty()) renderer->AddActor(mask_actor);

    if (follow_camera) {
        vtk_utils::track_camera_setting(renderer);
    }
    if (!name_camera.empty()) {
        vtk_utils::import_camera_settings(name_camera, renderer);
    }

    double dt = (tmax-tmin)/(nframes-1);
    std::cout << "dt = " << dt << '\n';
    for (size_t i=0; i<nframes ; ++i) {
        double t = tmin + i*dt;

        double days = t/24;
        std::string month;
        if (days < 30) month = "April";
        else if (days < 61) {
            month = "May";
            days -= 30;
        }
        else if (days < 91) {
            month = "June";
            days -= 61;
        }
        else if (days < 122) {
            month = "July";
            days -= 91;
        }
        else {
            month = "August";
            days -= 122;
        }

        int day = static_cast<int>(days)+1;
        int hour = static_cast<int>(t-std::floor(t/24.0)*24);
        std::ostringstream os;
        os << month << " " << std::setw(2) << day << ", " << std::setw(2) << hour << ":00h";

        std::cout << os.str() << '\n';

        VTK_CREATE(vtkTextActor, txt);
        txt->SetInput(os.str().c_str());
        txt->GetTextProperty()->SetFontFamilyToArial();
        txt->GetTextProperty()->SetFontSize(22);
        txt->GetTextProperty()->SetColor(1,1,1);
        txt->SetDisplayPosition(20, 20);
        txt->GetPosition2Coordinate()->SetCoordinateSystemToNormalizedViewport();
        txt->GetPosition2Coordinate()->SetValue(0.1, 0.1);

        VTK_SMART(vtkActor) actor = create_frame_actor(trajectories, t);

        std::cout << "generated frame actor for t=" << std::setw(12) << t << '\n';

        renderer->AddActor(actor);
        renderer->AddActor(txt);

        if (name_out.empty()) {
            // enter interactive mode
            VTK_CREATE(vtkRenderWindowInteractor, interactor);
            interactor->SetRenderWindow(window);
            interactor->Initialize();
            window->Render();
            interactor->Start();
        }
        else {
            if (extension == "") extension = ".jpeg";
            std::cout << "extension=" << extension << '\n';
            window->Render();
            std::ostringstream os;
            os << name_out << std::setfill('0') << std::setw(6) << i+start[1] << extension;
            vtk_utils::save_frame(window, os.str(), quality);
            std::cout << "just saved " << os.str() << '\n';
        }

        renderer->RemoveActor(actor);
        renderer->RemoveActor(txt);
    }

    return 0;
}
