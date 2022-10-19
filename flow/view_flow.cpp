#include <regex>
#include <string>

#include <vtk/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <misc/time_helper.hpp>
#include <math/stat.hpp>
#include <graphics/colors.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

typedef spurt::fixed_vector< double, 2 > vec2d;
typedef spurt::fixed_vector< double, 3 > vec3d;
typedef spurt::fixed_vector< int, 2 >    vec2i;
typedef spurt::fixed_vector< int, 3 >    vec3i;
typedef nvis::bounding_box< vec2d >    bbox2d;

std::string name_in, name_out, name_mask, name_cmap;
bool verbose;
bool use_pts = false;
bool use_spheres = false;
bool follow_camera = false;
std::string name_camera;
vec2i res(800, 800), data_res;
vec2d spc;
size_t npts, natt;
int down_factor=1;
bbox2d bounds, gbounds;
std::array< double, 4 > bounds_as_array;
Nrrd *mask=0, *fmap;
int nsteps = -1;
double tinit=0;
double radius=0.1;
double scale=0.001;

constexpr double invalid = -30000;

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
        parser.add_value("mask", name_mask, "Mask filename", optional_group);
        parser.add_tuple< 2 >("res", res, res, "Output image resolution", optional_group);
        parser.add_tuple< 4 >("bounds", bounds_as_array, bounds_as_array, "Sampling bounds (min longitude and latitude followed by max's)", optional_group);
        parser.add_value("down", down_factor, down_factor, "Downsampling factor", optional_group);
        parser.add_value("length", nsteps, nsteps, "Number of files to consider", optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbosity level", optional_group);
        parser.add_flag("spheres", use_spheres, "Use spheres instead of points", optional_group);
        parser.add_value("radius", radius, radius, "Tube radius", optional_group);
        parser.add_value("scale", scale, scale, "Time axis scaling factor", optional_group);
        parser.add_flag("track", follow_camera, "Track camera settings", optional_group);
        parser.add_value("camera", name_camera, "Import camera settings", optional_group);
        parser.add_value("cmap", name_cmap, "Import color map definition", optional_group);
        
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
    return spurt::nrrd_data_wrapper<T>(nin)[n];
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

vtkSmartPointer<vtkActor> create_points(std::vector<Nrrd*> datasets, const std::vector<int>& times) {
    size_t nskipped = 0;
    std::vector< vec3d > points;
    std::vector< std::vector<size_t> > orbits;
    std::vector<float> values;
    
    double threshold = 1./down_factor;
    
    // longitude: 262 - 279.5
    // latitude: 18 - 30.74
    
    for (int t=0; t<datasets.size(); ++t) {
        size_t orbit_id=0;
        srand48(123456);
        for (size_t i=0; i<npts; ++i) {
            if (mask != NULL && nrrd_value(mask, i)==0) {
                ++nskipped;
                continue;
            }
            else if (drand48() > threshold) {
                ++nskipped;
                continue;
            }
            
            vec3d x( nrrd_value( datasets[t], natt*i ), nrrd_value( datasets[t], natt*i+1), 
                     scale*(times[t]-tinit) );
            if (!bounds.empty() && !bounds.inside(where(i))) continue;
                     
            if (t==0) {
                orbits.push_back(std::vector<size_t>());
                orbits.back().push_back(points.size());
            }
            else orbits[orbit_id].push_back(points.size());
            points.push_back(x);
            values.push_back(nrrd_value( datasets[t], natt*i+2)/(0.1*static_cast<double>(times[t]-tinit)));
            
            ++orbit_id;
        }
        std::cout << "after " << t+1 << " iterations: " << points.size() << " points and " 
                  << orbits.size() << " orbits\n";
    }
    
    if (verbose) {
        std::cout << nskipped << " skipped entries\n";
    }
    
    vtkSmartPointer<vtkPolyData> poly = vtk_utils::make_points(points);
    std::cout << "points added\n";
    poly = vtk_utils::add_scalars(poly, values);
    std::cout << "values added\n";
    poly = vtk_utils::add_polylines(poly, orbits);
    std::cout << "polylines added\n";
    
    vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
    tubes->SetInputData(poly);
    tubes->SetVaryRadiusToVaryRadiusByScalar();
    tubes->SetRadius(radius);
    tubes->SetNumberOfSides(6);
    
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    
    mapper->SetInputConnection(tubes->GetOutputPort());
    mapper->ScalarVisibilityOn();
    
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    
    if (name_cmap.empty()) {
        std::vector<spurt::fvec3> scale;
        
        // heat map:
        std::vector<spurt::fvec3> heat_scale;
        heat_scale.push_back(spurt::black);
        heat_scale.push_back(spurt::red);
        // heat_scale.push_back(spurt::orange);
        heat_scale.push_back(spurt::yellow);
        heat_scale.push_back(spurt::white);
        // heat_scale.push_back(spurt::fvec3(0.5, 0.5, 1));
        // heat_scale.push_back(spurt::blue);
        
        std::cout << "heat_scale created\n";
        
        spurt::adaptive_color_map<float> cmap(values, heat_scale, true, 20);
        
        std::cout << "heat color map created\n";
        for (int i=0; i<cmap.t.size() ; ++i) {
            ctf->AddRGBPoint(cmap.t[i], cmap.colors[i][0], cmap.colors[i][1], cmap.colors[i][2]);
        }
    }
    else {
        vtk_utils::import_colormap(name_cmap, ctf);
    }
    
    std::cout << "color transfer function created\n";
    
    mapper->SetLookupTable(ctf);
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    
    actor->GetProperty()->SetLineWidth(1);
    // actor->GetProperty()->EdgeVisibilityOn();
    
    return actor;
}

int main(int argc, char* argv[]) {
    initialize(argc, (const char**)argv);
    
    if (!name_mask.empty()) {
        mask = spurt::readNrrd( name_mask );
    }
    
    std::regex time_regex("([0-9]+)h");
    std::smatch time_match;
    
    std::vector<Nrrd*> datasets;
    std::vector<int> times;
    std::fstream input(name_in, std::ios::in);
    size_t n=0;
    while (!input.eof() && input.good()) {
        if (nsteps > 0 && n>=nsteps) break;
        ++n;
        std::string name;
        input >> name;
        if (input.fail()) break;
        datasets.push_back(spurt::readNrrd(name));
        std::cout << "just imported: " << name << '\n';
        if (std::regex_search(name, time_match, time_regex)) {
            int t = std::stoi(time_match[1].str());
            times.push_back(t);
            std::cout << "corresponding time stamp: " << t << '\n';
        }
    }
    input.close();
    
    data_res[0] = datasets[0]->axis[1].size;
    data_res[1] = datasets[0]->axis[2].size;
    npts = data_res[0]*data_res[1];
    natt = datasets[0]->axis[0].size;
    spc[0] = datasets[0]->axis[1].spacing;
    spc[1] = datasets[0]->axis[2].spacing;
    gbounds.min() = vec2d( datasets[0]->axis[1].min, datasets[0]->axis[2].min );
    gbounds.max() = gbounds.min() + vec2d( data_res-vec2i(1) )*spc;
    
    std::cout << "gbounds.size()=" << gbounds.size() << '\n';
    size_t nskipped=0;
    
    vtkSmartPointer<vtkActor> actor = create_points(datasets, times);
    
/*
    std::vector<float> raster(data_res[0]*data_res[1]*3);
    for (size_t i=0; i<npts; ++i) {
        vec2d x( nrrd_value( fmap, natt*i ), nrrd_value( fmap, natt*i+1 ) );
        if (mask != NULL && nrrd_value(mask, i)==0) {
            raster[3*i] = invalid;
            // raster[i]=-1;
            ++nskipped;
            continue;
        }
        vec2d u = ( x-gbounds.min() )/gbounds.size()*vec2d( data_res );
        for (int j=0; j<2; ++j) {
            if (u[j]<0) u[j]=0;
            else if (u[j]>=data_res[j]) u[j]=data_res[j]-1;
        }
        int id = static_cast<int>(u[0]) + static_cast<int>(u[1])*data_res[0];
        int sx = i % data_res[0];
        int sy = i / data_res[1];
        raster[3*id]   += sx;
        raster[3*id+1] += sy;
        raster[3*id+2]++;
    }
    
    size_t n_nonzero=0;
    for (size_t i=0; i<npts ; ++i) {
        if ( raster[3*i+2] != 0 && raster[3*i+2] != -3000) {
            raster[3*i  ] /= raster[3*i+2];
            raster[3*i+1] /= raster[3*i+2]; 
            ++n_nonzero;
        }
    }
    
    if (verbose) {
        std::cout << n_nonzero << " non zero entries (" 
                  << 100.*static_cast<float>(n_nonzero)/
                     static_cast<float>(npts) << "%)\n";
        std::cout << nskipped << " skipped entries\n";
    }
    
    std::vector<double> _spc(2);
    std::vector<double> _min(2);
    std::vector<size_t> _dim(2);
    _spc[0] = spc[0];
    _spc[1] = spc[1];
    _min[0] = gbounds.min()[0];
    _min[1] = gbounds.min()[1];
    _dim[0] = data_res[0];
    _dim[1] = data_res[1];
    spurt::writeNrrd((int*)(&raster[0]), name_out, nrrdTypeInt, _dim, _spc, _min);
    
    vtkSmartPointer<vtkStructuredPoints> data = vtkSmartPointer<vtkStructuredPoints>::New();
    data->SetDimensions(data_res[0], data_res[1], 1);
    data->SetSpacing(spc[0], spc[1], 1.);
    data->SetOrigin(gbounds.min()[0], gbounds.min()[1], 0);
    
    vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
    colors->SetNumberOfComponents(3);
    colors->SetNumberOfTuples(npts);
    colors->SetName("Colors");
    spurt::vec3 black(0), red(255,0,0), white(255,255,255), blue(0,0,255), cf;
    spurt::fixed_vector<unsigned char, 3> c;
    for (size_t i=0; i<npts; ++i) {
        float x = raster[3*i];
        float y = raster[3*i+1];
        x /= (float)data_res[0];
        y /= (float)data_res[1];
        cf = (1-x)*(1-y)*black + x*(1-y)*blue + x*y*white + (1-x)*y*red;
        c[0] = static_cast<unsigned char>(cf[0]);
        c[1] = static_cast<unsigned char>(cf[1]);
        c[2] = static_cast<unsigned char>(cf[2]);
        colors->SetTypedTuple(i, (unsigned char*)&c[0]);
    }
    data->GetPointData()->SetScalars(colors);
    
    std::vector<float> values(npts);
    for (size_t i=0; i<npts; ++i) {
        if (raster[3*i+2] == -3000) values[i] = -3000;
        else if (raster[3*i+2] == 0) values[i] = 0;
        else values[i] = nrrd_value(fmap, natt*i+2);
    }
    vtk_utils::add_scalars(data, values);
    std::vector<float> copy(values.begin(), values.end());
    std::sort(copy.begin(), copy.end());
    auto upper = std::upper_bound(copy.begin(), copy.end(), 0);
    size_t id_0, n_ids;
    if (upper == copy.end()) {
        id_0 = copy.size()-1;
        n_ids = 0;
    } 
    else {
        id_0 = std::distance(copy.begin(), upper);
        n_ids = values.size() - id_0;
    }
    
    std::cout << "id_0 = " << id_0 << ", n_ids = " << n_ids << ", size = " << copy.size() << '\n';
    spurt::vec3 blue(0,0,1);
    spurt::vec3 yellow(1,1,0);
    
    vtkSmartPointer<vtkDataSetMapper> mapper = vtkSmartPointer<vtkDataSetMapper>::New();
    mapper->SetInputData(data);
    mapper->ScalarVisibilityOn();
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    ctf->AddRGBPoint(-3000, 0, 0.25, 0.25);
    ctf->AddRGBPoint(0, 0, 0, 0);
    for (int i=0; i<10; ++i) {
        float u = (float)i/10;
        spurt::vec3 c = (1.-u)*blue + u*yellow;
        float v = copy[id_0 + u*n_ids];
        ctf->AddRGBPoint(v, c[0], c[1], c[2]);
    }
    ctf->AddRGBPoint(copy.back(), 1, 1, 0);
    mapper->SetLookupTable(ctf);
    
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
*/
    
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->SetBackground(0, 0, 0);
    renderer->AddActor(actor);
    vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
    window->AddRenderer(renderer);
    window->SetSize(res[0], res[1]);
    vtk_utils::fill_window(renderer, gbounds);
    
    if (follow_camera) {
        vtk_utils::track_camera_setting(renderer);
    }
    if (!name_camera.empty()) {
        vtk_utils::import_camera_settings(name_camera, renderer);
    }
    
    double width=gbounds.size()[0];
    double height=gbounds.size()[1];
    double ratio = width/height;
    if (res[0]==-1 && res[1]==-1) {
        std::cout << "Invalid dimensions: " << res << "\n";
        res[0] = res[1] = 800;
    }
    else if (res[0]==-1) res[0]=ratio*res[1];
    else if (res[1]==-1) res[1]=res[0]/ratio;
    if (verbose) std::cout << "Setting resolution to " << res << '\n';
    window->SetSize(res[0], res[1]);
    
    if (name_out.empty()) {
        // enter interactive mode
        vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
        interactor->SetRenderWindow(window);
        interactor->Initialize();
        window->Render();
        interactor->Start();
    }
    else {
        window->Render();
        vtk_utils::save_frame(window, name_out);
    }
    
    return 0;
}



