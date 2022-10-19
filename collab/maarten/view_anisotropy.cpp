#include <VTK/vtk_utils.hpp>
#include <string>
#include <vector>
#include <iostream>
#include <set>
#include <list>
#include <image/nrrd_wrapper.hpp>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkWarpScalar.h>
#include <vtkDataSetReader.h>
#include <vtkImageReader.h>
#include <vtkImageData.h>
#include <vtkGlyph3D.h>
#include <vtkArrowSource.h>


// notational convenience
typedef VTK_SMART(PolyData)      smart_poly;
typedef VTK_SMART(Actor)         smart_actor;
typedef VTK_SMART(Renderer)      smart_renderer;
typedef VTK_SMART(RenderWindow)  smart_window;

template<typename T>
class anisotropy_helper {
public:
    typedef typename std::vector<T>::size_type size_type;
    
    anisotropy_helper(const std::vector<T>& data, size_type stride) 
        : _data(data), _stride(stride), _dtheta(2.*M_PI/(double)stride) {}
        
    double value(size_type n, size_type loc) const {
        size_type base_id = n*_stride;
        return (double)_data[base_id+loc];
    } 
        
    size_type max_id(size_type n) const {
        size_type begin = n*_stride;
        size_type end = begin + _stride;
        return std::distance(&_data[begin], 
                             std::max_element(&_data[begin], &_data[end]));
    }
    
    nvis::vec2 gradient(size_type n, size_type loc) const {
        double theta = _dtheta*((double)loc + 0.5);
        return value(n, loc) * nvis::vec2(cos(theta), sin(theta));    
    }
    
    nvis::vec2 max_gradient(size_t n) const {
        size_type id = max_id(n);
        return gradient(n, id); 
    }

private:
    const std::vector<T>& _data;
    size_t _stride;
    double _dtheta;
};

std::string me;
void usage(const std::string& message="") {
    if (!message.empty()) {
        std::cerr << "ERROR: " << message << '\n';
    }
    std::cout << '\n'
              << "DESCRIPTION: Visualize anisotropy of tomography data.\n"
              << '\n'
              << "USAGE: " << me << " [parameters] [options]\n"
              << '\n'
              << "PARAMETERS:\n"
              << " -i  | --input <string>         Anisotropy file name\n"
              << '\n'                             
              << "OPTIONS:\n"                     
              << " -h  | --help                   Print this information and exit\n"
              << " -m  | --magnitude <string>     Anisotropy magnitude\n"
              << " -hs | --height_scale <float>   Height field scaling\n"
              << " -vs | --vector_scale <float>   Vector glyph scaling\n"                             
              << " -p  | --path <string>          Path to prepend to all file names\n"
              << " -r  | --resolution <int> x2    Image resolution\n"
              << " -n  | --nopoly                 Turn off polygon depiction\n"
              << " -v  | --verbose                Turn on verbose mode\n"
              << std::endl;
    exit(!message.empty());
}

int main(int argc, char* argv[]) {
    std::string input_name="", mag_name="", path="";
    nvis::ivec2 resolution(800, 800);
    bool verbose = false;
    bool nopoly = false;
    double height_scale=1, vec_scale=1;
    
    me = argv[0];
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage();
        } else if (arg == "-i" || arg == "--input") {
            if (i == argc-1) {
                usage("missing input filename");
            }
            input_name = argv[++i];
        } else if (arg == "-m" || arg == "--magnitude") {
            if (i == argc-1) {
                usage("missing magnitude filename");
            }
            mag_name = argv[++i];
        } else if (arg == "-p" || arg == "--path") {
            if (i == argc-1) {
                usage("missing data path");
            }
            path = argv[++i];
        } else if (arg == "-r" || arg == "--resolution") {
            if (i >= argc-2) {
                usage("missing image resolution");
            }
            resolution[0] = atoi(argv[++i]);
            resolution[1] = atoi(argv[++i]);
        } else if (arg == "-hs" || arg == "--height_scale") {
            if (i == argc-1) {
                usage("missing height field scale");
            }
            height_scale = atof(argv[++i]);
        } else if (arg == "-vs" || arg == "--vector_scale") {
            if (i == argc-1) {
                usage("missing vector field scale");
            }
            vec_scale = atof(argv[++i]);
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else if (arg == "-n" || arg == "--nopoly") {
            nopoly = true;
        }else {
            usage("unrecognized argument: " + arg);
        }
    }
    if (input_name.empty()) usage("Missing input filename");
    
    VTK_CREATE(vtkRenderer, renderer);
    renderer->SetBackground(0,0,0);
    
    srand48(time(0));
    
    if (!mag_name.empty()) {
        VTK_CREATE(vtkDataSetReader, img_reader);
        img_reader->SetFileName(mag_name.c_str());
        img_reader->Update();
        
        if (verbose) {
            double* bds = img_reader->GetOutput()->GetBounds();
            std::cout << "bounds are "
            << bds[0] << " - " << bds[1] << ", "
            << bds[2] << " - " << bds[3] << ", " 
            << bds[4] << " - " << bds[5] << '\n';
        }
        
        VTK_CREATE(vtkStructuredPoints, img);
        img = img_reader->GetStructuredPointsOutput();
        VTK_CREATE(vtkWarpScalar, warp);
        warp->SetInput(img);
        warp->SetScaleFactor(height_scale);
        vtkDataSetMapper* warp_mapper = vtkDataSetMapper::New();
        warp_mapper->SetInputConnection(warp->GetOutputPort());
        warp_mapper->ScalarVisibilityOff();
        vtkActor *warp_actor = vtkActor::New();
        warp_actor->SetMapper(warp_mapper);
        warp_actor->GetProperty()->SetColor(0.7,0.5,0);
        renderer->AddActor(warp_actor);
    }
    
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, input_name.c_str(), NULL)) {
        std::cerr << "readNrrd: " << biffGetDone(NRRD) << std::endl;
        throw;
    }
    std::vector<float> __array;
    xavier::to_vector(__array, nin);
    size_t res_x = nin->axis[1].size;
    size_t res_y = nin->axis[2].size;
    size_t nbuckets = nin->axis[0].size;
    double dtheta = 2.*M_PI/(double)nbuckets;
    nvis::vec2 spc(nin->axis[1].spacing, nin->axis[2].spacing);
    nvis::vec2 origin(nin->axis[1].min, nin->axis[2].min);
    if (verbose) {
        std::cout << input_name << " has been imported and converted\n";
        std::cout << "dataset properties:\n"
                  << "number of points:   " << res_x*res_y << '\n'
                  << "resolution:         " << res_x << " x " << res_y << '\n'
                  << "angular resolution: " << nbuckets << " x " << dtheta << '\n'
                  << "bounding box:       " << origin << " -> " << origin + nvis::vec2(res_x, res_y)*spc << '\n';
    }
    anisotropy_helper<float> helper(__array, nbuckets);
    
    size_t npts = res_x * res_y;
    std::vector<nvis::vec2> points(2*npts);
    std::vector<nvis::vec2> vecs(2*npts);
    std::vector<float> cell_values;
    std::vector<nvis::vec2> glyph_pts;
    VTK_CREATE(vtkCellArray, glyph_polys);
    double min, max;
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::min();
    nvis::bbox2 glyph_bounds;
    for (size_t n=0 ; n<npts ; ++n) {
        size_t i = n % res_x;
        size_t j = n / res_x;
        points[2*n  ] = origin + nvis::vec2(i,j)*spc;
        points[2*n+1] = origin + nvis::vec2(i,j)*spc;
        
        std::vector<nvis::vec2> poly_points;
        if (!nopoly) {
            // draw anisotropy polygon
            for (size_t k=0 ; k<nbuckets ; ++k) {
                nvis::vec2 v = helper.gradient(n, k);
                if (nvis::norm(v) == 0) continue;
                poly_points.push_back(points[2*n] + vec_scale*v);
            }
            
            size_t start_at = glyph_pts.size();
            
            glyph_polys->InsertNextCell(poly_points.size());
            for (size_t k=0 ; k<poly_points.size() ; ++k) {
                glyph_pts.push_back(poly_points[k]);
                glyph_polys->InsertCellPoint(start_at + k);
            }
        }
        vecs[2*n] = helper.max_gradient(n);
        vecs[2*n+1] = -1.*vecs[2*n];
        if (nvis::norm(vecs[2*n]) < min) min = nvis::norm(vecs[2*n]);
        if (nvis::norm(vecs[2*n]) > max) max = nvis::norm(vecs[2*n]);
        if (!nopoly) {
            double val = nvis::norm(vecs[2*n]);
            for (size_t k=0 ; k<poly_points.size() ; ++k) {
                cell_values.push_back(val);
            }
        }
    }
    if (verbose) {
        std::cout << "magnitude ranges from " << min << " to " << max << '\n';
        std::cout << "there are " << glyph_pts.size() << " points and "
                  << glyph_polys->GetNumberOfCells() << " polygons\n";
        std::cout << "glyph bounding box is " << glyph_bounds << '\n';
    }
    
    VTK_CREATE(vtkColorTransferFunction, ctf);
    ctf->AddRGBPoint(min, 0, 0, 1);
    ctf->AddRGBPoint(0.5*(min+max), 0.5, 0.5, 0.5);
    ctf->AddRGBPoint(max, 1, 1, 0);
    
    if (!nopoly) {
        VTK_CREATE(vtkPolyData, polygons);
        polygons = vtk_utils::make_points(glyph_pts);
        vtk_utils::add_vertices(polygons);
        vtk_utils::add_scalars(polygons, cell_values);
        polygons->SetPolys(glyph_polys);
        VTK_MAKE_ACTOR(poly_actor, polygons);
        poly_actor->GetMapper()->ScalarVisibilityOn();
        poly_actor->GetMapper()->SetLookupTable(ctf);
        renderer->AddActor(poly_actor);
    }
    VTK_CREATE(vtkPolyData, vec_poly);
    vec_poly = vtk_utils::make_points(points);
    vtk_utils::add_vectors(vec_poly, vecs);
    
    VTK_CREATE(vtkArrowSource, arrow);
    arrow->SetTipRadius(0.2);
    arrow->SetShaftRadius(0.06);
    VTK_CREATE(vtkGlyph3D, arrows);
    VTK_CONNECT(arrows, vec_poly);
    arrows->SetSourceConnection(arrow->GetOutputPort());
    arrows->OrientOn();
    arrows->ScalingOn();
    arrows->SetScaleModeToScaleByVector();
    arrows->SetColorModeToColorByVector();
    arrows->SetScaleFactor(vec_scale);
    arrows->ClampingOff();
    VTK_CREATE(vtkPolyDataMapper, vec_mapper);
    VTK_PLUG(vec_mapper, arrows);
    // vec_mapper->SetLookupTable(ctf);
    vec_mapper->ScalarVisibilityOff();
    VTK_CREATE(vtkActor, vec_actor);
    vec_actor->SetMapper(vec_mapper);
    vec_actor->GetProperty()->SetColor(1, 0, 0);
    renderer->AddActor(vec_actor);
    
    VTK_CREATE(vtkRenderWindow, window);
    window->PointSmoothingOn();
    window->LineSmoothingOn();
    window->PolygonSmoothingOn();
    window->AddRenderer(renderer);
    window->SetSize(resolution[0], resolution[1]);

    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    interactor->Initialize();
    window->Render();
    interactor->Start();

    exit(0);
}