#include <string>

#include <vtkArrowSource.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtk/vtk_utils.hpp>

#include <format/filename.hpp>

std::string me;
void usage(const std::string& msg = "") {
    if (!msg.empty()) 
        std::cerr << "\nERROR: " << msg << '\n';
    std::cout 
        << '\n'
        << "USAGE: " << me << " [parameters] [options]\n"
        << '\n'
        << "PARAMETERS:\n"
        << " -i | --input <string>     Input file name\n"
        << '\n'
        << "OPTIONS:\n"
        << " -h | --help               Print this information\n"
        << " -o | --output <string>    Output image filename\n"
        << " -n | --nb_steps <int>     Number of integration steps\n"
        << " -e | --epsilon <float>    Integration step size\n"
        << " -r | --radius <float>     Radius of arrow glyphs\n"
        << " -s | --scale <float>      Scaling factor for glyphs\n"
        << " -v | --verbose            Turn on verbose mode\n"
        << "      --no-lic             Turn off LIC computation\n"
        << "      --no-glyphs          Turn off glyph visualization\n"
        << '\n';
        
    exit(!msg.empty());
}

int main(int argc, char* argv[]) {
    me                   = argv[0];
    std::string in_name  = "";
    std::string out_name = "";
    // computation parameters
    int nsteps           = 100;
    double eps           = 0.5;
    // visualization parameters
    double radius        = 0.01; // in %
    double scale         = 1;
    // control flags
    bool nolic           = false;
    bool noglyphs        = false;
    bool verbose         = false;
    
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage();
        }
        else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        }
        else if (arg == "--no-lic" || arg == "--no-LIC") {
            nolic = true;
        }
        else if (arg == "--no-glyphs" || arg == "--no-glyph") {
            noglyphs = true;
        }
        else if (arg == "-i" || arg == "--input") {
            if (i == argc-1) usage("Missing input filename");
            in_name = argv[++i];
        }
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) usage("Missing output filename");
            out_name = argv[++i];
        }
        else if (arg == "-n" || arg == "--nb_steps") {
            if (i == argc-1) usage("Missing number of steps");
            nsteps = atoi(argv[++i]);
        }
        else if (arg == "-e" || arg == "--epsilon") {
            if (i == argc-1) usage("Missing integration steps size");
            eps = atof(argv[++i]);
        }
        else if (arg == "-r" || arg == "--radius") {
            if (i == argc-1) usage("Missing radius value");
            radius = atof(argv[++i]);
        }
        else if (arg == "-s" || arg == "--scale") {
            if (i == argc-1) usage("Missing scale value");
            scale = atof(argv[++i]);
        }
        else usage("unrecognized parameters: " + arg);
    }
    
    if (in_name.empty()) 
        usage("Missing input filename");
    else if (out_name.empty() && !nolic) 
        usage("Missing output filename for LIC image");
    
    vtkImageData *rhs = vtk_utils::load_nrrd(in_name);
    if (verbose) {
        std::cout << in_name << " successfully imported\n";
        int res[3];
        rhs->GetDimensions(res);
        std::cout << "image resolution: " << res[0] << " x " << res[1] << '\n';
    }
    double* b = rhs->GetBounds();
    nvis::bbox2 bounds(nvis::vec2(b[0], b[2]), nvis::vec2(b[1], b[3]));
    double diameter = nvis::norm(bounds.size());
    if (verbose) std::cout << "domain bounds:\n" << bounds << '\n';
    radius *= diameter;
    
    VTK_CREATE(vtkRenderer, renderer);
    renderer->SetBackground(0,0,0);
    
    if (!nolic) {
        vtkImageData *lic = vtk_utils::do_lic(rhs, nsteps, eps, 4);
        if (verbose) std::cout << "LIC image successfully computed\n";
        vtk_utils::save_image(lic, out_name);
        if (verbose) std::cout << out_name << " successfully exported\n";
        lic->Delete();
    }
    
    if (!noglyphs) {
        VTK_CREATE(vtkArrowSource, arrow);
        arrow->SetTipRadius(3*radius);
        arrow->SetShaftRadius(radius);
        VTK_CREATE(vtkGlyph3D, arrows);
        VTK_CONNECT(arrows, rhs);
        arrows->SetSourceConnection(arrow->GetOutputPort());
        arrows->OrientOn();
        arrows->ScalingOn();
        arrows->SetScaleModeToScaleByVector();
        arrows->SetColorModeToColorByVector();
        arrows->SetScaleFactor(scale);
        arrows->ClampingOff();
        VTK_CREATE(vtkPolyDataMapper, vec_mapper);
        VTK_PLUG(vec_mapper, arrows);
        vec_mapper->ScalarVisibilityOff();
        VTK_CREATE(vtkActor, vec_actor);
        vec_actor->SetMapper(vec_mapper);
        vec_actor->GetProperty()->SetColor(1, 0, 0);
        renderer->AddActor(vec_actor);
        if (verbose) std::cout << "arrow glyphs created\n";
    }
    
    VTK_CREATE(vtkRenderWindow, window);
    window->PointSmoothingOn();
    window->LineSmoothingOn();
    window->PolygonSmoothingOn();
    window->AddRenderer(renderer);
    window->SetSize(800, 800); // spurt's laptop likes 800x800

    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    interactor->Initialize();
    window->Render();
    
    if (verbose) std::cout << "interactor starts now...\n"; 
    interactor->Start();
  
    rhs->Delete();
    return 0;
}