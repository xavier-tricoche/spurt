#include <VTK/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <misc/progress.hpp>
#include <math/stat.hpp>
#include <graphics/colors.hpp>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

typedef nvis::fixed_vector< double, 2 > vec2d;
typedef nvis::fixed_vector< double, 3 > vec3d;
typedef nvis::fixed_vector< int, 2 >    vec2i;
typedef nvis::fixed_vector< int, 3 >    vec3i;
typedef nvis::bounding_box< vec2d >    bbox2d;

std::string name_in, name_out, name_cmap_in;
bool verbose;
bool follow_camera = false;
std::string name_camera;
vec2i res(800, 800), data_res;
int down_factor=1;
float quality=100;

void initialize(int argc, const char* argv[])
{
    namespace xcl = xavier::command_line;
	
	srand48(123456);
        
    xcl::option_traits 
            required_group(true, false, "Required Options"), 
            optional_group(false, false, "Optional Group");
            
    xcl::option_parser parser(argv[0],
            "Visualize colored map version of an image in VTK or NRRD format");
    
    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input filename", required_group);
		parser.add_value("cmap", name_cmap_in, "Color map definition", required_group);
        parser.add_value("output", name_out, "Output filename", optional_group);
        parser.add_tuple< 2 >("res", res, res, "Output image resolution", optional_group);
        parser.add_value("down", down_factor, down_factor, "Downsampling factor", optional_group);
		parser.add_value("quality", quality, quality, "Image export quality", optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbosity level", optional_group);
		parser.add_flag("track", follow_camera, "Track camera settings", optional_group);
		parser.add_value("camera", name_camera, "Import camera settings", optional_group);
        
        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

int main(int argc, char* argv[]) {
    initialize(argc, (const char**)argv);
	
	VTK_SMART(vtkImageData) img = vtk_utils::load_image(name_in);
	
	double range[2];
	img->GetScalarRange(range);
	std::cout << "scalar range: " << range[0] << " -> " << range[1] << '\n';
	
	VTK_CREATE(vtkColorTransferFunction, ctf);
	vtk_utils::import_colormap(name_cmap_in, ctf); 
    
    VTK_CREATE(vtkDataSetMapper, mapper);
    mapper->SetInputData(img);
    mapper->ScalarVisibilityOn();
	mapper->SetLookupTable(ctf);
	
	VTK_CREATE(vtkActor, actor);
    actor->SetMapper(mapper);
    
    VTK_CREATE(vtkRenderer, renderer);
    renderer->SetBackground(0, 0, 0);
    renderer->AddActor(actor);
	
    VTK_CREATE(vtkRenderWindow, window);
    window->AddRenderer(renderer);
	window->SetSize(res[0], res[1]);
	
	if (follow_camera) {
		vtk_utils::track_camera_setting(renderer);
	}
	if (!name_camera.empty()) {
		vtk_utils::import_camera_settings(name_camera, renderer);
	}
    
    if (name_out.empty()) {
        // enter interactive mode
        VTK_CREATE(vtkRenderWindowInteractor, interactor);
        interactor->SetRenderWindow(window);
        interactor->Initialize();
        window->Render();
        interactor->Start();
    }
    else {
        window->Render();
        vtk_utils::save_frame(window, name_out, quality);
    }
    
    return 0;
}



