#include <string>
#include <vector>
#include <fstream>
#include <memory>
#include <ctime>

#include <vtk/vtk_utils.hpp>
#include <vtk/vtk_camera_helper.hpp>
#include <misc/option_parse.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <format/filename.hpp>
#include <math/stat.hpp>

#include <vtkColorTransferFunction.h>
#include <vtkDataSet.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkSmartVolumeMapper.h>
#include <vtkPiecewiseFunction.h>
#include <vtkVolume.h>
#include <vtkVolumeProperty.h>
#include <vtkCommand.h>
#include <vtkBoxWidget.h>
#include <vtkPlanes.h>
#include <vtkCallbackCommand.h>
#include <vtkCommand.h>
#include <vtkImageData.h>
#include <vtkNrrdReader.h>
#include <vtkClipVolume.h>
#include <vtkImageSliceMapper.h>
#include <vtkImageSlice.h>
#include <vtkImageProperty.h>
#include <vtkContourFilter.h>
#include <vtkCutter.h>
#include <vtkBox.h>
#include <vtkClipPolyData.h>

typedef nvis::fixed_vector<double, 3> color_type;
typedef nvis::fixed_vector<double, 3> alpha_region;

std::string input_name, alpha_name, color_name, img_name, cam_in, cam_out;
std::string log_name="log.txt";
int grain_id=-1;
bool verbose=false;
bool clip=false;
bool use_parallel=false;
bool save_and_quit=false;
color_type bg_col(0, 0, 0);
nvis::ivec2 res(1280, 800);
double dist=0.1;
double stddev=0;
double threshold_min = 5;
double threshold_max = 100.;
double specular = 0;
std::vector<double> alpha_tf;
std::vector<double> color_tf;
std::string input_alpha, input_color;
nvis::bbox3 bounds, global_bounds;
vtkSmartPointer<vtkImageData> image, filtered;
std::fstream log_file;
vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderWindow> window;
vtkSmartPointer<vtkColorTransferFunction> colors;
nvis::ivec3 minc, maxc;
std::vector< std::pair<double, double> > alpha_cps;
std::vector< std::pair<double, color_type > > color_cps;
std::string bounds_as_string;

// global domain / image contents information
double* domain;
vtkDoubleArray* scalars;
nvis::vec3 dmin, dmax;
nvis::ivec3 dims;
double volume;
nvis::vec3 center, up(-1, 0, 0);
double dfront, dside, dtop, dbottom;
double* valrange;

inline nvis::vec3 vec3(double* v) {
	return nvis::vec3(v[0], v[1], v[2]);
}

inline nvis::ivec3 ivec3(int* v) {
    return nvis::ivec3(v[0], v[1], v[2]);
}

inline nvis::ivec3 id2c(int id, const nvis::ivec3& dims) {
    nvis::ivec3 coord;
    coord[0] = id % dims[0];
    id /= dims[0];
    coord[1] = id % dims[1];
    coord[2] = id / dims[1];
    return coord;
}

inline std::string bounds_to_str(const nvis::bbox3& b) {
	std::ostringstream os;
	os << b.min() << " -> " << b.max();
	return os.str();
}

std::vector< VTK_SMART(vtkRenderer) >               slice_renderers(3);
std::vector< VTK_SMART(vtkRenderWindow) >           slice_windows(3);
std::vector< VTK_SMART(vtkRenderWindowInteractor) > slice_interactors(3);
std::vector< VTK_SMART(vtkPolyDataMapper) >         slice_image_mappers(3);
std::vector< VTK_SMART(vtkCutter) >                 slice_cutters(3);
std::vector< VTK_SMART(vtkContourFilter) >          slice_contour_filters(3);
std::vector< VTK_SMART(vtkClipPolyData) >           slice_box_clippers(3);

VTK_SMART(vtkBoxWidget) box_widget;
VTK_SMART(vtkSmartVolumeMapper) vr_mapper;

double isovalue;

// int which_slice(int dim, double u) {
//     double min = domain[2*dim];
//     double max = domain[2*dim+1];
//     double v = (u-min)/(max-min);
//     return 1 + static_cast<int>(v*dims[dim]);
// }

bool selected_region(nvis::ivec3& min, nvis::ivec3& max) {
	// make sampling region fit within dataset bounds
	nvis::vec3 _min = nvis::max(bounds.min(), global_bounds.min());
	nvis::vec3 _max = nvis::min(bounds.max(), global_bounds.max());
	if (nvis::any(_min > _max)) {
		std::cerr << "invalid sampling region: " << bounds_to_str(bounds) << '\n';
		return false;
	}
    int minid = image->FindPoint(_min[0], _min[1], _min[2]);
    int maxid = image->FindPoint(_max[0], _max[1], _max[2]);
	int tmp = minid;
	minid = std::min(minid, maxid);
	maxid = std::max(tmp, maxid);
	// std::cout << "minid=" << minid << ", maxid=" << maxid << '\n';
    min = id2c(minid, dims);
    max = id2c(maxid, dims);
	return true;
}

double measure() {
	double minv = std::max(isovalue, threshold_min);
    if (!selected_region(minc, maxc)) {
    	return 0;
    }
	if (verbose) {
		std::cout << "sampled region: " << bounds_to_str(bounds)  << '\n';
		std::cout << "measure interval: [" << minv << ", " << threshold_max << "]\n";
		std::cout << "minc=" << minc << ", maxc=" << maxc << '\n';
		std::cout << "overall domain: " << bounds_to_str(global_bounds) << '\n';
		std::cout << "voxel volume:   " << volume << '\n';
	}
	std::vector<double> all_values, selected_values;
    double sum=0;
	VTK_CREATE(vtkPlanes, planes);
	box_widget->GetPlanes(planes);
	size_t ntotal=0;
	size_t nbad=0;
    for (int k=minc[2]; k<=maxc[2]; ++k) {
        for (int j=minc[1]; j<=maxc[1]; ++j) {
            for (int i=minc[0]; i<=maxc[0]; ++i, ++ntotal) {
				vtkIdType id = i+dims[0]*(j+dims[1]*k);
                double value = scalars->GetTuple(id)[0];
				all_values.push_back(value);
                if (value >= minv && value <= threshold_max) {
					selected_values.push_back(value);
					sum += value*volume;
				}
            }
        }
    }
	if (verbose) std::cout << "sampling (" << all_values.size() << ") done\n";
	std::pair<double, double> mv = spurt::meanvariance(selected_values);
	std::pair<double, double> mv_total = spurt::meanvariance(all_values);
	if (verbose) {
		std::cout << all_values.size() << " total samples\n";
		std::cout << selected_values.size() << " selected samples ("
			      << 100.*selected_values.size() / all_values.size() << "%)\n";
		std::cout << "total mean / variance within selected region: "
				  << mv_total.first << " / " << mv_total.second << '\n';
		std::cout << "following results for selected value range x region selection:\n";
		std::cout << "sum / mean / variance / bounds:\n";
		std::cout << sum << '\n';
		std::cout << mv.first << '\n';
		std::cout << mv.second << '\n';
		std::cout << bounds_to_str(bounds) << '\n';
 		// std::cout << ntotal << " values sampled of which " << nbad << " where outside implicit region\n";
	}
    return sum;
}

void initialize_slices() {
	VTK_CREATE(vtkColorTransferFunction, ctf);
	if (!color_cps.empty()) {
		for (int i=0; i<color_cps.size() ; ++i) {
			ctf->AddRGBPoint(color_cps[i].first,
							 color_cps[i].second[0], color_cps[i].second[1], color_cps[i].second[2]);
		}
	}
	else {
		ctf->AddRGBPoint(valrange[0], 0, 0, 1);
		ctf->AddRGBPoint(valrange[1], 1, 1, 0);
	}
	for (int i=0; i<3; ++i) {
		slice_renderers[i] =        VTK_SMART(vtkRenderer)::New();
		slice_windows[i]   =        VTK_SMART(vtkRenderWindow)::New();
		slice_interactors[i] =      VTK_SMART(vtkRenderWindowInteractor)::New();
		slice_cutters[i] =          VTK_SMART(vtkCutter)::New();
		slice_contour_filters[i] =  VTK_SMART(vtkContourFilter)::New();
		slice_box_clippers[i] =     VTK_SMART(vtkClipPolyData)::New();

		VTK_CREATE(vtkBox, box);
		box->SetBounds(bounds.min()[0], bounds.max()[0],
					   bounds.min()[1], bounds.max()[1],
					   bounds.min()[2], bounds.max()[2]);
		slice_box_clippers[i]->InsideOutOn();
		slice_box_clippers[i]->SetClipFunction(box);
		slice_box_clippers[i]->SetInputConnection(slice_cutters[i]->GetOutputPort());

		slice_cutters[i]->SetInputData(image);
		VTK_CREATE(vtkPlane, plane);
		plane->SetOrigin(&bounds.center()[0]);
		nvis::vec3 normal(0);
		normal[i] = 1.;
		plane->SetNormal(&normal[0]);
		slice_cutters[i]->SetCutFunction(plane);
		VTK_CREATE(vtkPolyDataMapper, image_mapper);
		image_mapper->SetInputConnection(slice_box_clippers[i]->GetOutputPort());
		image_mapper->SetLookupTable(ctf);
		VTK_CREATE(vtkActor, image_actor);
		image_actor->SetMapper(image_mapper);
		slice_renderers[i]->AddActor(image_actor);

		slice_contour_filters[i]->SetInputConnection(slice_box_clippers[i]->GetOutputPort());
		VTK_CREATE(vtkPolyDataMapper, iso_mapper);
		iso_mapper->SetInputConnection(slice_contour_filters[i]->GetOutputPort());
		iso_mapper->ScalarVisibilityOff();
		VTK_CREATE(vtkActor, iso_actor);
		iso_actor->SetMapper(iso_mapper);
		iso_actor->GetProperty()->SetColor(0.5,0.5,0.5);
		iso_actor->GetProperty()->SetLineWidth(3);
		slice_renderers[i]->AddActor(iso_actor);

		// image -> cut -> color map
		// image -> cut -> contour filter

		slice_windows[i]->AddRenderer(slice_renderers[i]);
		slice_interactors[i]->SetRenderWindow(slice_windows[i]);
		slice_windows[i]->SetSize(800, 800);
        slice_interactors[i]->Initialize();
		std::ostringstream os;
		switch (i) {
			case 0: os << "X="; break;
			case 1: os << "Y="; break;
			case 2: os << "Z="; break;
		}
		os << bounds.center()[i];
		slice_windows[i]->SetWindowName(os.str().c_str());
        slice_renderers[i]->GetActiveCamera()->SetFocalPoint(&bounds.center()[0]);
	}
    slice_renderers[2]->GetActiveCamera()->SetViewUp(-1, 0, 0);
    slice_renderers[0]->GetActiveCamera()->SetViewUp(0, 0, -1);
    slice_renderers[1]->GetActiveCamera()->SetViewUp(-1, 0, 0);
    slice_renderers[0]->GetActiveCamera()->SetPosition(&(bounds.center()+100.*nvis::vec3(1, 0, 0))[0]);
    slice_renderers[1]->GetActiveCamera()->SetPosition(&(bounds.center()+100.*nvis::vec3(0, 1, 0))[0]);
    slice_renderers[2]->GetActiveCamera()->SetPosition(&(bounds.center()+100.*nvis::vec3(0, 0, 1))[0]);
    slice_windows[0]->Render();
    slice_windows[1]->Render();
    slice_windows[2]->Render();
}

void update_slices() {
	nvis::vec3 c = bounds.center();
	for (int i=0; i<3; ++i) {
        vtkPlane::SafeDownCast(slice_cutters[i]->GetCutFunction())->SetOrigin(&c[0]);
		slice_cutters[i]->Update();
		vtkBox::SafeDownCast(slice_box_clippers[i]->GetClipFunction())->SetBounds(
			bounds.min()[0], bounds.max()[0],
			bounds.min()[1], bounds.max()[1],
			bounds.min()[2], bounds.max()[2]);
		slice_box_clippers[i]->Update();
		slice_contour_filters[i]->Update();
		std::ostringstream os;
		switch (i) {
			case 0: os << "X="; break;
			case 1: os << "Y="; break;
			case 2: os << "Z="; break;
		}
		os << c[i];
		slice_windows[i]->SetWindowName(os.str().c_str());
		slice_windows[i]->Render();
	}
}

void update_isosurface() {
	if (verbose) std::cout << "Setting isovalue to " << isovalue << '\n';
	for (int i=0; i<3; ++i) {
		slice_contour_filters[i]->SetValue(0, isovalue);
		slice_contour_filters[i]->Update();
		slice_windows[i]->Render();
	}
}

void reset_camera_and_render(const nvis::vec3& focal, const nvis::vec3& dist, const nvis::vec3& up) {
	renderer->GetActiveCamera()->SetViewUp(&up[0]);
	renderer->GetActiveCamera()->SetFocalPoint(&focal[0]);
	renderer->GetActiveCamera()->SetPosition(&(focal+dist)[0]);
	window->Render();
}

void set_clipping_planes(int crop=0) { // 0: no crop, 1: inside out crop, 2: inside crop
	vr_mapper->SetCroppingRegionPlanes(bounds.min()[0], bounds.max()[0],
					  				   bounds.min()[1], bounds.max()[1],
					                   bounds.min()[2], bounds.max()[2]);
	if (crop == 0) {
		   vr_mapper->SetCroppingRegionFlags(134217727);
	}
	else if (crop == 1) {
		vr_mapper->SetCroppingRegionFlags(8192);
	}
	else if (crop == 2) {
		vr_mapper->SetCroppingRegionFlags(134209535);
	}
	else if (crop == 3) {
		vr_mapper->SetCroppingRegionFlagsToInvertedFence();
	}
	else if (crop == 4) {
		vr_mapper->SetCroppingRegionFlagsToCross();
	}
	else if (crop == 5) {
		vr_mapper->SetCroppingRegionFlagsToInvertedCross();
	}
}

void update_box(vtkBoxWidget*);

std::ostream& bold_on(std::ostream& os) {
	return os << "\e[1m";
}
std::ostream& bold_off(std::ostream& os) {
	return os << "\e[0m";
}

void handle_event(const std::string& what, bool quiet=false) {
	if (what == "question" || what == "h") {
		std::cout << "\nKeyboard codes:\n"
			<< bold_on << " * Widget control\n" << bold_off
			<< "\tb:   print bounds\n"
			<< "\tx:   (re)enable box widget\n"
			<< "\tz:   reset widget shape, orientation, and position\n"
			<< bold_on << " * Visibility control\n" << bold_off
			<< "\t-:   clip outside of box widget\n"
			<< "\t\\:   crop outside of box widget\n"
			<< "\t/:   crop inside of box widget\n"
			<< "\t=:   restore full visibility\n"
			<< bold_on << " * Level set control\n" << bold_off
			<< "\ti:   set isovalue\n"
			<< bold_on << " * Camera control\n" << bold_off
			<< "\tc:   print camera settings\n"
            << "\ts:   save current frame to file\n"
			<< "\tf:   place camera in front\n"
			<< "\tl:   place camera to the left\n"
			<< "\tr:   place camera to the right\n"
			<< "\tt:   place camera above\n"
			<< "\tu:   place camera underneath\n\n";
	}
	else if (what == "b") {
		std::cout << "bounds: " << bounds.min() << " -> " << bounds.max() << '\n';
	}
	else if (what == "m") {
		std::cout << "measure: " << measure() << '\n';
	}
	else if (what == "c") {
		vtk_utils::export_camera_settings("", renderer);
	}
    else if (what == "s") {
        std::srand(std::time(0));
        int r = std::rand();
        std::string outname = "snapshot" + std::to_string(r) + ".png";
        vtk_utils::save_frame(window, outname);
        if (verbose) {
            std::cout << "Current frame has been exported to " << outname << '\n';
        }
    }
	else if (what == "f") {
		reset_camera_and_render(center, dfront*nvis::vec3(0, 0, 1), up);
	}
	else if (what == "r") {
		reset_camera_and_render(center, dside*nvis::vec3(0, 1, 0), nvis::vec3(-1, 0, 0));
	}
	else if (what == "u") {
		reset_camera_and_render(center, dtop*nvis::vec3(1, 0, 0), nvis::vec3(0, 0, 1));
	}
	else if (what == "t") {
		reset_camera_and_render(center, dtop*nvis::vec3(-1, 0, 0), nvis::vec3(0, 0, -1));
	}
	else if (what == "l") {
		reset_camera_and_render(center, dside*nvis::vec3(0, -1, 0), nvis::vec3(-1, 0, 0));
	}
	else if (what == "z") {
		box_widget->PlaceWidget(); // restore axis aligned orientation
	}
	else if (what == "i") {
		std::cout << "Enter isovalue: ";
		std::cin >> isovalue;
		update_isosurface();
	}
	else if (what == "x") {
		box_widget->EnabledOn();
	}
	else if (what == "equal") {
		set_clipping_planes(0);
		std::cout << "cropping flags: " << vr_mapper->GetCroppingRegionFlags() << '\n';
		vr_mapper->CroppingOff();
		vr_mapper->Update();
		window->Render();
	}
	else if (what == "minus" || what == "backslash") {
		set_clipping_planes(1);
		std::cout << "cropping flags: " << vr_mapper->GetCroppingRegionFlags() << '\n';
		vr_mapper->CroppingOn();
		vr_mapper->Update();
		window->Render();
	}
	else if (what == "slash") {
		set_clipping_planes(2);
		std::cout << "cropping flags: " << vr_mapper->GetCroppingRegionFlags() << '\n';
		vr_mapper->CroppingOn();
		vr_mapper->Update();
		window->Render();
	}
	else if (what == "bracketleft") {
		set_clipping_planes(3);
		std::cout << "cropping flags: " << vr_mapper->GetCroppingRegionFlags() << '\n';
		vr_mapper->CroppingOn();
		vr_mapper->Update();
		window->Render();
	}
	else if (what == "bracketright") {
		set_clipping_planes(4);
		std::cout << "cropping flags: " << vr_mapper->GetCroppingRegionFlags() << '\n';
		vr_mapper->CroppingOn();
		vr_mapper->Update();
		window->Render();
	}
	else if (what == "semicolon") {
		set_clipping_planes(5);
		std::cout << "cropping flags: " << vr_mapper->GetCroppingRegionFlags() << '\n';
		vr_mapper->CroppingOn();
		vr_mapper->Update();
		window->Render();
	}
	else if (what == "Left") {
		double width = bounds.size()[0];
		VTK_CREATE(vtkTransform, t);
		box_widget->GetTransform(t);
		t->Translate(-0.05*width, 0, 0);
		box_widget->SetTransform(t);
		update_box(box_widget);
		window->Render();
	}
	else if (what == "Right") {
		double width = bounds.size()[0];
		VTK_CREATE(vtkTransform, t);
		box_widget->GetTransform(t);
		t->Translate(0.05*width, 0, 0);
		box_widget->SetTransform(t);
		update_box(box_widget);
		window->Render();
	}
	else if (what == "Up") {
		double width = bounds.size()[1];
		VTK_CREATE(vtkTransform, t);
		box_widget->GetTransform(t);
		t->Translate(0, 0.05*width, 0);
		box_widget->SetTransform(t);
		update_box(box_widget);
		window->Render();
	}
	else if (what == "Down") {
		double width = bounds.size()[1];
		VTK_CREATE(vtkTransform, t);
		box_widget->GetTransform(t);
		t->Translate(0, -0.05*width, 0);
		box_widget->SetTransform(t);
		update_box(box_widget);
		window->Render();
	}
	else if (what == "9") {
		double width = bounds.size()[2];
		VTK_CREATE(vtkTransform, t);
		box_widget->GetTransform(t);
		t->Translate(0, 0, -0.05*width);
		box_widget->SetTransform(t);
		update_box(box_widget);
		window->Render();
	}
	else if (what == "0") {
		double width = bounds.size()[2];
		VTK_CREATE(vtkTransform, t);
		box_widget->GetTransform(t);
		t->Translate(0, 0, 0.05*width);
		box_widget->SetTransform(t);
		update_box(box_widget);
		window->Render();
	}
	else {
		if (!quiet)
			std::cout << "WARNING: unrecognized event character code: " << what << '\n';
	}
}

void update_bounds(vtkBoxWidget *widget) {
    VTK_CREATE(vtkPlanes, planes);
    widget->GetPlanes(planes);

    double box[6];

    for (int i=0; i<planes->GetNumberOfPlanes(); ++i) {
        box[i] = planes->GetPlane(i)->GetOrigin()[i/2];
			// nvis::vec3 normal = vec3(a_plane->GetNormal());
    }
    bounds.min() = nvis::vec3(box[0], box[2], box[4]);
    bounds.max() = nvis::vec3(box[1], box[3], box[5]);
}

void update_box(vtkBoxWidget *widget) {
	update_bounds(widget);
    selected_region(minc, maxc);
	update_slices();
}

// Callback for moving the planes from the box widget to the mapper
class vtkBoxWidgetCallback : public vtkCommand {
public:
    static vtkBoxWidgetCallback *New()
        { return new vtkBoxWidgetCallback; }
    void Execute(vtkObject *caller, unsigned long, void*)
    {
        vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
		update_box(widget);
		window->Render();
    }
};

void KeypressCB ( vtkObject* caller, long unsigned int vtkNotUsed(eventId),
				  void* vtkNotUsed(clientData), void* vtkNotUsed(callData) ) {
	vtkRenderWindowInteractor *iren =
    	static_cast<vtkRenderWindowInteractor*>(caller);
	if (verbose) {
    	std::cout << "Pressed: " << iren->GetKeySym() << std::endl;
	}

	handle_event(iren->GetKeySym());
}

void MouseRightClickCB ( vtkObject* caller,
						 long unsigned int vtkNotUsed(eventId),
						 void* vtkNotUsed(clientData),
						 void* vtkNotUsed(callData) ) {
	vtkRenderWindowInteractor* iren =
		static_cast<vtkRenderWindowInteractor*>(caller);
	if (verbose) {
    	std::cout << "Pressed: right-click / " << iren->GetKeySym() << std::endl;
	}
	handle_event(iren->GetKeySym());
}

nvis::bbox3 str_to_bbox(const std::string& str) {
    std::string copy(str);
    for (int i=0; i<str.size(); ++i) {
        if (str[i] == '[' || str[i] == ']' || str[i] == ',') {
            copy[i] = ' ';
        }
        else if (i<str.size()-1 && str[i] == '-' && str[i+1] == '>') {
            copy[i++] = ' ';
            copy[i] = ' ';
        }
    }
    nvis::bbox3 box;
    std::istringstream iss(copy);
    for (int i=0; i<6; ++i) {
        if (i<3) iss >> box.min()[i];
        else iss >> box.max()[i-3];
        if (iss.fail()) {
            std::cerr << "unable to parse bounding box information.\n";
            exit(1);
        }
    }
    return box;
}

void initialize(int argc, char* argv[]) {
    namespace xcl = spurt::command_line;

    xcl::option_traits
            required_group(true, false, "Required Options"),
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "GPU-based volume rendering of scalar volume");

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("input", input_name, "Input filename",
		required_group);
        parser.add_value("bounds", bounds_as_string, "selected region", optional_group);
        parser.add_value("alpha", alpha_name,
		"Transfer function filename or inline definition", optional_group);
        parser.add_value("color", color_name,
						 "Transfer function filename or inline definition",
						  optional_group);
		parser.add_value("specular", specular, specular, "Specular coefficient", optional_group);
		parser.add_value("smooth", stddev, stddev, "Apply gaussian blurring (std dev)", optional_group);
        parser.add_value("min", threshold_min, threshold_min, "Min threshold for signal", optional_group);
        parser.add_value("max", threshold_max, threshold_max, "Max threshold for signal", optional_group);
        parser.add_value("log", log_name, log_name, "Log filename for measures", optional_group);
        parser.add_value("save", img_name, "Snapshot filename",
		optional_group);
        parser.add_value("quit", save_and_quit, "Quit after saving frame", optional_group);
        parser.add_value("step", dist, dist, "Sampling step size",
		optional_group);
        parser.add_value("camin", cam_in, "Camera input filename",
		optional_group);
        parser.add_value("camout", cam_out, "Camera output filename",
		optional_group);
        parser.add_value("parallel", use_parallel, "Use parallel projection", optional_group);
        parser.add_tuple<3>("bg", bg_col, bg_col, "Background color",
		optional_group);
        parser.add_tuple<2>("res", res, res, "Image resolution",
		optional_group);
				parser.add_value("clip", clip, clip, "Use clipping planes",
			optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbose output",
		optional_group);

        parser.parse(argc, const_cast<const char**>(argv));
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception: "
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

void read_alpha(const std::string& filename) {
	alpha_cps.clear();
    // check to see if alpha tf is defined inline
    size_t where = filename.find_first_of(" ;,");
    if (where != std::string::npos) {
        if (verbose) std::cout << "alpha input is not a filename" << std::endl;
        std::string keyword = filename.substr(0, where);
        std::string args = filename.substr(where+1);
        std::cout << "keyword=" << keyword << ", args=" << args << '\n';
        for (int i=0; i<args.size(); ++i) {
            if (args[i] == '(' || args[i] == ')' || args[i] == ','
            || args[i] == ';' || args[i] == '-') {
                args[i] = ' ';
            }
        }
        std::istringstream iss(args);
        if (keyword == "points") {
            // parse control points <value, alpha> pairs
            double val, alpha;
            while (!iss.eof()) {
                iss >> val >> alpha;
                if (!iss.fail()) {
                    alpha_cps.push_back(std::pair<double, double>(val, alpha));
                    if (verbose)
                    std::cout << "read (" << val << ", "
                    << alpha  << ")" << std::endl;
                }
                else break;
            }
        }
        else if (keyword == "tents") {
            // parse alpha tents
            double val, width, alpha;
            while (!iss.eof()) {
                iss >> val >> width >> alpha;
                if (!iss.fail()) {
                    alpha_cps.push_back(std::pair<double, double>(val-0.5*width, 0));
                    alpha_cps.push_back(std::pair<double, double>(val, alpha));
                    alpha_cps.push_back(std::pair<double, double>(val+0.5*width, 0));
                    if (verbose)
                    std::cout << "read (" << val << ", " << width << ", "
                    << alpha  << ")" << std::endl;
                }
                else {
                    std::cout << "unable to parse tent info\n";
                    break;
                }
            }
        }
        else {
            std::cerr << "unrecognized alpha xf keyword: " << keyword << '\n';
            exit(1);
        }
    }
    else {
        std::fstream input(filename.c_str(), std::ios::in);
        if (!input) {
            throw std::runtime_error("Unable to open " + filename +
            ": no such file or directory");
        }
        alpha_cps.clear();
        std::pair< double, double> point;
        while (!input.eof()) {
            input >> point.first >> point.second;
            alpha_cps.push_back(point);
        }
        input.close();
    }
    if (verbose) {
        std::cout << alpha_cps.size()
        << " alpha control points successfully created\n";
    }
}

void read_color(const std::string filename) {
    // check to see if alpha tf is defined inline
    size_t where = filename.find_first_of(" ;,");
	color_cps.clear();
    if (where != std::string::npos) {
        if (verbose) std::cout << "color input is not a filename" << std::endl;
        std::string args = filename;
        std::cout << "args=" << args << '\n';
        for (int i=0; i<args.size(); ++i) {
            if (args[i] == '(' || args[i] == ')' || args[i] == ','
                || args[i] == ';' || args[i] == '-') {
                args[i] = ' ';
            }
        }
        std::istringstream iss(args);
        // parse control points <value, alpha> pairs
        double val;
        color_type col;
        while (!iss.eof()) {
            iss >> val >> col[0] >> col[1] >> col[2];
            if (!iss.fail()) {
                color_cps.push_back(std::make_pair(val, col));
                if (verbose)
                std::cout << "read (" << val << ", "
                << col << ")" << std::endl;
            }
            else break;
        }
    }
    else {
        std::fstream input(filename.c_str(), std::ios::in);
        if (!input) {
            throw std::runtime_error("Unable to open " + filename +
            ": no such file or directory");
        }
        color_cps.clear();
        std::pair<double, color_type> point;
        while (!input.eof()) {
            input >> point.first
            >> point.second[0] >> point.second[1] >> point.second[2];
            color_cps.push_back(point);
        }
        input.close();
    }
    if (verbose) {
        std::cout << color_cps.size()
        << " color control points successfully created\n";
    }
}

void set_parameters() {
	scalars =
		vtkDoubleArray::SafeDownCast(image->GetPointData()->GetScalars());

	image->ComputeBounds();
	domain = image->GetBounds();
	dmin = nvis::vec3(domain[0], domain[2], domain[4]);
	dmax = nvis::vec3(domain[1], domain[3], domain[5]);
	global_bounds.min() = dmin;
	global_bounds.max() = dmax;
	if (verbose) {
		std::cout << "domain: " << bounds_to_str(global_bounds) << '\n';
	}
    dims = ivec3(image->GetDimensions());
	if (verbose) {
		std::cout << "image dimensions = " << dims << '\n';
	}
	double* spacing = image->GetSpacing();
	volume = spacing[0] * spacing[1] * spacing[2];
	center = 0.5*(dmin + dmax);
	dfront = dside = dtop = dbottom = 100.;
}

void raycast() {
	std::string ext = spurt::filename::extension(input_name);
	if (ext == "vtk") {
    	VTK_CREATE(vtkDataSetReader, reader);
    	reader->SetFileName(input_name.c_str());
    	reader->Update();
    	image = vtkImageData::SafeDownCast(reader->GetOutput());
		std::cout << "just imported:\n";
		// image->PrintSelf(std::cout, vtkIndent());
	}
	else if (ext == "nrrd" || ext == "nhdr") {
		VTK_CREATE(vtkNrrdReader, reader);
		reader->SetFileName(input_name.c_str());
		reader->Update();
    	image = vtkImageData::SafeDownCast(reader->GetOutput());
		std::cout << "just imported:\n";
		// image->PrintSelf(std::cout, vtkIndent());
	}
	set_parameters();

    double* b = image->GetBounds();

    if (!bounds_as_string.empty()) {
        bounds = str_to_bbox(bounds_as_string);
    }
    else {
        bounds.min() = nvis::vec3(b[0], b[2], b[4]);
        bounds.max() = nvis::vec3(b[1], b[3], b[5]);
        std::cout << "global bounds: " << bounds_to_str(bounds) << '\n';
	    nvis::vec3 save_min = bounds.min();
	    nvis::vec3 save_max = bounds.max();
	    bounds.min() += 0.01 * (save_max - save_min);
	    bounds.max() -= 0.01 * (save_max - save_min);
    }

	if (stddev > 0) {
		VTK_CREATE(vtkImageGaussianSmooth, smooth);
		smooth->SetInputData(image);
		smooth->SetStandardDeviation(stddev, stddev, stddev);
		smooth->SetDimensionality(3);
		smooth->Update();
		filtered = smooth->GetOutput();
	}
	else {
		filtered = vtkSmartPointer<vtkImageData>::New();
		filtered->DeepCopy(image);
	}
    filtered->GetDimensions(&dims[0]);

	isovalue = -1;

    valrange=filtered->GetPointData()->GetScalars()->GetRange();
	if (verbose)
		std::cout << "value range: "
			<< valrange[0] << " -> "
			<< valrange[1] << '\n';
    double deltaval = valrange[1]-valrange[0];

    // convert dataset to unsigned short if necessary
    int input_type = filtered->GetPointData()->GetScalars()->GetDataType();
    vtkSmartPointer<vtkDataArray> new_data;
    switch(input_type) {
        case VTK_CHAR:
        case VTK_UNSIGNED_CHAR:
        case VTK_SHORT:
        case VTK_UNSIGNED_SHORT: break;
        default: {
            new_data=vtk_utils::implicit_scale<unsigned short>(
                filtered->GetPointData()->GetScalars(),
                std::make_pair(valrange[0], valrange[1]),
                std::make_pair(0., static_cast<double>(1 << 16)));
            break;
        }
    }
    bool converted = new_data.Get() != NULL;
    if (converted) {
        filtered->GetPointData()->SetScalars(new_data);
    }
	double *newvalrange=filtered->GetPointData()->GetScalars()->GetRange();
	if (verbose)
		std::cout << "after conversion: value range: "
			<< newvalrange[0] << " -> "
			<< newvalrange[1] << '\n';

    VTK_CREATE(vtkPiecewiseFunction, alpha_tf);
    if (!alpha_name.empty()) {
        size_t maxushort = (1 << 16);
        double maxv = static_cast<double>(maxushort);

        read_alpha(alpha_name);
        if (converted) {
            for (int i=0; i<alpha_cps.size(); ++i) {
                double u = (alpha_cps[i].first-valrange[0])/deltaval;
                std::cout << "cp=" << alpha_cps[i].first << '\n';
                std::cout << "u=" << u << '\n';
                u*=maxv;
                std::cout << "newval=" << u << '\n';
                alpha_tf->AddPoint(u, alpha_cps[i].second);
				std::cout << "control point: " << u << ", " << alpha_cps[i].second << '\n';
            }
        }
        else {
            for (int i=0; i<alpha_cps.size(); ++i) {
                alpha_tf->AddPoint(alpha_cps[i].first, alpha_cps[i].second);
				std::cout << "control point: " << alpha_cps[i].first << ", " << alpha_cps[i].second << '\n';
            }
        }
    }
    else {
        std::cerr << "Warning: no alpha transfer function provided: "
			<< "using default ramp up one\n";
        if (converted) {
            alpha_tf->AddPoint(0, 0);
            alpha_tf->AddPoint(1 << 16, 1);
        }
        else {
            alpha_tf->AddPoint(valrange[0], 0);
            alpha_tf->AddPoint(valrange[1], 1);
        }
    }

    colors = vtkSmartPointer<vtkColorTransferFunction>::New();
    if (!color_name.empty()) {
        read_color(color_name);
        if (converted) {
            for (int i=0; i<color_cps.size(); ++i) {
                double newval=(color_cps[i].first-valrange[0])/
					deltaval*(1 << 16);
                std::cout << "color cp: " << color_cps[i].first << '\n';
                std::cout << "color newval: " << newval << '\n';
                colors->AddRGBPoint(newval,
				color_cps[i].second[0], color_cps[i].second[1], color_cps[i].second[2]);
            }
        }
        else {
            for (int i=0; i<color_cps.size(); ++i) {
                colors->AddRGBPoint(color_cps[i].first,
				color_cps[i].second[0], color_cps[i].second[1], color_cps[i].second[2]);
            }
        }
    }
    else {
        std::cerr << "Warning: no color transfer function provided: "
			<< "using default blue to yellow one\n";
        if (converted) {
            colors->AddRGBPoint(0, 0, 0, 1);
            colors->AddRGBPoint(1<<16, 1, 1, 0);
        }
        else {
            colors->AddRGBPoint(valrange[0], 0, 0, 1);
            colors->AddRGBPoint(valrange[1], 1, 1, 0);
        }
    }

	initialize_slices();

    VTK_CREATE(vtkVolume, volume);
    vr_mapper = VTK_SMART(vtkSmartVolumeMapper)::New();
	vr_mapper->SetBlendModeToComposite();
    vr_mapper->SetSampleDistance(dist);
    vr_mapper->SetAutoAdjustSampleDistances(0);
    vr_mapper->SetRequestedRenderModeToOSPRay();
    renderer = vtkRenderer::New();
	window = vtkRenderWindow::New();
    VTK_CREATE(vtkRenderWindowInteractor, interactor);

    // Set connections
    vr_mapper->SetInputData(filtered);
    volume->SetMapper(vr_mapper);
    renderer->AddViewProp(volume);
    window->AddRenderer(renderer);
	window->SetSize(res[0], res[1]);
    interactor->SetRenderWindow(window);

	// renderer->AddActor(isoactor);

    if (use_parallel) {
        renderer->GetActiveCamera()->ParallelProjectionOn();
    }

    // Now set the opacity and the color
    vtkVolumeProperty* vol_prop = volume->GetProperty();
    vol_prop->SetIndependentComponents(1);
 // 	vol_prop->SetInterpolationTypeToLinear();
    vol_prop->SetScalarOpacity(0, alpha_tf);
    vol_prop->SetColor(0, colors);
	if (specular > 0) vol_prop->SetSpecular(specular);
	vol_prop->ShadeOn();
    vol_prop->SetAmbient(0.1);
    vol_prop->SetDiffuse(0.9);
    vol_prop->SetSpecular(0.2);
    vol_prop->SetSpecularPower(10.0);
	volume->SetProperty(vol_prop);

	// Add a box widget if the clip option was selected
	box_widget = VTK_SMART(vtkBoxWidget)::New();
	if (clip) {
        log_file.open(log_name.c_str(), std::ios::out | std::ios::app );
		box_widget->SetInteractor(interactor);
		box_widget->SetPlaceFactor(0.99);
		box_widget->SetInputData(filtered);

		box_widget->SetDefaultRenderer(renderer);
		box_widget->InsideOutOn();
		box_widget->PlaceWidget();
		vtkBoxWidgetCallback *callback = vtkBoxWidgetCallback::New();
		// callback->SetMapper(vr_mapper);
		box_widget->AddObserver(vtkCommand::InteractionEvent, callback);
		callback->Delete();
		box_widget->EnabledOn();
		box_widget->GetSelectedFaceProperty()->SetOpacity(0.0);
	}
    VTK_CREATE(vtkCallbackCommand, keycmd);
    keycmd->SetCallback(KeypressCB);
    interactor->AddObserver( vtkCommand::KeyPressEvent, keycmd);

    interactor->Initialize();
    renderer->SetBackground(bg_col[0], bg_col[1], bg_col[2]);
    if (!cam_in.empty()) {
        vtk_utils::import_camera_settings(cam_in, renderer);
    }
    else {
		renderer->ResetCamera();
		renderer->GetActiveCamera()->SetViewUp(&up[0]);
	}
    window->Render();
    if (!img_name.empty()) {
        vtk_utils::save_frame(window, img_name);
        if (verbose) {
            std::cout << "Current frame has been saved to " << img_name << '\n';
        }
        if (save_and_quit) exit(0);
    }
    interactor->Start();

    if (!cam_out.empty()) {
        vtk_utils::export_camera_settings(cam_out, renderer);
    }
}

int main(int argc, char* argv[]) {
    initialize(argc, argv);

    raycast();

    return 0;
}
