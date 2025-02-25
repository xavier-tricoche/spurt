#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include <vtk/vtk_utils.hpp>
#include <vtk/vtk_camera_helper.hpp>
#include <misc/option_parse.hpp>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <format/filename.hpp>

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

typedef nvis::fixed_vector<double, 3> color_type;
typedef nvis::fixed_vector<double, 3> alpha_region;

std::string input_name, alpha_name, color_name, img_name, cam_in, cam_out;
std::string log_name="log.txt";
int grain_id=-1;
bool verbose=false;
bool clip=false;
color_type bg_col(0, 0, 0);
nvis::ivec2 res(1280, 800);
double dist=0.1;
double stddev=0;
double threshold = 10;
std::vector<double> alpha_tf;
std::vector<double> color_tf;
std::string input_alpha, input_color;
nvis::bbox3 bounds;
vtkSmartPointer<vtkImageData> image, filtered;
std::fstream log_file;
vtkSmartPointer<vtkRenderer> renderer;
vtkSmartPointer<vtkRenderWindow> window;
vtkSmartPointer<vtkColorTransferFunction> colors;
nvis::ivec3 minc, maxc;

// global domain / image contents information
double* domain;
vtkDoubleArray* scalars;
nvis::vec3 dmin, dmax;
nvis::ivec3 dims;
double volume;
nvis::vec3 center, up(-1, 0, 0);
double dfront, dside, dtop, dbottom;

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

void selected_region(nvis::ivec3& min, nvis::ivec3& max) {
	const nvis::vec3& minp = bounds.min();
	const nvis::vec3& maxp = bounds.max();
    int minid = image->FindPoint(minp[0], minp[1], minp[2]);
    int maxid = image->FindPoint(maxp[0], maxp[1], maxp[2]);
    min = id2c(minid, dims);
    max = id2c(maxid, dims);
}

double measure() {
    selected_region(minc, maxc);
	if (verbose) {
		std::cout << "minc=" << minc << ", maxc=" << maxc << '\n';
	}
    double sum=0;
    for (int k=minc[2]; k<=maxc[2]; ++k) {
        for (int j=minc[1]; j<=maxc[1]; ++j) {
            for (int i=minc[0]; i<=maxc[0]; ++i) {
                double value = scalars->GetTuple(i+dims[0]*(j+dims[1]*k))[0];
                if (value >= threshold) sum += value*volume;
            }
        }
    }
    return sum;
}


std::vector< VTK_SMART(Renderer) >               slice_renderers(3);
std::vector< VTK_SMART(RenderWindow) >           slice_windows(3);
std::vector< VTK_SMART(RenderWindowInteractor) > slice_interactors(3);
std::vector< VTK_SMART(ImageSliceMapper) >       slice_mappers(3);
std::vector< VTK_SMART(ImageSlice) >             slices(3);

int which_slice(int dim, double u) {
    double min = domain[2*dim];
    double max = domain[2*dim+1];
    double v = (u-min)/(max-min);
    return 1 + static_cast<int>(v*dims[dim]);
}

void initialize_slices() {
	for (int i=0; i<3; ++i) {
		slice_renderers[i] =   VTK_SMART(Renderer)::New();
		slice_windows[i]   =   VTK_SMART(RenderWindow)::New();
		slice_interactors[i] = VTK_SMART(RenderWindowInteractor)::New();
		slice_mappers[i] =     VTK_SMART(ImageSliceMapper)::New();
		slices[i] =            VTK_SMART(ImageSlice)::New();
		slice_windows[i]->AddRenderer(slice_renderers[i]);
		slice_interactors[i]->SetRenderWindow(slice_windows[i]);
		slice_windows[i]->SetSize(800, 800);
        slice_mappers[i]->SetInputData(filtered);
        slice_mappers[i]->SetOrientation(i);
        slice_mappers[i]->SetSliceNumber(which_slice(i, bounds.center()[i]));
        slice_interactors[i]->Initialize();
        slices[i]->SetMapper(slice_mappers[i]);
        slices[i]->GetProperty()->SetLookupTable(colors);
		slice_renderers[i]->AddViewProp(slices[i]);
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
        vtkImageSliceMapper::SafeDownCast(slices[i]->GetMapper())->SetSliceNumber(which_slice(i, bounds.center()[i]));
        slice_mappers[i]->SetCroppingRegion(minc[0], maxc[0], minc[1], maxc[1], minc[2], maxc[2]);
        slice_mappers[i]->CroppingOn();
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

// Callback for moving the planes from the box widget to the mapper
class vtkBoxWidgetCallback : public vtkCommand
{
public:
    static vtkBoxWidgetCallback *New()
        { return new vtkBoxWidgetCallback; }
    void Execute(vtkObject *caller, unsigned long, void*) VTK_OVERRIDE
    {
        vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
        if (this->Mapper)
        {
            VTK_CREATE(vtkPlanes, planes);
            widget->GetPlanes(planes);

            double box[6];

            for (int i=0; i<planes->GetNumberOfPlanes(); ++i) {
                box[i] = planes->GetPlane(i)->GetOrigin()[i/2];
    				// nvis::vec3 normal = vec3(a_plane->GetNormal());
            }
    	    bounds.min() = nvis::vec3(box[0], box[2], box[4]);
    	    bounds.max() = nvis::vec3(box[1], box[3], box[5]);

			std::cout << "region: " << bounds.min() << " -> " << bounds.max() << '\n';

            selected_region(minc, maxc);

            this->Mapper->SetClippingPlanes(planes);
			this->Mapper->Update();

			update_slices();
			window->Render();
        }
    }
  void SetMapper(vtkSmartVolumeMapper* m)
    { this->Mapper = m; }

protected:
  vtkBoxWidgetCallback()
    { this->Mapper = 0; }

  vtkSmartVolumeMapper *Mapper;
};

void KeypressCB ( vtkObject* caller, long unsigned int vtkNotUsed(eventId), void* vtkNotUsed(clientData), void* vtkNotUsed(callData) )
{
    vtkRenderWindowInteractor *iren =
    static_cast<vtkRenderWindowInteractor*>(caller);
    if (*iren->GetKeySym() == 'b') {
        std::cout << "bounds: " << bounds.min() << " -> " << bounds.max() << '\n';
    }
    else if (*iren->GetKeySym() == 'm') {
        std::cout << "measure: " << measure() << '\n';
    }
    else if (*iren->GetKeySym() == 's') {
        std::cout << "Enter name for measure: ";
        std::string name;
        std::cin >> name;
        double m = measure();
        if (verbose) {
            std::cout <<  "file: " << input_name << ", grain: " << grain_id
				<< ", bounds: " << bounds.min() << "-" << bounds.max()
				<< ", what :" << name << ", intensity :" << m << std::endl;
        }
        else std::cout << name << " measure: " << m << '\n';
        if (!log_file.fail()) {
            // save to csv file...
            log_file << "file: " << input_name << ", grain: " << grain_id
				<< ", bounds: " << bounds.min() << "-" << bounds.max()
				<< ", what :" << name << ", intensity :" << m << std::endl;
        }
    }
    else if (*iren->GetKeySym() == 'g') {
        std::cout << "Enter current grain ID: ";
        std::cin >> grain_id;
    }
    else if (*iren->GetKeySym() == 'c') {
		vtk_utils::export_camera_settings("", renderer);
    }
	else if (*iren->GetKeySym() == 'f') {
		renderer->GetActiveCamera()->SetViewUp(&up[0]);
		renderer->GetActiveCamera()->SetFocalPoint(&center[0]);
		renderer->GetActiveCamera()->SetPosition(&(center + dfront*nvis::vec3(0, 0, 1))[0]);
		window->Render();
	}
	else if (*iren->GetKeySym() == 'r') {
		nvis::vec3 _up(-1, 0, 0);
		renderer->GetActiveCamera()->SetViewUp(&_up[0]);
		renderer->GetActiveCamera()->SetFocalPoint(&center[0]);
		renderer->GetActiveCamera()->SetPosition(&(center + dside*nvis::vec3(0, 1, 0))[0]);
		window->Render();
	}
	else if (*iren->GetKeySym() == 'b') {
		renderer->GetActiveCamera()->SetViewUp(&up[0]);
		renderer->GetActiveCamera()->SetFocalPoint(&center[0]);
		renderer->GetActiveCamera()->SetPosition(&(center + dside*nvis::vec3(1, 0, 0))[0]);
		window->Render();
	}
	else if (*iren->GetKeySym() == 't') {
		nvis::vec3 _up(0, 0, -1);
		renderer->GetActiveCamera()->SetViewUp(&_up[0]);
		renderer->GetActiveCamera()->SetFocalPoint(&center[0]);
		renderer->GetActiveCamera()->SetPosition(&(center + dtop*nvis::vec3(-1, 0, 0))[0]);
		window->Render();
	}
	else if (*iren->GetKeySym() == 'l') {
		nvis::vec3 _up(-1, 0, 0);
		renderer->GetActiveCamera()->SetViewUp(&_up[0]);
		renderer->GetActiveCamera()->SetFocalPoint(&center[0]);
		renderer->GetActiveCamera()->SetPosition(&(center + dtop*nvis::vec3(0, -1, 0))[0]);
		window->Render();
	}

    std::cout << "Pressed: " << iren->GetKeySym() << std::endl;
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
        parser.add_value("alpha", alpha_name,
		"Transfer function filename or inline definition", optional_group);
        parser.add_value("color", color_name,
						 "Transfer function filename or inline definition",
						  optional_group);
		parser.add_value("smooth", stddev, stddev, "Apply gaussian blurring (std dev)", optional_group);
        parser.add_value("min", threshold, threshold, "Min threshold for signal", optional_group);
        parser.add_value("log", log_name, log_name, "Log filename for measures", optional_group);
        parser.add_value("save", img_name, "Snapshot filename",
		optional_group);
        parser.add_value("step", dist, dist, "Sampling step size",
		optional_group);
        parser.add_value("camin", cam_in, "Camera input filename",
		optional_group);
        parser.add_value("camout", cam_out, "Camera output filename",
		optional_group);
        parser.add_value("bg", bg_col, bg_col, "Background color",
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

void read_alpha(std::vector< std::pair<double, double> >& cps, const std::string& filename) {
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
                    cps.push_back(std::pair<double, double>(val, alpha));
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
                    cps.push_back(std::pair<double, double>(val-0.5*width, 0));
                    cps.push_back(std::pair<double, double>(val, alpha));
                    cps.push_back(std::pair<double, double>(val+0.5*width, 0));
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
        cps.clear();
        std::pair< double, double> point;
        while (!input.eof()) {
            input >> point.first >> point.second;
            cps.push_back(point);
        }
        input.close();
    }
    if (verbose) {
        std::cout << cps.size()
        << " alpha control points successfully created\n";
    }
}

void read_color(std::vector< std::pair<double, color_type> >& cps, const std::string filename) {
    // check to see if alpha tf is defined inline
    size_t where = filename.find_first_of(" ;,");
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
                cps.push_back(std::make_pair(val, col));
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
        cps.clear();
        std::pair<double, color_type> point;
        while (!input.eof()) {
            input >> point.first
            >> point.second[0] >> point.second[1] >> point.second[2];
            cps.push_back(point);
        }
        input.close();
    }
    if (verbose) {
        std::cout << cps.size()
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
	if (verbose) {
		std::cout << "domain: " << dmin << " -> " << dmax << '\n';
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
		image->PrintSelf(std::cout, vtkIndent());
	}
	else if (ext == "nrrd" || ext == "nhdr") {
		VTK_CREATE(vtkNrrdReader, reader);
		reader->SetFileName(input_name.c_str());
		reader->Update();
    	image = vtkImageData::SafeDownCast(reader->GetOutput());
		std::cout << "just imported:\n";
		image->PrintSelf(std::cout, vtkIndent());
	}
	set_parameters();

    double* b = image->GetBounds();
    bounds.min() = nvis::vec3(b[0], b[2], b[4]);
    bounds.max() = nvis::vec3(b[1], b[3], b[5]);
    std::cout << "global bounds: " << bounds.min() << " -> " << bounds.max() << '\n';
	nvis::vec3 save_min = bounds.min();
	nvis::vec3 save_max = bounds.max();
	bounds.min() += 0.01 * (save_max - save_min);
	bounds.max() -= 0.01 * (save_max - save_min);

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

    double* valrange=filtered->GetPointData()->GetScalars()->GetRange();
	if (verbose)
		std::cout << "value range: "
			<< valrange[0] << " -> "
			<< valrange[1] << '\n';

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

    VTK_CREATE(vtkPiecewiseFunction, alpha_tf);
    if (!alpha_name.empty()) {
        std::vector< std::pair<double, double> > cps;
        read_alpha(cps, alpha_name);
        if (converted) {
            for (int i=0; i<cps.size(); ++i) {
                double newval=(cps[i].first-
					valrange[0])/(valrange[1]-valrange[0])*(1 << 16);
                alpha_tf->AddPoint(newval, cps[i].second);
				std::cout << "control point: " << newval << ", " << cps[i].second << '\n';
            }
        }
        else {
            for (int i=0; i<cps.size(); ++i) {
                alpha_tf->AddPoint(cps[i].first, cps[i].second);
				std::cout << "control point: " << cps[i].first << ", " << cps[i].second << '\n';
            }
        }
    }
    else {
        std::cerr << "Warning: no alpha transfer function provided: "
			<< "using default one\n";
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
        std::vector< std::pair<double, color_type> > cps;
        read_color(cps, color_name);
        if (converted) {
            for (int i=0; i<cps.size(); ++i) {
                double newval=(cps[i].first-valrange[0])/
					(valrange[1]-valrange[0])*(1 << 16);
                colors->AddRGBPoint(newval,
				cps[i].second[0], cps[i].second[1], cps[i].second[2]);
            }
        }
        else {
            for (int i=0; i<cps.size(); ++i) {
                colors->AddRGBPoint(cps[i].first,
				cps[i].second[0], cps[i].second[1], cps[i].second[2]);
            }
        }
    }
    else {
        std::cerr << "Warning: no color transfer function provided: "
			<< "using default one\n";
        if (converted) {
            colors->AddRGBPoint(0, 1, 1, 1);
            colors->AddRGBPoint(1<<16, 1, 1, 1);
        }
        else {
            colors->AddRGBPoint(valrange[0], 1, 1, 1);
            colors->AddRGBPoint(valrange[1], 1, 1, 1);
        }
    }

	initialize_slices();

    VTK_CREATE(vtkVolume, volume);
    VTK_CREATE(vtkSmartVolumeMapper, mapper);
	mapper->SetBlendModeToComposite();
    mapper->SetSampleDistance(dist);
    mapper->SetAutoAdjustSampleDistances(1);
    renderer = vtkRenderer::New();
	window = vtkRenderWindow::New();
    VTK_CREATE(vtkRenderWindowInteractor, interactor);

    // Set connections
    mapper->SetInputData(filtered);
    volume->SetMapper(mapper);
    renderer->AddViewProp(volume);
    window->AddRenderer(renderer);
	window->SetSize(res[0], res[1]);
    interactor->SetRenderWindow(window);

    renderer->GetActiveCamera()->ParallelProjectionOn();

    // Now set the opacity and the color
    vtkVolumeProperty* vol_prop = volume->GetProperty();
    vol_prop->SetIndependentComponents(1);
    vol_prop->SetScalarOpacity(0, alpha_tf);
    vol_prop->SetColor(0, colors);
	vol_prop->ShadeOn();
	volume->SetProperty(vol_prop);

	// Add a box widget if the clip option was selected
	vtkBoxWidget *box = vtkBoxWidget::New();
	if (clip) {
        log_file.open(log_name.c_str(), std::ios::out | std::ios::app );
		box->SetInteractor(interactor);
		box->SetPlaceFactor(0.99);
		box->SetInputData(filtered);

		box->SetDefaultRenderer(renderer);
		box->InsideOutOn();
		box->PlaceWidget();
		vtkBoxWidgetCallback *callback = vtkBoxWidgetCallback::New();
		callback->SetMapper(mapper);
		box->AddObserver(vtkCommand::InteractionEvent, callback);
		callback->Delete();
		box->EnabledOn();
		box->GetSelectedFaceProperty()->SetOpacity(0.0);
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
