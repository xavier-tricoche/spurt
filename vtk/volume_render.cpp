/*=========================================================================

  Program:   Visualization Toolkit
  Module:    GPURenderDemo.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// VTK includes
#include "vtkBoxWidget.h"
#include "vtkCamera.h"
#include "vtkCommand.h"
#include "vtkColorTransferFunction.h"
#include "vtkDataArray.h"
#include "vtkDICOMImageReader.h"
#include "vtkImageData.h"
#include "vtkImageResample.h"
#include "vtkImageShiftScale.h"
#include "vtkIndent.h"
#include "vtkMetaImageReader.h"
#include "vtkNrrdReader.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPlanes.h"
#include "vtkPointData.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkXMLImageDataReader.h"
#include "vtkSmartVolumeMapper.h"

#include <format/filename.hpp>
#include <math/fixed_vector.hpp>

#include <sstream>
#include <iterator>

// Callback for moving the planes from the box widget to the mapper
class vtkBoxWidgetCallback : public vtkCommand
{
public:
  static vtkBoxWidgetCallback *New()
    { return new vtkBoxWidgetCallback; }
  void Execute(vtkObject *caller, unsigned long, void*) override
  {
      vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
      if (this->Mapper)
      {
        vtkPlanes *planes = vtkPlanes::New();
        widget->GetPlanes(planes);
        this->Mapper->SetClippingPlanes(planes);
        planes->Delete();
      }
  }
  void SetMapper(vtkSmartVolumeMapper* m)
    { this->Mapper = m; }

protected:
  vtkBoxWidgetCallback()
    { this->Mapper = nullptr; }

  vtkSmartVolumeMapper *Mapper;
};

std::string me;
void printUsage(const std::string& msg = "") {
    if (!msg.empty())
        std::cerr << "ERROR: " << msg << '\n';
    std::cout <<  "USAGE: " << me << " [parameters] [options]\n";
    std::cout <<  "DESCRIPTION: Volume rendering of scalar dataset\n";
    std::cout <<  "PARAMETERS:\n";
    std::cout <<  " -i | --input <string> Input filename\n";
    std::cout <<  "OPTIONS:\n";
    std::cout <<  " -h | --help                    Print this information\n";
    std::cout <<  " -g | --gradient <string>       Filename to export gradient magnitude\n";
    std::cout <<  " -d | --distance <float>        Sampling distance along ray (default: 0->5)\n";
    std::cout <<  " -r | --res <int> (x2)          Render window resolution (default: 800 800)\n";
    std::cout <<  " -b | --background <float> (x3) Set background color (default: black)\n";
    std::cout <<  " -t | --tent <int> (x3)         Add scalar opacity tent (center, width, height)\n";
    std::cout <<  " -o | --opacity <float> (x2)    Add opacity control point (val, alpha)\n";
    std::cout <<  "      --opacity-file <string>   Opacity transfer function file\n";
    std::cout <<  " -c | --color <float> (x4)      Add RGB color control point (val, R, G, B)\n";
    std::cout <<  "      --color-file <string>     Color transfer function file\n";
	std::cout <<  "    | --clip                    Toggle clipping box\n";
	std::cout <<  " -q | --quantize <int>          Convert dataset to 8 or 16 bits before rendering\n";
    std::cout <<  " -v | --verbose                 Turn on verbose mode\n";
    std::cout <<  "      --no-gpu                  Disable GPU rendering\n";
    std::cout <<  "      --no-shade                Disable shading\n";
    std::cout <<  "      --no-color-bar            Do not display color bar\n";
    std::cout <<  "      --no-linear               Do not use trilinear interpolation\n";
    exit(0);
}

typedef nvis::fvec2 opacity_control_point;
typedef nvis::fvec4 color_control_point;

nvis::fvec3 bg_color(0,0,0);
nvis::ivec2 res(800, 800);
double distance=0;
bool verbose = false;
bool gpu = true;
bool shade = true;
bool clip = false;
int maxval = 0;

double convert_value(double value, double range[]) {
	if (maxval > 0)
		return (value-range[0])/(range[1]-range[0])*static_cast<double>(maxval);
	else return value;
}

double convert_window(double window, double range[]) {
	if (maxval > 0)
		return window/(range[1]-range[0])*static_cast<double>(maxval);
	else return window;
}

int main(int argc, char *argv[])
{
	std::string filename;
	std::vector<opacity_control_point> opacity_tf;
	std::vector<color_control_point> color_tf;
	std::string arg;
	me = argv[0];
	bool distance_is_relative = false;
	for (int i=1; i<argc; ++i) {
		arg = argv[i];
		if (arg == "-i" || arg == "--input") {
			filename = argv[++i];
		}
		else if (arg == "-d" || arg == "--distance") {
			arg = argv[++i];
			if (arg.back() == '%') {
				distance_is_relative = true;
				distance = std::stof(arg.substr(0, arg.size()-1));
			}
			else {
				distance_is_relative = false;
				distance = std::stof(arg);
			}
		}
		else if (arg == "-r" || arg == "--resolution") {
			res[0] = std::stoi(argv[++i]);
			res[1] = std::stoi(argv[++i]);
		}
		else if (arg == "-b" || arg == "--background") {
			bg_color[0] = std::stof(argv[++i]);
			bg_color[1] = std::stof(argv[++i]);
			bg_color[2] = std::stof(argv[++i]);
		}
		else if (arg == "-t" || arg == "--tent") {
			float center = std::stof(argv[++i]);
			float width = std::stof(argv[++i]);
			float alpha = std::stof(argv[++i]);
			opacity_tf.push_back(opacity_control_point(center-width, 0));
			opacity_tf.push_back(opacity_control_point(center, alpha));
			opacity_tf.push_back(opacity_control_point(center+width, 0));
		}
		else if (arg == "-o" || arg == "--opacity") {
			float val = std::stof(argv[++i]);
			float alpha = std::stof(argv[++i]);
			opacity_tf.push_back(opacity_control_point(val, alpha));
		}
		else if (arg == "-c" || arg == "--color") {
			color_tf.push_back(color_control_point());
			color_control_point& col = color_tf.back();
			col[0] = std::stof(argv[++i]);
			col[1] = std::stof(argv[++i]);
			col[2] = std::stof(argv[++i]);
			col[3] = std::stof(argv[++i]);
		}
		else if (arg == "-q" || arg == "--quantize") {
			int nbits = std::stoi(argv[++i]);
			if (nbits > 16 || nbits<=0) nbits = 16;
			maxval = 2 << nbits;
		}
		else if (arg == "-v" || arg == "--verbose") {
			verbose = true;
		}
		else if (arg == "--no-gpu") {
			gpu = false;
		}
		else if (arg == "--no-shade") {
			shade = false;
		}
		else if (arg == "-h" || arg == "--help") {
			printUsage();
		}
		else if (arg == "--clip") {
			clip = true;
		}
		else {
			printUsage("Invalid option");
		}
	}

	if (filename.empty()) {
		printUsage("Missing input filename");
	}
	else if (verbose) {
		std::cout << "filename = " << filename << '\n';
		std::cout << "background color = " << bg_color << '\n';
		std::cout << "opacity transfer function: ";
		std::copy(opacity_tf.begin(), opacity_tf.end(), std::ostream_iterator<opacity_control_point>(std::cout, ", "));
		std::cout << "\ncolor transfer function: ";
		std::copy(color_tf.begin(), color_tf.end(), std::ostream_iterator<color_control_point>(std::cout, ", "));
		std::cout << "\ndistance=" << distance << '\n';
		std::cout << "resolution: " << res << '\n';
	}

	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->AddRenderer(renderer);
	renderer->SetBackground(bg_color[0], bg_color[1], bg_color[2]);

	vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renWin);

	iren->GetInteractorStyle()->SetDefaultRenderer(renderer);

	// Read the data
	vtkSmartPointer<vtkImageData> input = vtkSmartPointer<vtkImageData>::New();

	std::string ext = spurt::filename::extension(filename);

	if (ext == "nrrd" || ext == "NRRD" || ext == "nhdr")
	{
		vtkSmartPointer<vtkNrrdReader> reader = vtkSmartPointer<vtkNrrdReader>::New();
		reader->SetFileName(filename.c_str());
		reader->Update();
		input->ShallowCopy(reader->GetOutput());
	}
	else if (ext == "vti" || ext == "VTI")
	{
		vtkSmartPointer<vtkXMLImageDataReader> reader = vtkSmartPointer<vtkXMLImageDataReader>::New();
		reader->SetFileName(filename.c_str());
		reader->Update();
		input->ShallowCopy(reader->GetOutput());
	}
	else if (ext == "vtk" || ext == "VTK")
	{
	    vtkSmartPointer<vtkStructuredPointsReader> reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
	    reader->SetFileName(filename.c_str());
	    reader->Update();
	    input->ShallowCopy(static_cast<vtkImageData*>(reader->GetOutput()));
	}
	else
	{
		printUsage("Unrecognized file type in input");
	}

	if (verbose) {
		std::cout << "input image: ";
		input->PrintSelf(std::cout, vtkIndent(4));
	}

	double range[2];
	input->GetPointData()->GetScalars()->GetRange(range);
	if (verbose) {
		std::cout << "Value range: " << range[0] << ", " << range[1] << '\n';
	}

	// Create our volume and mapper
	vtkSmartPointer<vtkVolume> volume = vtkSmartPointer<vtkVolume>::New();
	vtkSmartPointer<vtkSmartVolumeMapper> mapper = vtkSmartPointer<vtkSmartVolumeMapper>::New();

	if (maxval > 0) {
		vtkSmartPointer<vtkImageShiftScale> shift = vtkSmartPointer<vtkImageShiftScale>::New();
		shift->SetInputData(input);
 		shift->SetShift(-range[0]);
		double magnitude = range[1]-range[0];
		if(magnitude==0.0){  magnitude=1.0; }
		shift->SetScale(static_cast<double>(maxval)/magnitude);
		if (maxval < 256)
			shift->SetOutputScalarTypeToUnsignedChar();
		else
			shift->SetOutputScalarTypeToUnsignedShort();
		shift->Update();
		mapper->SetInputConnection(shift->GetOutputPort());
	}
	else mapper->SetInputData(input);

	// Add a box widget if the clip option was selected
	vtkSmartPointer<vtkBoxWidget> box = vtkSmartPointer<vtkBoxWidget>::New();
	if (clip)
	{
		box->SetInteractor(iren);
		box->SetPlaceFactor(1.01);
		box->SetInputData(input);
		box->SetDefaultRenderer(renderer);
		box->InsideOutOn();
		box->PlaceWidget();
		vtkSmartPointer<vtkBoxWidgetCallback> callback = vtkSmartPointer<vtkBoxWidgetCallback>::New();
		callback->SetMapper(mapper);
		box->AddObserver(vtkCommand::InteractionEvent, callback);
		box->EnabledOn();
		box->GetSelectedFaceProperty()->SetOpacity(0.0);
	}

	// Set the sample distance on the ray to be 1/2 the average spacing
	if (distance_is_relative || distance == 0) {
		double spacing[3];
		input->GetSpacing(spacing);
		double mind = *std::min_element(spacing, spacing+3);
		if (distance == 0) distance = 0.25*mind;
		else distance *= mind/100.;
	}

	mapper->SetSampleDistance(distance);
	if (gpu) mapper->SetRequestedRenderModeToGPU();
	else mapper->SetRequestedRenderModeToRayCast();

	// Create our transfer function
	vtkSmartPointer<vtkColorTransferFunction> colorFun = vtkSmartPointer<vtkColorTransferFunction>::New();
	vtkSmartPointer<vtkPiecewiseFunction> opacityFun = vtkSmartPointer<vtkPiecewiseFunction>::New();

	for (size_t i=0; i<color_tf.size(); ++i) {
		const color_control_point& col = color_tf[i];
		colorFun->AddRGBPoint(convert_value(col[0], range), col[1], col[2], col[3]);
	}

	for (size_t i=0; i<opacity_tf.size(); ++i) {
		const opacity_control_point& op = opacity_tf[i];
		opacityFun->AddPoint(convert_value(op[0], range), op[1]);
	}

	if (verbose) {
		std::cout << "color tf: ";
		colorFun->PrintSelf(std::cout, vtkIndent(4));
		std::cout << "\nopacity tf: ";
		opacityFun->PrintSelf(std::cout, vtkIndent(4));
	}

	// Create the property and attach the transfer functions
	vtkSmartPointer<vtkVolumeProperty> property = vtkSmartPointer<vtkVolumeProperty>::New();
	// property->SetIndependentComponents(independentComponents);
	property->SetColor( colorFun );
	property->SetScalarOpacity( opacityFun );
	property->SetInterpolationTypeToLinear();
	if (!shade) {
		property->ShadeOff();
	}
	else {
		property->ShadeOn();
	}

	// connect up the volume to the property and the mapper
	volume->SetProperty( property );
	volume->SetMapper( mapper );

	// Set the default window size
	renWin->SetSize(res[0],res[1]);
	renWin->Render();

	// Add the volume to the scene
	renderer->AddVolume( volume );

	renderer->ResetCamera();

	// interact with data
	renWin->Render();

	iren->Start();

	return 0;
}
