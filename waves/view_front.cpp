#include "vtkPolyDataReader.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkPolyDataNormals.h"
#include "vtkCommand.h"
#include "vtkInteractorObserver.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkTransform.h"
#include "vtkCylinderSource.h"
#include "vtkWindowToImageFilter.h"
#include "vtkTIFFWriter.h"
#include "vtkDataSetReader.h"
#include "vtkTubeFilter.h"
#include "vtkUnstructuredGrid.h"
#include "vtkCellData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkCubeSource.h"
#include "vtkContourFilter.h"
#include "vtkStructuredPointsReader.h"
#include "vtkScalarBarActor.h"
#include "vtkColorTransferFunction.h"
#include "vtkTextProperty.h"
#include "vtkStructuredPoints.h"
#include "vtkConeSource.h"
#include "vtkCylinderSource.h"

#include <string>
#include <math/fixed_vector.hpp>
#include <VTK/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>
#include "boost/shared_ptr.hpp"
#include "data/field_wrapper.hpp"


void wait(int s)
{
    nvis::timer t;
    while (t.elapsed() < s) {}
}

char* name_in;
float length;

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1, 1, &name_in,             NULL,       "input file (NRRD)");
    hestOptAdd(&hopt, "l",      "length",           airTypeFloat,   0, 1, &length,              "0.1",      "vector length");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize location and direction of gaussian wave packets",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
	initialize(argc, argv);
	
	Nrrd* nin = xavier::readNrrd(name_in);
	std::vector<float> data;
	xavier::to_vector<float>(data, nin);
	
	size_t npts = data.size()/3;
	
	vtkSmartPointer<vtkFloatArray> coords = vtkSmartPointer<vtkFloatArray>::New();
	coords->SetNumberOfComponents(3);
	coords->SetNumberOfTuples(npts);
	for (int i = 0 ; i < npts ; ++i) {
		nvis::fvec3 x(data[3*i], data[3*i+1], 0);
		coords->SetTuple(i, &x[0]);
	}
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	points->SetData(coords);
	
	vtkSmartPointer<vtkFloatArray> vals = vtkSmartPointer<vtkFloatArray>::New();
	vals->SetNumberOfComponents(3);
	vals->SetNumberOfTuples(npts);
	for (int i = 0 ; i < npts ; ++i) {
		float theta = data[3*i+2];
		nvis::fvec3 d(cos(theta), sin(theta), 0);
		coords->SetTuple(i, &d[0]);
	}

	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata->SetPoints(points);
	polydata->GetPointData()->SetVectors(vals);
	
	vtkSmartPointer<vtkConeSource> cone = vtkSmartPointer<vtkConeSource>::New();
	cone->SetHeight(0.3);
	cone->SetRadius(0.2);
	cone->SetResolution(6);
	vtkSmartPointer<vtkCylinderSource> cylinder = vtkSmartPointer<vtkCylinderSource>::New();
	cylinder->SetHeight(1.0);
	cylinder->SetRadius(0.05);
	cylinder->SetResolution(6);
	
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
#if VTK_MAJOR_VERSION >= 6
	mapper->SetInputData(polydata);
#else
    mapper->SetInput(polydata);
#endif
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1,0,0);
    
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    ren->SetUseDepthPeeling(1);
    ren->SetMaximumNumberOfPeels(100);
    ren->SetOcclusionRatio(0.1);
    ren->SetBackground(0, 0, 0);
    
    ren->AddActor(actor);
    ren->ResetCamera();
    
    vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    ren_win->SetAlphaBitPlanes(1);
    ren_win->SetMultiSamples(0);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(640, 480);
    ren_win->Render();
    
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(ren_win);
    iren->Start();
}































































