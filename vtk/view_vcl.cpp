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

#include <string>
#include <math/fixed_vector.hpp>
#include <VTK/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>

char* vcl_name;
float maxz, r, sp;
int minsize;

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
    hestOptAdd(&hopt, "i",      "file name",        airTypeString,  1, 1, &vcl_name,            NULL,       "Input file name (VCL format)");
    hestOptAdd(&hopt, "maxz",   "max z",            airTypeFloat,   0, 1, &maxz,                "0",        "max z threshold");
    hestOptAdd(&hopt, "minl",   "min length",       airTypeInt,     0, 1, &minsize,             "0",        "min curve length threshold");
    hestOptAdd(&hopt, "r",      "radius",           airTypeFloat,   0, 1, &r,                   "0.001",    "tube radius");
    hestOptAdd(&hopt, "sp",     "min spc",          airTypeFloat,   0, 1, &sp,                  "0",        "min distance between points");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize VCL curves in DPL dataset",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef std::vector<nvis::vec3> curve_type;

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    vtkCellArray* selected_curves = vtkCellArray::New();
    vtkDoubleArray* coords = vtkDoubleArray::New();
    coords->SetNumberOfComponents(3);
    
    int npts = 0;
    int ncurves = 0;
    
    std::fstream in(vcl_name, std::ios::in);
    nvis::vec3 last_read;
    curve_type curve;
    while (!in.eof()) {
        double x, y, z;
        char ch;
        
        in >> ch;
        
        if (in.eof()) {
            break;
        }
        
        switch (ch) {
            case 'p':
                in >> x >> y >> z;
                if (z > maxz) {
                    break;
                }
                last_read = nvis::vec3(x, y, z);
                if (!curve.size() || nvis::norm(last_read - curve.back()) > sp) {
                    curve.push_back(last_read);
                }
                break;
            case 'n':
                if (curve.size() && nvis::norm(last_read - curve.back()) > sp) {
                    curve.push_back(last_read);
                }
                if (curve.size() >= minsize) {
                    vtkIdType ids[curve.size()];
                    for (int i = 0 ; i < curve.size() ; ++i) {
                        coords->InsertNextTuple(curve[i].begin());
                        ids[i] = npts;
                        ++npts;
                    }
                    selected_curves->InsertNextCell(curve.size(), ids);
                }
                curve.clear();
                break;
        }
        
        std::string skipped;
        std::getline(in, skipped);
    }
    in.close();
    vtkPoints* pts = vtkPoints::New();
    pts->SetData(coords);
    
    std::cerr << selected_curves->GetNumberOfCells() << " lines passed the test\n";
    std::cerr << "longest curve has cardinal " << selected_curves->GetMaxCellSize() << '\n';
    
    vtkPolyData* polydata = vtkPolyData::New();
    polydata->SetPoints(pts);
    polydata->SetLines(selected_curves);
    
    vtkTubeFilter* tubes = vtkTubeFilter::New();
    tubes->SetInput(polydata);
    tubes->SetRadius(r);
    tubes->SetNumberOfSides(6);
    vtkPolyDataMapper* edge_mapper = vtkPolyDataMapper::New();
    edge_mapper->SetInputConnection(tubes->GetOutputPort());
    vtkActor* edge_actor = vtkActor::New();
    edge_actor->SetMapper(edge_mapper);
    edge_actor->GetProperty()->SetColor(1, 0, 0);
    
    vtkRenderer* ren = vtkRenderer::New();
    ren->SetUseDepthPeeling(1);
    ren->SetMaximumNumberOfPeels(100);
    ren->SetOcclusionRatio(0.1);
    ren->AddActor(edge_actor);
    ren->SetBackground(0, 0, 0);
    ren->ResetCamera();
    
    vtkRenderWindow* ren_win = vtkRenderWindow::New();
    ren_win->PointSmoothingOn();
    ren_win->LineSmoothingOn();
    ren_win->PolygonSmoothingOn();
    ren_win->SetAlphaBitPlanes(1);
    ren_win->SetMultiSamples(10);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1280, 800);
    
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(ren_win);
    ren_win->Render();
    iren->Initialize();
    iren->Start();
}





