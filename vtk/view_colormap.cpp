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
#include "vtkPNGWriter.h"
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
#include "vtkSphereSource.h"
#include "vtkArrowSource.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkClipPolyData.h"
#include "vtkPlane.h"

#include <string>
#include <math/fixed_vector.hpp>
#include <VTK/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>
#include "Garcia_vis_helper.hpp"

bool    param_v;
float   param_g;
char*   param_f;
int     param_k;
float   param_b;

void wait(int s)
{
    nvis::timer t;
    while (t.elapsed() < s) {}
}

template<typename T>
inline T sign(T a)
{
    return (a < 0 ? -1 : 1);
}

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
    hestOptAdd(&hopt, "k",      "kind",             airTypeInt,     0, 1, &param_k,     "0",    "color scale kind (0: RGS, 1:RGB, 2:GRS, 3: GRB, 4:BYS, 5:YBS");
    hestOptAdd(&hopt, "g",      "gamma",            airTypeFloat,   0, 1, &param_g,     "1",    "color scale gamma factor");
    hestOptAdd(&hopt, "b",      "brightness",       airTypeFloat,   0, 1, &param_b,     "1",    "color scale brightness factor");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &param_v,     "0",    "verbose mode (debugging)");
    hestOptAdd(&hopt, "f",      "image file",       airTypeString,  1, 1, &param_f,     NULL,   "output screenshot");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize color map",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    using namespace Garcia_vis_helper;
    
    vtkSmartPointer<vtkSphereSource> sphere = vtkSmartPointer<vtkSphereSource>::New();
    sphere->SetThetaResolution(100);
    sphere->SetPhiResolution(100);
    
    vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
    plane->SetOrigin(0, 0, 0);
    plane->SetNormal(0, 0, 1);
    
    vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
    clipper->SetInputConnection(sphere->GetOutputPort());
    clipper->SetClipFunction(plane);
    clipper->GenerateClippedOutputOn();
    clipper->SetValue(0);
    clipper->Update();
    
    vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::New();
    mesh->DeepCopy(clipper->GetOutput());
    
    vtkSmartPointer<vtkUnsignedCharArray> color = vtkSmartPointer<vtkUnsignedCharArray>::New();
    color->SetNumberOfComponents(3);
    color->SetName("Colors");
    
    for (int i = 0 ; i < mesh->GetNumberOfPoints() ; ++i) {
        nvis::vec3 x;
        mesh->GetPoint(i, x.begin());
        uchar_color_type c = to_color(upcol(x, param_b, param_g));
        color->InsertNextTypedTuple(c.begin());
    }
    mesh->GetPointData()->SetScalars(color);
    
    vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
    normals->SetInputData(mesh);
    
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(normals->GetOutputPort());
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    ren->SetUseDepthPeeling(1);
    ren->SetMaximumNumberOfPeels(100);
    ren->SetOcclusionRatio(0.1);
    ren->SetBackground(1, 1, 1);
    ren->ResetCamera();
    
    ren->GetActiveCamera()->SetPosition(2.71106, -2.91483, 2.24705);
    ren->GetActiveCamera()->SetFocalPoint(0.0155195, -0.0228576, -0.0470498);
    ren->GetActiveCamera()->SetViewUp(-0.348574, 0.361191, 0.864892);
    ren->GetActiveCamera()->SetClippingRange(1.09368, 8.95884);
    
    if (param_v) {
        vtk_utils::camera_setting_callback* cb = vtk_utils::camera_setting_callback::New();
        ren->AddObserver(vtkCommand::StartEvent, cb);
        cb->Delete();
    }
    
    vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    ren_win->SetAlphaBitPlanes(1);
    ren_win->SetMultiSamples(0);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1600, 1200);
    
    vtkSmartPointer<vtkArrowSource> arrow = vtkSmartPointer<vtkArrowSource>::New();
    arrow->SetTipResolution(20);
    arrow->SetShaftResolution(20);
    arrow->SetTipRadius(0.06);
    arrow->SetTipLength(0.2);
    arrow->SetShaftRadius(0.02);
    
    vtkSmartPointer<vtkPolyDataNormals> normals2 = vtkSmartPointer<vtkPolyDataNormals>::New();
    normals2->SetInputConnection(arrow->GetOutputPort());
    
    nvis::vec3 col;
    
    vtkSmartPointer<vtkTransform> transformX = vtkSmartPointer<vtkTransform>::New();
    // transformX->Scale(1.01, 1.01, 1.01);
    vtkSmartPointer<vtkTransformPolyDataFilter> tpdX = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    tpdX->SetInputConnection(normals2->GetOutputPort());
    tpdX->SetTransform(transformX);
    
    vtkSmartPointer<vtkPolyDataMapper> mapperX = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperX->SetInputConnection(tpdX->GetOutputPort());
    vtkSmartPointer<vtkActor> actorX = vtkSmartPointer<vtkActor>::New();
    col = upcol(nvis::vec3(1, 0, 0), param_b, param_g);
    actorX->GetProperty()->SetColor(col[0], col[1], col[2]);
    actorX->SetMapper(mapperX);
    
    vtkSmartPointer<vtkTransform> transformY = vtkSmartPointer<vtkTransform>::New();
    // transformY->Scale(1.25, 1.25, 1.25);
    transformY->RotateZ(90);
    vtkSmartPointer<vtkTransformPolyDataFilter> tpdY = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    tpdY->SetInputConnection(normals2->GetOutputPort());
    tpdY->SetTransform(transformY);
    
    vtkSmartPointer<vtkPolyDataMapper> mapperY = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperY->SetInputConnection(tpdY->GetOutputPort());
    vtkSmartPointer<vtkActor> actorY = vtkSmartPointer<vtkActor>::New();
    col = upcol(nvis::vec3(0, 1, 0), param_b, param_g);
    actorY->GetProperty()->SetColor(col[0], col[1], col[2]);
    actorY->SetMapper(mapperY);
    
    vtkSmartPointer<vtkTransform> transformZ = vtkSmartPointer<vtkTransform>::New();
    // transformZ->Scale(1.25, 1.25, 1.25);
    transformZ->RotateY(-90);
    vtkSmartPointer<vtkTransformPolyDataFilter> tpdZ = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    tpdZ->SetInputConnection(normals2->GetOutputPort());
    tpdZ->SetTransform(transformZ);
    
    vtkSmartPointer<vtkPolyDataMapper> mapperZ = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperZ->SetInputConnection(tpdZ->GetOutputPort());
    vtkSmartPointer<vtkActor> actorZ = vtkSmartPointer<vtkActor>::New();
    col = upcol(nvis::vec3(0, 0, 1), param_b, param_g);
    actorZ->GetProperty()->SetColor(col[0], col[1], col[2]);
    actorZ->SetMapper(mapperZ);
    
    //
    
    vtkSmartPointer<vtkTransform> transformXX = vtkSmartPointer<vtkTransform>::New();
    // transformX->Scale(1.01, 1.01, 1.01);
    transformXX->RotateZ(180);
    vtkSmartPointer<vtkTransformPolyDataFilter> tpdXX = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    tpdXX->SetInputConnection(normals2->GetOutputPort());
    tpdXX->SetTransform(transformXX);
    
    vtkSmartPointer<vtkPolyDataMapper> mapperXX = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperXX->SetInputConnection(tpdXX->GetOutputPort());
    vtkSmartPointer<vtkActor> actorXX = vtkSmartPointer<vtkActor>::New();
    col = upcol(nvis::vec3(-1, 0, 0), param_b, param_g);
    actorXX->GetProperty()->SetColor(col[0], col[1], col[2]);
    actorXX->SetMapper(mapperXX);
    
    vtkSmartPointer<vtkTransform> transformYY = vtkSmartPointer<vtkTransform>::New();
    // transformY->Scale(1.25, 1.25, 1.25);
    transformYY->RotateZ(-90);
    vtkSmartPointer<vtkTransformPolyDataFilter> tpdYY = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    tpdYY->SetInputConnection(normals2->GetOutputPort());
    tpdYY->SetTransform(transformYY);
    
    vtkSmartPointer<vtkPolyDataMapper> mapperYY = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapperYY->SetInputConnection(tpdYY->GetOutputPort());
    vtkSmartPointer<vtkActor> actorYY = vtkSmartPointer<vtkActor>::New();
    col = upcol(nvis::vec3(0, -1, 0), param_b, param_g);
    actorYY->GetProperty()->SetColor(col[0], col[1], col[2]);
    actorYY->SetMapper(mapperYY);
    
    // vtkSmartPointer<vtkTransform> transformZZ = vtkSmartPointer<vtkTransform>::New();
    // // transformZ->Scale(1.25, 1.25, 1.25);
    // transformZZ->RotateY(90);
    // vtkSmartPointer<vtkTransformPolyDataFilter> tpdZZ = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
    // tpdZZ->SetInputConnection(normals2->GetOutputPort());
    // tpdZZ->SetTransform(transformZZ);
    //
    // vtkSmartPointer<vtkPolyDataMapper> mapperZZ = vtkSmartPointer<vtkPolyDataMapper>::New();
    // mapperZZ->SetInputConnection(tpdZZ->GetOutputPort());
    // vtkSmartPointer<vtkActor> actorZZ = vtkSmartPointer<vtkActor>::New();
    // col = cmap(nvis::vec3(0, 0, -1));
    // actorZZ->GetProperty()->SetColor(col[0], col[1], col[2]);
    // actorZZ->SetMapper(mapperZZ);
    
    ren->AddActor(actor);
    ren->AddActor(actorX);
    ren->AddActor(actorY);
    ren->AddActor(actorZ);
    ren->AddActor(actorXX);
    ren->AddActor(actorYY);
    // ren->AddActor(actorZZ);
    
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);
    
    ren_win->Render();
    
    vtkWindowToImageFilter* capture = vtkWindowToImageFilter::New();
    capture->SetInput(ren_win);
    
    vtkPNGWriter* writer = vtkPNGWriter::New();
    writer->SetInputConnection(capture->GetOutputPort());
    
    writer->SetFileName(param_f);
    writer->Write();
    
    
    iren->Initialize();
    iren->Start();
    
    
    return 0;
}








































