#include <vector>
#include <map>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <teem/hest.h>
#include <teem/nrrd.h>
#include <fstream>
#include <iostream>

#include <vtk/vtk_utils.hpp>


char* ftle_in, *coord_in, *out;
float ratio;
bool v, c;


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
    hestOptAdd(&hopt, "i",      "input coordinates",    airTypeString,  1, 1, &coord_in,            NULL,       "input coordinates (NRRD)");
    hestOptAdd(&hopt, "f",      "ftle values",          airTypeString,  1, 1, &ftle_in,             NULL,       "FTLE file");
    hestOptAdd(&hopt, "r",      "threshold ratio",      airTypeFloat,   1, 1, &ratio,               NULL,       "FTLE cutoff ratio");
    hestOptAdd(&hopt, "c",      "clip",                 airTypeBool,    0, 1, &c,                   "0",        "clip scene along X and Y axes");
    hestOptAdd(&hopt, "v",      "verbose",              airTypeBool,    0, 1, &v,                   "0",        "display camera setting");
    hestOptAdd(&hopt, "s",      "snapshot",             airTypeString,  0, 1, &out,                 "none",     "output image");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize particle systems filtered by FTLE value",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    Nrrd* nin_coord = nrrdNew();
    if (nrrdLoad(nin_coord, coord_in, NULL)) {
        std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
        return -1;
    }
    float* coord = (float*)nin_coord->data;
    
    Nrrd* nin_ftle = nrrdNew();
    if (nrrdLoad(nin_ftle, ftle_in, NULL)) {
        std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
        return -1;
    }
    float* ftle = (float*)nin_ftle->data;
    
    int N = nin_coord->axis[1].size;
    
    std::vector<double> __tmp(ftle, &ftle[2*N]);
    std::sort(__tmp.begin(), __tmp.end());
    
    float threshold = ratio * __tmp[floor(ratio*2.*(float)N)];
    
    std::vector<nvis::fvec3>    pos;
    std::vector<float>          val;
    
    for (int i = 0 ; i < N ; ++i) {
        float ftle_p = ftle[2*i  ];
        float ftle_m = ftle[2*i+1];
        // if (coord[3*i] < 0 || coord[3*i+1] < 0) continue;
        
        if (ftle_p > ftle_m && ftle_p > threshold) {
            pos.push_back(nvis::fvec3(coord[3*i], coord[3*i+1], coord[3*i+2]));
            val.push_back(ftle_p);
        }
        
        if (ftle_m > ftle_p && ftle_m > threshold) {
            pos.push_back(nvis::fvec3(coord[3*i], coord[3*i+1], coord[3*i+2]));
            val.push_back(-ftle_m);
        }
    }
    
    std::cerr << pos.size() << " (/" << N << ") particles selected\n";
    
#if 1
    
    int nbpts = pos.size();
    
    std::vector<float> copy_val(val.begin(), val.end());
    std::sort(copy_val.begin(), copy_val.end());
    int zero_pos = std::distance(copy_val.begin(), std::lower_bound(copy_val.begin(), copy_val.end(), 0));
    float control_points[] = { copy_val[0],
                               copy_val[zero_pos/2],
                               0,
                               copy_val[(zero_pos + nbpts)/2],
                               copy_val.back()
                             };
    std::cerr << "control values for color mapping are " << control_points[0]
              << ", " << control_points[1] << ", " << control_points[2] << ", "
              << control_points[3] << ", " << control_points[4] << '\n';
              
    vtkColorTransferFunction* colors = vtkColorTransferFunction::New();
    colors->AddRGBPoint(control_points[0], 0, 0, 1);
    colors->AddRGBPoint(control_points[1], 0.5, 0.5, 1);
    colors->AddRGBPoint(control_points[2], 1, 1, 1);
    colors->AddRGBPoint(control_points[3], 1, 0.5, 0.5);
    colors->AddRGBPoint(control_points[4], 1, 0, 0);
    
    // creat a vtkPolyData
    vtkFloatArray* parray = vtkFloatArray::New();
    vtkFloatArray* varray = vtkFloatArray::New();
    parray->SetNumberOfComponents(3);
    parray->SetNumberOfTuples(nbpts);
    varray->SetNumberOfComponents(1);
    varray->SetNumberOfTuples(nbpts);
    for (int i = 0 ; i < nbpts ; ++i) {
        parray->SetTuple(i, pos[i].begin());
        varray->SetTuple1(i, val[i]);
    }
    vtkPolyData* particles = vtkPolyData::New();
    vtkPoints* pts = vtkPoints::New();
    pts->SetData(parray);
    particles->SetPoints(pts);
    particles->GetPointData()->SetScalars(varray);
    vtkIdTypeArray* ids = vtkIdTypeArray::New();
    ids->SetNumberOfValues(nbpts);
    for (int i = 0 ; i < nbpts ; ++i) {
        ids->InsertValue(i, i);
    }
    vtkCellArray* verts = vtkCellArray::New();
    verts->SetCells(nbpts, ids);
    particles->SetVerts(verts);
    
    std::cerr << "polydata created with " << particles->GetNumberOfPoints() << '\n';
    
    vtkPolyData* clipped;
    if (c) {
        vtkPlane* plane = vtk_utils::make_plane(nvis::vec3(1, 0, 0), nvis::vec3(0, 0, 0));
        clipped = vtk_utils::clip_polydata(particles, plane);
        
        std::cerr << "polydata clipped\n";
        std::cerr << "clipped polydata contains " << clipped->GetNumberOfPoints() << " particles\n";
    }
    
    // vtkPlane *plane = vtkPlane::New();
    // plane->SetOrigin(0, 0, 0);
    // plane->SetNormal(1, 0, 0);
    //
    // vtkClipPolyData *clip = vtkClipPolyData::New();
    // clip->SetInput(particles);
    // clip->SetClipFunction(plane);
    // clip->GenerateClipScalarsOff();
    // clip->GenerateClippedOutputOn();
    // clip->SetValue(0);
    // clip->Update();
    
    vtkSphereSource* ball = vtkSphereSource::New();
    ball->SetRadius(0.001);
    ball->SetThetaResolution(6);
    ball->SetPhiResolution(6);
    
    vtkGlyph3D* spheres = vtkGlyph3D::New();
    spheres->SetSourceConnection(ball->GetOutputPort());
    if (!c)
#if VTK_MAJOR_VERSION >= 6
        spheres->SetInputData(particles);
#else
        spheres->SetInput(particles);
#endif
    else
#if VTK_MAJOR_VERSION >= 6
        spheres->SetInputData(clipped);
#else
        spheres->SetInput(clipped);
#endif
    spheres->SetScaling(0);
    spheres->Update();
    
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
    mapper->SetInputConnection(spheres->GetOutputPort());
    mapper->SetLookupTable(colors);
    
    vtkActor* actor = vtkActor::New();
    actor->SetMapper(mapper);
    
    vtkRenderWindow* renWin = vtkRenderWindow::New();
    renWin->PointSmoothingOn();
    renWin->LineSmoothingOn();
    renWin->PolygonSmoothingOn();
    
    // Create graphics stuff
    vtkRenderer* ren = vtkRenderer::New();
    
    if (v) {
        vtk_utils::camera_setting_callback* cb = vtk_utils::camera_setting_callback::New();
        ren->AddObserver(vtkCommand::StartEvent, cb);
        cb->Delete();
    }
    
    vtkCylinderSource* cylinder = vtkCylinderSource::New();
    cylinder->SetRadius(0.15);
    cylinder->SetResolution(100);
    cylinder->SetCenter(0, 0, 0.016);
    
    renWin->AddRenderer(ren);
    renWin->SetSize(1200, 800);
    
    ren->AddActor(actor);
    
    ren->SetBackground(0.0, 0.0, 0.0);
    ren->ResetCamera();
    ren->ResetCameraClippingRange();
    ren->GetActiveCamera()->SetPosition(-0.15957, 0.106947, 0.0862912);
    ren->GetActiveCamera()->SetFocalPoint(0.0737662, 0.0997238, 0.0339788);
    ren->GetActiveCamera()->SetViewUp(0.220237, 0.0664091, 0.973183);
    ren->GetActiveCamera()->SetClippingRange(0.074651, 0.446994);
    
    if (strcmp(out, "none")) {
        vtkWindowToImageFilter* capture = vtkWindowToImageFilter::New();
        capture->SetInput(renWin);
        
        vtkPNGWriter* writer = vtkPNGWriter::New();
        writer->SetInputConnection(capture->GetOutputPort());
        
        writer->SetFileName(out);
        writer->Write();
    }
    
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);
    renWin->Render();
    iren->Initialize();
    iren->Start();
    
#else
    
    std::fstream tmp(out, std::ios::out);
    tmp << "# vtk DataFile Version 2.0\n"
        << "particles with FTLE values for " << coord_in << "\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << pos.size() << " float\n";
    for (int i = 0 ; i < pos.size() ; ++i) {
        tmp << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << '\n';
    }
    tmp << "POINT_DATA " << pos.size() << '\n'
        << "SCALARS ftle float 1\n"
        << "LOOKUP_TABLE default\n";
    for (int i = 0 ; i < val.size() ; ++i) {
        tmp << val[i] << '\n';
    }
    tmp.close();
    std::cerr << out << " was exported\n";
    
#endif
    
    return -1;
}

















