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

#include <string>
#include <math/fixed_vector.hpp>
#include <vtk/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>
#include "Garcia_vis_helper.hpp"

int     param_id;
double  param_mins;
double  param_maxs;
double  param_ds;
char*   param_ssf;
bool    param_v;
bool    param_grains = true;
float   param_g;

void wait(int s)
{
    nvis::timer t;
    while (t.elapsed() < s) {}
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
    hestOptAdd(&hopt, "id",     "dataset ID",       airTypeInt,     1, 1, &param_id,            NULL,       "dataset ID: 0: textured, 1: untextured, 2: MC_r00b09, 3: MC_r06b09, 4: MC_r10b09");
    hestOptAdd(&hopt, "mins",   "min stress",       airTypeDouble,  0, 1, &param_mins,          "-1e+10",   "min stress value");
    hestOptAdd(&hopt, "maxs",   "max stress",       airTypeDouble,  0, 1, &param_maxs,          "1e+10",    "max stress value");
    hestOptAdd(&hopt, "ds",     "delta value",      airTypeDouble,  0, 1, &param_ds,            "0.1",      "stress increment (in %)");
    hestOptAdd(&hopt, "g",      "gamma",            airTypeFloat,   0, 1, &param_g,             "1",        "color scale gamma factor");
    hestOptAdd(&hopt, "grain",  "show grains",      airTypeBool,    0, 0, &param_grains,        "0",        "show intersected grains");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &param_v,             "0",        "verbose mode (debugging)");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize stress field in granular microstructure",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    using namespace Garcia_vis_helper;
    set_paths();
    
    dataset_info* __info;
    switch (param_id) {
        case 0:
            __info = &textured_info;
            break;
        case 1:
            __info = &untextured_info;
            break;
        case 2:
            __info = &mc_info_00;
            break;
        case 3:
            __info = &mc_info_06;
            break;
        case 4:
            __info = &mc_info_10_09;
            break;
        case 5:
            __info = &mc_info_10_05;
            break;
        default:
            std::cerr << "unknown dataset\n";
            return 1;
    }
    const dataset_info& info = *__info;
    
    std::cerr << "mesh base = " << info.mesh_base << '\n';
    std::cerr << "base dir = " << info.base_dir << '\n';
    std::cerr << "grain = " << (param_grains ? "true" : "false") << '\n';
    
    std::string mesh_name(info.mesh_base), edge_name(info.mesh_base), corner_name(info.mesh_base),
        vertex_tag_name(info.mesh_base);
    mesh_name.append("-mesh.vtk");
    edge_name.append("-edges.vtk");
    corner_name.append("-corners.vtk");
    vertex_tag_name.append("-point-attributes.txt");
    
    vtkSmartPointer<vtkDataSetReader> mesh_reader = vtkSmartPointer<vtkDataSetReader>::New();
    mesh_reader->SetFileName(mesh_name.c_str());
    mesh_reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    mesh->DeepCopy(mesh_reader->GetOutput());
    std::cerr << "mesh contains " << std::flush << mesh->GetNumberOfPoints() << " points\n";
    std::cerr << "mesh contains " << std::flush << mesh->GetNumberOfCells() << " cells\n";
    mesh_reader->Delete();
    std::cerr << mesh_name << " loaded.\n";
    
    vtkSmartPointer<vtkPolyDataReader> edge_reader = vtkSmartPointer<vtkPolyDataReader>::New();
    edge_reader->SetFileName(edge_name.c_str());
    edge_reader->Update();
    vtkSmartPointer<vtkPolyData> edges = vtkSmartPointer<vtkPolyData>::New();
    edges->DeepCopy(edge_reader->GetOutput());
    vtkCellArray* lines = edges->GetLines();
    std::cerr << "there are " << lines->GetNumberOfCells() << " edges in input\n";
    edge_reader->Delete();
    std::cerr << edge_name << " loaded.\n";
    
    vtkSmartPointer<vtkPolyDataReader> corner_reader = vtkSmartPointer<vtkPolyDataReader>::New();
    corner_reader->SetFileName(corner_name.c_str());
    corner_reader->Update();
    vtkSmartPointer<vtkPolyData> corners = vtkSmartPointer<vtkPolyData>::New();
    corners->DeepCopy(corner_reader->GetOutput());
    corner_reader->Delete();
    std::cerr << corner_name << " loaded.\n";
    
    vtkSmartPointer<vtkStructuredPointsReader> field_reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    field_reader->SetFileName(info.dfield_norm.c_str());
    field_reader->Update();
    vtkStructuredPoints* field = field_reader->GetOutput();
    std::cerr << info.dfield_norm << " loaded.\n";
    
    Nrrd* __ids = spurt::nrrd_utils::readNrrd(info.id_to_tags);
    spurt::nrrd_utils::nrrd_data_wrapper<int> ids(__ids);
    std::cerr << info.id_to_tags << " loaded.\n";
    
    Nrrd* __span = spurt::nrrd_utils::readNrrd(info.dfield_span);
    spurt::nrrd_utils::nrrd_data_wrapper<float> grain_field(__span);
    std::cerr << info.dfield_span << " loaded.\n";
    int nb_grains = __span->axis[1].size;
    
    typedef color_map<double>::color_type   color_type;
    
    double mins, maxs, ds;
    std::vector<double> vals;
    {
        if (param_mins > param_maxs) {
            maxs = param_mins;
            mins = param_maxs;
        } else {
            maxs = param_maxs;
            mins = param_mins;
        }
        for (int i = 0 ; i < 2*nb_grains ; ++i) {
            vals.push_back(grain_field[i]);
        }
        mins = std::max(mins, *std::min_element(vals.begin(), vals.end()));
        maxs = std::min(maxs, *std::max_element(vals.begin(), vals.end()));
        
        ds = param_ds * (maxs - mins);
    }
    std::cerr << "min value = " << mins << ", max value = " << maxs << ", delta val = " << ds << '\n';
    color_map<double> cmap(vals, param_g, true);
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    {
        double span = maxs - mins;
        double dv = span / 100.;
        for (double v = mins ; v <= maxs ; v += dv) {
            color_type c = cmap(v, color_map<double>::REDUNDANT_RAINBOW);
            ctf->AddRGBPoint(v, c[0], c[1], c[2]);
        }
    }
    vtkSmartPointer<vtkScalarBarActor> color_bar = vtkSmartPointer<vtkScalarBarActor>::New();
    color_bar->SetLookupTable(ctf);
    color_bar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    color_bar->GetPositionCoordinate()->SetValue(0.6, 0.95);
    color_bar->SetOrientationToHorizontal();
    color_bar->SetWidth(0.35);
    color_bar->SetHeight(0.05);
    color_bar->SetNumberOfLabels(7);
    
    std::map<int, tag_type> vertex_tags;
    std::fstream attributes(vertex_tag_name.c_str());
    while (!attributes.eof()) {
        int i, n, id;
        attributes >> i >> n;
        tag_type tags;
        for (int j = 0 ; j < n ; ++j) {
            attributes >> id;
            tags.insert(id);
        }
        vertex_tags[i] = tags;
    }
    attributes.close();
    std::cerr << vertex_tag_name << " loaded.\n";
    
    // selected triangles
    vtkSmartPointer<vtkDoubleArray> coords = vtkSmartPointer<vtkDoubleArray>::New();
    coords->SetNumberOfComponents(3);
    for (int i = 0 ; i < mesh->GetNumberOfPoints() ; ++i) {
        double pt[3];
        mesh->GetPoint(i, pt);
        coords->InsertNextTuple(pt);
    }
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetData(coords);
    std::cerr << "mesh vertices duplicated.\n";
    
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    ren->SetUseDepthPeeling(1);
    ren->SetMaximumNumberOfPeels(100);
    ren->SetOcclusionRatio(0.1);
    ren->SetBackground(0, 0, 0);
    ren->ResetCamera();
    
    nvis::bbox3 bounds;
    {
        double tmp[6];
        field->GetBounds(tmp);
        bounds.min() = nvis::vec3(tmp[0], tmp[2], tmp[4]);
        bounds.max() = nvis::vec3(tmp[1], tmp[3], tmp[5]);
    }
    vtkSmartPointer<vtkActor> frame_actor = draw_frame(bounds, 0.1);
    ren->AddActor(frame_actor);
    std::cerr << "cube created\n";
    
    if (param_v) {
        vtk_utils::camera_setting_callback* cb = vtk_utils::camera_setting_callback::New();
        ren->AddObserver(vtkCommand::StartEvent, cb);
        cb->Delete();
    }
    
    ren->GetActiveCamera()->SetPosition(info.position.begin());
    ren->GetActiveCamera()->SetFocalPoint(info.focal_point.begin());
    ren->GetActiveCamera()->SetViewUp(info.up.begin());
    ren->GetActiveCamera()->SetClippingRange(info.near, info.far);
    std::cerr << "Camera set\n";
    
    vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    ren_win->SetAlphaBitPlanes(1);
    ren_win->SetMultiSamples(0);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1600, 1200);
    
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);
    
    vtkSmartPointer<vtkContourFilter> contour = vtkSmartPointer<vtkContourFilter>::New();
    contour->SetNumberOfContours(1);
    contour->ComputeNormalsOn();
    contour->SetInputConnection(field_reader->GetOutputPort());
    
    std::cerr << "contour filter initialized\n";
    std::cerr << "min val = " << mins << ", max val = " << maxs << '\n';
    
    int frame_counter = 0;
    for (double s = mins ; s < maxs + 0.1*ds ; s += ds) {
        std::cerr << "isovalue = " << s << '\n';
        contour->SetValue(0, s);
        contour->Update();
        
        std::cerr << "current isosurface comprises "
                  << contour->GetOutput()->GetPolys()->GetNumberOfCells()
                  << " cells\n";
                  
        if (!contour->GetOutput()->GetPolys()->GetNumberOfCells()) {
            std::cerr << "skipping\n";
            continue;
        }
        vtkSmartPointer<vtkPolyDataNormals> normals = vtkSmartPointer<vtkPolyDataNormals>::New();
        normals->SetInputConnection(contour->GetOutputPort());
        normals->SplittingOff();
        
        vtkSmartPointer<vtkSmoothPolyDataFilter> smooth = vtkSmartPointer<vtkSmoothPolyDataFilter>::New();
        smooth->SetInputConnection(normals->GetOutputPort());
        smooth->SetRelaxationFactor(0.25);
        smooth->SetNumberOfIterations(10);
        
        vtkSmartPointer<vtkPolyDataMapper> iso_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        iso_mapper->SetInputConnection(smooth->GetOutputPort());
        iso_mapper->ScalarVisibilityOff();
        
        vtkSmartPointer<vtkActor> iso_actor = vtkSmartPointer<vtkActor>::New();
        iso_actor->SetMapper(iso_mapper);
        nvis::vec3 col = cmap(s, color_map<double>::REDUNDANT_RAINBOW);
        iso_actor->GetProperty()->SetColor(col[0], col[1], col[2]);
        iso_actor->GetProperty()->SetOpacity(1);
        iso_actor->GetProperty()->SetInterpolationToPhong();
        vtkSmartPointer<vtkActor> mesh_actor, edge_actor;
        
        if (param_grains) {
            std::cerr << "we are showing grains\n";
            // included grains
            std::set<int> selected_grains;
            for (int i = 0 ; i < nb_grains ; ++i) {
                if (ids[i] >= 0 &&
                        grain_field[2*i] <= s && grain_field[2*i+1] >= s) {
                    selected_grains.insert(ids[i]);
                }
            }
            std::cerr << selected_grains.size() << " from " << nb_grains << " grains passed the test.\n";
            
            // included vertices
            std::set<int> selected_vertices;
            for (int i = 0 ; i < vertex_tags.size() ; ++i) {
                const tag_type& tag = vertex_tags[i];
                std::vector<int> inter;
                std::set_intersection(tag.begin(), tag.end(),
                                      selected_grains.begin(), selected_grains.end(),
                                      std::back_inserter(inter));
                if (!inter.empty()) {
                    selected_vertices.insert(i);
                }
            }
            std::cerr << selected_vertices.size() << " vertices passed the test.\n";
            
            // included triangles
            vtkSmartPointer<vtkCellArray> selected_triangles = vtkSmartPointer<vtkCellArray>::New();
            for (int n = 0 ; n < mesh->GetNumberOfCells() ; ++n) {
                vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
                mesh->GetCellPoints(n, ids);
                if (included(selected_vertices, ids->GetPointer(0), 3)) {
                    selected_triangles->InsertNextCell(3, ids->GetPointer(0));
                }
            }
            std::cerr << selected_triangles->GetNumberOfCells() << " triangles passed the test\n";
            
            vtkSmartPointer<vtkPolyData> selected_mesh = vtkSmartPointer<vtkPolyData>::New();
            selected_mesh->SetPoints(pts);
            selected_mesh->SetPolys(selected_triangles);
            
            vtkSmartPointer<vtkCellArray> selected_lines = vtkSmartPointer<vtkCellArray>::New();
            lines->InitTraversal();
            while (true) {
                vtkIdType _n, *ids;
                if (!lines->GetNextCell(_n, ids)) {
                    break;
                }
                if (included(selected_vertices, ids, 2)) {
                    selected_lines->InsertNextCell(_n, ids);
                }
            }
            std::cerr << selected_lines->GetNumberOfCells() << " edges passed the test\n";
            
            //
            vtkSmartPointer<vtkPolyDataMapper> mesh_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            mesh_mapper->SetInput(selected_mesh);
            mesh_actor = vtkSmartPointer<vtkActor>::New();
            mesh_actor->SetMapper(mesh_mapper);
            mesh_actor->GetProperty()->SetColor(1, 1, 1);
            mesh_actor->GetProperty()->SetOpacity(0.5);
            
            vtkSmartPointer<vtkPolyData> selected_edges = vtkSmartPointer<vtkPolyData>::New();
            selected_edges->SetPoints(edges->GetPoints());
            selected_edges->SetLines(selected_lines);
            //
            vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
            tubes->SetInput(selected_edges);
            
            double radius = 0.05 * nvis::norm(info.step);
            tubes->SetRadius(radius);
            tubes->SetNumberOfSides(6);
            vtkSmartPointer<vtkPolyDataMapper> edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
            edge_mapper->SetInputConnection(tubes->GetOutputPort());
            edge_actor = vtkSmartPointer<vtkActor>::New();
            edge_actor->SetMapper(edge_mapper);
            edge_actor->GetProperty()->SetColor(0.5, 0, 0);
            
            ren->AddActor(mesh_actor);
            ren->AddActor(edge_actor);
        }
        ren->AddActor2D(color_bar);
        ren->AddActor(iso_actor);
        ren_win->Render();
        
        vtkSmartPointer<vtkWindowToImageFilter> capture = vtkSmartPointer<vtkWindowToImageFilter>::New();
        capture->SetInput(ren_win);
        
        vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
        writer->SetInputConnection(capture->GetOutputPort());
        
        std::ostringstream os;
        os << info.movie_path << (param_grains ? "grains/" : "no-grains/")
           << "frame_with_color_bar" << std::setw(3) << std::setfill('0') << frame_counter++ << ".tiff";
           
        writer->SetFileName(os.str().c_str());
        std::cerr << "about to write to file " << os.str() << '\n';
        writer->Write();
        
        ren->RemoveActor(iso_actor);
        ren->RemoveActor(color_bar);
        if (param_grains) {
            ren->RemoveActor(mesh_actor);
            ren->RemoveActor(edge_actor);
        }
    }
}



































