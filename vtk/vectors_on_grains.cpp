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
#include <VTK/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>
#include "Garcia_vis_helper.hpp"
#include "boost/shared_ptr.hpp"
#include "data/field_wrapper.hpp"

int     param_id;
double  param_d;
char*   param_ssf;
char*   param_val;
char*   param_span;
char*   param_val_name;
bool    param_v;
float   param_g;
bool    param_signed;


using namespace Garcia_vis_helper;


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
    hestOptAdd(&hopt, "s",      "scalar field",     airTypeString,  1, 1, &param_val,           NULL,       "input scalar field");
    hestOptAdd(&hopt, "n",      "scalar name",      airTypeString,  1, 1, &param_val_name,      NULL,       "name of scalar quantity");
    hestOptAdd(&hopt, "span",   "grain span",       airTypeString,  1, 1, &param_span,          NULL,       "scalar value span per grain");
    hestOptAdd(&hopt, "d",      "delta value",      airTypeDouble,  0, 1, &param_d,             "0.1",      "value increment (in %)");
    hestOptAdd(&hopt, "g",      "gamma",            airTypeFloat,   0, 1, &param_g,             "1",        "color scale gamma factor");
    hestOptAdd(&hopt, "signed", "signed values",    airTypeBool,    0, 1, &param_signed,        "1",        "are values signed?");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &param_v,             "0",        "verbose mode (debugging)");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize value distribution over granular microstructure",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
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
    
    std::string mesh_name(info.mesh_base), edge_name(info.mesh_base),
        vertex_tag_name(info.mesh_base), point_coord_name(info.mesh_base);
        
    mesh_name.append("-mesh.vtk");
    edge_name.append("-edges.vtk");
    vertex_tag_name.append("-point-attributes.txt");
    point_coord_name.append("-init_coords.nrrd");
    
    Nrrd* __coords = xavier::nrrd_utils::readNrrd(point_coord_name);
    xavier::nrrd_utils::nrrd_data_wrapper<float> coords(__coords);
    
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
    
    typedef xavier::nrrd_data_traits<Nrrd*>  nrrd_data_traits;
    Nrrd* __val = xavier::nrrd_utils::readNrrd(param_val);
    boost::shared_ptr<nrrd_data_traits> scalar_field =
        boost::shared_ptr<nrrd_data_traits>(new nrrd_data_traits(__val));
        
    Nrrd* __span = xavier::nrrd_utils::readNrrd(param_span);
    std::map<int, std::pair<float, float> > span;
    std::vector<float> all_values;
    {
        xavier::nrrd_utils::nrrd_data_wrapper<float> tmp(__span);
        for (int i = 0 ; i < __span->axis[1].size ; ++i) {
            span[tmp[3*i]] = std::pair<float, float>(tmp[3*i+1], tmp[3*i+2]);
            all_values.push_back(tmp[3*i+1]);
            all_values.push_back(tmp[3*i+2]);
        }
    }
    
    int nb_grains = span.size();
    
    std::cerr << "there are " << nb_grains << " grains\n";
    
    color_map<float> cmap(all_values, param_g, true);
    color_map<float>::map_kind map_type = (param_signed ?
                                           color_map<float>::DOUBLE_ENDED :
                                           color_map<float>::REDUNDANT_RAINBOW);
                                           
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    double min = *std::min_element(all_values.begin(), all_values.end());
    double max = *std::max_element(all_values.begin(), all_values.end());
    double d = param_d * (max - min);
    for (double v = min ; v <= max ; v += 0.01 * (max - min)) {
        nvis::vec3 c = cmap(v, map_type);
        ctf->AddRGBPoint(v, c[0], c[1], c[2]);
    }
    vtkSmartPointer<vtkScalarBarActor> color_bar;
    color_bar = vtkSmartPointer<vtkScalarBarActor>::New();
    color_bar->SetLookupTable(ctf);
    color_bar->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    color_bar->GetPositionCoordinate()->SetValue(0.6, 0.95);
    color_bar->SetOrientationToHorizontal();
    color_bar->SetWidth(0.35);
    color_bar->SetHeight(0.05);
    color_bar->SetNumberOfLabels(7);
    std::cerr << "colors set.\n";
    std::cerr << "min value = " << min << ", max value = " << max << ", delta val = " << d << '\n';
    
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
    
    vtkSmartPointer<vtkDoubleArray> coord = vtkSmartPointer<vtkDoubleArray>::New();
    coord->SetNumberOfComponents(3);
    vtkSmartPointer<vtkDoubleArray> mesh_values = vtkSmartPointer<vtkDoubleArray>::New();
    mesh_values->SetNumberOfComponents(1);
    for (int i = 0 ; i < mesh->GetNumberOfPoints() ; ++i) {
        double pt[3];
        mesh->GetPoint(i, pt);
        coord->InsertNextTuple(pt);
        double f = 0;
        scalar_field->get_value(nvis::vec3(pt[0], pt[1], pt[2]), f);
        mesh_values->InsertNextTuple(&f);
    }
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetData(coord);
    std::cerr << "mesh vertices duplicated.\n";
    
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    ren->SetUseDepthPeeling(1);
    ren->SetMaximumNumberOfPeels(100);
    ren->SetOcclusionRatio(0.1);
    ren->SetBackground(0, 0, 0);
    ren->ResetCamera();
    
    nvis::bbox3 bound;
    {
        bound.min() = info.valid_bounds.min() * info.step;
        bound.max() = info.valid_bounds.max() * info.step;
    }
    vtkSmartPointer<vtkActor> frame_actor = draw_frame(bound);
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
    
    ren->AddActor2D(color_bar);
    
    vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    ren_win->SetAlphaBitPlanes(1);
    ren_win->SetMultiSamples(0);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1600, 1200);
    
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);
    
    int frame_counter = 0;
    for (double s = min ; s < max + 0.1*d ; s += d) {
        std::cerr << "value = " << s << '\n';
        
        double __min = s;
        double __max = s + d;
        
        // included grains
        std::set<int> selected_grains;
        for (std::map< int, std::pair<float, float> >::const_iterator it = span.begin() ;
                it != span.end() ; ++it) {
            int id = it->first;
            if (it->second.first <= s && it->second.second >= s) {
                selected_grains.insert(id);
            }
        }
        std::cerr << selected_grains.size() << " from " << nb_grains << " grains passed the test.\n";
        vtkSmartPointer<vtkActor> mesh_actor, edge_actor;
        
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
            vtkSmartPointer<vtkIdList> id = vtkSmartPointer<vtkIdList>::New();
            mesh->GetCellPoints(n, id);
            if (included(selected_vertices, id->GetPointer(0), 3)) {
                selected_triangles->InsertNextCell(3, id->GetPointer(0));
            }
        }
        std::cerr << selected_triangles->GetNumberOfCells() << " triangles passed the test\n";
        
        vtkSmartPointer<vtkCellArray> selected_lines = vtkSmartPointer<vtkCellArray>::New();
        lines->InitTraversal();
        while (true) {
            vtkIdType _n, *id;
            if (!lines->GetNextCell(_n, id)) {
                break;
            }
            if (included(selected_vertices, id, 2)) {
                selected_lines->InsertNextCell(_n, id);
            }
        }
        std::cerr << selected_lines->GetNumberOfCells() << " edges passed the test\n";
        
        vtkSmartPointer<vtkPolyData> selected_mesh = vtkSmartPointer<vtkPolyData>::New();
        selected_mesh->SetPoints(pts);
        selected_mesh->SetPolys(selected_triangles);
        selected_mesh->GetPointData()->SetScalars(mesh_values);
        
        vtkSmartPointer<vtkContourFilter> contour = vtkSmartPointer<vtkContourFilter>::New();
        contour->SetInput(selected_mesh);
        contour->SetNumberOfContours(1);
        contour->SetValue(0, s);
        contour->Update();
        
        vtkSmartPointer<vtkTubeFilter> contour_tubes = vtkSmartPointer<vtkTubeFilter>::New();
        contour_tubes->SetInputConnection(contour->GetOutputPort());
        contour_tubes->SetRadius(0.1*nvis::norm(info.step));
        contour_tubes->SetNumberOfSides(8);
        
        vtkSmartPointer<vtkPolyDataMapper> contour_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        contour_mapper->SetInputConnection(contour_tubes->GetOutputPort());
        contour_mapper->SetLookupTable(ctf);
        contour_mapper->ScalarVisibilityOn();
        vtkSmartPointer<vtkActor> contour_actor = vtkSmartPointer<vtkActor>::New();
        contour_actor->SetMapper(contour_mapper);
        contour_actor->GetProperty()->SetOpacity(1);
        
        //
        vtkSmartPointer<vtkPolyDataMapper> mesh_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mesh_mapper->SetInput(selected_mesh);
        mesh_mapper->ScalarVisibilityOn();
        mesh_mapper->SetLookupTable(ctf);
        mesh_actor = vtkSmartPointer<vtkActor>::New();
        mesh_actor->SetMapper(mesh_mapper);
        // mesh_actor->GetProperty()->SetColor(1, 1, 1);
        mesh_actor->GetProperty()->SetOpacity(0.75);
        
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
        ren->AddActor(contour_actor);
        ren_win->Render();
        
        vtkSmartPointer<vtkWindowToImageFilter> capture = vtkSmartPointer<vtkWindowToImageFilter>::New();
        capture->SetInput(ren_win);
        
        vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
        writer->SetInputConnection(capture->GetOutputPort());
        
        std::ostringstream os;
        os << info.movie_path << param_val_name << "_g=" << param_g << "_frame_" << std::setw(3) << std::setfill('0') << frame_counter++ << "_of_"
           << (int)floor((max - min) / d) + 1 << ".tiff";
           
        writer->SetFileName(os.str().c_str());
        std::cerr << "about to write to file " << os.str() << '\n';
        writer->Write();
        
        ren->RemoveActor(mesh_actor);
        ren->RemoveActor(edge_actor);
        ren->RemoveActor(contour_actor);
    }
}































































