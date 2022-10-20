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
double  param_min;
double  param_max;
double  param_d;
char*   param_ssf;
bool    param_v;
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
    hestOptAdd(&hopt, "min",    "min angle",        airTypeDouble,  0, 1, &param_min,           "90",       "min angle");
    hestOptAdd(&hopt, "g",      "gamma",            airTypeFloat,   0, 1, &param_g,             "1",        "color scale gamma factor");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &param_v,             "0",        "verbose mode (debugging)");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize granular microstructure",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef std::set<int>   tag_type;

template<typename T>
inline bool included(const tag_type& tags, const T values, int size)
{
    for (int i = 0 ; i < size ; ++i) {
        if (tags.find(values[i]) == tags.end()) {
            return false;
        }
    }
    return true;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s)
{
    os << "[ ";
    for (typename std::set<T>::const_iterator i = s.begin() ; i != s.end() ; ++i) {
        os << *i << " ";
    }
    os << "]";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& s)
{
    os << "[ ";
    for (typename std::vector<T>::const_iterator i = s.begin() ; i != s.end() ; ++i) {
        os << *i << " ";
    }
    os << "]";
    return os;
}

template<typename T>
std::ostream& operator<<(std::ostream& os, const std::list<T>& s)
{
    os << "[ ";
    for (typename std::list<T>::const_iterator i = s.begin() ; i != s.end() ; ++i) {
        os << *i << " ";
    }
    os << "]";
    return os;
}

typedef unsigned char                   uchar;
typedef nvis::fixed_vector<uchar, 4>    uchar_color_type;

inline nvis::vec3 zdirection(const nvis::vec3& dir)
{
    return nvis::vec3(sin(dir[2])*sin(dir[1]),
                      cos(dir[2])*sin(dir[1]),
                      cos(dir[1]));
}

inline uchar_color_type a2col(double a)
{
    double u = pow((fabs(a) - param_min)/(180.-param_min), param_g);
    nvis::vec4 c = (1. - u) * nvis::vec4(1, 1, 1, 0.5) + u * nvis::vec4(1, 0, 0, 1);
    return uchar_color_type(255.*c);
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
    
    std::string mesh_name(info.mesh_base), edge_name(info.mesh_base), vertex_tag_name(info.mesh_base);
    mesh_name.append("-mesh.vtk");
    edge_name.append("-edges.vtk");
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
    
    Nrrd* __id = spurt::nrrd_utils::readNrrd(info.id_to_tags);
    spurt::nrrd_utils::nrrd_data_wrapper<int> id(__id);
    std::cerr << info.id_to_tags << " loaded.\n";
    
    std::string stat_base(info.stat_base);
    std::string id_name(stat_base);
    id_name.append("-ids_to_tags.nrrd");
    int nb_grains = __id->axis[0].size;
    
    std::cerr << "there are " << nb_grains << " grains\n";
    
    Nrrd* __orient = spurt::nrrd_utils::readNrrd(info.orientation);
    spurt::nrrd_utils::nrrd_data_wrapper<float> orient(__orient);
    std::cerr << info.orientation << " loaded.\n";
    
    
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    {
        for (double a = param_min ; a <= 180 ; a += 1) {
            nvis::vec4 c = a2col(a);
            c /= 255.;
            ctf->AddRGBPoint(a, c[0], c[1], c[2]);
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
    vtkSmartPointer<vtkDoubleArray> coord = vtkSmartPointer<vtkDoubleArray>::New();
    coord->SetNumberOfComponents(3);
    for (int i = 0 ; i < mesh->GetNumberOfPoints() ; ++i) {
        double pt[3];
        mesh->GetPoint(i, pt);
        coord->InsertNextTuple(pt);
    }
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetData(coord);
    std::cerr << "mesh vertices duplicated.\n";
    
    // included triangles
    vtkSmartPointer<vtkUnsignedCharArray> color = vtkSmartPointer<vtkUnsignedCharArray>::New();
    color->SetNumberOfComponents(4);
    color->SetName("Colors");
    vtkSmartPointer<vtkCellArray> selected_triangles = vtkSmartPointer<vtkCellArray>::New();
    for (int n = 0 ; n < mesh->GetNumberOfCells() ; ++n) {
        vtkSmartPointer<vtkIdList> ids = vtkSmartPointer<vtkIdList>::New();
        mesh->GetCellPoints(n, ids);
        
        std::vector<int> tmp1, tmp2;
        vtkIdType id0, id1, id2;
        id0 = ids->GetId(0);
        id1 = ids->GetId(1);
        id2 = ids->GetId(2);
        std::set_intersection(vertex_tags[id0].begin(),
                              vertex_tags[id0].end(),
                              vertex_tags[id1].begin(),
                              vertex_tags[id1].end(),
                              std::back_inserter(tmp1));
        std::set_intersection(vertex_tags[id2].begin(),
                              vertex_tags[id2].end(),
                              tmp1.begin(), tmp1.end(),
                              std::back_inserter(tmp2));
                              
        if (tmp2.size() == 2 && *std::max_element(tmp2.begin(), tmp2.end()) > 0) {
            int i1 = tmp2[0], i2 = tmp2[1];
            nvis::vec3 dir1 = zdirection(nvis::vec3(orient[3*i1], orient[3*i1+1], orient[3*i1+2]));
            nvis::vec3 dir2 = zdirection(nvis::vec3(orient[3*i2], orient[3*i2+1], orient[3*i2+2]));
            dir1 /= nvis::norm(dir1);
            dir2 /= nvis::norm(dir2);
            double alpha = acos(nvis::inner(dir1, dir2)) * 180. / M_PI;
            if (alpha > param_min) {
                // add this triangle
                selected_triangles->InsertNextCell(3, ids->GetPointer(0));
                // add corresponding color
                color->InsertNextTypedTuple(a2col(alpha).begin());
            }
        }
    }
    
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    ren->SetUseDepthPeeling(1);
    ren->SetMaximumNumberOfPeels(100);
    ren->SetOcclusionRatio(0.1);
    ren->SetBackground(0.5, 0.5, 0.5);
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
    
    // ren->AddActor2D(color_bar);
    
    vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    ren_win->SetAlphaBitPlanes(1);
    ren_win->SetMultiSamples(0);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1600, 1200);
    
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);
    
    vtkSmartPointer<vtkActor> mesh_actor, edge_actor;
    
    // vtkSmartPointer<vtkCellArray> selected_lines = vtkSmartPointer<vtkCellArray>::New();
    // lines->InitTraversal();
    // while (true) {
    //  vtkIdType _n, *id;
    //  if (!lines->GetNextCell(_n, id)) break;
    //  if (included(selected_vertices, id, 2))
    //      selected_lines->InsertNextCell(_n, id);
    // }
    // std::cerr << selected_lines->GetNumberOfCells() << " edges passed the test\n";
    //
    // // assign colors to triangles based on the orientation of their grain(s)
    // selected_triangles->InitTraversal();
    // vtkSmartPointer<vtkUnsignedCharArray> color = vtkSmartPointer<vtkUnsignedCharArray>::New();
    // color->SetNumberOfComponents(3);
    // color->SetName("Colors");
    // while (true) {
    //  vtkIdType *ids;
    //  vtkIdType npts;
    //  if (!selected_triangles->GetNextCell(npts, ids)) break;
    //  std::vector<int> tmp1, tmp2, tmp3;
    //
    //  std::set_intersection(vertex_tags[ids[0]].begin(),
    //                        vertex_tags[ids[0]].end(),
    //                        vertex_tags[ids[1]].begin(),
    //                        vertex_tags[ids[1]].end(),
    //                        std::back_inserter(tmp1));
    //  std::set_intersection(vertex_tags[ids[2]].begin(),
    //                        vertex_tags[ids[2]].end(),
    //                        tmp1.begin(), tmp1.end(),
    //                        std::back_inserter(tmp2));
    //  std::set_intersection(tmp2.begin(), tmp2.end(),
    //                        selected_grains.begin(),
    //                        selected_grains.end(),
    //                        std::back_inserter(tmp3));
    //  if (!tmp3.size()) {
    //      std::cerr << "WARNING: current triangle has invalid corner tags\n";
    //  }
    //  else {
    //      int ref_id = tmp3.back(); // largest included index common to all 3 vertices
    //      // std::cerr << "ref_id is " << ref_id << " with color " << nvis::ivec3(colors[ref_id]) << '\n';
    //      // std::cerr << "inserting color " << colors[ref_id] << " = " << nvis::vec3(colors[ref_id]) / nvis::vec3(255, 255, 255) << '\n';
    //      color->InsertNextTypedTuple(colors[ref_id].begin());
    //  }
    // }
    
    vtkSmartPointer<vtkPolyData> selected_mesh = vtkSmartPointer<vtkPolyData>::New();
    selected_mesh->SetPoints(pts);
    selected_mesh->SetPolys(selected_triangles);
    selected_mesh->GetCellData()->SetScalars(color);
    
    //
    vtkSmartPointer<vtkPolyDataMapper> mesh_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mesh_mapper->SetInput(selected_mesh);
    mesh_actor = vtkSmartPointer<vtkActor>::New();
    mesh_actor->SetMapper(mesh_mapper);
    // mesh_actor->GetProperty()->SetColor(1, 1, 1);
    // mesh_actor->GetProperty()->SetOpacity(0.5);
    
    // vtkSmartPointer<vtkPolyData> selected_edges = vtkSmartPointer<vtkPolyData>::New();
    // selected_edges->SetPoints(edges->GetPoints());
    // selected_edges->SetLines(selected_lines);
    // //
    // vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
    // tubes->SetInput(selected_edges);
    //
    // double radius = 0.05 * nvis::norm(info.step);
    // tubes->SetRadius(radius);
    // tubes->SetNumberOfSides(6);
    // vtkSmartPointer<vtkPolyDataMapper> edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    // edge_mapper->SetInputConnection(tubes->GetOutputPort());
    // edge_actor = vtkSmartPointer<vtkActor>::New();
    // edge_actor->SetMapper(edge_mapper);
    // edge_actor->GetProperty()->SetColor(0.5, 0, 0);
    
    ren->AddActor(mesh_actor);
    // ren->AddActor(edge_actor);
    ren->AddActor2D(color_bar);
    ren_win->Render();
    
    vtkSmartPointer<vtkWindowToImageFilter> capture = vtkSmartPointer<vtkWindowToImageFilter>::New();
    capture->SetInput(ren_win);
    
    if (true) {
        iren->Initialize();
        iren->Start();
    } else {
        vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
        writer->SetInputConnection(capture->GetOutputPort());
        
        std::ostringstream os;
        os << info.movie_path << "grains_orientation_segmentation.tiff";
        
        writer->SetFileName(os.str().c_str());
        std::cerr << "about to write to file " << os.str() << '\n';
        writer->Write();
    }
    return 0;
}































