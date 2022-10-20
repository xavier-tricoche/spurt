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

#include <string>
#include <math/fixed_vector.hpp>
#include <vtk/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>

struct dataset_info {
    std::string mesh_base;
    std::string stress;
    std::string id_to_tags;
    std::string stress_span;
    nvis::vec3 up, position, focal_point;
    double near, far;
    nvis::vec3 step;
    nvis::bbox3 valid_bounds;
};

int     param_id;
double  param_s;
char*   param_ssf;
bool    param_v;
bool    param_grains;

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
    hestOptAdd(&hopt, "s",      "delta stress",     airTypeDouble,  0, 1, &param_s,             "0",        "stress isovalue");
    hestOptAdd(&hopt, "grain",  "show grains",      airTypeBool,    0, 0, &param_grains,        "0",        "show intersected grains");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &param_v,             "0",        "verbose mode (debugging)");

    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize stress field in granular microstructure",
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

typedef unsigned char                   uchar;
typedef nvis::fixed_vector<uchar, 3>    color_type;

inline color_type vec2col(const nvis::vec3& dir, double cl)
{
    nvis::vec3 e(dir);
    double weight = pow(cl, 0.4);
    double norm = nvis::norm(e);
    if (!norm || cl == 0) {
        return color_type(127, 127, 127);
    }
    e = nvis::abs(e * weight / norm);
    color_type c(127 + uchar(floor(128.*e[0])),
                 127 + uchar(floor(128.*e[1])),
                 127 + uchar(floor(128.*e[2])));
    if (*std::min_element(&c[0], &c[3]) < 127) {
        std::cerr << "ERROR!!!!! wrong color : " << c << " generated for dir = " << dir
                  << " (e = " << e << ") and cl = " << cl << '\n';
    }

    return c;
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

nvis::vec3 color(float v, float min, float max)
{
    nvis::vec3 __colors[] = {
        nvis::vec3(0, 0, 0.5),  // 0:  dark blue
        nvis::vec3(0, 0, 1),    // 1:  blue
        nvis::vec3(0, 0.5, 1),  // 2:  sky blue
        nvis::vec3(0, 1, 1),    // 3:  cyan
        nvis::vec3(0, 1, 0.5),  // 4:
        nvis::vec3(0, 1, 0),    // 5:  green
        nvis::vec3(0.5, 1, 0),  // 6:
        nvis::vec3(1, 1, 0),    // 7:  yellow
        nvis::vec3(1, 0.5, 0),  // 8:  orange
        nvis::vec3(1, 0, 0),    // 9:  red
        nvis::vec3(1, 0, 0.5),  // 10:
        nvis::vec3(1, 0, 1),    // 11: magenta
        nvis::vec3(1, 0.5, 1),  // 12:
        nvis::vec3(1, 1, 1)     // 13: white
    };

    int i = std::min((int)floor((v - min) / (max - min) * 14), 13);
    float __min = i * (max - min) / 14.;
    float __max = (i + 1) * (max - min) / 14.;
    float u = (v - __min) / (__max - __min);
    return (1. - u)*__colors[i] + u*__colors[i+1];
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);

    // hard-coded settings for available datasets
    dataset_info textured_info, untextured_info, mc_info_00, mc_info_06, mc_info_10_09;
    std::string mesh_base = "/scratch4/data/Garcia/Microstructure/geometry/";
    std::string stat_base = "/scratch4/data/Garcia/Stat/";
    std::string stress_base = "/scratch4/data/Garcia/Stress/vtk/";

    textured_info.mesh_base = mesh_base + "TexBNKT_190_190_37";
    textured_info.stress = stress_base + "trace-Estress_TexturedBNKT_b09.vtk";
    textured_info.id_to_tags = stat_base + "Textured_b09-ids_to_tags.nrrd";
    textured_info.stress_span = stat_base + "Textured_b09-grain-span.nrrd";
    textured_info.position = nvis::vec3(-36.3331, -8.47804, 28.695);
    textured_info.focal_point = nvis::vec3(16.9142, 19.7463, 3.89301);
    textured_info.up = nvis::vec3(0.291739, 0.258784, 0.920825);
    textured_info.near = 18.0374;
    textured_info.far = 124.656;
    textured_info.step = nvis::vec3(0.2, 0.2, 0.4);
    textured_info.valid_bounds.min() = nvis::vec3(10., 10., 5.);
    textured_info.valid_bounds.max() = nvis::vec3(179., 179., 31.);


    untextured_info.mesh_base = mesh_base + "UnTexBNKT_01_140_140_70";
    untextured_info.stress = stress_base + "trace-Estress_UntexturedBNKT_b09.vtk";
    untextured_info.id_to_tags = stat_base + "Untextured_b09-ids_to_tags.nrrd";
    untextured_info.stress_span = stat_base + "Untextured_b09-grain-span.nrrd";
    untextured_info.position = nvis::vec3(-12.6458, -3.45589, 10.6859);
    untextured_info.focal_point = nvis::vec3(17.4109, 12.476, -3.31413);
    untextured_info.up = nvis::vec3(0.291739, 0.258784, 0.920825);
    untextured_info.near = 8.08427;
    untextured_info.far = 36.5545;
    untextured_info.step = nvis::vec3(0.07, 0.07, 0.1);
    untextured_info.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    untextured_info.valid_bounds.max() = nvis::vec3(129., 129., 59.);


    mc_info_00.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info_00.stress = stress_base + "trace-Estress_CG_588Grains_r00b09.vtk";
    mc_info_00.id_to_tags = stat_base + "CG_588Grains_r00b09-ids_to_tags.nrrd";
    mc_info_00.stress_span = stat_base + "CG_588Grains_r00b09-grain-span.nrrd";
    mc_info_00.position = nvis::vec3(-318.363, -111.902, 286.273);
    mc_info_00.focal_point = nvis::vec3(39.8581, 77.9774, 119.417);
    mc_info_00.up = nvis::vec3(0.291739, 0.258784, 0.920825);
    mc_info_00.near = 214.331;
    mc_info_00.far = 888.836;
    mc_info_00.step = nvis::vec3(2, 2, 2);
    mc_info_00.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_00.valid_bounds.max() = nvis::vec3(99., 99., 99.);


    mc_info_06.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info_06.stress = stress_base + "trace-Estress_CG_588Grains_r06b09.vtk";
    mc_info_06.id_to_tags = stat_base + "CG_588Grains_r06b09-ids_to_tags.nrrd";
    mc_info_06.stress_span = stat_base + "CG_588Grains_r06b09-grain-span.nrrd";
    mc_info_06.position = nvis::vec3(-318.363, -111.902, 286.273);
    mc_info_06.focal_point = nvis::vec3(39.8581, 77.9774, 119.417);
    mc_info_06.up = nvis::vec3(0.291739, 0.258784, 0.920825);
    mc_info_06.near = 214.331;
    mc_info_06.far = 888.836;
    mc_info_06.step = nvis::vec3(2, 2, 2);
    mc_info_06.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_06.valid_bounds.max() = nvis::vec3(99., 99., 99.);


    mc_info_10_09.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info_10_09.stress = stress_base + "trace-Estress_CG_588Grains_r10b09.vtk";
    mc_info_10_09.id_to_tags = stat_base + "CG_588Grains_r10b09-ids_to_tags.nrrd";
    mc_info_10_09.stress_span = stat_base + "CG_588Grains_r10b09-grain-span.nrrd";
    mc_info_10_09.position = nvis::vec3(-318.363, -111.902, 286.273);
    mc_info_10_09.focal_point = nvis::vec3(39.8581, 77.9774, 119.417);
    mc_info_10_09.up = nvis::vec3(0.291739, 0.258784, 0.920825);
    mc_info_10_09.near = 214.331;
    mc_info_10_09.far = 888.836;
    mc_info_10_09.step = nvis::vec3(2, 2, 2);
    mc_info_10_09.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_10_09.valid_bounds.max() = nvis::vec3(99., 99., 99.);

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
        default:
            std::cerr << "unknown dataset\n";
            return 1;
    }
    const dataset_info& info = *__info;

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

    vtkSmartPointer<vtkStructuredPointsReader> stress_field_reader = vtkSmartPointer<vtkStructuredPointsReader>::New();
    stress_field_reader->SetFileName(info.stress.c_str());
    stress_field_reader->Update();
    vtkStructuredPoints* stress_field = stress_field_reader->GetOutput();

    Nrrd* __ids = spurt::nrrd_utils::readNrrd(info.id_to_tags);
    spurt::nrrd_utils::nrrd_data_wrapper<int> ids(__ids);
    std::cerr << info.id_to_tags << " loaded.\n";

    Nrrd* __gstress = spurt::nrrd_utils::readNrrd(info.stress_span);
    spurt::nrrd_utils::nrrd_data_wrapper<float> grain_stress(__gstress);
    std::cerr << info.stress_span << " loaded.\n";
    int nb_grains = __gstress->axis[1].size;

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

    // border mesh
    nvis::bbox3 frame;
    frame.min() = info.step * info.valid_bounds.min();
    frame.max() = info.step * info.valid_bounds.max();
    nvis::vec3 diag = frame.size();
    std::cerr << "frame = " << frame << '\n';
    nvis::vec3 logo[] = {
        nvis::vec3(0, 0, 0), nvis::vec3(diag[0], 0, 0),
        nvis::vec3(0, diag[1], 0), nvis::vec3(-diag[0], 0, 0),
        nvis::vec3(0, -diag[1], diag[2]), nvis::vec3(diag[0], 0, 0),
        nvis::vec3(0, diag[1], 0), nvis::vec3(-diag[0], 0, 0)
    };
    vtkIdType __edges[][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 4},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    };

    // create cube outline
    vtkSmartPointer<vtkDoubleArray> __coords = vtkSmartPointer<vtkDoubleArray>::New();
    __coords->SetNumberOfComponents(3);
    nvis::vec3 p = frame.min();
    __coords->InsertNextTuple(p.begin());
    for (int i = 1 ; i < 8 ; ++i) {
        p += logo[i];
        __coords->InsertNextTuple(p.begin());
    }
    vtkSmartPointer<vtkPoints> cube_pts = vtkSmartPointer<vtkPoints>::New();
    cube_pts->SetData(__coords);

    vtkSmartPointer<vtkCellArray> cube_lines = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType n = 2;
    for (int i = 0 ; i < 12 ; ++i) {
        cube_lines->InsertNextCell(n, __edges[i]);
    }

    vtkSmartPointer<vtkPolyData> cube_edges = vtkSmartPointer<vtkPolyData>::New();
    cube_edges->SetPoints(cube_pts);
    cube_edges->SetLines(cube_lines);

    vtkSmartPointer<vtkUnsignedCharArray> color = vtkSmartPointer<vtkUnsignedCharArray>::New();
    color->SetNumberOfComponents(3);
    color->SetName("Colors");
    unsigned char col[][3] = {
        {0, 0, 0}, {255, 0, 0}, {255, 255, 0}, {0, 255, 0},
        {0, 0, 255}, {255, 0, 255}, {255, 255, 255}, {0, 255, 255}
    };
    color->InsertNextTypedTuple(col[0]);
    color->InsertNextTypedTuple(col[1]);
    color->InsertNextTypedTuple(col[2]);
    color->InsertNextTypedTuple(col[3]);
    color->InsertNextTypedTuple(col[4]);
    color->InsertNextTypedTuple(col[5]);
    color->InsertNextTypedTuple(col[6]);
    color->InsertNextTypedTuple(col[7]);
    cube_edges->GetPointData()->SetScalars(color);

    vtkSmartPointer<vtkTubeFilter> edge_tubes = vtkSmartPointer<vtkTubeFilter>::New();
    edge_tubes->SetInputData(cube_edges);
    edge_tubes->SetRadius(0.3 * nvis::norm(info.step));
    edge_tubes->SetNumberOfSides(6);

    std::cerr << "cube created\n";

    vtkSmartPointer<vtkPolyDataMapper> cube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cube_mapper->SetInputConnection(edge_tubes->GetOutputPort());
    vtkSmartPointer<vtkActor> cube_actor = vtkSmartPointer<vtkActor>::New();
    cube_actor->SetMapper(cube_mapper);
    cube_actor->GetProperty()->SetColor(1, 1, 1);
    cube_actor->GetProperty()->SetOpacity(1);

    std::cerr << "cube actor created\n";

    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    ren->SetUseDepthPeeling(1);
    ren->SetMaximumNumberOfPeels(100);
    ren->SetOcclusionRatio(0.1);
    ren->SetBackground(0, 0, 0);
    ren->ResetCamera();

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
    // ren_win->PointSmoothingOn();
    // ren_win->LineSmoothingOn();
    // ren_win->PolygonSmoothingOn();
    ren_win->SetAlphaBitPlanes(1);
    // ren_win->SetMultiSamples(10);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1600, 1200);

    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);

    vtkSmartPointer<vtkContourFilter> contour = vtkSmartPointer<vtkContourFilter>::New();
    contour->SetNumberOfContours(1);
    contour->ComputeNormalsOn();
    contour->SetInputConnection(stress_field_reader->GetOutputPort());

    std::cerr << "contour filter initialized\n";
    std::cerr << "stress isovalue = " << param_s << '\n';
    contour->SetValue(0, param_s);

    std::cerr << "isosurface comprises "
              << contour->GetOutput()->GetPolys()->GetNumberOfCells()
              << " cells\n";

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
    iso_actor->GetProperty()->SetColor(1, 0, 1);
    iso_actor->GetProperty()->SetOpacity(0.5);
    iso_actor->GetProperty()->SetInterpolationToPhong();
    // iso_actor->GetProperty()->EdgeVisibilityOn();
    // iso_actor->GetProperty()->BackfaceCullingOff();
    // iso_actor->GetProperty()->SetAmbientColor(1, 0, 1);
    // iso_actor->GetProperty()->SetEdgeColor(1, 0, 1);

    vtkSmartPointer<vtkActor> mesh_actor, edge_actor;

    if (param_grains) {
        std::cerr << "we are showing grains\n";
        // included grains
        std::set<int> selected_grains;
        for (int i = 0 ; i < nb_grains ; ++i) {
            if (ids[i] >= 0 &&
                    grain_stress[2*i] <= param_s && grain_stress[2*i+1] >= param_s) {
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
            vtkIdType _n;
            const vtkIdType* ids;
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
        mesh_mapper->SetInputData(selected_mesh);
        mesh_actor = vtkSmartPointer<vtkActor>::New();
        mesh_actor->SetMapper(mesh_mapper);
        mesh_actor->GetProperty()->SetColor(1, 1, 1);
        mesh_actor->GetProperty()->SetOpacity(0.3);

        vtkSmartPointer<vtkPolyData> selected_edges = vtkSmartPointer<vtkPolyData>::New();
        selected_edges->SetPoints(edges->GetPoints());
        selected_edges->SetLines(selected_lines);
        //
        vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
        tubes->SetInputData(selected_edges);

        double radius = 0.1 * nvis::norm(info.step);
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
    ren->AddActor(cube_actor);
    ren->AddActor(iso_actor);

    ren_win->Render();
    iren->Initialize();
    iren->Start();
}
