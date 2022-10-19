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
#include <VTK/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>

struct dataset_info {
    std::string mesh_base;
    std::string stress;
    std::string fa;
    std::string size;
    std::string westin;
    std::string dir;
    std::string shot;
    std::string id_to_tags;
    std::string stress_span;
    nvis::vec3 up, position, focal_point;
    double near, far;
    nvis::vec3 step;
    nvis::bbox3 valid_bounds;
};

int     param_id;
double  param_mina;
double  param_maxa;
double  param_minv;
double  param_maxv;
char*   param_ssf;
bool    param_v;
double  param_g;
bool    param_shot;

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
    hestOptAdd(&hopt, "mina",   "min anisotropy",   airTypeDouble,  0, 1, &param_mina,          "0",        "min grain anisotropy");
    hestOptAdd(&hopt, "maxa",   "max anisotropy",   airTypeDouble,  0, 1, &param_maxa,          "1",        "max grain anisotropy");
    hestOptAdd(&hopt, "minv",   "min volume",       airTypeDouble,  0, 1, &param_minv,          "0",        "min grain volume");
    hestOptAdd(&hopt, "maxv",   "max volume",       airTypeDouble,  0, 1, &param_maxv,          "-1",       "max grain volume");
    hestOptAdd(&hopt, "g",      "gamma",            airTypeDouble,  0, 1, &param_g,             "1",        "color gamma");
    hestOptAdd(&hopt, "shot",   "screenshot",       airTypeBool,    0, 0, &param_shot,          "0",        "screenshot mode");
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
    double weight = pow(cl, param_g);
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

    if (param_maxv < 0) {
        param_maxv = std::numeric_limits<double>::max();
    }

    // hard-coded settings for available datasets
    dataset_info textured_info, untextured_info, mc_info_00, mc_info_06, mc_info_10_09;
    std::string mesh_base = "/scratch4/data/Garcia/Microstructure/geometry/";
    std::string stat_base = "/scratch4/data/Garcia/Stat/";
    std::string stress_base = "/scratch4/data/Garcia/Stress/vtk/";
    std::string shot_base = "/scratch4/data/Garcia/Screenshots/";

    textured_info.mesh_base = mesh_base + "TexBNKT_190_190_37";
    textured_info.stress = stress_base + "trace-Estress_TexturedBNKT_b09.vtk";
    textured_info.fa = stat_base + "Textured_b09-fa.nrrd";
    textured_info.size = stat_base + "Textured_b09-size.nrrd";
    textured_info.westin = stat_base + "Textured_b09-westin.nrrd";
    textured_info.dir = stat_base + "Textured_b09-direction.nrrd";
    textured_info.shot = shot_base + "Textured_b09.tiff";
    textured_info.id_to_tags = stat_base + "Textured_b09-ids_to_tags.nrrd";
    textured_info.stress_span = stat_base + "Textured_b09-grain-span.nrrd";
    textured_info.position = nvis::vec3(-38.2382, -20.2557, 27.1584);
    textured_info.focal_point = nvis::vec3(17.1675, 19.0469, 3.32833);
    textured_info.up = nvis::vec3(0.254529, 0.213122, 0.943289);
    textured_info.near = 17.6447;
    textured_info.far = 140.982;
    textured_info.step = nvis::vec3(0.2, 0.2, 0.4);
    textured_info.valid_bounds.min() = nvis::vec3(10., 10., 5.);
    textured_info.valid_bounds.max() = nvis::vec3(179., 179., 31.);


    untextured_info.mesh_base = mesh_base + "UnTexBNKT_01_140_140_70";
    untextured_info.stress = stress_base + "trace-Estress_UntexturedBNKT_b09.vtk";
    untextured_info.fa = stat_base + "Untextured_b09-fa.nrrd";
    untextured_info.size = stat_base + "Untextured_b09-size.nrrd";
    untextured_info.westin = stat_base + "Untextured_b09-westin.nrrd";
    untextured_info.dir = stat_base + "Untextured_b09-direction.nrrd";
    untextured_info.shot = shot_base + "Untextured_b09.tiff";
    untextured_info.id_to_tags = stat_base + "Untextured_b09-ids_to_tags.nrrd";
    untextured_info.stress_span = stat_base + "Untextured_b09-grain-span.nrrd";
    untextured_info.position = nvis::vec3(-11.7274, -4.0038, 11.1954);
    untextured_info.focal_point = nvis::vec3(15.4841, 11.0543, -3.40036);
    untextured_info.up = nvis::vec3(0.352451, 0.239877, 0.904565);
    untextured_info.near = 5.28852;
    untextured_info.far = 39.2914;
    untextured_info.step = nvis::vec3(0.07, 0.07, 0.1);
    untextured_info.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    untextured_info.valid_bounds.max() = nvis::vec3(129., 129., 59.);


    mc_info_00.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info_00.stress = stress_base + "trace-Estress_CG_588Grains_r00b09.vtk";
    mc_info_00.fa = stat_base + "CG_588Grains_r00b09-fa.nrrd";
    mc_info_00.size = stat_base + "CG_588Grains_r00b09-size.nrrd";
    mc_info_00.westin = stat_base + "CG_588Grains_r00b09-westin.nrrd";
    mc_info_00.dir = stat_base + "CG_588Grains_r00b09-direction.nrrd";
    mc_info_00.shot = shot_base + "CG_588Grains.tiff";
    mc_info_00.id_to_tags = stat_base + "CG_588Grains_r00b09-ids_to_tags.nrrd";
    mc_info_00.stress_span = stat_base + "CG_588Grains_r00b09-grain-span.nrrd";
    mc_info_00.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_00.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_00.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_00.near = 157.65;
    mc_info_00.far = 961.655;
    mc_info_00.step = nvis::vec3(2, 2, 2);
    mc_info_00.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_00.valid_bounds.max() = nvis::vec3(99., 99., 99.);


    mc_info_06.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info_06.stress = stress_base + "trace-Estress_CG_588Grains_r06b09.vtk";
    mc_info_06.fa = stat_base + "CG_588Grains_r00b09-fa.nrrd";
    mc_info_06.size = stat_base + "CG_588Grains_r00b09-size.nrrd";
    mc_info_06.westin = stat_base + "CG_588Grains_r00b09-westin.nrrd";
    mc_info_06.dir = stat_base + "CG_588Grains_r00b09-direction.nrrd";
    mc_info_06.shot = shot_base + "CG_588Grains.tiff";
    mc_info_06.id_to_tags = stat_base + "CG_588Grains_r06b09-ids_to_tags.nrrd";
    mc_info_06.stress_span = stat_base + "CG_588Grains_r06b09-grain-span.nrrd";
    mc_info_06.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_06.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_06.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_06.near = 157.65;
    mc_info_06.far = 961.655;
    mc_info_06.step = nvis::vec3(2, 2, 2);
    mc_info_06.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info_06.valid_bounds.max() = nvis::vec3(99., 99., 99.);


    mc_info_10_09.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info_10_09.stress = stress_base + "trace-Estress_CG_588Grains_r10b09.vtk";
    mc_info_10_09.fa = stat_base + "CG_588Grains_r00b09-fa.nrrd";
    mc_info_10_09.size = stat_base + "CG_588Grains_r00b09-size.nrrd";
    mc_info_10_09.westin = stat_base + "CG_588Grains_r00b09-westin.nrrd";
    mc_info_10_09.dir = stat_base + "CG_588Grains_r00b09-direction.nrrd";
    mc_info_10_09.shot = shot_base + "CG_588Grains.tiff";
    mc_info_10_09.id_to_tags = stat_base + "CG_588Grains_r10b09-ids_to_tags.nrrd";
    mc_info_10_09.stress_span = stat_base + "CG_588Grains_r10b09-grain-span.nrrd";
    mc_info_10_09.position = nvis::vec3(-324.175, -118.824, 264.273);
    mc_info_10_09.focal_point = nvis::vec3(38.2538, 78.0644, 115.621);
    mc_info_10_09.up = nvis::vec3(0.286275, 0.183103, 0.940489);
    mc_info_10_09.near = 157.65;
    mc_info_10_09.far = 961.655;
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
//    vtkStructuredPoints* stress_field = stress_field_reader->GetOutput();

    Nrrd* __ids = xavier::nrrd_utils::readNrrd(info.id_to_tags);
    xavier::nrrd_utils::nrrd_data_wrapper<int> ids(__ids);
    std::cerr << info.id_to_tags << " loaded.\n";

    Nrrd* __fa = xavier::nrrd_utils::readNrrd(info.fa);
    xavier::nrrd_utils::nrrd_data_wrapper<float> fa(__fa);
    std::cerr << info.fa << " loaded.\n";
    int nb_grains = __fa->axis[0].size;

    Nrrd* __size = xavier::nrrd_utils::readNrrd(info.size);
    xavier::nrrd_utils::nrrd_data_wrapper<float> size(__size);
    std::cerr << info.size << " loaded.\n";

    Nrrd* __gstress = xavier::nrrd_utils::readNrrd(info.stress_span);
    xavier::nrrd_utils::nrrd_data_wrapper<float> grain_stress(__gstress);
    std::cerr << info.stress_span << " loaded.\n";

    Nrrd* __westin = xavier::nrrd_utils::readNrrd(info.westin);
    xavier::nrrd_utils::nrrd_data_wrapper<float> westin(__westin);
    std::cerr << info.westin << " loaded.\n";

    Nrrd* __dir = xavier::nrrd_utils::readNrrd(info.dir);
    xavier::nrrd_utils::nrrd_data_wrapper<float> dir(__dir);
    std::cerr << info.dir << " loaded.\n";

    std::map<int, color_type> colors;
    for (int i = 0 ; i < nb_grains ; ++i) {
        nvis::vec3 _d(dir[3*i], dir[3*i+1], dir[3*i+2]);
        colors[ids[i]] = vec2col(_d, westin[3*i]);
        // std::cerr << "color(" << ids[i] << ") is " << nvis::ivec3(colors[ids[i]]) << '\n';
    }
    std::cerr << "colors set.\n";

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

    // vtkSmartPointer<vtkUnsignedCharArray> color = vtkSmartPointer<vtkUnsignedCharArray>::New();
    // color->SetNumberOfComponents(3);
    // color->SetName("Colors");
    // unsigned char col[][3] = {
    //  {0, 0, 0}, {255, 0, 0}, {255, 255, 0}, {0, 255, 0},
    //  {0, 0, 255}, {255, 0, 255}, {255, 255, 255}, {0, 255, 255}
    // };
    // color->InsertNextTypedTuple(col[0]);
    // color->InsertNextTypedTuple(col[1]);
    // color->InsertNextTypedTuple(col[2]);
    // color->InsertNextTypedTuple(col[3]);
    // color->InsertNextTypedTuple(col[4]);
    // color->InsertNextTypedTuple(col[5]);
    // color->InsertNextTypedTuple(col[6]);
    // color->InsertNextTypedTuple(col[7]);
    // cube_edges->GetPointData()->SetScalars(color);

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

    vtkSmartPointer<vtkActor> mesh_actor, edge_actor;

    if (true) {
        std::cerr << "we are showing grains\n";
        // included grains
        std::set<int> selected_grains;
        for (int i = 0 ; i < nb_grains ; ++i) {
            if (ids[i] >= 0 &&
                    size[i] <= param_maxv && size[i] >= param_minv &&
                    fa[i] >= param_mina && fa[i] <= param_maxa) {
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

        // assign colors to triangles based on the orientation of their grain(s)
        selected_triangles->InitTraversal();
        vtkSmartPointer<vtkUnsignedCharArray> color = vtkSmartPointer<vtkUnsignedCharArray>::New();
        color->SetNumberOfComponents(3);
        color->SetName("Colors");
        while (true) {
            const vtkIdType* ids;
            vtkIdType npts;
            if (!selected_triangles->GetNextCell(npts, ids)) {
                break;
            }
            std::vector<int> tmp1, tmp2, tmp3;

            std::set_intersection(vertex_tags[ids[0]].begin(),
                                  vertex_tags[ids[0]].end(),
                                  vertex_tags[ids[1]].begin(),
                                  vertex_tags[ids[1]].end(),
                                  std::back_inserter(tmp1));
            std::set_intersection(vertex_tags[ids[2]].begin(),
                                  vertex_tags[ids[2]].end(),
                                  tmp1.begin(), tmp1.end(),
                                  std::back_inserter(tmp2));
            std::set_intersection(tmp2.begin(), tmp2.end(),
                                  selected_grains.begin(),
                                  selected_grains.end(),
                                  std::back_inserter(tmp3));
            if (!tmp3.size()) {
                std::cerr << "WARNING: current triangle has invalid corner tags\n";
            } else {
                int ref_id = tmp3.back(); // largest included index common to all 3 vertices
                // std::cerr << "ref_id is " << ref_id << " with color " << nvis::ivec3(colors[ref_id]) << '\n';
                // std::cerr << "inserting color " << colors[ref_id] << " = " << nvis::vec3(colors[ref_id]) / nvis::vec3(255, 255, 255) << '\n';
                color->InsertNextTypedTuple(colors[ref_id].begin());
            }
        }


        vtkSmartPointer<vtkPolyData> selected_mesh = vtkSmartPointer<vtkPolyData>::New();
        selected_mesh->SetPoints(pts);
        selected_mesh->SetPolys(selected_triangles);
        selected_mesh->GetCellData()->SetScalars(color);

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
        // mesh_actor->GetProperty()->SetColor(1, 1, 1);
        mesh_actor->GetProperty()->SetOpacity(0.8);

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

    ren_win->Render();

    if (param_shot) {
        vtkSmartPointer<vtkWindowToImageFilter> capture = vtkSmartPointer<vtkWindowToImageFilter>::New();
        capture->SetInput(ren_win);

        vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
        writer->SetInputConnection(capture->GetOutputPort());

        writer->SetFileName(info.shot.c_str());
        std::cerr << "about to write to file " << info.shot << '\n';
        writer->Write();
    } else {
        iren->Initialize();
        iren->Start();
    }
}
