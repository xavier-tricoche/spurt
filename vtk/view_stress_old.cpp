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
#include <teem/hest_helper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>

char* mesh, *stress, *idf, *gstress, *img;
bool verbose;
double mins, maxs, ds;

struct dataset_info {
    std::string mesh_base;
    std::string id_to_tags;
    std::string stress_span;
    nvis::vec3 up, direction, focal_point;
    double near, far;
    nvis::vec3 d;
    nvis::bbox3 valid_bounds;
};

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "m",      "mesh",             airTypeString,  1, 1, &mesh,            NULL,       "mesh base name");
    hestOptAdd(&hopt, "s",      "stress",           airTypeString,  1, 1, &stress,          NULL,       "stress field");
    hestOptAdd(&hopt, "id",     "id to tags",       airTypeString,  1, 1, &idf,             NULL,       "id to tag mapping");
    hestOptAdd(&hopt, "gs",     "grain stress",     airTypeString,  1, 1, &gstress,         NULL,       "grain stress");
    hestOptAdd(&hopt, "mins",   "min stress",       airTypeDouble,  1, 1, &mins,            NULL,       "min stress value");
    hestOptAdd(&hopt, "maxs",   "max stress",       airTypeDouble,  1, 1, &maxs,            NULL,       "max stress value");
    hestOptAdd(&hopt, "ds",     "delta stress",     airTypeDouble,  1, 1, &ds,              NULL,       "stress increment (in %)");
    hestOptAdd(&hopt, "ssf",    "screen shot file", airTypeString,  1, 1, &img,             NULL,       "screenshot base name (batch mode)");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &verbose,         "0",        "verbose mode (debugging)");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Visualize stress field in granular microstructure",
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
    std::cerr << "input mins = " << mins << ", maxs = " << maxs << '\n';
    
    struct dataset_info {
        std::string mesh_base;
        std::string id_to_tags;
        std::string stress_span;
        nvis::vec3 up, position, focal_point;
        double near, far;
        nvis::vec3 d;
        nvis::bbox3 valid_bounds;
    };
    
    data_set_info textured_info, untextured_info, mc_info;
    std::string mesh_base = "/scratch4/data/Garcia/Microstructure/geometry/";
    std::string stat_base = "/scratch4/data/Garcia/Stat/";
    
    textured_info.mesh_base = mesh_base + "TexBNKT_190_190_37";
    textured_info.id_to_tags = stat_base + "Textured_b09-ids_to_tags.nrrd";
    textured_info.stress_span = stat_base + "Textured_b09-grain-span.nrrd";
    textured_info.position = nvis::vec3(-36.3331, -8.47804, 28.695);
    textured_info.focal_point = nvis::vec3(16.9142, 19.7463, 3.89301);
    textured_info.up = nvis::vec3(0.291739, 0.258784, 0.920825);
    textured_info.near = 18.0374;
    textured_info.far = 124.656;
    textured_info.d = nvis::vec3(0.2, 0.2, 0.4);
    textured_info.valid_bounds.min() = nvis::vec3(10., 10., 5.);
    textured_info.valid_bounds.max() = nvis::vec3(179., 179., 31.);
    
    
    untextured_info.mesh_base = mesh_base + "UnTexBNKT_01_140_140_70";
    untextured_info.id_to_tags = stat_base + "Untextured_b09-ids_to_tags.nrrd";
    untextured_info.stress_span = stat_base + "Untextured_b09-grain-span.nrrd";
    untextured_info.position = nvis::vec3(-12.6458, -3.45589, 10.6859);
    untextured_info.focal_point = nvis::vec3(17.4109, 12.476, -3.31413);
    untextured_info.up = nvis::vec3(0.291739, 0.258784, 0.920825);
    untextured_info.near = 8.08427;
    untextured_info.far = 36.5545;
    untextured_info.d = nvis::vec3(0.07, 0.07, 0.1);
    untextured_info.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    untextured_info.valid_bounds.max() = nvis::vec3(129., 129., 59.);
    
    
    mc_info.mesh_base = mesh_base + "MC_588grains_110cubed";
    mc_info.id_to_tags = stat_base + "CG_588Grains_r00b09-ids_to_tags.nrrd";
    mc_info.stress_span = stat_base + "CG_588Grains_r00b09-grain-span.nrrd";
    mc_info.position = nvis::vec3(-318.363, -111.902, 286.273);
    mc_info.focal_point = nvis::vec3(39.8581, 77.9774, 119.417);
    mc_info.up = nvis::vec3(0.291739, 0.258784, 0.920825);
    mc_info.near = 214.331;
    mc_info.far = 888.836;
    mc_info.d = nvis::vec3(2, 2, 2);
    mc_info.valid_bounds.min() = nvis::vec3(10., 10., 10.);
    mc_info.valid_bounds.max() = nvis::vec3(99., 99., 99.);
    
    
    std::string mesh_base(mesh);
    std::string mesh_name(mesh_base), edge_name(mesh_base), corner_name(mesh_base),
        vertex_tag_name(mesh_base);
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
    stress_field_reader->SetFileName(stress);
    stress_field_reader->Update();
    vtkStructuredPoints* stress_field = stress_field_reader->GetOutput();
    
    Nrrd* __ids = xavier::nrrd_utils::readNrrd(idf);
    xavier::nrrd_utils::nrrd_data_wrapper<int> ids(__ids);
    std::cerr << std::string(idf) << " loaded.\n";
    
    Nrrd* __gstress = xavier::nrrd_utils::readNrrd(gstress);
    xavier::nrrd_utils::nrrd_data_wrapper<float> grain_stress(__gstress);
    std::cerr << gstress << " loaded.\n";
    int nb_grains = __gstress->axis[1].size;
    
    {
        if (mins > maxs) {
            double tmp = mins;
            mins = maxs;
            maxs = tmp;
        }
        std::vector<double> gstress_vals;
        for (int i = 0 ; i < 2*nb_grains ; ++i) {
            gstress_vals.push_back(grain_stress[i]);
        }
        mins = std::max(mins, *std::min_element(gstress_vals.begin(), gstress_vals.end()));
        maxs = std::min(maxs, *std::max_element(gstress_vals.begin(), gstress_vals.end()));
        
        ds *= (maxs - mins);
    }
    
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
    vtkSmartPointer<vtkCubeSource> cube = vtkSmartPointer<vtkCubeSource>::New();
    
    // BNTK_01
    // double dx = 7.0e-8;
    // double dy = 7.0e-8;
    // double dz = 1.0e-7;
    // cube->SetBounds(10.*dx, 129.*dx, 10.*dy, 129.*dy, 10.*dz, 59.*dz);
    
    // // Textured_09
    // double dx = 0.2;
    // double dy = 0.2;
    // double dz = 0.4;
    // cube->SetBounds(10.*dx, 179.*dx, 10.*dy, 179.*dy, 5.*dz, 31.*dz);
    
    // Untextured_09
    double dx = 0.07;
    double dy = 0.07;
    double dz = 0.1;
    cube->SetBounds(10.*dx, 129.*dx, 10.*dy, 129.*dy, 10.*dz, 59.*dz);
    
    // // MC
    // double dx = 2;
    // double dy = 2;
    // double dz = 2;
    // cube->SetBounds(10.*dx, 99.*dx, 10.*dy, 99.*dy, 10.*dz, 99.*dz);
    
    nvis::bbox3 frame(nvis::vec3(10.*dx, 10.*dy, 10.*dz),
                      nvis::vec3(99.*dx, 99.*dy, 99.*dz));
                      
    // create cube outline
    vtkSmartPointer<vtkDoubleArray> __coords = vtkSmartPointer<vtkDoubleArray>::New();
    __coords->SetNumberOfComponents(3);
    nvis::vec3 p = frame.min();
    __coords->InsertNextTuple(p.begin());
    p[0] = frame.max()[0];
    __coords->InsertNextTuple(p.begin());
    p[1] = frame.max()[1];
    __coords->InsertNextTuple(p.begin());
    p[0] = frame.min()[0];
    __coords->InsertNextTuple(p.begin());
    p[1] = frame.min()[1];
    p[2] = frame.max()[2];
    __coords->InsertNextTuple(p.begin());
    p[0] = frame.max()[0];
    __coords->InsertNextTuple(p.begin());
    p[1] = frame.max()[1];
    __coords->InsertNextTuple(p.begin());
    p[0] = frame.min()[0];
    __coords->InsertNextTuple(p.begin());
    
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetData(__coords);
    
    vtkSmartPointer<vtkCellArray> __lines = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType n;
    vtkIdType edges[][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0},
        {4, 5}, {5, 6}, {6, 7}, {7, 8},
        {0, 4}, {1, 5}, {2, 6}, {3, 7}
    };
    for (int i = 0 ; i < 12 ; ++i) {
        __lines->InsertNextCell(n, edges[i]);
    }
    
    vtkSmartPointer<vtkPolyData> cube_edges = vtkSmartPointer<vtkPolyData>::New();
    cube_edges->SetPoints(pts);
    cube_edges->SetLines(__lines);
    
    vtkSmartPointer<vtkTubeFilter> edge_tubes = vtkSmartPointer<vtkTubeFilter>::New();
    edge_tubes->SetInput(cube_edges);
    edge_tubes->SetRadius(0.3 * sqrt(dx * dx + dy * dy + dz * dz));
    edge_tubes->SetNumberOfSides(6);
    
    
    
    std::cerr << "cube created\n";
    
    vtkSmartPointer<vtkPolyDataMapper> cube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    // cube_mapper->SetInputConnection(cube->GetOutputPort());
    cube_mapper->SetInputConnection(edge_tubes);
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
    
    if (verbose) {
        vtk_utils::camera_setting_callback* cb = vtk_utils::camera_setting_callback::New();
        ren->AddObserver(vtkCommand::StartEvent, cb);
        cb->Delete();
    }
    
// // BNTK_01
// ren->GetActiveCamera()->SetPosition(-1.17812e-05, 1.46389e-05, -3.96052e-06);
// ren->GetActiveCamera()->SetFocalPoint(4.45208e-06, 4.7677e-06, 4.33919e-06);
// ren->GetActiveCamera()->SetViewUp(0.35207, -0.19136, -0.916203);
// ren->GetActiveCamera()->SetClippingRange(6.51136e-06, 3.85059e-05);

// // Textured_09
// ren->GetActiveCamera()->SetPosition(-36.3331, -8.47804, 28.695);
// ren->GetActiveCamera()->SetFocalPoint(16.9142, 19.7463, 3.89301);
// ren->GetActiveCamera()->SetViewUp(0.291739, 0.258784, 0.920825);
// ren->GetActiveCamera()->SetClippingRange(18.0374, 124.656);

// Untextured_09
    ren->GetActiveCamera()->SetPosition(-12.6458, -3.45589, 10.6859);
    ren->GetActiveCamera()->SetFocalPoint(17.4109, 12.476, -3.31413);
    ren->GetActiveCamera()->SetViewUp(0.291739, 0.258784, 0.920825);
    ren->GetActiveCamera()->SetClippingRange(8.08427, 36.5545);
    
// // MC
// ren->GetActiveCamera()->SetPosition(-318.363, -111.902, 286.273);
// ren->GetActiveCamera()->SetFocalPoint(39.8581, 77.9774, 119.417);
// ren->GetActiveCamera()->SetViewUp(0.291739, 0.258784, 0.920825);
// ren->GetActiveCamera()->SetClippingRange(214.331, 888.836);

    std::cerr << "Camera set\n";
    
    vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    ren_win->PointSmoothingOn();
    ren_win->LineSmoothingOn();
    ren_win->PolygonSmoothingOn();
    ren_win->SetAlphaBitPlanes(1);
    ren_win->SetMultiSamples(10);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1600, 1200);
    
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);
    
    vtkSmartPointer<vtkContourFilter> contour = vtkSmartPointer<vtkContourFilter>::New();
    contour->SetNumberOfContours(1);
    contour->ComputeNormalsOn();
    contour->SetInputConnection(stress_field_reader->GetOutputPort());
    
    std::cerr << "contour filter initialized\n";
    std::cerr << "min stress = " << mins << ", max stress = " << maxs << '\n';
    
    for (double s = mins ; s < maxs + 0.1*ds ; s += ds) {
        std::cerr << "stress isovalue = " << s << '\n';
        contour->SetValue(0, s);
        contour->Update();
        
        vtkSmartPointer<vtkPolyDataMapper> iso_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        iso_mapper->SetInputConnection(contour->GetOutputPort());
        vtkSmartPointer<vtkActor> iso_actor = vtkSmartPointer<vtkActor>::New();
        iso_actor->SetMapper(iso_mapper);
        nvis::vec3 col = color(s, mins, maxs);
        iso_actor->GetProperty()->SetColor(col[0], col[1], col[2]);
        iso_actor->GetProperty()->SetOpacity(0.8);
        
        // included grains
        std::set<int> selected_grains;
        for (int i = 0 ; i < nb_grains ; ++i) {
            if (ids[i] >= 0 &&
                    grain_stress[2*i] <= s && grain_stress[2*i+1] >= s) {
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
        // mesh_mapper->SetColorModeToMapScalars();
        // mesh_mapper->ScalarVisibilityOn();
        // mesh_mapper->SelectColorArray("Colors");
        // mesh_mapper->SetScalarModeToUseCellData();
        vtkSmartPointer<vtkActor> mesh_actor = vtkSmartPointer<vtkActor>::New();
        mesh_actor->SetMapper(mesh_mapper);
        mesh_actor->GetProperty()->SetColor(1, 1, 1);
        mesh_actor->GetProperty()->SetOpacity(0.8);
        
        vtkSmartPointer<vtkPolyData> selected_edges = vtkSmartPointer<vtkPolyData>::New();
        selected_edges->SetPoints(edges->GetPoints());
        selected_edges->SetLines(selected_lines);
        //
        vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
        tubes->SetInput(selected_edges);
        
        double radius = 0.1 * sqrt(dx * dx + dy * dy + dz * dz);
        tubes->SetRadius(radius);
        tubes->SetNumberOfSides(6);
        vtkSmartPointer<vtkPolyDataMapper> edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        edge_mapper->SetInputConnection(tubes->GetOutputPort());
        vtkSmartPointer<vtkActor> edge_actor = vtkSmartPointer<vtkActor>::New();
        edge_actor->SetMapper(edge_mapper);
        edge_actor->GetProperty()->SetColor(0.5, 0, 0);
        
        ren->AddActor(mesh_actor);
        ren->AddActor(edge_actor);
        ren->AddActor(cube_actor);
        ren->AddActor(iso_actor);
        ren_win->Render();
        
        vtkSmartPointer<vtkWindowToImageFilter> capture = vtkSmartPointer<vtkWindowToImageFilter>::New();
        capture->SetInput(ren_win);
        
        vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
        writer->SetInputConnection(capture->GetOutputPort());
        
        std::ostringstream os;
        os << img << "-stress=" << s << ".tiff";
        
        writer->SetFileName(os.str().c_str());
        writer->Write();
        
        // std::cerr << "before clearing the screen:\n"
        //           << "@iso_actor  = " << (int)&(*iso_actor) << '\n'
        //           << "@mesh_actor = " << (int)&(*mesh_actor) << '\n'
        //           << "@edge_actor = " << (int)&(*edge_actor) << '\n'
        //           << std::flush;
        //
        // ren->RemoveActor(iso_actor);
        // ren->RemoveActor(mesh_actor);
        // ren->RemoveActor(edge_actor);
        //
        // std::cerr << "after clearing the screen:\n"
        //           << "@iso_actor  = " << (int)&(*iso_actor) << '\n'
        //           << "@mesh_actor = " << (int)&(*mesh_actor) << '\n'
        //           << "@edge_actor = " << (int)&(*edge_actor) << '\n'
        //           << std::flush;
    }
}















