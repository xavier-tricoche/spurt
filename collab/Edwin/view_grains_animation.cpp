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

#include <string>
#include <math/small_vector.hpp>
#include <vtk/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>

char* mesh, *stat, *img;
bool verbose;
double mina, maxa, da;
int mins, maxs, ds;
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
    hestOptAdd(&hopt, "s",      "stat",             airTypeString,  1, 1, &stat,            NULL,       "grain statistics base name");
    hestOptAdd(&hopt, "mina",   "min aniso",        airTypeDouble,  0, 1, &mina,            "0",        "min grain anisotropy");
    hestOptAdd(&hopt, "maxa",   "min aniso",        airTypeDouble,  0, 1, &maxa,            "1",        "max grain anisotropy");
    hestOptAdd(&hopt, "da",     "delta aniso",      airTypeDouble,  0, 1, &da,              "0",        "anisotropy increment");
    hestOptAdd(&hopt, "mins",   "min size",         airTypeInt,     0, 1, &mins,            "0",        "min grain size");
    hestOptAdd(&hopt, "maxs",   "min size",         airTypeInt,     0, 1, &maxs,            "-1",       "max grain size");
    hestOptAdd(&hopt, "ds",     "delta size",       airTypeInt,     0, 1, &ds,              "0",        "size increment");
    hestOptAdd(&hopt, "ssf",    "screen shot file", airTypeString,  1, 1, &img,             NULL,       "screenshot base name (batch mode)");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &verbose,         "0",        "verbose mode (debugging)");

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


int main(int argc, char* argv[])
{
    initialize(argc, argv);
    if (maxs == -1) {
        maxs = std::numeric_limits<int>::max();
    }

    std::string mesh_base(mesh);
    std::string mesh_name(mesh_base), edge_name(mesh_base), corner_name(mesh_base),
        vertex_tag_name(mesh_base);
    mesh_name.append("-mesh.vtk");
    edge_name.append("-edges.vtk");
    corner_name.append("-corners.vtk");
    vertex_tag_name.append("-point-attributes.txt");

    vtkDataSetReader* mesh_reader = vtkDataSetReader::New();
    mesh_reader->SetFileName(mesh_name.c_str());
    mesh_reader->Update();
    vtkUnstructuredGrid* mesh = vtkUnstructuredGrid::New();
    mesh->DeepCopy(mesh_reader->GetOutput());
    std::cerr << "mesh contains " << std::flush << mesh->GetNumberOfPoints() << " points\n";
    std::cerr << "mesh contains " << std::flush << mesh->GetNumberOfCells() << " cells\n";
    mesh_reader->Delete();
    std::cerr << mesh_name << " loaded.\n";

    vtkPolyDataReader* edge_reader = vtkPolyDataReader::New();
    edge_reader->SetFileName(edge_name.c_str());
    edge_reader->Update();
    vtkPolyData* edges = vtkPolyData::New();
    edges->DeepCopy(edge_reader->GetOutput());
    vtkCellArray* lines = edges->GetLines();
    std::cerr << "there are " << lines->GetNumberOfCells() << " edges in input\n";
    edge_reader->Delete();
    std::cerr << edge_name << " loaded.\n";

    vtkPolyDataReader* corner_reader = vtkPolyDataReader::New();
    corner_reader->SetFileName(corner_name.c_str());
    corner_reader->Update();
    vtkPolyData* corners = vtkPolyData::New();
    corners->DeepCopy(corner_reader->GetOutput());
    corner_reader->Delete();
    std::cerr << corner_name << " loaded.\n";

    std::string stat_base(stat);
    std::string fa_name(stat_base), size_name(stat_base), id_name(stat_base),
        westin_name(stat_base), dir_name(stat_base);
    fa_name.append("-fa.nrrd");
    size_name.append("-size.nrrd");
    id_name.append("-ids_to_tags.nrrd");
    westin_name.append("-westin.nrrd");
    dir_name.append("-direction.nrrd");

    Nrrd* __ids = spurt::nrrd_utils::readNrrd(id_name);
    spurt::nrrd_utils::nrrd_data_wrapper<int> ids(__ids);
    std::cerr << id_name << " loaded.\n";

    Nrrd* __fa = spurt::nrrd_utils::readNrrd(fa_name);
    spurt::nrrd_utils::nrrd_data_wrapper<float> fa(__fa);
    std::cerr << fa_name << " loaded.\n";
    int nb_grains = __fa->axis[0].size;

    Nrrd* __size = spurt::nrrd_utils::readNrrd(size_name);
    spurt::nrrd_utils::nrrd_data_wrapper<int> size(__size);
    std::cerr << size_name << " loaded.\n";

    Nrrd* __westin = spurt::nrrd_utils::readNrrd(westin_name);
    spurt::nrrd_utils::nrrd_data_wrapper<float> westin(__westin);
    std::cerr << westin_name << " loaded.\n";

    Nrrd* __dir = spurt::nrrd_utils::readNrrd(dir_name);
    spurt::nrrd_utils::nrrd_data_wrapper<float> dir(__dir);
    std::cerr << dir_name << " loaded.\n";

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

    if (true) {
        int* sizes = (int*)__size->data;
        std::vector<int> __s;
        std::copy(sizes, sizes + nb_grains, std::back_inserter(__s));
        std::sort(__s.begin(), __s.end());
        std::cerr << "size ranges from " << __s[0] << " and " << __s[__s.size()-3]
                  << ", median value is " << __s[__s.size()/2] << '\n';
        maxs = std::min(maxs, __s[__s.size()-3]);

        float* aniso = (float*)__fa->data;
        std::vector<float> __a;
        std::copy(aniso, aniso + nb_grains, std::back_inserter(__a));
        std::sort(__a.begin(), __a.end());
        std::cerr << "FA ranges from " << __a[0] << " and " << __a.back()
                  << ", median value is " << __a[__a.size()/2] << '\n';
    }

    // border mesh
    vtkSmartPointer<vtkCubeSource> cube = vtkSmartPointer<vtkCubeSource>::New();
    // MC
    double dx = 2;
    double dy = 2;
    double dz = 2;
    cube->SetBounds(10.*dx, 99.*dx, 10.*dy, 99.*dy, 10.*dz, 99.*dz);

    vtkSmartPointer<vtkPolyDataMapper> cube_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    cube_mapper->SetInputConnection(cube->GetOutputPort());
    vtkSmartPointer<vtkActor> cube_actor = vtkSmartPointer<vtkActor>::New();
    cube_actor->SetMapper(cube_mapper);
    cube_actor->GetProperty()->SetColor(1, 1, 1);
    cube_actor->GetProperty()->SetOpacity(0.2);

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

    // // Untextured_09
    // ren->GetActiveCamera()->SetPosition(-12.6458, -3.45589, 10.6859 );
    // ren->GetActiveCamera()->SetFocalPoint(17.4109, 12.476, -3.31413);
    // ren->GetActiveCamera()->SetViewUp(0.291739, 0.258784, 0.920825);
    // ren->GetActiveCamera()->SetClippingRange(8.08427, 36.5545);

    // MC
    ren->GetActiveCamera()->SetPosition(-318.363, -111.902, 286.273);
    ren->GetActiveCamera()->SetFocalPoint(39.8581, 77.9774, 119.417);
    ren->GetActiveCamera()->SetViewUp(0.291739, 0.258784, 0.920825);
    ren->GetActiveCamera()->SetClippingRange(214.331, 888.836);

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

    if (da <= 0 || da > (maxa - mina)) {
        da = maxa - mina;
    }
    if (ds <= 0 || ds > (maxs - mins)) {
        ds = maxs - mins;
    }

    for (double _mina = mina ; _mina + 0.001*da < maxa ; _mina += da) {
        double _maxa = _mina + da;
        for (int _mins = mins ; _mins + 0.001*ds < maxs ; _mins += ds) {
            int _maxs = _mins + ds;

            std::cerr << "anisotropy bounds are: " << _mina << " - " << _maxa << '\n';
            std::cerr << "size bounds are: " << _mins << " - " << _maxs << '\n';

            // included grains
            std::set<int> selected_grains;
            for (int i = 0 ; i < nb_grains ; ++i) {
                if (ids[i] >= 0 &&
                        fa[i] >= _mina && fa[i] <= _maxa &&
                        size[i] >= _mins && size[i] <= _maxs) {
                    selected_grains.insert(ids[i]);
                }
            }
            std::cerr << selected_grains.size() << " from " << nb_grains << " grains passed the test.\n";
            // std::cerr << "these grains are " << selected_grains << "'\n";

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

                // std::cerr << "triangle ids are " << ids[0] << ", " << ids[1] << ", " << ids[2] << '\n';
                // std::cerr << "associated tags are:\n";
                // std::cerr << "0: " << vertex_tags[ids[0]] << '\n'
                //           << "1: " << vertex_tags[ids[1]] << '\n'
                //           << "2: " << vertex_tags[ids[2]] << '\n';

                std::set_intersection(vertex_tags[ids[0]].begin(),
                                      vertex_tags[ids[0]].end(),
                                      vertex_tags[ids[1]].begin(),
                                      vertex_tags[ids[1]].end(),
                                      std::back_inserter(tmp1));
                // std::cerr << "first intersection yields: " << tmp1 << '\n';
                std::set_intersection(vertex_tags[ids[2]].begin(),
                                      vertex_tags[ids[2]].end(),
                                      tmp1.begin(), tmp1.end(),
                                      std::back_inserter(tmp2));
                // std::cerr << "second intersection yields: " << tmp2 << '\n';
                std::set_intersection(tmp2.begin(), tmp2.end(),
                                      selected_grains.begin(),
                                      selected_grains.end(),
                                      std::back_inserter(tmp3));
                // std::cerr << "final intersection yields: " << tmp3 << '\n';
                if (!tmp3.size()) {
                    std::cerr << "WARNING: current triangle has invalid corner tags\n";
                } else {
                    int ref_id = tmp3.back(); // largest included index common to all 3 vertices
                    // std::cerr << "ref_id is " << ref_id << " with color " << nvis::ivec3(colors[ref_id]) << '\n';
                    // std::cerr << "inserting color " << colors[ref_id] << " = " << nvis::vec3(colors[ref_id]) / nvis::vec3(255, 255, 255) << '\n';
                    color->InsertNextTypedTuple(colors[ref_id].begin());
                }
            }
            // std::cerr << "after setting colors we have:\n";
            // for (int i=0 ; i<color->GetNumberOfTuples() ; ++i) {
            //  uchar c[3];
            //  color->GetTypedTuple(i, c);
            //  std::cerr << i << ": " << nvis::ivec3(c[0], c[1], c[2]) << '\n';
            // }

            vtkSmartPointer<vtkPolyData> selected_mesh = vtkSmartPointer<vtkPolyData>::New();
            selected_mesh->SetPoints(pts);
            selected_mesh->SetPolys(selected_triangles);
            selected_mesh->GetCellData()->SetScalars(color);

            if (_mina == mina) {
                vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
                writer->SetInputData(selected_mesh);
                writer->SetFileName("/Users/xmt/Desktop/blah.vtk");
                writer->Update();
            }

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
            // mesh_mapper->SetColorModeToMapScalars();
            // mesh_mapper->ScalarVisibilityOn();
            // mesh_mapper->SelectColorArray("Colors");
            // mesh_mapper->SetScalarModeToUseCellData();
            vtkSmartPointer<vtkActor> mesh_actor = vtkSmartPointer<vtkActor>::New();
            mesh_actor->SetMapper(mesh_mapper);
            // mesh_actor->GetProperty()->SetColor(1, 1, 1);
            mesh_actor->GetProperty()->SetOpacity(0.8);

            vtkSmartPointer<vtkPolyData> selected_edges = vtkSmartPointer<vtkPolyData>::New();
            selected_edges->SetPoints(edges->GetPoints());
            selected_edges->SetLines(selected_lines);
            //
            vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
            tubes->SetInputData(selected_edges);

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
            ren_win->Render();

            vtkSmartPointer<vtkWindowToImageFilter> capture = vtkSmartPointer<vtkWindowToImageFilter>::New();
            capture->SetInput(ren_win);

            vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
            writer->SetInputConnection(capture->GetOutputPort());

            std::ostringstream os;
            os << img << "-a-" << _mina*100 << "-" << _maxa*100 << "-s-" << _mins << "-" << _maxs << ".tiff";

            writer->SetFileName(os.str().c_str());
            writer->Write();

            ren->RemoveActor(mesh_actor);
            ren->RemoveActor(edge_actor);
        }
    }
}
