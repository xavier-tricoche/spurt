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
#include "vtkSmartPointer.h"

#include <string>
#include <math/fixed_vector.hpp>
#include <vtk/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <set>
#include <sstream>
#include "Garcia_vis_helper.hpp"

char* mesh, *stat, *img;
bool verbose, do_screenshot;
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
    hestOptAdd(&hopt, "maxa",   "max aniso",        airTypeDouble,  0, 1, &maxa,            "1",        "max grain anisotropy");
    hestOptAdd(&hopt, "mins",   "min size",         airTypeInt,     0, 1, &mins,            "0",        "min grain size");
    hestOptAdd(&hopt, "maxs",   "max size",         airTypeInt,     0, 1, &maxs,            "-1",       "max grain size");
    hestOptAdd(&hopt, "ss",     "screen shot",      airTypeBool,    0, 0, &do_screenshot,   "0",        "save screen shot in hard-coded position (batch mode)");
    hestOptAdd(&hopt, "v",      "verbose",          airTypeBool,    0, 0, &verbose,         "0",        "verbose mode (debugging)");
    hestOptAdd(&hopt, "ssf",    "screen shot file", airTypeString,  0, 1, &img,             "/tmp/img", "screenshot base name (batch mode)");

    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize granular microstructure",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

struct vertex_info {
    int             id;
    std::set<int>   tags;
    float           aniso;
    int             size;
};

struct Lt_vertex_aniso {
    int operator()(const vertex_info& i0, const vertex_info& i1) {
        return i0.aniso < i1.aniso;
    }
};

struct Lt_vertex_size {
    int operator()(const vertex_info& i0, const vertex_info& i1) {
        return i0.size < i1.size;
    }
};

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    if (maxs == -1) {
        maxs = std::numeric_limits<int>::max();
    }
    std::cerr << "verbose is " << (verbose ? "true" : "false") << '\n';
    std::cerr << "do_screenshot is " << (do_screenshot ? "true" : "false") << '\n';

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
    std::cerr << mesh_name << " loaded.\n";

    vtkPolyDataReader* edge_reader = vtkPolyDataReader::New();
    edge_reader->SetFileName(edge_name.c_str());
    edge_reader->Update();
    std::cerr << edge_name << " loaded.\n";

    vtkPolyDataReader* corner_reader = vtkPolyDataReader::New();
    corner_reader->SetFileName(corner_name.c_str());
    corner_reader->Update();
    std::cerr << corner_name << " loaded.\n";

    std::string stat_base(stat);
    std::string fa_name(stat_base), size_name(stat_base), id_name(stat_base);
    fa_name.append("-fa.nrrd");
    size_name.append("-size.nrrd");
    id_name.append("-ids_to_tags.nrrd");

    Nrrd* __ids = spurt::nrrd_utils::readNrrd(id_name);
    spurt::nrrd_utils::nrrd_data_wrapper<int> ids(__ids);
    std::cerr << id_name << " loaded.\n";

    Nrrd* __fa = spurt::nrrd_utils::readNrrd(fa_name);
    spurt::nrrd_utils::nrrd_data_wrapper<float> fa(__fa);
    std::cerr << fa_name << " loaded.\n";
    int nb_grains = __fa->axis[0].size;

    Nrrd* __size = spurt::nrrd_utils::readNrrd(size_name);
    spurt::nrrd_utils::nrrd_data_wrapper<float> size(__size);
    std::cerr << size_name << " loaded.\n";

    if (false) {
        int* sizes = (int*)__size->data;
        std::vector<int> __s;
        std::copy(sizes, sizes + nb_grains, std::back_inserter(__s));
        std::sort(__s.begin(), __s.end());
        std::cerr << "size ranges from " << __s[0] << " and " << __s.back()
                  << ", median value is " << __s[__s.size()/2] << '\n';

        float* aniso = (float*)__fa->data;
        std::vector<float> __a;
        std::copy(aniso, aniso + nb_grains, std::back_inserter(__a));
        std::sort(__a.begin(), __a.end());
        std::cerr << "FA ranges from " << __a[0] << " and " << __a.back()
                  << ", median value is " << __a[__a.size()/2] << '\n';
    }

    std::cerr << "anisotropy bounds are: " << mina << " - " << maxa << '\n';
    std::cerr << "size bounds are: " << mins << " - " << maxs << '\n';

    std::set<int> selected_grains;
    for (int i = 0 ; i < nb_grains ; ++i) {
        if (ids[i] < 0) {
            continue;
        }
        if (fa[i] >= mina && fa[i] <= maxa &&
                size[i] >= mins && size[i] <= maxs) {
            selected_grains.insert(ids[i]);
        }
    }
    std::cerr << selected_grains.size() << " from " << nb_grains << " grains have passed the test.\n";

    std::set<int> selected_vertices;
    std::fstream attributes(vertex_tag_name.c_str());
    std::cerr << vertex_tag_name << " file opened\n";
    while (!attributes.eof()) {
        int i, n, id;
        attributes >> i >> n;
        for (int j = 0 ; j < n ; ++j) {
            attributes >> id;
            if (selected_grains.find(id) != selected_grains.end()) {
                selected_vertices.insert(i);
            }
        }
        // std::cerr << "read vertex #" << i << '\n';
    }
    attributes.close();
    std::cerr << selected_vertices.size() << " vertices passed the test.\n";

    // selected triangles
    std::cerr << "mesh contains " << std::flush << mesh_reader->GetOutput()->GetNumberOfPoints() << " points\n";
    std::cerr << "mesh contains " << std::flush << mesh_reader->GetOutput()->GetNumberOfCells() << " cells\n";
    vtkCellArray* selected_triangles = vtkCellArray::New();
    vtkDoubleArray* coords = vtkDoubleArray::New();
    coords->SetNumberOfComponents(3);
    for (int i = 0 ; i < mesh_reader->GetOutput()->GetNumberOfPoints() ; ++i) {
        double pt[3];
        mesh_reader->GetOutput()->GetPoint(i, pt);
        coords->InsertNextTuple(pt);
    }
    vtkPoints* pts = vtkPoints::New();
    pts->SetData(coords);

    for (int n = 0 ; n < mesh_reader->GetOutput()->GetNumberOfCells() ; ++n) {
        vtkIdList* ids = vtkIdList::New();
        mesh_reader->GetOutput()->GetCellPoints(n, ids);
        bool included = true;
        vtkIdType _ids[ids->GetNumberOfIds()];
        for (int i = 0 ; i < ids->GetNumberOfIds() && included ; ++i) {
            vtkIdType id = ids->GetId(i);
            _ids[i] = id;
            included = (selected_vertices.find(id) != selected_vertices.end());
        }
        if (included) {
            selected_triangles->InsertNextCell(ids->GetNumberOfIds(), _ids);
        }
    }
    std::cerr << selected_triangles->GetNumberOfCells() << " triangles passed the test\n";

    vtkCellArray* selected_lines = vtkCellArray::New();
    vtkCellArray* lines = edge_reader->GetOutput()->GetLines();
    std::cerr << "there are " << lines->GetNumberOfCells() << " edges in input\n";
    while (true) {
        vtkIdType _n;
        const vtkIdType* _ids;
        if (!lines->GetNextCell(_n, _ids)) {
            break;
        }
        bool included = true;
        for (int i = 0 ; i < _n && included ; ++i) {
            vtkIdType id = _ids[i];
            included = (selected_vertices.find(id) != selected_vertices.end());
        }
        if (included) {
            selected_lines->InsertNextCell(_n, _ids);
        } else {
            // std::cerr << "edge " << _ids[0] << " - " << _ids[1] << " is not contained\n";
        }
    }
    std::cerr << "there are " << selected_lines->GetNumberOfCells() << " selected edges"
              << " of which the largest has cardinal " << selected_lines->GetMaxCellSize() << "\n";

    vtkPolyData* selected_mesh = vtkPolyData::New();
    selected_mesh->SetPoints(pts);
    selected_mesh->SetPolys(selected_triangles);
    //
    vtkPolyDataMapper* mesh_mapper = vtkPolyDataMapper::New();
    mesh_mapper->SetInputData(selected_mesh);
    vtkActor* mesh_actor = vtkActor::New();
    mesh_actor->SetMapper(mesh_mapper);
    mesh_actor->GetProperty()->SetColor(1, 1, 1);
    mesh_actor->GetProperty()->SetOpacity(0.5);

    vtkPolyData* selected_edges = vtkPolyData::New();
    selected_edges->SetPoints(edge_reader->GetOutput()->GetPoints());
    selected_edges->SetLines(selected_lines);
    //
    vtkTubeFilter* tubes = vtkTubeFilter::New();
    tubes->SetInputData(selected_edges);
    tubes->SetRadius(2.0e-8);
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
    ren->AddActor(mesh_actor);
    ren->AddActor(edge_actor);
    ren->SetBackground(0, 0, 0);
    ren->ResetCamera();

    if (verbose) {
        std::cerr << "verbose mode: camera setting will be exported to stdio\n";
        vtk_utils::camera_setting_callback* cb = vtk_utils::camera_setting_callback::New();
        ren->AddObserver(vtkCommand::StartEvent, cb);
        cb->Delete();
    }

    // BNTK_01
    // ren->GetActiveCamera()->SetPosition(-9.09581e-06, 1.4439e-05, -8.51636e-06);
    // ren->GetActiveCamera()->SetFocalPoint(4.87001e-06, 4.86777e-06, 3.45e-06);
    // ren->GetActiveCamera()->SetViewUp(0.384824, -0.447644, -0.807171);
    // ren->GetActiveCamera()->SetClippingRange(8.07201e-06, 3.672e-05);

    // Textured_09
    ren->GetActiveCamera()->SetPosition(-36.3331, -8.47804, 28.695);
    ren->GetActiveCamera()->SetFocalPoint(16.9142, 19.7463, 3.89301);
    ren->GetActiveCamera()->SetViewUp(0.291739, 0.258784, 0.920825);
    ren->GetActiveCamera()->SetClippingRange(18.0374, 124.656);

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

    if (do_screenshot) {
        vtkWindowToImageFilter* capture = vtkWindowToImageFilter::New();
        capture->SetInput(ren_win);

        vtkTIFFWriter* writer = vtkTIFFWriter::New();
        writer->SetInputConnection(capture->GetOutputPort());

        std::ostringstream os;
        os << img << "-a-" << mina*100 << "-" << maxa*100 << "-s-" << mins << "-" << maxs << ".tiff";

        writer->SetFileName(os.str().c_str());
        writer->Write();
    } else {
        iren->Initialize();
        iren->Start();
    }

    mesh_actor->Delete();
    edge_actor->Delete();
    ren_win->Delete();
    ren->Delete();

    return 0;
}
