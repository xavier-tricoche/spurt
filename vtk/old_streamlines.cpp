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
#include "vtkStructuredPoints.h"
#include "vtkStructuredPointsReader.h"
#include "vtkScalarBarActor.h"
#include "vtkColorTransferFunction.h"
#include "vtkTextProperty.h"
#include "vtkGeometryFilter.h"

#include <string>
#include <math/fixed_vector.hpp>
#include <vtk/vtk_utils.hpp>
#include <math/bounding_box.hpp>
#include <image/nrrd_wrapper.hpp>
#include <teem/hest_helper.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>
#include <vis/streamline.hpp>
#include <boost/shared_ptr.hpp>
#include <data/field_wrapper.hpp>
#include "Garcia_vis_helper.hpp"
#include "math/inverse_transform.hpp"

inline void wait(int s)
{
    nvis::timer t;
    while (t.elapsed() < s) {}
}

int     dataset;
char*    scalar;
int     nblines;
float   length;
float   center[3];
float   radius[3];
float   eps;
float   h;
float   lmax;
float   min_step;
float   col_gamma;
int     discretization;

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
    hestOptAdd(&hopt, "id",     "dataset ID",       airTypeInt,     1, 1, &dataset,             NULL,       "dataset ID: 0: textured, 1: untextured, 2: MC_r00b09, 3: MC_r06b09, 4: MC_r10b09");
    hestOptAdd(&hopt, "s",      "scalar",           airTypeString,  0, 1, &scalar,              "none",     "scalar field used for color coding (NRRD)");
    hestOptAdd(&hopt, "n",      "nb lines",         airTypeInt,     0, 1, &nblines,             "100",      "number of lines");
    hestOptAdd(&hopt, "l",      "length",           airTypeFloat,   1, 1, &length,              NULL,       "integration length");
    hestOptAdd(&hopt, "c",      "center",           airTypeFloat,   3, 3, &center,              NULL,       "seed center");
    hestOptAdd(&hopt, "r",      "radius",           airTypeFloat,   3, 3, &radius,              NULL,       "seed radius");
    hestOptAdd(&hopt, "e",      "eps",              airTypeFloat,   0, 1, &eps,                 "1.0e-6",   "integration precision");
    hestOptAdd(&hopt, "max",    "max step",         airTypeFloat,   0, 1, &lmax,                "0.5",      "max integration step length");
    hestOptAdd(&hopt, "h",      "hinit",            airTypeFloat,   0, 1, &h,                   "1.",       "initial step size guess");
    hestOptAdd(&hopt, "min",    "min step",         airTypeFloat,   0, 1, &min_step,            "1.0e-8",   "step size underflow threshold");
    hestOptAdd(&hopt, "d",      "discretization",   airTypeInt,     0, 1, &discretization,      "20",       "number of lines");
    hestOptAdd(&hopt, "g",      "gamma",            airTypeFloat,   0, 1, &col_gamma,           "1.",       "gamma factor for color scale");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Compute and visualize streamlines in 3D vector field",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline void add_line(std::list<nvis::vec3>& all_verts,
                     std::list<double>& all_vals,
                     std::list<int>& start_at,
                     const std::list<nvis::vec3>& line,
                     const std::list<double>& vals)
{
    start_at.push_back(all_verts.size());
    all_verts.insert(all_verts.end(), line.begin(), line.end());
    all_vals.insert(all_vals.end(), vals.begin(), vals.end());
    
    // std::cerr << "current line starts at " << start_at.back() << '\n';
    // std::cerr << "there are " << all_verts.size() << " total vertices and "
    //           << all_vals.size() << " total values stored\n";
}

template<typename T>
inline T sign(const T& t)
{
    return (t < 0 ? -1 : 1);
}

struct check {

    check(size_t max_count, const nvis::bbox3& box)
        : _counter(0), _max_count(max_count), _bounds(box) {}
        
    void reset() {
        _counter = 0;
    }
    
    bool operator()(const nvis::streamline::int_step& step) {
        // std::cerr << _counter << nvis::subv<0, 3>(step.y1()) << '\n';
        
        if (++_counter >= _max_count) {
            return true;
        } else if (!_bounds.inside(nvis::subv<0, 3>(step.y1()))) {
            return true;
        }
        
        // std::cerr << "continue\n";
        return false;
    }
    
    size_t _counter, _max_count;
    nvis::bbox3 _bounds;
};

struct i2x {
    i2x(const Nrrd* nrrd) {
        for (int i = 0 ; i < 3 ; ++i) {
            step[i] = nrrd->axis[nrrd->dim-3+i].spacing;
            size[i] = nrrd->axis[nrrd->dim-3+i].size;
        }
        bounds = spurt::nrrd_utils::get_bounds<3>(nrrd);
    }
    
    nvis::vec3 operator()(int id) const {
        int i = id % size[0];
        id /= size[0];
        int j = id % size[1];
        int k = id / size[1];
        return bounds.min() + nvis::vec3(i, j, k)*step;
    }
    
    nvis::vec3  step;
    nvis::ivec3 size;
    nvis::bbox3 bounds;
};

template<typename T>
struct right_hand_side {
    typedef T   field_type;
    
    right_hand_side(const field_type& field)
        : _field(field) {}
        
    bool operator()(const nvis::vec3& x, nvis::vec3& f) const {
        return _field.get_value(x, f);
    }
    
    const field_type& _field;
};

template<typename RHS>
struct euler {

    typedef RHS     rhs_type;
    
    enum state {
        OK = 0,
        LEFT_DOMAIN,
        STOPPED,
    };
    
    euler(const rhs_type& rhs, double h, double lmax, double eps = 1.0e-6) : _rhs(rhs), _h(h), _lmax(lmax), _eps(eps) {}
    
    state advance(const nvis::vec3& in, std::list<nvis::vec3>& out, double length, double& actual_length) const {
        double h = (length < 0 ? -_h : _h);
        actual_length = 0;
        nvis::vec3 f, x = in;
        while (actual_length < fabs(length)) {
            if (_rhs(x, f)) {
                out.push_back(x);
                if (nvis::norm(f) < _eps) {
                    return STOPPED;
                } else if (fabs(h)*nvis::norm(f) > _lmax) {
                    h = sign(length) * 0.5 * _lmax / nvis::norm(f);
                    x += h * f;
                } else {
                    x += h * f;
                    if (fabs(h)*nvis::norm(f) < 10.*_lmax) {
                        h *= 10;
                    }
                }
                actual_length += fabs(h);
            } else {
                return LEFT_DOMAIN;
            }
        }
        return OK;
    }
    
    const rhs_type& _rhs;
    double _h, _lmax, _eps;
};


typedef spurt::nrrd_data_traits<Nrrd*>  field_type;
typedef right_hand_side<field_type> rhs_type;
typedef euler<rhs_type>             euler_type;

inline vtkSmartPointer<vtkPolyData> to_polydata(const std::list<nvis::vec3>& line)
{
    vtkSmartPointer<vtkDoubleArray> coords = vtkSmartPointer<vtkDoubleArray>::New();
    coords->SetNumberOfComponents(3);
    for (std::list<nvis::vec3>::const_iterator it = line.begin() ;
            it != line.end() ; ++it) {
        coords->InsertNextTuple(it->begin());
    }
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetData(coords);
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType n = 2;
    for (int i = 0 ; i < line.size() - 1 ; ++i) {
        vtkIdType ids[] = { i, i + 1 };
        lines->InsertNextCell(n, ids);
    }
    
    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(pts);
    pd->SetLines(lines);
    return pd;
}

inline vtkSmartPointer<vtkPolyData>
to_tubes(const std::list<nvis::vec3>& all_verts,
         const std::list<double>& all_vals,
         const std::list<int>& start_at, double radius)
{
    vtkSmartPointer<vtkDoubleArray> coords = vtkSmartPointer<vtkDoubleArray>::New();
    coords->SetNumberOfComponents(3);
    for (std::list<nvis::vec3>::const_iterator it = all_verts.begin() ;
            it != all_verts.end() ; ++it) {
        coords->InsertNextTuple(it->begin());
    }
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetData(coords);
    
    std::vector<int> ids(start_at.begin(), start_at.end());
    ids.push_back(all_verts.size());
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    for (int i = 0 ; i < ids.size() - 1 ; ++i) {
        vtkIdType n = ids[i+1] - ids[i];
        vtkIdType v[n];
        for (int j = ids[i] ; j < ids[i+1] ; ++j) {
            v[j-ids[i]] = j;
        }
        lines->InsertNextCell(n, v);
    }
    
    vtkSmartPointer<vtkDoubleArray> values = vtkSmartPointer<vtkDoubleArray>::New();
    values->SetNumberOfComponents(1);
    values->SetName("Scalars");
    for (std::list<double>::const_iterator it = all_vals.begin() ;
            it != all_vals.end() ; ++it) {
        values->InsertNextTypedTuple(&(*it));
    }
    
    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(pts);
    pd->SetLines(lines);
    pd->GetPointData()->SetScalars(values);
    
    vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
    tubes->SetInput(pd);
    tubes->SetRadius(radius);
    tubes->SetNumberOfSides(12);
    // tubes->SetVaryRadiusToVaryRadiusByScalar();
    // tubes->SetRadiusFactor(10);
    
    return tubes->GetOutput();
}

typedef nvis::vec3                              value_type;
typedef spurt::nrrd_data_traits<Nrrd*>              nrrd_data_traits;

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    using namespace Garcia_vis_helper;
    
    set_paths();
    
    dataset_info* __info;
    switch (dataset) {
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
    
    Nrrd* __coords = spurt::nrrd_utils::readNrrd(point_coord_name);
    spurt::nrrd_utils::nrrd_data_wrapper<float> coords(__coords);
    
    vtkSmartPointer<vtkDataSetReader> mesh_reader = vtkSmartPointer<vtkDataSetReader>::New();
    mesh_reader->SetFileName(mesh_name.c_str());
    mesh_reader->Update();
    vtkSmartPointer<vtkUnstructuredGrid> mesh = vtkSmartPointer<vtkUnstructuredGrid>::New();
    mesh->DeepCopy(mesh_reader->GetOutput());
    std::cerr << "mesh contains " << std::flush << mesh->GetNumberOfPoints() << " points\n";
    std::cerr << "mesh contains " << std::flush << mesh->GetNumberOfCells() << " cells\n";
    mesh_reader->Delete();
    std::cerr << mesh_name << " loaded.\n";
    
    vtkSmartPointer<vtkDoubleArray> coord = vtkSmartPointer<vtkDoubleArray>::New();
    coord->SetNumberOfComponents(3);
    for (int i = 0 ; i < mesh->GetNumberOfPoints() ; ++i) {
        double pt[3];
        mesh->GetPoint(i, pt);
        coord->InsertNextTuple(pt);
    }
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetData(coord);
    
    vtkSmartPointer<vtkPolyDataReader> edge_reader = vtkSmartPointer<vtkPolyDataReader>::New();
    edge_reader->SetFileName(edge_name.c_str());
    edge_reader->Update();
    vtkSmartPointer<vtkPolyData> edges = vtkSmartPointer<vtkPolyData>::New();
    edges->DeepCopy(edge_reader->GetOutput());
    vtkCellArray* lines = edges->GetLines();
    std::cerr << "there are " << lines->GetNumberOfCells() << " edges in input\n";
    edge_reader->Delete();
    std::cerr << edge_name << " loaded.\n";
    
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
    
    // compute bounding box of each grain
    std::map<int, nvis::bbox3> grain_bounds;
    for (int i = 0 ; i < mesh->GetNumberOfPoints() ; ++i) {
        double pt[3];
        mesh->GetPoint(i, pt);
        nvis::vec3 x(pt[0], pt[1], pt[2]);
        for (tag_type::const_iterator it = vertex_tags[i].begin() ; it != vertex_tags[i].end() ; ++it) {
            int id = *it;
            if (grain_bounds.find(id) == grain_bounds.end()) {
                grain_bounds[id] = nvis::bbox3(x, x);
            } else {
                grain_bounds[id].add(x);
            }
        }
    }
    int nb_grains = grain_bounds.size();
    
    std::string name = info.base_dir + "Efield.nrrd";
    Nrrd* nin_vec = spurt::nrrd_utils::readNrrd(name);
    if (nin_vec->dim != 4) {
        throw;
    }
    field_type      vf(nin_vec);
    rhs_type        rhs(vf);
    euler_type      integrator(rhs, h, lmax);
    
    std::vector<double> scalars;
    std::vector<double> tmp_scl;
    spurt::to_vector(tmp_scl, nin_vec);
    std::cerr << tmp_scl.size() / 3 << " vector values\n";
    scalars.resize(tmp_scl.size() / 3);
    for (int i = 0 ; i < tmp_scl.size() / 3 ; ++i) {
        nvis::vec3 v(tmp_scl[3*i], tmp_scl[3*i+1], tmp_scl[3*i+2]);
        scalars[i] = nvis::norm(v);
    }
    std::vector<double> tmp(scalars.begin(), scalars.end());
    std::sort(tmp.begin(), tmp.end());
    std::cerr << "min norm = " << tmp.front() << ", max norm = " << tmp.back() << ", median norm = " << tmp[tmp.size()/2] << '\n';
    
    spurt::inverse_transform_sampling<double> itf(scalars);
    
    typedef Garcia_vis_helper::color_map<double>    color_map;
    typedef color_map::color_type                   color_type;
    color_map cmap(scalars, col_gamma, false);
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    {
        double min = *std::min_element(scalars.begin(), scalars.end());
        double max = *std::max_element(scalars.begin(), scalars.end());
        double span = max - min;
        
        std::cerr << "min scalar = " << min << ", max scalar = " << max << '\n';
        
        double dv = span / 100.;
        for (double v = min ; v <= max ; v += dv) {
            color_type c = cmap(v, color_map::REDUNDANT_RAINBOW);
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
    
    vtkSmartPointer<vtkRenderer> ren = vtkSmartPointer<vtkRenderer>::New();
    ren->SetUseDepthPeeling(1);
    ren->SetMaximumNumberOfPeels(100);
    ren->SetOcclusionRatio(0.1);
    ren->SetBackground(0, 0, 0);
    ren->ResetCamera();
    
    ren->GetActiveCamera()->SetPosition(info.position.begin());
    ren->GetActiveCamera()->SetFocalPoint(info.focal_point.begin());
    ren->GetActiveCamera()->SetViewUp(info.up.begin());
    ren->GetActiveCamera()->SetClippingRange(info.near, info.far);
    std::cerr << "Camera set\n";
    
    ren->AddActor2D(color_bar);
    
    vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    ren_win->SetAlphaBitPlanes(1);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1600, 1200);
    
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);
    
    nvis::bbox3 bounds = spurt::nrrd_utils::get_bounds<3>(nin_vec);
    vtkSmartPointer<vtkActor> frame_actor = Garcia_vis_helper::draw_frame(bounds);
    ren->AddActor(frame_actor);
    
    nvis::vec3 _center(center[0], center[1], center[2]);
    nvis::vec3 _radius(radius[0], radius[1], radius[2]);
    _center *= bounds.size();
    _radius *= bounds.size();
    nvis::vec3 min = _center - _radius;
    nvis::vec3 diameter = 2 * _radius;
    
    std::list<nvis::vec3>   all_verts;
    std::list<double>       all_vals;
    std::list<nvis::bbox3>  all_bounds;
    std::list<int>          start_at;
    
    nvis::vec3 diag = bounds.size();
    double radius = 0.01* *std::max_element(diag.begin(), diag.end());
    
    srand48(time(0));
    double intg_time = 0;
    for (int i = 0 ; i < nblines ; ++i) {
    
        nvis::vec3 seed = min + nvis::vec3(drand48(), drand48(), drand48()) * diameter;
        // int seed_id = itf.sample();
        // nvis::vec3 seed = to_pos(seed_id);
        
        std::cerr << "seeding streamline #" << i << "/" << nblines << " at " << seed << '\n';
        
        std::list<nvis::vec3>   line, aux;
        double                  actual_l;
        nvis::timer             one_timer;
        
        nvis::bbox3 box(seed, seed);
        
        integrator.advance(seed, line, length, actual_l);
        integrator.advance(seed, aux, -length, actual_l);
        intg_time += one_timer.elapsed();
        aux.pop_front();
        line.insert(line.begin(), aux.rbegin(), aux.rend());
        std::cerr << line.size() << " points in line\n";
        
        std::list<double>   vals;
        for (std::list<nvis::vec3>::const_iterator it = line.begin() ; it != line.end() ; ++it) {
            nvis::vec3 x = *it;
            box.add(x);
            double scl;
            nvis::vec3 f;
            rhs(x, f);
            scl = nvis::norm(f);
            vals.push_back(scl);
        }
        all_bounds.push_back(box);
        add_line(all_verts, all_vals, start_at, line, vals);
    }
    vtkSmartPointer<vtkPolyData> tubes =
        to_tubes(all_verts, all_vals, start_at, radius);
        
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput(tubes);
    mapper->SetLookupTable(ctf);
    mapper->ScalarVisibilityOn();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    ren->AddActor(actor);
    
    int frame_counter = 0;
    for (double y = bounds.min()[1] ; y <= bounds.max()[1] ; y += 10) {
        std::cerr << " y = " << y << '\n';
        
        // included grains
        std::set<int> selected_grains;
        for (std::map< int, nvis::bbox3 >::const_iterator it = grain_bounds.begin() ;
                it != grain_bounds.end() ; ++it) {
            if (it->first < 0 || it->second.min()[1] < y) {
                continue;
            }
            selected_grains.insert(it->first);
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
        
        //
        vtkSmartPointer<vtkPolyDataMapper> mesh_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mesh_mapper->SetInput(selected_mesh);
        mesh_mapper->ScalarVisibilityOff();
        vtkSmartPointer<vtkActor> mesh_actor = vtkSmartPointer<vtkActor>::New();
        mesh_actor->SetMapper(mesh_mapper);
        mesh_actor->GetProperty()->SetColor(1, 1, 1);
        mesh_actor->GetProperty()->SetOpacity(0.75);
        
        vtkSmartPointer<vtkPolyData> selected_edges = vtkSmartPointer<vtkPolyData>::New();
        selected_edges->SetPoints(edges->GetPoints());
        selected_edges->SetLines(selected_lines);
        vtkSmartPointer<vtkTubeFilter> etubes = vtkSmartPointer<vtkTubeFilter>::New();
        etubes->SetInput(selected_edges);
        
        double eradius = 0.05 * nvis::norm(info.step);
        etubes->SetRadius(eradius);
        etubes->SetNumberOfSides(6);
        vtkSmartPointer<vtkPolyDataMapper> edge_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        edge_mapper->SetInputConnection(etubes->GetOutputPort());
        vtkSmartPointer<vtkActor> edge_actor = vtkSmartPointer<vtkActor>::New();
        edge_actor->SetMapper(edge_mapper);
        edge_actor->GetProperty()->SetColor(0.5, 0, 0);
        
        ren->AddActor(mesh_actor);
        ren->AddActor(edge_actor);
        ren_win->Render();
        
        vtkSmartPointer<vtkWindowToImageFilter> capture = vtkSmartPointer<vtkWindowToImageFilter>::New();
        capture->SetInput(ren_win);
        
        vtkSmartPointer<vtkTIFFWriter> writer = vtkSmartPointer<vtkTIFFWriter>::New();
        writer->SetInputConnection(capture->GetOutputPort());
        
        std::ostringstream os;
        os << info.movie_path << "streamlines_Y_frame_" << std::setw(3) << std::setfill('0') << frame_counter++ << "_of_"
           << (int)floor(bounds.size()[1] / 10.) + 1 << ".tiff";
           
        writer->SetFileName(os.str().c_str());
        std::cerr << "about to write to file " << os.str() << '\n';
        writer->Write();
        
        ren->RemoveActor(mesh_actor);
        ren->RemoveActor(edge_actor);
    }
    
    if (nin_vec) {
        nrrdNuke(nin_vec);
    }
    
    return 0;
}












































































































