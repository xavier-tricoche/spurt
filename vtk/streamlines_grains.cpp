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
#include <VTK/vtk_utils.hpp>
#include <math/bounding_box.hpp>
#include <image/nrrd_wrapper.hpp>
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
char*    density;
char*    movie_base;
char*    grain_norm;
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
int     nf;

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
    hestOptAdd(&hopt, "d",      "density",          airTypeString,  1, 1, &density,             NULL,       "per-grain scalar field used to filter grains (NRRD)");
    hestOptAdd(&hopt, "gn",     "norm",             airTypeString,  1, 1, &grain_norm,          NULL,       "per-grain norm of vector field (NRRD)");
    hestOptAdd(&hopt, "nf",     "nb frames",        airTypeInt,     0, 1, &nf,                  "100",      "number of animation frames");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &movie_base,          NULL,       "frame base name (including path)");
    hestOptAdd(&hopt, "n",      "nb lines",         airTypeInt,     0, 1, &nblines,             "100",      "number of lines");
    hestOptAdd(&hopt, "l",      "length",           airTypeFloat,   1, 1, &length,              NULL,       "integration length");
    hestOptAdd(&hopt, "e",      "eps",              airTypeFloat,   0, 1, &eps,                 "1.0e-6",   "integration precision");
    hestOptAdd(&hopt, "max",    "max step",         airTypeFloat,   0, 1, &lmax,                "0.5",      "max integration step length");
    hestOptAdd(&hopt, "h",      "hinit",            airTypeFloat,   0, 1, &h,                   "1.",       "initial step size guess");
    hestOptAdd(&hopt, "min",    "min step",         airTypeFloat,   0, 1, &min_step,            "1.0e-8",   "step size underflow threshold");
    hestOptAdd(&hopt, "ds",     "discretization",   airTypeInt,     0, 1, &discretization,      "20",       "number of lines");
    hestOptAdd(&hopt, "g",      "gamma",            airTypeFloat,   0, 1, &col_gamma,           "1.",       "gamma factor for color scale");

    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute and visualize streamlines in 3D vector field",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
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
        bounds = xavier::nrrd_utils::get_bounds<3>(nrrd);
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


typedef xavier::nrrd_data_traits<Nrrd*>  field_type;
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

typedef std::pair<nvis::vec3, double>   curve_point;
typedef std::list<curve_point>          curve_type;

inline vtkSmartPointer<vtkPolyData>
to_tubes(const std::list<curve_type>& curves,
         double radius)
{
    vtkSmartPointer<vtkDoubleArray> coords = vtkSmartPointer<vtkDoubleArray>::New();
    coords->SetNumberOfComponents(3);
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkSmartPointer<vtkDoubleArray> values = vtkSmartPointer<vtkDoubleArray>::New();
    values->SetNumberOfComponents(1);
    values->SetName("Scalars");
    std::list<curve_type>::const_iterator cit;
    int n = 0;
    for (cit = curves.begin() ; cit != curves.end() ; ++cit) {
        const curve_type& curve = *cit;
        vtkIdType ids[curve.size()];
        int m = 0;
        for (curve_type::const_iterator it = curve.begin() ; it != curve.end() ; ++it) {
            nvis::vec3 x = it->first;
            x += nvis::vec3(20, 20, 20);
            double val = it->second;
            coords->InsertNextTuple(x.begin());
            ids[m++] = n++;
            values->InsertNextTypedTuple(&val);
        }
        lines->InsertNextCell(curve.size(), ids);
    }
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetData(coords);

    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(pts);
    pd->SetLines(lines);
    pd->GetPointData()->SetScalars(values);

    vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
    tubes->SetInputData(pd);
    tubes->SetRadius(radius);
    tubes->SetNumberOfSides(12);
    // tubes->SetVaryRadiusToVaryRadiusByScalar();
    // tubes->SetRadiusFactor(10);

    return tubes->GetOutput();
}

typedef nvis::vec3                              value_type;
typedef xavier::nrrd_data_traits<Nrrd*>              nrrd_data_traits;
typedef Garcia_vis_helper::color_map<double>    color_map_type;
typedef color_map_type::color_type              color_type;

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

    std::map<int, float>    grain_vals;
    std::vector<double>     all_grain_vals;
    {
        std::cerr << "opening " << density << '\n';
        Nrrd* nin_val = xavier::nrrd_utils::readNrrd(density);
        std::vector<float> vals;
        xavier::nrrd_utils::to_vector(vals, nin_val);
        std::cerr << vals.size() << " values in file\n";
        int n = vals.size() / 2;
        for (int i = 0 ; i < n ; ++i) {
            int id = vals[2*i];
            float v = vals[2*i+1];
            if (id < 0) {
                continue;
            }
            grain_vals[id] = v;
            all_grain_vals.push_back(v);
        }
        std::cerr << all_grain_vals.size() << " values recorded\n";
    }
    std::sort(all_grain_vals.begin(), all_grain_vals.end());
    std::vector<double> thresholds(nf+1);
    int delta_i = all_grain_vals.size()/(nf-1);
    for (int i=0 ; i<nf ; ++i) {
        thresholds[i] = all_grain_vals[i*delta_i];
    }
    thresholds.push_back(all_grain_vals.back() + 1);

    std::map<int, float>    grain_norms;
    {
        std::cerr << "opening " << grain_norm << '\n';
        Nrrd* nin_val = xavier::nrrd_utils::readNrrd(grain_norm);
        std::vector<float> vals;
        xavier::nrrd_utils::to_vector(vals, nin_val);
        std::cerr << vals.size() << " values in file\n";
        int n = vals.size() / 2;
        for (int i = 0 ; i < n ; ++i) {
            int id = vals[2*i];
            float v = vals[2*i+1];
            if (id < 0) {
                continue;
            }
            grain_norms[id] = v;
        }
    }

    std::string name = info.base_dir + "Efield.nrrd";
    Nrrd* nin_vec = xavier::nrrd_utils::readNrrd(name);
    if (nin_vec->dim != 4) {
        throw;
    }
    field_type      vf(nin_vec);
    rhs_type        rhs(vf);
    euler_type      integrator(rhs, h, lmax);

    std::vector<double> scalars;
    std::vector<double> tmp_scl;
    xavier::nrrd_utils::to_vector(tmp_scl, nin_vec);
    std::cerr << tmp_scl.size() / 3 << " vector values\n";
    scalars.resize(tmp_scl.size() / 3);
    for (int i = 0 ; i < tmp_scl.size() / 3 ; ++i) {
        nvis::vec3 v(tmp_scl[3*i], tmp_scl[3*i+1], tmp_scl[3*i+2]);
        scalars[i] = nvis::norm(v);
    }
    std::vector<double> tmp(scalars.begin(), scalars.end());
    std::sort(tmp.begin(), tmp.end());
    std::cerr << "min norm = " << tmp.front() << ", max norm = " << tmp.back() << ", median norm = " << tmp[tmp.size()/2] << '\n';

    xavier::inverse_transform_sampling<double> itf(scalars);

    color_map_type cmap(tmp, col_gamma, false);
    vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
    {
        double min = *std::min_element(tmp.begin(), tmp.end());
        double max = *std::max_element(tmp.begin(), tmp.end());
        double span = max - min;

        std::cerr << "min scalar = " << min << ", max scalar = " << max << '\n';

        double dv = span / 100.;
        for (double v = min ; v <= max ; v += dv) {
            color_type c = cmap(v, color_map_type::REDUNDANT_RAINBOW);
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
    // ren->AddActor2D(color_bar2);

    vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    ren_win->SetAlphaBitPlanes(1);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1600, 1200);

    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);

    nvis::bbox3 bounds = xavier::nrrd_utils::get_bounds<3>(nin_vec);
    bounds.min() += nvis::vec3(20, 20, 20);
    bounds.max() += nvis::vec3(20, 20, 20);
    std::cerr << "frame bounds are " << bounds << '\n';
    vtkSmartPointer<vtkActor> frame_actor = Garcia_vis_helper::draw_frame(bounds);
    bounds.min() -= nvis::vec3(20, 20, 20);
    bounds.max() -= nvis::vec3(20, 20, 20);
    ren->AddActor(frame_actor);

    nvis::vec3 diag = bounds.size();
    double radius = 0.005* *std::max_element(diag.begin(), diag.end());

    srand48(time(0));
    double intg_time = 0;

    std::list<curve_type>   all_curves;
    std::list<nvis::bbox3>  all_bounds;

    int nx, ny;
    nx = ny = sqrt(nblines);
    nblines = nx * ny;
    double dx = diag[0] / (double)(nx + 1);
    double dy = diag[1] / (double)(ny + 1);

    for (int i = 0 ; i < nblines ; ++i) {

        int _i = i % nx;
        int _j = i / nx;

        double x = bounds.min()[0] + dx * (_i + 0.5);
        double y = bounds.min()[1] + dy * (_j + 0.5);

        nvis::vec3 seed(x, y, 0.001);

        // nvis::vec3 seed = min + nvis::vec3(drand48(), drand48(), drand48()) * diameter;
        // int seed_id = itf.sample();
        // nvis::vec3 seed = to_pos(seed_id);

        // nvis::vec3 seed = bounds.min() + nvis::vec3(drand48(), drand48(), 0.001) * bounds.size();

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

        curve_type      curve;
        for (std::list<nvis::vec3>::const_iterator it = line.begin() ; it != line.end() ; ++it) {
            nvis::vec3 x = *it;
            box.add(x);
            double scl;
            nvis::vec3 f;
            rhs(x, f);
            scl = nvis::norm(f);
            curve.push_back(curve_point(x, scl));
            box.add(x);
        }
        all_curves.push_back(curve);
        all_bounds.push_back(box);
    }

    for (int i=0 ; i<nf ; ++i) {
        double dmin = thresholds[i];
        double dmax = thresholds[i+1];

        vtkSmartPointer<vtkPolyData> tubes =
            to_tubes(all_curves, radius);

        vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mapper->SetInputData(tubes);
        mapper->SetLookupTable(ctf);
        mapper->ScalarVisibilityOn();
        vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
        actor->SetMapper(mapper);
        ren->AddActor(actor);

        // included grains
        std::set<int> selected_grains;
        for (std::map< int, float >::const_iterator it = grain_vals.begin() ;
                it != grain_vals.end() ; ++it) {
            if (it->first < 0 || it->second < dmin || it->second >= dmax) {
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
            vtkIdType _n;
            const vtkIdType* id;
            if (!lines->GetNextCell(_n, id)) {
                break;
            }
            if (included(selected_vertices, id, 2)) {
                selected_lines->InsertNextCell(_n, id);
            }
        }
        std::cerr << selected_lines->GetNumberOfCells() << " edges passed the test\n";

        // assign colors to triangles based on the value of their grain(s)
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
                float v = grain_norms[ref_id];
                uchar_color_type col = to_color(cmap(v, color_map_type::REDUNDANT_RAINBOW));
                color->InsertNextTypedTuple(col.begin());
            }
        }


        vtkSmartPointer<vtkPolyData> selected_mesh = vtkSmartPointer<vtkPolyData>::New();
        selected_mesh->SetPoints(pts);
        selected_mesh->SetPolys(selected_triangles);
        selected_mesh->GetCellData()->SetScalars(color);

        //
        vtkSmartPointer<vtkPolyDataMapper> mesh_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
        mesh_mapper->SetInputData(selected_mesh);
        // mesh_mapper->ScalarVisibilityOn();
        vtkSmartPointer<vtkActor> mesh_actor = vtkSmartPointer<vtkActor>::New();
        mesh_actor->SetMapper(mesh_mapper);
        // mesh_actor->GetProperty()->SetColor(1, 1, 1);
        // mesh_actor->GetProperty()->SetOpacity(0.75);
        mesh_actor->GetProperty()->SetOpacity(0.75);

        vtkSmartPointer<vtkPolyData> selected_edges = vtkSmartPointer<vtkPolyData>::New();
        selected_edges->SetPoints(edges->GetPoints());
        selected_edges->SetLines(selected_lines);
        vtkSmartPointer<vtkTubeFilter> etubes = vtkSmartPointer<vtkTubeFilter>::New();
        etubes->SetInputData(selected_edges);

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
        os << movie_base << "_frame_" << std::setw(3) << std::setfill('0') << i+1 << "_of_"
           << nf << ".tiff";

        writer->SetFileName(os.str().c_str());
        std::cerr << "about to write to file " << os.str() << '\n';
        writer->Write();

        ren->RemoveActor(mesh_actor);
        ren->RemoveActor(edge_actor);
        ren->RemoveActor(actor);
    }

    if (nin_vec) {
        nrrdNuke(nin_vec);
    }

    return 0;
}
