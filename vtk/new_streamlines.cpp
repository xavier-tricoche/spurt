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

#include <string>
#include <math/fixed_vector.hpp>
#include <VTK/vtk_utils.hpp>
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

char*    file;
char*    scalar;
int     nblines;
float   length;
float   center[3];
float   radius[3];
float   eps;
float   h0;
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
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1, 1, &file,                NULL,       "input vector field (NRRD or VTK)");
    hestOptAdd(&hopt, "s",      "scalar",           airTypeString,  0, 1, &scalar,              "none",     "scalar field used for color coding (NRRD)");
    hestOptAdd(&hopt, "n",      "nb lines",         airTypeInt,     0, 1, &nblines,             "100",      "number of lines");
    hestOptAdd(&hopt, "l",      "length",           airTypeFloat,   1, 1, &length,              NULL,       "integration length");
    hestOptAdd(&hopt, "c",      "center",           airTypeFloat,   3, 3, &center,              NULL,       "seed center");
    hestOptAdd(&hopt, "r",      "radius",           airTypeFloat,   3, 3, &radius,              NULL,       "seed radius");
    hestOptAdd(&hopt, "e",      "eps",              airTypeFloat,   0, 1, &eps,                 "1.0e-6",   "integration precision");
    hestOptAdd(&hopt, "h",      "hinit",            airTypeFloat,   0, 1, &h0,                  "1.",       "initial step size guess");
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
    tubes->SetNumberOfSides(6);
    tubes->SetVaryRadiusToVaryRadiusByScalar();
    tubes->SetRadiusFactor(30);
    
    return tubes->GetOutput();
}

typedef nvis::vec3                              value_type;
typedef xavier::nrrd_data_traits<Nrrd*>              nrrd_data_traits;

struct my_rhs {
    my_rhs(const nrrd_data_traits& field) : _field(field) {}
    
    bool operator()(double, const nvis::vec3& x, nvis::vec3& f) const {
        bool r = _field.get_value(x, f);
        f *= 1.0e+6;
        return r;
    }
    
    const nrrd_data_traits& _field;
};

template<typename RHS>
struct euler {

    typedef RHS     rhs_type;
    
    enum state {
        OK = 0,
        LEFT_DOMAIN,
    };
    
    euler(const rhs_type& rhs, double h) : _rhs(rhs), _h(h) {}
    
    state integrate(const nvis::vec3& in, std::list<nvis::vec3>& out,
                    double length, double& actual_length) const {
        double h = (length < 0 ? -_h : _h);
        actual_length = 0;
        out.clear();
        out.push_back(in);
        nvis::vec3 x = in;
        nvis::vec3 f;
        int n = floor(fabs(length) / _h);
        for (int i = 0 ; i < n ; ++i) {
            if (_rhs(x, f)) {
                x += h * f;
                std::cerr << "at " << x << '\n';
                actual_length += fabs(h);
                out.push_back(x);
            } else {
                return LEFT_DOMAIN;
            }
        }
        double dt = length - n * h;
        if (_rhs(x, f)) {
            x += dt * f;
            out.push_back(x);
            actual_length += dt;
        } else {
            return LEFT_DOMAIN;
        }
        return OK;
    }
    
    const rhs_type& _rhs;
    double _h;
};

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    boost::shared_ptr<nrrd_data_traits> nrrd_vf, nrrd_sf;
    
    std::string name(file);
    Nrrd* nin_vec = 0, *nin_scl = 0;
    if (name.substr(name.find_last_of(".") + 1) == "nrrd") {
        std::cerr << "NRRD file extension recognized.\n";
        nin_vec = xavier::nrrd_utils::readNrrd(file);
        if (nin_vec->dim != 4) {
            throw;
        }
        
        nrrd_vf = boost::shared_ptr<nrrd_data_traits>(new nrrd_data_traits(nin_vec));
    } else {
        std::cerr << "file extension in " << name << " was not recognized.\n";
        return 1;
    }
    i2x to_pos(nin_vec);
    
    bool use_norm = !strcmp(scalar, "none");
    if (!use_norm) {
        nin_scl = xavier::nrrd_utils::readNrrd(scalar);
        if (nin_scl->dim != 3) {
            throw;
        }
        nrrd_sf = boost::shared_ptr<nrrd_data_traits>(new nrrd_data_traits(nin_scl));
    }
    
    std::vector<double> scalars;
    if (use_norm) {
        std::vector<double> tmp_scl;
        xavier::to_vector(tmp_scl, nin_vec);
        scalars.resize(tmp_scl.size() / 3);
        for (int i = 0 ; i < tmp_scl.size() / 3 ; ++i) {
            nvis::vec3 v(tmp_scl[3*i], tmp_scl[3*i+1], tmp_scl[3*i+2]);
            scalars[i] = nvis::norm(v);
        }
    } else {
        xavier::to_vector(scalars, nin_scl);
    }
    
    xavier::inverse_transform_sampling<double> itf(scalars);
    
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
            color_type c = cmap(v, (use_norm ? color_map::REDUNDANT_RAINBOW : color_map::DOUBLE_ENDED));
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
    
    std::cerr << "Camera set\n";
    
    vtkSmartPointer<vtkRenderWindow> ren_win = vtkSmartPointer<vtkRenderWindow>::New();
    ren_win->SetAlphaBitPlanes(1);
    ren_win->AddRenderer(ren);
    ren_win->SetSize(1600, 1200);
    
    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(ren_win);
    
    nvis::bbox3 bounds = xavier::nrrd_utils::get_bounds<3>(nin_vec);
    vtkSmartPointer<vtkActor> frame_actor = Garcia_vis_helper::draw_frame(bounds);
    ren->AddActor(frame_actor);
    
    nvis::vec3 min(center[0] - radius[0],
                   center[1] - radius[1],
                   center[2] - radius[2]);
    nvis::vec3 r(2.*radius[0], 2.*radius[1], 2.*radius[2]);
    
    std::list<nvis::vec3>   all_verts;
    std::list<double>       all_vals;
    std::list<int>          start_at;
    
    srand48(time(0));
    
    my_rhs rhs(*nrrd_vf);
#if 0
    euler<my_rhs> intg(rhs, h0);
    double intg_time = 0;
    for (int i = 0 ; i < nblines ; ++i) {
        std::list<double>       vals;
        
        nvis::vec3 seed = min + nvis::vec3(drand48(), drand48(), drand48()) * r;
        // int seed_id = itf.sample();
        // nvis::vec3 seed = to_pos(seed_id);
        
        std::cerr << "seeding streamline #" << i << "/" << nblines << " at " << seed << '\n';
        
        nvis::timer one_timer;
        std::list<nvis::vec3> forward, backward;
        double tmp;
        int fwd = intg.integrate(seed, forward, length, tmp);
        int bwd = intg.integrate(seed, backward, -length, tmp);
        intg_time += one_timer.elapsed();
        
        backward.pop_front();
        forward.insert(forward.begin(), backward.rbegin(), backward.rend());
        for (std::list<nvis::vec3>::const_iterator it = forward.begin() ; it != forward.end() ; ++it) {
            double scl;
            if (!use_norm) {
                nrrd_sf->get_value(*it, scl);
            } else {
                nvis::vec3 f;
                nrrd_vf->get_value(*it, f);
                scl = nvis::norm(f);
            }
            vals.push_back(scl);
        }
        add_line(all_verts, all_vals, start_at, forward, vals);
    }
#else
    check stop(1000, bounds);
    double intg_time = 0;
    for (int i = 0 ; i < nblines ; ++i) {
        std::list<nvis::vec3>   line;
        std::list<double>       vals;
    
        nvis::vec3 seed = min + nvis::vec3(drand48(), drand48(), drand48()) * r;
        // int seed_id = itf.sample();
        // nvis::vec3 seed = to_pos(seed_id);
    
        std::cerr << "seeding streamline #" << i << "/" << nblines << " at " << seed << '\n';
    
        nvis::streamline sl(seed);
        sl.record = true;
        sl.reltol = sl.abstol = eps;
        sl.stepsz = h0;
    
        nvis::timer one_timer;
        stop.reset();
        int fwd = sl.advance(rhs, length, stop);
        stop.reset();
        int bwd = sl.advance(rhs, -length, stop);
        intg_time += one_timer.elapsed();
    
        double dt = (sl.t_max() - sl.t_min()) / (float)discretization;
        if (dt == 0) {
            continue;
        }
        for (double t = sl.t_min() ; t <= sl.t_max() ; t += dt) {
            nvis::vec3 x = sl(t);
            // std::cerr << "x(" << t << ") = " << x << '\n';
            double scl;
            if (!use_norm) {
                nrrd_sf->get_value(x, scl);
            } else {
                nvis::vec3 f;
                nrrd_vf->get_value(x, f);
                scl = nvis::norm(f);
            }
            line.push_back(x);
            vals.push_back(scl);
        }
        add_line(all_verts, all_vals, start_at, line, vals);
    }
#endif
    vtkSmartPointer<vtkPolyData> tubes =
        to_tubes(all_verts, all_vals, start_at, 0.0001 * nvis::norm(nrrd_vf->bounds().size()));
        
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInput(tubes);
    mapper->SetLookupTable(ctf);
    mapper->ScalarVisibilityOn();
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    ren->AddActor(actor);
    
    ren->AddActor2D(color_bar);
    
    std::cerr << "done.\n";
    std::cerr << "integration took " << intg_time << " seconds, i.e. "
              << intg_time / (double)nblines << " seconds per streamline.\n";
              
    ren->ResetCamera();
    
    ren_win->Render();
    iren->Initialize();
    iren->Start();
    
    if (nin_vec) {
        nrrdNuke(nin_vec);
    }
    if (nin_scl) {
        nrrdNuke(nin_scl);
    }
    
    return 0;
}
























































































