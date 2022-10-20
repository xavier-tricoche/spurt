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
#include "vtkTubeFilter.h"
#include "vtkCellData.h"
#include "vtkSmartPointer.h"
#include "vtkPolyDataWriter.h"
#include "vtkScalarBarActor.h"
#include "vtkColorTransferFunction.h"
#include "vtkTextProperty.h"
#include "vtkGeometryFilter.h"
#include "vtkHedgeHog.h"
#include "vtkGlyph3D.h"
#include "vtkSphereSource.h"
#include "vtkPolyDataNormals.h"
#include "vtkInteractorStyleUser.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkInteractorStyleJoystickCamera.h"

#include <string>
#include <math/fixed_vector.hpp>
#include <vtk/vtk_utils.hpp>
#include <math/bounding_box.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>
#include <vis/integral_curve.hpp>
#include <boost/shared_ptr.hpp>
#include <exception>
#include "formulae.hpp"

inline void wait(int s)
{
    nvis::timer t;
    while (t.elapsed() < s) {}
}

double d=0.02, f=7.5, rho=1.2, e=0.92, K1=2.8e+5, a_over_d=1, eps=1.0e-7, T=0.5;
int nb_particles=20;
nvis::bbox2 bounds;
float seed_radius = 0.1;

void printUsageAndExit( const std::string& argv0, std::string what = "", bool doExit = true )
{
    if (what != "") {
        std::cerr << what << '\n';
    }
    std::cerr
            << "Usage  : " << argv0 << " [options]\n"
            << "Options:\n"
            << "    -h   | --help                      Print this information\n"
            << "    -d   | --diameter <float>          Particle diameter (m)\n"
            << "    -f   | --frequency <float>         Tap frequency (Hz)\n"
            << "    -rho | --density <float>           Material density (kg/m^3)\n"
            << "    -e   | --restitution <float>       Coefficient of restitution\n"
            << "    -K   | --stiffness <float>         Loading stiffness (N/m^2)\n"
            << "    -a   | --amplitude <float>         Normalized tap amplitude\n"
            << "    -np  | --particles <int>           Number of particles in column\n"
            << "    -n   | --seeds <int>               Number of integration seeds\n"
            << "    -p   | --precision <float>         Integration precision\n"
            << "    -T   | --period <float>            Time interval between taps (s)\n"
            << "    -b   | --bounds 4 x <float>        Velocity / Acceleration bounds\n"
            << "    -sr  | --seedradius <float>        Seed radius\n"
            << std::endl;
            
    if (doExit) {
        exit(1);
    }
}

std::vector<nvis::vec2> large_steps;
nvis::bbox2 large_bounds;
struct ODE_wrapper {
    ODE_wrapper(const njit::tapping_parameters& param) : _param(param) {}
    
    bool operator()(double t, const nvis::vec2& x, nvis::vec2& f) const {
        f = njit::ode_BRTU_31(x, t, _param);
        if (nvis::norm(f) > 50) {
            large_steps.push_back(x);
            large_bounds.add(x);
        }
        return true;
    }
    
    const njit::tapping_parameters& _param;
};

struct Euler_solver {
    Euler_solver() : _seed(0), abstol(1.0e-6), reltol(1.0e-6) {}
    Euler_solver(const nvis::vec2& x) : _seed(x), abstol(1.0e-6), reltol(1.0e-6) {}
    
    template<typename RHS>
    void advance(const RHS& rhs, double t) {
        double dt;
        if (t > 0) {
            dt = abstol;
        } else {
            dt = -abstol;
        }
        
        _trajectory[0.] = _seed;
        nvis::vec2 x = _seed;
        for (double u=0 ; u<t ; u+=dt) {
            nvis::vec2 f;
            rhs(u, x, f);
            x += f;
            _trajectory[u+dt] = x;
        }
    }
    
    double t_max() const {
        return _trajectory.rbegin()->first;
    }
    
    double t_min() const {
        return _trajectory.begin()->first;
    }
    
    nvis::vec2 operator()(double t) const {
        if (t < t_min() || t > t_max()) {
            throw std::runtime_error("invalid time coordinate");
        }
        std::map<double, nvis::vec2>::const_iterator prev, next;
        next = _trajectory.lower_bound(t);
        prev = next;
        --prev;
        double tlo = prev->first;
        double thi = next->first;
        double u = (t-tlo)/(thi-tlo);
        return (1.-u)*prev->second + u*next->second;
    }
    
    double abstol, reltol;
    
    nvis::vec2 _seed;
    double _t;
    std::map<double, nvis::vec2> _trajectory;
};


typedef ODE_wrapper             rhs_type;
typedef nvis::integral_curve<2> sl_type;
// typedef Euler_solver sl_type;

template<typename Iterator>
inline vtkSmartPointer<vtkPolyData> to_polydata(const Iterator& _begin, const Iterator& _end)
{
    vtkSmartPointer<vtkDoubleArray> coords = vtkSmartPointer<vtkDoubleArray>::New();
    coords->SetNumberOfComponents(3);
    int _size = 0;
    for (Iterator it = _begin ; it != _end ; ++it, ++_size) {
        double x[] = {(*it)[0], (*it)[1], 0.};
        coords->InsertNextTuple(x);
    }
    vtkSmartPointer<vtkPoints> pts = vtkSmartPointer<vtkPoints>::New();
    pts->SetData(coords);
    vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType n = 2;
    for (int i = 0 ; i < _size-1 ; ++i) {
        vtkIdType ids[] = { i, i + 1 };
        lines->InsertNextCell(n, ids);
    }
    
    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(pts);
    pd->SetLines(lines);
    return pd;
}

void draw_balls(const std::vector<nvis::vec2>& pts, const nvis::fvec3& col, float radius,
                vtkRenderer* ren)
{

    vtkSmartPointer<vtkPoints> _pts = vtkSmartPointer<vtkPoints>::New();
    for (int i=0 ; i<pts.size() ; ++i) {
        _pts->InsertPoint(i, pts[i][0], pts[i][1], 0);
    }
    vtkSmartPointer<vtkPolyData> pd = vtkSmartPointer<vtkPolyData>::New();
    pd->SetPoints(_pts);
    vtkSmartPointer<vtkSphereSource> balls = vtkSmartPointer<vtkSphereSource>::New();
    balls->SetRadius(radius);
    balls->SetPhiResolution(20);
    balls->SetThetaResolution(20);
    vtkSmartPointer<vtkGlyph3D> glyphPoints = vtkSmartPointer<vtkGlyph3D>::New();
#if VTK_MAJOR_VERSION >= 6
    glyphPoints->SetInputData(pd);
#else
    glyphPoints->SetInput(pd);
#endif
    glyphPoints->SetSourceConnection(balls->GetOutputPort());
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(glyphPoints->GetOutputPort());
    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(col[0], col[1], col[2]);
    ren->AddActor(actor);
}

inline nvis::vec3 tovec(double* a)
{
    nvis::vec3 v(a[0], a[1], a[2]);
    return v;
}

class cameraCB : public vtkCommand {
public:
    static cameraCB* New() {
        return new cameraCB;
    }
    virtual void Execute(vtkObject* caller, unsigned long, void*) {
        vtkRenderer* ren = reinterpret_cast<vtkRenderer*>(caller);
        double foc[3], pos[3], up[3];
        double clip[2];
        ren->GetActiveCamera()->GetPosition(pos);
        ren->GetActiveCamera()->GetFocalPoint(foc);
        ren->GetActiveCamera()->GetViewUp(up);
        ren->GetActiveCamera()->GetClippingRange(clip);
        
        std::cout << "camera position:       " << tovec(pos) << '\n';
        std::cout << "camera focal point:    " << tovec(foc) << '\n';
        std::cout << "camera up vector:      " << tovec(up) << '\n';
        std::cout << "camera clipping range: " << "(" << clip[0] << ", " << clip[1] << ")\n";
    }
};

class passiveCB : public vtkCommand {
public:
    static passiveCB* New() {
        return new passiveCB;
    }
    
    virtual void Execute(vtkObject*, unsigned long, void*) {
    }
};

struct mystop {
    mystop() : max_dist(1) {}
    
    bool operator()(const sl_type::int_step& step) {
        double d = nvis::norm(step.y0() - step.y1());
        if (d > max_dist) {
            std::cerr << "y(" << step.t0() << ")=" << step.y0() << ", y(" << step.t1() << ")=" << step.y1() << '\n';
            max_dist = d;
        }
        return false;
    }
    
    double max_dist;
};

int main(int argc, char* argv[])
{
    bool bounds_set = false;
    int nseeds = 10;
    for (int i=1; i<argc ; ++i) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            printUsageAndExit(argv[0]);
        } else if (arg == "-np" || arg == "--particles") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            nb_particles = atoi(argv[++i]);
        } else if (arg == "-n" || arg == "--seeds") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            nseeds = atoi(argv[++i]);
        } else if (arg == "-d" || arg == "--diameter") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            d = atof(argv[++i]);
        } else if (arg == "-f" || arg == "--frequency") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            f = atof(argv[++i]);
        } else if (arg == "-rho" || arg == "--density") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            rho = atof(argv[++i]);
        } else if (arg == "-K" || arg == "--stiffness") {
            if (i <= argc-2) {
                printUsageAndExit(argv[0]);
            }
            K1 = atof(argv[++i]);
        } else if (arg == "-e" || arg == "--restitution") {
            if (i <= argc-2) {
                printUsageAndExit(argv[0]);
            }
            e = atof(argv[++i]);
        } else if (arg == "-a" || arg == "--amplitude") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            a_over_d = atof(argv[++i]);
        } else if (arg == "-p" || arg == "--precision") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            eps = atof(argv[++i]);
        } else if (arg == "-T" || arg == "--period") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            T = atof(argv[++i]);
        } else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-4) {
                printUsageAndExit(argv[0]);
            }
            float minx = atof(argv[++i]);
            float miny = atof(argv[++i]);
            float maxx = atof(argv[++i]);
            float maxy = atof(argv[++i]);
            bounds.min() = nvis::vec2(minx, miny);
            bounds.max() = nvis::vec2(maxx, maxy);
            bounds_set = true;
        } else if (arg == "-sr" || arg == "--seedradius") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            seed_radius = atof(argv[++i]);
        }
    }
    if (!bounds_set) {
        printUsageAndExit(argv[0], "Missing bounds");
    } else {
        std::cerr << "bounds are " << bounds << '\n';
    }
    
    vtkRenderWindow* renWin = vtkRenderWindow::New();
    renWin->PointSmoothingOn();
    renWin->LineSmoothingOn();
    renWin->PolygonSmoothingOn();
    vtkRenderer* ren = vtkRenderer::New();
    renWin->AddRenderer(ren);
    renWin->SetSize(1200, 800);
    
    njit::tapping_parameters param;
    param.N = nb_particles;
    param.e = e;
    param.K1 = K1;
    param.rho = rho;
    param.f = f;
    param.d = d;
    param.a = d*a_over_d;
    param.T = T;
    
    std::cerr << "total mass is " << param.total_mass() << '\n';
    
    srand48(time(0));
    // double intg_time = 0;
    
    typedef std::list<nvis::vec2>   curve_type;
    typedef std::list<nvis::vec2>   plot_type;
    
    std::list<vtkSmartPointer<vtkPolyData> > all_curves;
    
    std::vector<nvis::vec3> pts, vecs;
    
    // draw the seeding window
    vtkSmartPointer<vtkPoints> fpts = vtkSmartPointer<vtkPoints>::New();
    fpts->InsertPoint(0, bounds.min()[0], bounds.min()[1], 0);
    fpts->InsertPoint(1, bounds.max()[0], bounds.min()[1], 0);
    fpts->InsertPoint(2, bounds.max()[0], bounds.max()[1], 0);
    fpts->InsertPoint(3, bounds.min()[0], bounds.max()[1], 0);
    vtkSmartPointer<vtkCellArray> flines = vtkSmartPointer<vtkCellArray>::New();
    vtkIdType ids[] = { 0, 1, 2, 3, 0 };
    flines->InsertNextCell(2, ids  );
    flines->InsertNextCell(2, ids+1);
    flines->InsertNextCell(2, ids+2);
    flines->InsertNextCell(2, ids+3);
    vtkSmartPointer<vtkPolyData> fpd = vtkSmartPointer<vtkPolyData>::New();
    fpd->SetPoints(fpts);
    fpd->SetLines(flines);
    vtkSmartPointer<vtkTubeFilter> ftube = vtkSmartPointer<vtkTubeFilter>::New();
#if VTK_MAJOR_VERSION >= 6
    ftube->SetInputData(fpd);
#else
    ftube->SetInput(fpd);
#endif
    ftube->SetRadius(0.1);
    ftube->SetCapping(1);
    ftube->SetNumberOfSides(12);
    vtkSmartPointer<vtkPolyDataMapper> fmapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    fmapper->SetInputConnection(ftube->GetOutputPort());
    vtkSmartPointer<vtkActor> factor = vtkSmartPointer<vtkActor>::New();
    factor->SetMapper(fmapper);
    factor->GetProperty()->SetColor(1,0,0);
    // ren->AddActor(factor);
    
    // compute integral curves started within seeding window
    nvis::bbox2 actual_bounds;
    std::vector<nvis::vec2> seeds;
    mystop stop;
    stop.max_dist = 0.5;
    for (int i = 0 ; i < nseeds ; ++i) {
        std::ostringstream os;
        os << "seed_" << i << ".csv";
        std::fstream out(os.str().c_str(), std::ios::out);
        
        nvis::vec2 seed = bounds.min() + nvis::vec2(drand48(), drand48()) * bounds.size();
        std::cerr << "seeding streamline #" << i << "/" << nseeds << " at " << seed << '\n';
        sl_type sl(seed);
        sl.abstol = sl.reltol = eps;
        rhs_type rhs(param);
        sl.advance<rhs_type, mystop>(rhs, 10.*T, stop);
        curve_type curve;
        nvis::vec2 last = sl(0);
        curve.push_back(last);
        double dt = 0.01;
        double t = 0;
        int counter = 0;
        std::cerr  << "tmax = " << sl.t_max() << '\n';
        while (t<sl.t_max()) {
            nvis::vec2 cur = sl(t+dt);
            if (nvis::norm(cur-last) > 0.1 && dt>1.0e-6) {
                dt /= 2;
                continue;
            } else if (nvis::norm(cur-last) < 1.0e-3) {
                dt *= 2;
            }
            curve.push_back(cur);
            out << cur[0] << "," << cur[1] << "," << t+dt << '\n';
            
            std::cout << counter++ << ": " << t+dt << ": " << dt << " : " << curve.back() << '\n';
            t += dt;
            last = cur;
            actual_bounds.add(curve.back());
        }
        out.close();
        
        all_curves.push_back(to_polydata(curve.begin(), curve.end()));
        seeds.push_back(seed);
    }
    std::cerr << "bounds = " << actual_bounds << '\n';
    std::cerr << "large bounds = " << large_bounds << '\n';
    
    // visualize seed points
    draw_balls(seeds, nvis::fvec3(1,1,0), seed_radius, ren);
    
    // visualize problematic locations
    // draw_balls(large_steps, nvis::fvec3(1,0,0), seed_radius, ren);
    
    // visualize integral curves as randomly colored curves
    for (std::list<vtkSmartPointer<vtkPolyData> >::const_iterator it=all_curves.begin() ;
            it!=all_curves.end() ; ++it) {
        vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
#if VTK_MAJOR_VERSION >= 6
        mapper->SetInputData(*it);
#else
        mapper->SetInput(*it);
#endif
        vtkActor* actor = vtkActor::New();
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(drand48(), drand48(), drand48());
        ren->AddActor(actor);
    }
    
    // finalize visual setting
    ren->SetBackground(0.0, 0.0, 0.0);
    ren->ResetCamera();
    // display camera setting to standard output
    // cameraCB *cb = cameraCB::New();
    // ren->AddObserver(vtkCommand::StartEvent, cb);
    // cb->Delete();
    // initialize camera to 2D friendly setting
    nvis::vec2 center = actual_bounds.center();
    ren->GetActiveCamera()->SetPosition(center[0], center[1], 10);
    ren->GetActiveCamera()->SetFocalPoint(center[0], center[1], 0);
    ren->GetActiveCamera()->SetViewUp(0,1,0);
    // start interactor
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin);
    renWin->Render();
    
    // select trackball interactor style
    // vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
    //  vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();
    // select joystick interactor style
    vtkSmartPointer<vtkInteractorStyleJoystickCamera> style =
        vtkSmartPointer<vtkInteractorStyleJoystickCamera>::New();
    // passiveCB* pcb1 = passiveCB::New();
    // style->AddObserver(vtkCommand::MouseWheelForwardEvent, pcb1);
    // passiveCB* pcb2 = passiveCB::New();
    // style->AddObserver(vtkCommand::MouseWheelBackwardEvent, pcb2);
    iren->SetInteractorStyle(style);
    iren->Initialize();
    
    iren->Start();
    
    return 0;
}
