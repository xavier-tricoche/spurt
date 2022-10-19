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
#include "vtkInteractorStyleImage.h"
#include "vtkMaskPoints.h"
#include "vtkSliderWidget.h"
#include "vtkSliderRepresentation.h"
#include "vtkSliderRepresentation2D.h"
#include "vtkCommand.h"
#include "vtkProperty2D.h"

#include <string>
#include <sstream>
#include <math/fixed_vector.hpp>
#include <VTK/vtk_utils.hpp>
#include <math/bounding_box.hpp>
#include <set>
#include <sstream>
#include <map>
#include <list>
#include <util/timer.hpp>
#include <vis/integral_curve.hpp>
#include <boost/shared_ptr.hpp>
#include <exception>

const double g = 9.80665;

double __w, __T, __g, __r, __gstar, __a;
int __N, nseeds;
nvis::bbox2 bounds;

inline double modulo_omegaT(double t)
{
    static const double omegaT = __w*__T;
    return fmod(t, omegaT);
}

inline double W(double s)
{
    if (s <= M_PI) {
        return cos(s);
    } else {
        return 0.;
    }
}

inline nvis::vec2 Phi(const nvis::vec2& x)
{
    const double& theta = x[0];
    const double& v     = x[1];
    return nvis::vec2(modulo_omegaT(theta), __r*v + __g*W(theta+v));
}

vtkSmartPointer<vtkRenderer>     renderer;
vtkSmartPointer<vtkRenderWindow> window;

std::list<vtkSmartPointer<vtkPolyData> > all_orbits;
std::vector<nvis::vec2> seeds;

std::list<vtkSmartPointer<vtkActor> >          my_actors;
std::list<vtkSmartPointer<vtkPolyDataMapper> > my_mappers;
void redraw()
{
    // erase the contents of the renderer
    typedef std::list<vtkSmartPointer<vtkActor> > actor_list;
    for (actor_list::iterator it=my_actors.begin() ; it!=my_actors.end() ; ++it) {
        renderer->RemoveActor(*it);
    }
    my_actors.clear();
    my_mappers.clear();
    
    std::cerr << "redraw: " << all_orbits.size() << " orbits to draw\n";
    
    typedef std::list<vtkSmartPointer<vtkPolyData> > poly_list;
    typedef poly_list::const_iterator iterator_type;
    for (iterator_type it=all_orbits.begin() ; it!=all_orbits.end() ; ++it) {
        // use vtkMaskPoints to turn points into vertices
        vtkSmartPointer<vtkMaskPoints> mask = vtkSmartPointer<vtkMaskPoints>::New();
#if VTK_MAJOR_REVISION == 6
        mask->SetInputData(*it);
#else
        mask->SetInput(*it);
#endif
        mask->SetOnRatio(1);
        mask->GenerateVerticesOn();
        mask->Update();
        
        std::cerr << "current orbit contains " << mask->GetOutput()->GetNumberOfPoints()
                  << " points and " << mask->GetOutput()->GetNumberOfCells() << " cells\n";
                  
                  
        my_mappers.push_back(vtkSmartPointer<vtkPolyDataMapper>::New());
        my_mappers.back()->SetInputConnection(mask->GetOutputPort());
        
        my_actors.push_back(vtkSmartPointer<vtkActor>::New());
        my_actors.back()->SetMapper(my_mappers.back());
        my_actors.back()->GetProperty()->SetColor(drand48(), drand48(), drand48());
        
        double* _dx, *_dy, *_dz, *_c;
        _dx = my_actors.back()->GetXRange();
        _dy = my_actors.back()->GetYRange();
        _dz = my_actors.back()->GetZRange();
        _c  = my_actors.back()->GetCenter();
        nvis::vec2 dx(_dx[0], _dx[1]);
        nvis::vec2 dy(_dy[0], _dy[1]);
        nvis::vec2 dz(_dz[0], _dz[1]);
        nvis::vec2  c( _c[0],  _c[1]);
        
        std::cerr << "created actor has xrange " << dx << ", yrange " << dy << ", zrange" << dz << ", center " << c << '\n';
        
        renderer->AddActor(my_actors.back());
    }
    
    // force rerender
    renderer->ResetCamera();
    window->Render();
    std::cerr << "rendering complete\n";
}

void recompute()
{
    // recompute model parameters
    __gstar = g/(double)__N;
    __w     = sqrt(__g*__gstar/(2.*__a*(1.+__r)));
    
    all_orbits.clear();
    for (int i=0 ; i<seeds.size() ; ++i) {
        std::list<nvis::vec2> curve;
        curve.push_back(seeds[i]);
        for (int n=0 ; n<200 ; ++n) {
            curve.push_back(Phi(curve.back()));
        }
        std::cerr << "curve #" << i << " contains " << curve.size() << " points\n";
        all_orbits.push_back(vtk_utils::create_points<double, 2>(curve));
    }
    std::cerr << "(re)computed " << seeds.size() << " orbits\n";
    redraw();
}

void reseed()
{
    std::cerr << "creating " << nseeds << " seeds\n";
    seeds.resize(nseeds);
    for (int i = 0 ; i < nseeds ; ++i) {
        seeds[i] = bounds.min() + nvis::vec2(drand48(), drand48()) * bounds.size();
    }
    recompute();
}


inline void wait(int s)
{
    nvis::timer t;
    while (t.elapsed() < s) {}
}

void printUsageAndExit( const std::string& argv0, std::string what = "", bool doExit = true )
{
    if (what != "") {
        std::cerr << what << '\n';
    }
    std::cerr
            << "Description: Visualization of discrete dynamical system described\n"
            << "             in Equation 34 of Chaos submission.\n"
            << "Usage  : " << argv0 << " [options]\n"
            << "Options:\n"
            << "    -h   | --help                      Print this information\n"
            << "    -r   | --rho <float>               Coefficient of restitution\n"
            << "    -g   | --gamma <float>             Gamma coefficient\n"
            << "    -N   | --particles <int>           Number of particles in column\n"
            << "    -n   | --seeds <int>               Number of integration seeds\n"
            << "    -T   | --period <float>            Time interval between taps (s)\n"
            << std::endl;
    std::cerr
            << "Note: Following recommendations are made in the paper\n"
            << "    Pi/omega << T\n"
            << "    rho = 0.8\n"
            << "    gamma initially \"small\"\n"
            << "    system is time periodic over [0, omega*T] with tap occurring in interval [0, Pi]\n"
            << "Default values:\n"
            << "    rho     = " << __r << '\n'
            << "    gamma   = " << __g << '\n'
            << "    T       = " << __T << '\n'
            << "    N       = " << __N << '\n'
            << "    nseeds  = " << nseeds << '\n'
            << std::endl;
            
    if (doExit) {
        exit(1);
    }
}

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

    vtkSmartPointer<vtkPolyData> points = vtk_utils::create_points<double, 2>(pts);
    vtkSmartPointer<vtkSphereSource> balls = vtkSmartPointer<vtkSphereSource>::New();
    balls->SetRadius(radius);
    balls->SetPhiResolution(12);
    balls->SetThetaResolution(12);
    vtkSmartPointer<vtkGlyph3D> glyphPoints = vtkSmartPointer<vtkGlyph3D>::New();
#if VTK_MAJOR_VERSION == 6
    glyphPoints->SetInputData(points);
#else
    glyphPoints->SetInput(points);
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


class Callback_Seeds : public vtkCommand {
    Callback_Seeds() {}
public:
    static Callback_Seeds* New() {
        return new Callback_Seeds();
    }
    
    virtual void Execute(vtkObject* caller, unsigned long, void*) {
        vtkSliderWidget* slider = reinterpret_cast<vtkSliderWidget*>(caller);
        nseeds = (int)static_cast<vtkSliderRepresentation*>(slider->GetRepresentation())->GetValue();
        reseed();
    }
};
class Callback_N : public vtkCommand {
    Callback_N() {}
public:
    static Callback_N* New() {
        return new Callback_N();
    }
    
    virtual void Execute(vtkObject* caller, unsigned long, void*) {
        vtkSliderWidget* slider = reinterpret_cast<vtkSliderWidget*>(caller);
        __N = (int)static_cast<vtkSliderRepresentation*>(slider->GetRepresentation())->GetValue();
        recompute();
    }
};
class Callback_Gamma : public vtkCommand {
    Callback_Gamma() {}
public:
    static Callback_Gamma* New() {
        return new Callback_Gamma();
    }
    
    virtual void Execute(vtkObject* caller, unsigned long, void*) {
        vtkSliderWidget* slider = reinterpret_cast<vtkSliderWidget*>(caller);
        __g = static_cast<vtkSliderRepresentation*>(slider->GetRepresentation())->GetValue();
        recompute();
    }
};
class Callback_Amplitude : public vtkCommand {
    Callback_Amplitude() {}
public:
    static Callback_Amplitude* New() {
        return new Callback_Amplitude();
    }
    
    virtual void Execute(vtkObject* caller, unsigned long, void*) {
        vtkSliderWidget* slider = reinterpret_cast<vtkSliderWidget*>(caller);
        __a = static_cast<vtkSliderRepresentation*>(slider->GetRepresentation())->GetValue();
        recompute();
    }
};

vtkSmartPointer<vtkSliderRepresentation2D>
make_rep(const std::string& text,
         double v, double x, double y, double w,
         double min, double max, int prec=2)
{
    std::ostringstream os;
    os << "%0." << prec << "f";
    vtkSmartPointer<vtkSliderRepresentation2D> rep = vtkSmartPointer<vtkSliderRepresentation2D>::New();
    rep->SetMinimumValue(min);
    rep->SetMaximumValue(max);
    rep->SetValue(v);
    rep->SetTitleText(text.c_str());
    rep->SetLabelFormat(os.str().c_str());
    rep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    rep->GetPoint1Coordinate()->SetValue(x, y);
    rep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    rep->GetPoint2Coordinate()->SetValue(x+w, y);
    rep->SetSliderLength(0.005);
    rep->SetSliderWidth(0.02);
    rep->SetEndCapLength(0.002);
    rep->SetEndCapWidth(0.03);
    rep->SetTubeWidth(0.002);
    rep->GetTubeProperty()->SetColor(0.9, 0.9, 0.9);
    rep->GetCapProperty()->SetColor(0.9, 0.9, 0.9);
    rep->GetSliderProperty()->SetColor(1, 1, 0);
    rep->GetTitleProperty()->SetColor(1, 1, 0);
    rep->GetTitleProperty()->ShadowOff();
    rep->GetTitleProperty()->SetFontFamilyToTimes();
    rep->SetTitleHeight(0.02);
    rep->GetTitleProperty()->BoldOff();
    rep->GetLabelProperty()->SetColor(1, 1, 1);
    rep->SetLabelHeight(0.02);
    rep->GetLabelProperty()->SetFontFamilyToTimes();
    rep->GetLabelProperty()->BoldOff();
    rep->GetLabelProperty()->ShadowOff();
    return rep;
}

vtkSmartPointer<vtkRenderWindowInteractor> interactor;
template<typename Callback>
vtkSmartPointer<vtkSliderWidget>
make_slider_widget(const std::string& text, double v, int id, int ntotal, double min, double max, int prec)
{
    double spc = 0.05;
    double len = (1.-(ntotal+1)*spc)/(double)ntotal;
    double x = spc + id*(spc + len);
    
    std::cerr << "slider anchored at " << x << ", length = " << len << '\n';
    
    vtkSmartPointer<vtkSliderRepresentation2D> rep = make_rep(text, v, x, 0.07, len, min, max, prec);
    vtkSmartPointer<vtkSliderWidget> widget = vtkSmartPointer<vtkSliderWidget>::New();
    widget->SetInteractor(interactor);
    widget->SetRepresentation(rep);
    widget->SetAnimationModeToJump();
    vtkSmartPointer<Callback> cb = vtkSmartPointer<Callback>::New();
    widget->AddObserver(vtkCommand::EndInteractionEvent, cb);
    return widget;
}

int main(int argc, char* argv[])
{
    bool bounds_set = false;
    
    nseeds  = 200;
    __N     = 20;
    __T     = 0.5;
    __a     = 0.1;
    __r     = 0.8;
    __g     = 0.1;
    
    for (int i=1; i<argc ; ++i) {
        std::string arg(argv[i]);
        if (arg == "-h" || arg == "--help") {
            printUsageAndExit(argv[0]);
        } else if (arg == "-np" || arg == "--particles") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            __N = atoi(argv[++i]);
        } else if (arg == "-n" || arg == "--seeds") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            nseeds = atoi(argv[++i]);
        } else if (arg == "-r" || arg == "--rho") {
            if (i <= argc-2) {
                printUsageAndExit(argv[0]);
            }
            __r = atof(argv[++i]);
        } else if (arg == "-g" || arg == "--gamma") {
            if (i <= argc-2) {
                printUsageAndExit(argv[0]);
            }
            __g = atof(argv[++i]);
        } else if (arg == "-T" || arg == "--period") {
            if (i == argc-1) {
                printUsageAndExit(argv[0]);
            }
            __T = atof(argv[++i]);
        }
    }
    
    renderer = vtkRenderer::New();
    renderer->SetBackground(0, 0, 0);
    
    window = vtkRenderWindow::New();
    window->PointSmoothingOn();
    window->LineSmoothingOn();
    window->PolygonSmoothingOn();
    window->AddRenderer(renderer);
    window->SetSize(1200, 800);
    
    srand48(time(0));
    nvis::bbox2 bounds;
    bounds.min() = nvis::vec2(-1,-1);
    bounds.max() = nvis::vec2(1,1);
    reseed();
    
    // set interactor first
    vtkSmartPointer<vtkInteractorStyleImage> style = vtkSmartPointer<vtkInteractorStyleImage>::New();
    interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(window);
    interactor->SetInteractorStyle(style);
    interactor->Initialize();
    
    // add widgets
    vtkSmartPointer<vtkSliderWidget> w1 = make_slider_widget<Callback_Seeds>("# seeds", nseeds, 0, 3, 100, 1000, 0);
    vtkSmartPointer<vtkSliderWidget> w2 = make_slider_widget<Callback_Gamma>("Gamma", __g, 1, 3, 0, 1, 5);
    vtkSmartPointer<vtkSliderWidget> w3 = make_slider_widget<Callback_Amplitude>("Amplitude", __a, 2, 3, 0, 1, 3);
    w1->EnabledOn();
    w2->EnabledOn();
    w3->EnabledOn();
    
    // start interactor
    window->Render();
    interactor->Start();
    
    return 0;
}
