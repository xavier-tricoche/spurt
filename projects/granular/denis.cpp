// stdlib
#include <exception>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <sstream>
#include <string>
// boost
#define BPO_WRAPPER_IS_BROKEN
#ifndef BPO_WRAPPER_IS_BROKEN
#   include <boost/shared_ptr.hpp>
#   include <boost/program_options.hpp>
#   include <boost/lexical_cast.hpp>
#endif
// nvis
#include <math/bounding_box.hpp>
#include <math/fixed_vector.hpp>
#include <util/timer.hpp>
#include <vis/integral_curve.hpp>
// VTK
#include "vtkActor.h"
#include "vtkAxisActor.h"
#include "vtkCamera.h"
#include "vtkCardinalSpline.h"
#include "vtkCellData.h"
#include "vtkColorTransferFunction.h"
#include "vtkCommand.h"
#include "vtkCylinderSource.h"
#include "vtkGeometryFilter.h"
#include "vtkGlyph3D.h"
#include "vtkHedgeHog.h"
#include "vtkInteractorObserver.h"
#include "vtkInteractorStyleImage.h"
#include "vtkLight.h"
#include "vtkMaskPoints.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataWriter.h"
#include "vtkProperty.h"
#include "vtkProperty2D.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkScalarBarActor.h"
#include "vtkSliderRepresentation.h"
#include "vtkSliderRepresentation2D.h"
#include "vtkSliderWidget.h"
#include "vtkSmartPointer.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkSphereSource.h"
#include "vtkTextProperty.h"
#include "vtkTIFFWriter.h"
#include "vtkTransform.h"
#include "vtkTubeFilter.h"
#include "vtkWindowToImageFilter.h"
// spurt
#include <vtk/vtk_utils.hpp>
#include <poincare/metric.hpp>
#ifndef BPO_WRAPPER_IS_BROKEN
#   include <misc/option_parse.hpp>
#endif
#include "denis.hpp"

// notational convenience
typedef VTK_SMART(PolyData)      smart_poly;
typedef VTK_SMART(Actor)         smart_actor;
typedef VTK_SMART(Renderer)      smart_renderer;
typedef VTK_SMART(RenderWindow)  smart_window;

typedef nvis::fvec3 Color;
typedef std::vector<nvis::vec2> orbit_type;

enum curve_type {
    CURVE_POINTS,
    CURVE_POLYLINES,
    CURVE_SPLINES
};

curve_type _which_rep;

using namespace spurt::denis;

// Equation parameters
double  _omega, _T, _gamma, _rho, _g_star, _a;
size_t  _N;

// Computation parameters
size_t _nseeds, _niter;
nvis::vec2 _vbounds, _ic;
double _height;
Color _bg;
bool _use_ic, _zoom;

// Visualization parameters
double _kappa;
nvis::ivec2 _res;
bool _verbose;

spurt::map_metric _metric;

inline nvis::vec2 normalize(const nvis::vec2& x) {
    nvis::vec2 y;
    y[0] = x[0] / (_omega*_T);
    y[1] = (x[1] - _vbounds[0])/(_vbounds[1] - _vbounds[0]);
    return y;
}

inline bool inside(const nvis::vec2& x) {
    static const nvis::bbox2 bounds(nvis::vec2(0, _vbounds[0]), 
                                    nvis::vec2(_omega*_T, _vbounds[1]));
    return bounds.inside(x);
}

double winding_number(const orbit_type& orb) {
    size_t n = orb.size();
    if (!n) return 0;
    
    nvis::vec2 x0(M_PI/2, 0);
    nvis::vec2 r = orb[0] - x0;
    double last_theta = atan2(r[1], r[0]);
    double sum_theta=0;
    
    for (int i=1 ; i<n ; ++i) {
        r = orb[i]-x0;
        double theta = atan2(r[1], r[0]);
        double dtheta = theta - last_theta;
        if (dtheta > 2*M_PI) dtheta -= 2*M_PI;
        else if (dtheta < -2*M_PI) dtheta += 2*M_PI;
        sum_theta += dtheta;
        last_theta = theta;
    }
    
    return sum_theta/(2*M_PI)/(n-1);
}

smart_poly points_to_curves(const std::vector<std::vector<nvis::vec2> >& polylines) {
    VTK_CREATE(vtkPolyData, polydata);
    
    switch (_which_rep) {
        case CURVE_POINTS:
        case CURVE_POLYLINES: {
            polydata = vtk_utils::make_polylines(polylines);
            break;
        }
        case CURVE_SPLINES: {
            const double max_dist = 0.01*(_vbounds[1]-_vbounds[0]);
            const size_t nbr_of_curves = polylines.size();
            std::vector<std::vector<nvis::vec2> > curves;
            for (size_t n=0 ; n<nbr_of_curves ; ++n) {
                const std::vector<nvis::vec2>& polyline = polylines[n];
                if (polyline.size() < 2) continue;
                VTK_CREATE(vtkCardinalSpline, sx);
                VTK_CREATE(vtkCardinalSpline, sy);
                for (size_t i=0 ; i<polyline.size() ; ++i) {
                    sx->AddPoint(i, polyline[i][0]);
                    sy->AddPoint(i, polyline[i][1]);
                }
                
                curves.push_back(std::vector<nvis::vec2>());
                const size_t nbr_of_ctrl_pts = polyline.size();
                std::vector<nvis::vec2>& curve = curves.back();
                curve.push_back(nvis::vec2(sx->Evaluate(0), sy->Evaluate(0)));
                // adaptive sampling of spline curve
                double dt = 0.005*(nbr_of_ctrl_pts-1); // initial step size
                double last_t=0, t=dt;
                while (t <= nbr_of_ctrl_pts-1) {
                    nvis::vec2 next(sx->Evaluate(t), sy->Evaluate(t));
                    double dist = nvis::norm(next-curve.back());
                    if (dist < max_dist || dt <= 1.0e-5) {
                        curve.push_back(next);
                        last_t = t;
                        if (dist < 0.1*max_dist) dt *= 2.;
                        t += dt;
                    }
                    else {
                        dt /= 5.;
                        t = last_t + dt;
                    }
                }
                
                polydata = vtk_utils::make_polylines(curves);
            }
            break;
        }
    }
    
    return polydata;
}

typedef std::vector<nvis::vec2>           segment;
typedef std::pair<nvis::vec2, nvis::vec2> point_pair;

inline void decompose(const nvis::vec2& from, const nvis::vec2& to, 
                      std::vector<point_pair>& parts) {
    nvis::vec2 current, next;
    std::vector<point_pair> clips;
    parts.clear();
    
    _metric.clip_segment(clips, from, to);
    
    for (int i=0 ; i<clips.size() ; ++i) {
        if (!_metric.bounds().inside(current)) break;
        parts.push_back(clips[i]);
    }
}

inline void half_orbit(const nvis::vec2& x0, size_t niter, 
                       const discrete_holmes& phi, bool fwd,
                       std::vector<segment>& orbit) 
{
    nvis::vec2 current, next;
    orbit.clear();
    
    current = _metric.modulo(x0);
    orbit.push_back(segment());
    orbit.back().push_back(current);
    for (int i=0 ; i<niter ; ++i) {
        next = phi(current, fwd, false);
        std::vector<point_pair> steps;
        decompose(current, next, steps);
        if (steps.empty()) {
            break;
        }
        else if (steps.size() == 1) {
            orbit.back().push_back(next);
            current = next;
        }
        else {
            // complete current segment
            orbit.back().push_back(steps[0].second);
            orbit.push_back(segment());
            // add intermediate segments
            for (int j=1 ; j<steps.size()-1 ; ++j) {
                orbit.back().push_back(steps[j].first);
                orbit.back().push_back(steps[j].second);
                orbit.push_back(segment());
            }
            // turn final portion into new segment
            orbit.back().push_back(steps.back().first);
            orbit.back().push_back(steps.back().second);
            current = orbit.back().back();
        }
    }
    if (orbit.back().size() < 2) orbit.pop_back();
}

inline void half_orbit2(const nvis::vec2& x0, size_t niter, 
                        const discrete_holmes& phi, bool fwd,
                        std::vector<segment>& orbit) 
{
    nvis::vec2 current, next;
    orbit.clear();
    
    current = x0;
    orbit.push_back(segment());
    orbit.back().push_back(current);
    for (int i=0 ; i<niter ; ++i) {
        next = phi(current, fwd, false);
        if (next[1] < _vbounds[0] || next[1] > _vbounds[1]) break;
        orbit.back().push_back(next);
        current = next;
    }
}

inline smart_poly
compute_one_orbit(const nvis::vec2& x0, size_t niter, const discrete_holmes& phi)
{
    nvis::vec2 current, next;
    std::vector<segment> bwd, fwd, both;
    
    half_orbit(x0, niter, phi, false, bwd);
    half_orbit(x0, niter, phi, true, fwd);
    
    // reverse bwd portion and preprend to fwd section
    both.push_back(segment());
    for (int i=bwd.size()-1 ; i>0 ; --i ) {
        const segment& s = bwd[i];
        for (int j=s.size()-1 ; j>=0 ; --j) {
            both.back().push_back(s[j]);
        }
        both.push_back(segment());
    }
    // combine fwd and bwd segments at seed
    // first point is already present in fwd branch, skip it
    for (int i=bwd[0].size()-1 ; i>0 ; --i) {
        both.back().push_back(bwd[0][i]);
    }
    std::copy(fwd[0].begin(), fwd[0].end(), std::back_inserter(both.back()));
    for (int i=1 ; i<fwd.size() ; ++i) {
        both.push_back(segment());
        std::copy(fwd[i].begin(), fwd[i].end(), std::back_inserter(both.back()));
    }
    if (both.back().size() < 2) both.pop_back();
    
    VTK_CREATE(vtkPolyData, polydata);
    polydata = points_to_curves(both);
    
    // std::cout << "winding number = " << winding_number(both.front()) << '\n';
    
    return polydata;
}

#if 0
inline smart_poly
compute_one_orbit(const nvis::vec2& x0, size_t niter, const discrete_holmes& phi)
{
    typedef std::vector<nvis::vec2>           segment;
    typedef std::pair<nvis::vec2, nvis::vec2> point_pair;
    std::vector<segment> points;
    std::vector<point_pair> clips;
    points.push_back(segment());
    points.back().push_back(_metric.modulo(x0));
    
    
    if (_verbose) std::cout << "\nforward iteration...\n";
    for (int i=0 ; i<niter ; ++i) {
        const nvis::vec2 current = points.back().back();
        std::pair<nvis::vec2,nvis::vec2> r = phi.step(current, true);
        nvis::vec2 next = r.second;
        nvis::vec2 vec = r.first;
        if (_verbose) {
            std::cout << "phi(" << current << ")=" << next << '\n';
            std::cout << "\t\tphi(" << normalize(current) << ")=" << normalize(next) << ")\n"; 
            std::cout << "\t\tstep=" << vec << '\n';
        }
        _metric.clip_segment(clips, current, next);
        if (clips.size() > 1) {
            points.back().push_back(clips[0].second);
            points.push_back(segment());
            points.back().push_back(clips[1].first);
            if (_verbose) {
                std::cout << "after clipping: " << current << "->" << clips[0].second << " + "
                          << clips[1].first << "->" << next << '\n';
                std::cout << "\t\t" << normalize(current) << "->" << normalize(clips[0].second) << " + "
                          << normalize(clips[1].first) << "->" << normalize(next) << '\n';
            }
        }
        if (!_metric.bounds().inside(next)) { points.pop_back(); break; }
        points.back().push_back(_metric.modulo(next));
    }
    if (_verbose) std::cout << "\nBackward iteration...\n";
    points.push_back(segment());
    points.back().push_back(_metric.modulo(x0));
    for (int i=0 ; i<niter ; ++i) {
        const nvis::vec2 current = points.back().back();
        nvis::vec2 next = phi(current, false, false);
        if (!inside(next)) {
            _metric.clip_segment(clips, current, next);
            points.back().push_back(clips[0].second);
            
            while (inside(clips[1].first)) {
                points.push_back(segment());
                points.back().push_back(clips[1].first);
                points.back().push_back(_metric.modulo(next))
        }
        
        
        std::pair<nvis::vec2,nvis::vec2> r = phi.step(current, false);
        nvis::vec2 next = r.second;
        nvis::vec2 vec = r.first;
        if (_verbose) {
            std::cout << "phi(" << current << ")=" << next << '\n';
            std::cout << "\t\tphi(" << normalize(current) << ")=" << normalize(next) << ")\n";
            std::cout << "\t\tstep=" << vec << '\n';
        }
        _metric.clip_segment(clips, current, next);
        if (clips.size() > 1) {
            points.back().push_back(clips[0].second);
            points.push_back(segment());
            points.back().push_back(clips[1].first);        
            if (_verbose) {
                std::cout << "after clipping: " << current << "->" << clips[0].second << " + "
                          << clips[1].first << "->" << next << '\n';
                std::cout << "\t\t" << normalize(current) << "->" << normalize(clips[0].second) << " + "
                          << normalize(clips[1].first) << "->" << normalize(next) << '\n';
            }
        }
        if (!_metric.bounds().inside(next)) { points.pop_back(); break; }
        points.back().push_back(_metric.modulo(next));
    }
    VTK_CREATE(vtkPolyData, polydata);
    polydata = points_to_curves(points);
    return polydata;
}
#endif

typedef std::pair<smart_poly, nvis::fvec3> color_poly;

std::list<color_poly> all_orbits;
std::list<nvis::vec2> all_seeds;
void compute_all_orbits(size_t niter, const discrete_holmes& phi)
{
    srand48(13081975);
    
    static const double _OMEGAT = _omega*_T;
    static const double _VSPAN  = _vbounds[1]-_vbounds[0];
    static const double WIDTH = _zoom ? M_PI : _OMEGAT;
    if (!_use_ic) {
        for (int i=0 ; i<_nseeds ; ++i) {
            nvis::vec2 seed;
            double x = drand48();
            double y = drand48();
            seed[0]  = x*WIDTH;
            seed[1]  = _vbounds[0] + y*_VSPAN;
            
            // color coding for the orbit
            nvis::fvec3 col(drand48(), drand48(), drand48());
            col /= nvis::max(col); // saturate one channel
            all_orbits.push_back(color_poly(compute_one_orbit(seed, niter, phi), col));
            all_seeds.push_back(seed);
        }
    }
    else {
        nvis::vec2 seed;
        seed[0] = _ic[0]*_OMEGAT;
        seed[1] = _vbounds[0] + _ic[1]*_VSPAN;
        nvis::fvec3 col(1,0,0);
        all_orbits.push_back(color_poly(compute_one_orbit(seed, niter, phi), col));
        all_seeds.push_back(seed);
    }
}

std::list<smart_actor> all_actors;
void update_actors()
{
    for (std::list<color_poly>::iterator it=all_orbits.begin() ; it!=all_orbits.end() ; ++it) {
        VTK_CREATE(vtkActor, one_actor);
        if (_which_rep == CURVE_POINTS) {
            VTK_CREATE(vtkMaskPoints, mask);
            VTK_CONNECT(mask, it->first);
            mask->RandomModeOff();
            mask->SetOnRatio(1);
            mask->GenerateVerticesOn();
            mask->Update();
            vtkPolyData* pd = mask->GetOutput();
            VTK_MAKE_ACTOR(actor_ptr, pd);
            one_actor = actor_ptr;
        }
        else {
            VTK_MAKE_ACTOR(actor_ptr, it->first);
            one_actor = actor_ptr;
        }
        one_actor->GetMapper()->ScalarVisibilityOff();
        const nvis::fvec3& col = it->second;
        one_actor->GetProperty()->SetColor(col[0], col[1], col[2]);
        one_actor->GetProperty()->SetRepresentationToWireframe();
        if (_which_rep == CURVE_POINTS)
            one_actor->GetProperty()->SetPointSize(2);
        all_actors.push_back(one_actor);
    }
}

smart_renderer renderer;
smart_window   window;
void draw()
{
    for (std::list<smart_actor>::iterator it=all_actors.begin() ; it!=all_actors.end() ; ++it) {
        renderer->AddActor(*it);
    }
    // window->GetRenderers()->GetFirstRenderer()->ResetCamera();
    // window->GetRenderers()->GetFirstRenderer()->GetActiveCamera()->Zoom(_zoom);
    
    // set up camera
    nvis::vec2 center;
    center[0] = 0.5*(_zoom ? M_PI : _omega*_T);
    center[1] = 0.5*(_vbounds[0] + _vbounds[1]);
    double scale = 0.5*(_vbounds[1]-_vbounds[0]);
    
    renderer->GetActiveCamera()->SetFocalPoint(center[0], center[1], 0);
    renderer->GetActiveCamera()->SetPosition(center[0], center[1], 10);
    renderer->GetActiveCamera()->SetViewUp(0, 1, 0);
    renderer->GetActiveCamera()->ParallelProjectionOn();
    renderer->GetActiveCamera()->SetParallelScale(scale);
    
    window->Render();
}

std::string me;
void usage(const std::string& what = "")
{
    if (what != "") {
        std::cerr << "Error: " << what << '\n';
    }
    std::cerr
        << "DESCRIPTION: Visualize discrete dynamical system describing\n"
        << "             the motion of the mass center of a column of\n"
        << "             particles undergoing tapping.\n"
        << '\n'
        << "USAGE:" << me << " [options]\n"
        << '\n'
        << "OPTIONS:\n"
        << " -h | --help               Print this information\n"
        << " -o | --output <string>    Save image to file\n"
        << " -a | --amplitude <float>  Normalized tap amplitude\n"
        << " -r | --rho <float>        Coefficient of restitution\n"
        << " -g | --gamma <float>      Gamma constant\n"
        << " -w | --omega <float>      Omega constant\n"
        << " -N | --particles <int>    Number of particles in column\n"
        << " -n | --seeds <int>        Number of integration seeds\n"
        << " -x | --ic <float> x2      Initial condition\n"
        << " -T | --period <float>     Overall tap period\n"
        << " -b | --bounds <float> x2  Velocity bounds\n"
        << " -s | --size <int> x2      Image size\n"
        << " -c | --color <float> x3   Background color\n"
        << " -v | --verbose            Turn on verbose mode\n"
        << " -l | --line <int>         Line type (0: points, 1: polyline, 2: spline)\n"
        << " -z | --zoom               Focus computation on tap region [0, \\pi]\n"
        << '\n';

    exit(!what.empty());
}

int main(int argc, char* argv[])
{    
    _vbounds[0] = -1;
    _vbounds[1] = 1.;
    _N          = 20;
    _nseeds     = 200;
    _niter      = 200;
    _T          = 10.;
    _rho        = 0.8;
    _a          = -1;    // invalid value
    _gamma      = -1;    // invalid value
    _omega      = -1;    // invalid value
    _verbose    = false;
    _res[0]     = 1400;
    _res[1]     = 850;
    _bg[0]      = 0;
    _bg[1]      = 0;
    _bg[2]      = 0;
    _ic[0]      = 0;
    _ic[1]      = 0;
    _use_ic     = false;
    _which_rep  = CURVE_SPLINES;
    _zoom       = false;
    std::string output;
  
    me = argv[0];
    for (int i=0 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") usage();
        else if (arg == "-a" || arg == "--amplitude") {
            if (i == argc-1) usage("Missing amplitude value");
            _a = atof(argv[++i]);
        }
        else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) usage("Missing output file name");
            output = argv[++i];
        }
        else if (arg == "-r" || arg == "--rho") {
            if (i == argc-1) usage("Missing rho value");
            _rho = atof(argv[++i]);
        }
        else if (arg == "-g" || arg == "--gamma") {
            if (i == argc-1) usage("Missing gamma value");
            _gamma = atof(argv[++i]);
        }
        else if (arg == "-w" || arg == "--omega") {
            if (i == argc-1) usage("Missing omega value");
            _omega = atof(argv[++i]);
        }
        else if (arg == "-N" || arg == "--particles") {
            if (i == argc-1) usage("Missing number of particles");
            _N = atoi(argv[++i]);
        }
        else if (arg == "-n" || arg == "--seeds") {
            if (i == argc-1) usage("Missing number of seeds");
            _nseeds = atoi(argv[++i]);
        }
        else if (arg == "-x" || arg == "--ic") {
            if (i >= argc-2) usage("Missing initial condition");
            _ic[0] = atof(argv[++i]);
            _ic[1] = atof(argv[++i]);
            _use_ic = true;
        }
        else if (arg == "-T" || arg == "--period") {
            if (i == argc-1) usage("Missing period value");
            _T = atof(argv[++i]);
        }
        else if (arg == "-b" || arg == "--bounds") {
            if (i >= argc-2) usage("Missing bounds");
            _vbounds[0] = atof(argv[++i]);
            _vbounds[1] = atof(argv[++i]);
        }
        else if (arg == "-s" || arg == "--size") {
            if (i >= argc-2) usage("Missing image size");
            std::string next;
            next = argv[++i];
            if (next == "-") _res[0] = -1;
            else _res[0] = atoi(argv[i]);
            next = argv[++i];
            if (next == "-") _res[0] = -1;
            else _res[1] = atoi(argv[i]);
            if (_res[0] < 0 && _res[1] < 0) {
                usage("Invalid image resolution");
            }
        }
        else if (arg == "-c" || arg == "--color") {
            if (i >= argc-3) usage("Missing background color");
            _bg[0] = atof(argv[++i]);
            _bg[1] = atof(argv[++i]);
            _bg[2] = atof(argv[++i]);
        }
        else if (arg == "-l" || arg == "--line") {
            if (i == argc-1) usage("Missing line type");
            _which_rep = curve_type(atoi(argv[++i]));
        }
        else if (arg == "-v" || arg == "--verbose") {
            _verbose = true;
        }
        else if (arg == "-z" || arg == "--zoom") {
            _zoom = true;
        }
    }
    
    _g_star = discrete_holmes::g/(double)_N;
    
    if (_a < 0 && _gamma > 0 && _omega > 0) {
        // compute amplitude as a function of gamma and omega
        // \gamma := 2*\omega^2*a*(1 + \rho)/g
        // a = g*\gamma/(2*\omega^2*(1 + \rho))
        _a = _gamma*_g_star/(2*_omega*_omega*(1+_rho));
    }
    else if (_gamma < 0 && _a > 0 && _omega > 0) {
        _gamma = 2*_omega*_omega*_a*(1+_rho)/_g_star;
    }
    else if (_omega < 0 && _a > 0 && _gamma > 0) {
    }
    else {
        if (_a < 0) _a = 0.1;
        if (_gamma < 0) _a = 0.1;
    }
    
    discrete_holmes phi(_T, _a, _rho, _gamma, _N);
    _omega = phi.omega();
    
    double M = 2*_a*_omega/(_g_star*_T)*(1+_rho)/(1-_rho);
    
    if (_verbose) {
        std::cout << "parameters:\n"
                  << "omega = " << _omega << ", pi/omega = " << M_PI/_omega << '\n'
                  << "gamma = " << _gamma << '\n'
                  << "amp   = " << _a << '\n'
                  << "T     = " << _T << '\n'
                  << "N     = " << _N << '\n'
                  << "M     = " << M << '\n';
                  
        if (M >= 1) {
            for (int m=1 ; m<M ; ++m) {
                double xm = acos(m*_omega*_T/_gamma*(1-_rho));
                double ym = m*_omega*_T;
                std::cout << "singularity #" << m << "/" << (int)floor(M) 
                          << ": " << nvis::vec2(xm, ym) << '\n';
            }
        }
    }
    
    // set metric for display
    _metric = phi.metric();
    _metric.bounds().min()[1] = _vbounds[0];
    _metric.bounds().max()[1] = _vbounds[1];
    if (_verbose) std::cout << "computation bounds are " << _metric.bounds() << '\n';
    
    srand48(time(0));
    
    renderer = smart_renderer::New();
    renderer->SetBackground(_bg[0], _bg[1], _bg[2]);
    
    // vtk_utils::camera_setting_callback *cb = vtk_utils::camera_setting_callback::New();
    // renderer->AddObserver(vtkCommand::StartEvent, cb);
    // cb->Delete();
    
    window = smart_window::New();
    window->PointSmoothingOn();
    window->LineSmoothingOn();
    window->PolygonSmoothingOn();
    window->AddRenderer(renderer);
    
    _height = _vbounds[1] - _vbounds[0];
    double zoomed_width = _zoom ? M_PI : _omega*_T;
    double ratio = _height/zoomed_width;
    double newheight = ratio*_res[0];
    double newwidth = _res[1]/ratio;
    if (_verbose) {
        std::cout << "Ratio = " << ratio << '\n'
        << "new height = " << newheight << '\n'
        << "new width = " << newwidth << '\n';
    }
    if (_res[0]<0) _res[0] = newwidth;
    else if (_res[1]<0) _res[1] = newheight;
    else if (newheight <= _res[1]) _res[1] = newheight;
    else if (newwidth <= _res[0]) _res[0] = newwidth;
    else std::cerr << "unable to resize window to fit!\n";
    window->SetSize(_res[0], _res[1]);
    
    if (_verbose) {
        std::cout << "Resolution is " << _res[0] << " x " << _res[1] << '\n';
    }
    
    compute_all_orbits(_niter, phi);
    update_actors();
    
    // double _hh = _vbounds[1]-_vbounds[0];
    // std::vector<nvis::vec2> frame_pts(4);
    // frame_pts[0] = nvis::vec2(0,         _vbounds[0]);
    // frame_pts[1] = nvis::vec2(_omega*_T, _vbounds[0]);
    // frame_pts[2] = nvis::vec2(_omega*_T, _vbounds[1]);
    // frame_pts[3] = nvis::vec2(0,         _vbounds[1]);
    // vtkPolyData* frame = vtk_utils::make_points(frame_pts);
    // vtk_utils::add_vertices(frame);
    // std::vector<int> frame_lines(8);
    // frame_lines[0] = 0;
    // frame_lines[1] = 1;
    // frame_lines[2] = 1;
    // frame_lines[3] = 2;
    // frame_lines[4] = 2;
    // frame_lines[5] = 3;
    // frame_lines[6] = 3;
    // frame_lines[7] = 0;
    // vtk_utils::add_lines(frame, frame_lines);
    // VTK_CREATE(vtkExtractEdges, frame_edges);
    // frame_edges->SetInput(frame);
    // VTK_CREATE(vtkTubeFilter, tubes);
    // // tubes->SetInput(frame);
    // VTK_PLUG(tubes, frame_edges);
    // tubes->SetRadius(0.01);
    // tubes->SetNumberOfSides(6);
    // VTK_CREATE(vtkPolyDataMapper, frame_mapper);
    // VTK_PLUG(frame_mapper, tubes);
    // VTK_CREATE(vtkActor, frame_actor);
    // frame_actor->SetMapper(frame_mapper);
    // frame_actor->GetProperty()->SetColor(1, 1, 1);
    // frame_actor->GetProperty()->SetLineWidth(4);
    // renderer->AddActor(frame_actor);
    
    VTK_CREATE(vtkInteractorStyleImage, style);
    VTK_CREATE(vtkRenderWindowInteractor, interactor);
    interactor->SetRenderWindow(window);
    interactor->SetInteractorStyle(style);
    interactor->Initialize();
    
    // First render must occur after initialization of interactor.
    // Insight copyright Steve Plite 2013, Purdue Computer Science
    draw();
    
    if (!output.empty()) {
        if (_verbose) {
            std::cout << "exporting image to \"" << output << "\"\n";
        }
        vtkWindowToImageFilter *capture = vtkWindowToImageFilter::New();
        capture->SetInput(window);
        capture->Update();
        vtkPNGWriter *writer = vtkPNGWriter::New();
        writer->SetInputConnection(capture->GetOutputPort());
        writer->SetFileName(output.c_str());
        writer->Write();
        writer->Delete();
        capture->Delete();
    }
    else 
        // start interactor
        interactor->Start();
    
    return 0;
}
