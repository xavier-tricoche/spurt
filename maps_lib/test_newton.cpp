#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <iomanip>
#include <functional>

#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <math/rational.hpp>

#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>

// NRRD interface
#include <image/nrrd_wrapper.hpp>
// data structure
#include <data/grid.hpp>
#include <data/edge.hpp>
#include <data/raster_data.hpp>
// display
#include <graphics/colors.hpp>
#include <graphics/GLUT_helper.hpp>
#include <graphics/GUI/GLUI_Wrapper.h>
#include "logical2physical.hpp"
// poincare map API
#include "xmt_poincare_map.hpp"
#include "map_field_wrapper.hpp"
#include "map_approximation.hpp"
// topological analysis
#include "map_analysis.hpp"
#include "newton.hpp"
#include "fixpoints.hpp"
#include "invariant_manifold.hpp"
// math
#include <math/rational.hpp>

#if _OPENMP
#include <omp.h>
#endif


using namespace spurt;

const double invalid_double = std::numeric_limits<double>::max();

int nbthreads;

// -------------------------
//
//      Data Structure
//
// -------------------------
typedef grid<double, 3>                                     volume_type;
typedef grid<double, 2>                                     plane_type;
typedef raster_data<nvis::vec3, double, 3>                  field_type;
typedef raster_data<orbit_data, double, 2>                  dataset_type;
typedef xmt_poincare_map<spurt::map::wrapper<field_type> > map_type;

map_metric  orbit_data::metric;
int         orbit_data::max_period;

field_type*                 _field;
std::vector<nvis::vec3>     _vectors;
volume_type*                _volume;
nvis::ivec3                 _volume_res;
nvis::bbox3                 _volume_bounds;
nvis::vec3                  _volume_spacing;
plane_type*                 _plane;
nvis::ivec2                 _plane_res;
nvis::bbox2                 _plane_bounds;
spurt::map_metric          _plane_metric;
nvis::vec2                  _plane_spacing;
dataset_type*               _dataset;

// -------------------------
//
//              UI
//
// -------------------------
char*    in, *phys, *topo;
int     niter, width, height, per, show_orb, show_vec, show_edges, cell_analysis;
int     nogfx, fast_angle, cellid[2];
double  eps;
float   pt_sz, ups;
bool    logical;

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
    hestOptAdd(&hopt, "i",      "input",                airTypeString,  1,  1,  &in,        NULL,           "input file name (NRRD)");
    hestOptAdd(&hopt, "p",      "period",               airTypeInt,     1,  1,  &per,       NULL,           "considered period");
    hestOptAdd(&hopt, "x",              "upsample factor",              airTypeFloat,   0,      1,      &ups,           "1",                    "plane upsampling factor");
    hestOptAdd(&hopt, "e",      "eps",                  airTypeDouble,  0,  1,  &eps,       "1.0e-8",       "integration precision");
    hestOptAdd(&hopt, "id",     "cell index",           airTypeInt,     2,  2,  &cellid,    NULL,           "cell index");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const  char*)me, "Test Newton solver",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline bool equal(bool a, bool b)
{
    return ((a && b) || (!a && !b));
}

template<typename InputIterator>
InputIterator epsilon_find(const InputIterator first, const InputIterator last,
                           const nvis::vec2& x, const map_metric& metric, double eps)
{
    double _dist;
    for (InputIterator i=first ; i!=last ; ++i) {
        const nvis::vec2& y = i->pos;
        _dist = metric.distance(x, y);
        if (_dist < eps) {
            return i;
        }
    }
    return last;
}

// -------------------------
//
//          Analysis
//
// -------------------------
typedef rational_surface_found::edge_type               segment_type;
typedef rational_surface_found::edge_point              edge_point;
typedef boost::rational<int>                            rational_type;
typedef edge<nvis::ivec2, nvis::lexicographical_order>  edge_type;
typedef std::vector<nvis::vec2>                         orbit_type;
typedef std::pair<orbit_type, nvis::fvec3>              color_orbit_type;

template<typename T, typename Compare = std::less<T> >
struct tagged {
    typedef tagged<T, Compare>  self_type;
    
    tagged() : t(), p(0) {}
    tagged(T _t, int _p) : t(_t), p(_p) {}
    
    int& tag() {
        return p;
    }
    int tag() const {
        return p;
    }
    T& object() {
        return t;
    }
    const T& object() const {
        return t;
    }
    
    bool operator<(const self_type& hp) const {
        if (p < hp.p) {
            return true;
        } else if (p > hp.p) {
            return false;
        }
        Compare Lt;
        return Lt(t, hp.t);
    }
    
    T   t;
    int p;
};
typedef tagged<edge_type>                                   p_edge_type;
typedef tagged<nvis::vec2, nvis::lexicographical_order>     p_cell_type;

std::map<p_edge_type, double>                   _edge_angles;
std::vector<tagged<p_cell_type> >               _cells_indices;
std::vector<p_edge_type>                        _failed_edges;
std::vector<std::vector<fixpoint> >             _chains;
std::vector<std::vector<separatrix> >           _separatrices;

// -------------------------
//
//          Display
//
// -------------------------
std::vector<std::vector<nvis::vec2> >           _orbits;
spurt::discrete_color_map<int>*                _cmap;
spurt::logical2physical*                       _converter;
spurt::map_analysis_param                      _params;
std::vector<nvis::ivec2>                        _saddle_cells, _center_cells;
nvis::vec2                                      _last;
std::vector<nvis::vec2>                         _problematic_seeds;
std::vector<color_orbit_type>                   _rational_surfaces;
std::vector<nvis::vec2>                         _saddles;
std::vector<nvis::vec2>                         _centers;

struct rational_segment {
    edge_point pt[2];
    rational_type sf;
};

std::vector<rational_segment>           _rational_segments;
std::vector<fixpoint>                   _fixpoints;
std::vector<int>                        _tagged_cells;
double                                  _index_dx;
std::vector<rational_type>              _valid_rationals;


std::vector<std::pair<nvis::vec2, nvis::vec2> > _vecs;
std::vector<std::vector<nvis::vec2> >           _jac_sls;
std::vector<std::vector<nvis::vec2> >           _newton_steps;

enum {
    FIXED_POINT = 0,
    CONVERGED,
    DEGENERATE
};

static void init()
{
    nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    
    Nrrd* nin = nrrdNew();
    nin = spurt::readNrrd(in);
    
    // verify data type
    if(nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        exit(-1);
    }
    
    std::vector<double> _array;
    spurt::to_vector(_array, nin);
    _volume_res = nvis::ivec3(nin->axis[1].size,
                              nin->axis[2].size,
                              nin->axis[3].size);
    _volume_spacing = nvis::vec3(nin->axis[1].spacing,
                                 nin->axis[2].spacing,
                                 nin->axis[3].spacing);
    _volume = new volume_type(_volume_res, _volume_spacing,
                              nvis::fixed_vector<bool, 3>(false, true, true));
    _vectors.resize(_array.size() / 3);
    for(int i = 0 ; i < _array.size() / 3 ; ++i) {
        _vectors[i][0] = _array[3*i  ];
        _vectors[i][1] = _array[3*i+1];
        _vectors[i][2] = _array[3*i+2];
    }
    _field = new field_type(*_volume, _vectors);
    _field->verbose(false);
    
    // check_data_structure(*_field);
    
    double h = eps;
    
    _volume_bounds = _volume->bounds();
    _plane_bounds = nvis::bbox2(nvis::vec2(_volume_bounds.min()[1], _volume_bounds.min()[0]),
                                nvis::vec2(_volume_bounds.max()[1], _volume_bounds.max()[0]));
                                
    int __nx = floor(ups*(_volume_res[1]-1))+1;
    int __ny = floor(ups*(_volume_res[0]-1))+1;
    _plane_res = nvis::ivec2(__nx, __ny);
    _plane_metric.bounds() = _plane_bounds;
    _plane_metric.periodic(0) = true;
    _plane_metric.periodic(1) = false;
    int npoints = _plane_res[0] * _plane_res[1];
    
    std::cerr << "plane resolution = " << _plane_res << std::endl;
    
    orbit_data::metric = _plane_metric;
    _params.metric = _plane_metric;
    
    map_type pmap(*_field);
    pmap.precision(eps);
    _plane = new plane_type(_plane_res, _plane_bounds);
    _dataset = new dataset_type(*_plane, orbit_data());
    
    _plane_spacing = _plane->spacing();
    std::cerr << "plane spacing = " << _plane_spacing << std::endl;
    
    nvis::ivec2 cell(cellid[0], cellid[1]);
    nvis::vec2 min = (*_plane)(cell);
    nvis::vec2 max = min + _plane_spacing;
    min -= 2*_plane_spacing;
    max += 2*_plane_spacing;
    GLUT_helper::box = _plane_bounds;
}

double __edge_rotation(const edge_type& e, int period, double lmin)
{

    map_type pmap(*_field);
    pmap.precision(eps);
    rhs_wrapper<map_type> rhs(pmap, _plane_metric, period);
    _params.verbose = true;
    
    nvis::ivec2 i0 = e[0];
    nvis::ivec2 i1 = e[1];
    nvis::vec2 v0 = rhs((*_plane)(i0));
    nvis::vec2 v1 = rhs((*_plane)(i1));
    
    double theta = 0;
    if (lmin > 0) {
        nvis::vec2 x0 = (*_plane)(i0);
        nvis::vec2 x1 = (*_plane)(i1);
        try {
            return adaptive_rotation_angle(x0, v0, x1, v1, pmap, period, lmin, _params);
        } catch(std::runtime_error& err) {
            std::cerr << "__edge_rotation caught: " << err.what() << std::endl;
            throw std::runtime_error("unable to compute rotation");
        }
    } else {
        return signed_angle(v0, v1);
    }
}

// per edge vector rotation
double compute_edge_rotation(const edge_type& e, int period)
{
    double theta;
    // easy
    // theta = __edge_rotation(e, period, 0);
    // if (fabs(theta) < MAX_ANGLE_VARIATION) {
    //     std::cerr << "direct computation successful\n";
    //     return theta;
    // }
    // difficult
    double lmin = 0.99*std::min(_plane_spacing[0], _plane_spacing[1])/32.;
    try {
        theta = __edge_rotation(e, period, lmin);
        std::cerr << "coarse adaptive computation successful\n";
        return theta;
    } catch(...) {}
    lmin = 0.99*std::min(_plane_spacing[0], _plane_spacing[1])/256;
    return __edge_rotation(e, period, lmin);
}

int compute_cell_index(const nvis::ivec2& cell_id, int period)
{
    bool verbose = false;
    nvis::ivec2 pt[5];
    pt[0] = cell_id;
    pt[1] = cell_id + nvis::ivec2(1,0);
    pt[2] = cell_id + nvis::ivec2(1,1);
    pt[3] = cell_id + nvis::ivec2(0,1);
    pt[4] = cell_id;
    
    double theta = 0;
    bool valid = true;
    std::ostringstream os;
    for (int i=0 ; i<4 ; ++i) {
        edge_type e(pt[i], pt[i+1]);
        try {
            double dtheta = compute_edge_rotation(e, period);
            std::cerr << "angle(edge " << i << ") = " << dtheta << std::endl;
            theta += (i<2 ? 1 : -1) * dtheta;
        } catch(...) {
            std::cerr << "unable to compute cell index\n";
            return 0;
        }
    }
    std::cerr << "total angle is " << theta << std::endl;
    theta /= 2.*M_PI;
    int index = lround(theta);
    std::cerr << "number of rotations = " << theta << ", index = " << index << '\n';
    return index;
}
int main_window;
nvis::vec2 wc;

void idle(void)
{
    // switch context back to main window after GLUI activity
    glutSetWindow(main_window);
    glutPostRedisplay();
}

void draw(void)
{
    using namespace GLUT_helper;
    
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    const plane_type& mesh = *_plane;
    nvis::ivec2 cell(cellid[0], cellid[1]);
    nvis::bbox2 box;
    box.min() = mesh(cell);
    box.max() = box.min() + _plane_spacing;
    
    // draw cell
    draw_quad(box, nvis::fvec3(1, 0, 1), 2);
    
    // draw vectors
    draw_vectors(_vecs, nvis::fvec3(0,0,1), 1);
    
    // draw streamlines
    for (int n=0 ; n<_jac_sls.size() ; ++n) {
        const std::vector<nvis::vec2>& line = _jac_sls[n];
        // std::cerr << "drawing streamline with " << line.size() << " points\n";
        draw_curve(line, nvis::fvec3(0.7,0.7,0.7), 1);
    }
    
    // draw interactively computed orbits
    for (int i=0 ; i<_orbits.size() ; ++i) {
        const std::vector<nvis::vec2>& o = _orbits[i];
        draw_dots(o, nvis::fvec3(1, 0, 0));
    }
    
    // draw newton steps
    // std::cerr << spurt::newton_steps.size() << " newton steps\n";
    draw_curve(spurt::newton_steps, nvis::fvec3(0,1,1), 1);
    
    // draw search steps
    // std::cerr << spurt::search_steps.size() << " search steps\n";
    draw_vectors(spurt::search_steps, nvis::fvec3(1,1,0), 1);
}

void display(void)
{
    GLUT_helper::setup_display(draw);
    glutSwapBuffers();
}

void guiCallback(int)
{
    display();
}

void check_orbit(const nvis::vec2& x)
{
    map_type pmap(*_field);
    pmap.precision(eps);
    
    std::cerr << "period is " << per << std::endl;
    
    std::vector<nvis::vec2> tmp;
    iterate<map_type>(x, tmp, pmap, 25);
    nvis::vec2 y = pmap.map(x, per);
    
    nvis::vec2 v = _plane_metric.displacement(x, y);
    
    _vecs.push_back(std::pair<nvis::vec2, nvis::vec2>(x, v));
    
    std::vector<nvis::vec2> steps;
    steps.push_back(x);
    for (int i=0 ; i<tmp.size() ; ++i) {
        steps.push_back(tmp[i]);
    }
    
    
    double sf = safety_factor(steps, _plane_metric);
    std::cerr << "safety factor at " << x << " is " << sf << '\n';
    for (int i=0 ; i<steps.size() ; ++i) {
        steps[i] = _plane_metric.modulo(steps[i]);
    }
    _orbits.push_back(steps);
}

struct linear_rhs {
    linear_rhs(const nvis::mat2& J, const nvis::vec2& x0, const nvis::vec2& v0)
        : _J(J), _x0(x0), _v0(v0) {}
        
    nvis::vec2 operator()(const nvis::vec2& x) const {
        return _v0 + _J*(x-_x0);
    }
    
    nvis::mat2 _J;
    nvis::vec2 _x0, _v0;
};

template<typename RHS>
void __euler(const RHS& rhs, std::vector<nvis::vec2>& sl, const nvis::vec2& seed,
             const nvis::bbox2& domain)
{

    double dx = 0.01*nvis::norm(domain.size());
    std::list<nvis::vec2> __sl;
    nvis::vec2 x = seed;
    while (domain.inside(x)) {
        __sl.push_back(x);
        nvis::vec2 v = rhs(x);
        double n = nvis::norm(v);
        if (n < 1.0e-6) {
            break;
        }
        v /= n;
        x += dx*v;
    }
    x = seed;
    while (domain.inside(x)) {
        nvis::vec2 v = rhs(x);
        double n = nvis::norm(v);
        if (n < 1.0e-6) {
            break;
        }
        v /= n;
        x -= dx*v;
        __sl.push_front(x);
    }
    sl.clear();
    std::copy(__sl.begin(), __sl.end(), std::back_inserter(sl));
}

void check_cubic(const nvis::vec2& x)
{
    map_type pmap(*_field);
    pmap.precision(eps);
    rhs_wrapper<map_type> rhs(pmap, _plane_metric, per, 0);
    
    nvis::ivec2 cell_id = _plane->local_coordinates(x).first;
    nvis::bbox2 bounds;
    bounds.min() = (*_plane)(cell_id);
    bounds.max() = bounds.min() + _plane_spacing;
    
    _jac_sls.clear();
    map_approximation<rhs_wrapper<map_type> > quartic_rhs(rhs, bounds, 4, nvis::ivec2(10,10));
    
    std::cerr << "entering check_cubic in cell " << cell_id << std::endl;
    
    nvis::vec2 step = bounds.size() / nvis::vec2(9, 9);
    for (int i=0 ; i<9 ; ++i) {
        for (int j=0 ; j<9 ; ++j) {
            nvis::vec2 p = bounds.min() + step*nvis::vec2(i,j);
            nvis::vec2 v = rhs(p);
            _vecs.push_back(std::pair<nvis::vec2, nvis::vec2>(p, v));
        }
    }
    
    double minx = bounds.min()[0];
    double miny = bounds.min()[1];
    double dx = bounds.size()[0]/5;
    double dy = bounds.size()[1]/5;
    
    for (int i=0 ; i<6 ; ++i) {
        double _x = minx + i*dx;
        for (int j=0 ; j<6 ; ++j) {
            double _y = miny + j*dy;
            nvis::vec2 seed(_x, _y);
            std::vector<nvis::vec2> sl;
            std::cerr << "integrating from " << seed << " (" << i << ", " << j << ")..." << std::flush;
            __euler(quartic_rhs, sl, seed, bounds);
            _jac_sls.push_back(sl);
            std::cerr << "done\n";
        }
    }
}

void check_jacobian(const nvis::vec2& x)
{
    double J_eps = _volume_bounds.size()[2]/(5*_volume_res[2]);
    map_type pmap(*_field);
    pmap.precision(eps);
    rhs_wrapper<map_type> rhs(pmap, _plane_metric, per, J_eps);
    
    _jac_sls.clear();
    nvis::vec2 v0 = rhs(x);
    nvis::mat2 J = rhs.jacobian(x);
    linear_rhs linrhs(J, x, v0);
    
    std::cerr << "entering check_jacobian at " << x << ", J = " << J << std::endl;
    
    nvis::bbox2 domain(x-0.25*_plane_spacing, x+0.25*_plane_spacing);
    
    double minx = domain.min()[0];
    double miny = domain.min()[1];
    double dx = domain.size()[0]/5;
    double dy = domain.size()[1]/5;
    
    for (int i=0 ; i<6 ; ++i) {
        double _x = minx + i*dx;
        for (int j=0 ; j<6 ; ++j) {
            double _y = miny + j*dy;
            nvis::vec2 seed(_x, _y);
            std::vector<nvis::vec2> sl;
            __euler(linrhs, sl, seed, domain);
            _jac_sls.push_back(sl);
        }
    }
}


bool save_to_file;
void mouse(int button, int state, int x, int y)
{
    map_type pmap(*_field);
    pmap.precision(eps);
    
    nvis::vec3 _wc = GLUT_helper::world_coordinates(x, y);
    wc = nvis::vec2(_wc[0], _wc[1]);
    
    if (nvis::all(wc == _last)) {
        return;
    } else {
        _last = wc;
    }
    display();
    
    std::ostringstream os;
    os << "Cursor at " << wc;
    glutSetWindowTitle(os.str().c_str());
    GLUT_helper::glut_helper_mouse(button, state, x, y);
}

void keyboard(unsigned char key, int x, int y)
{
    if(key == 'r') {
        GLUT_helper::resetCamera();
    } else if (key == 'i') {
        // iterate
        map_type pmap(*_field);
        pmap.precision(eps);
        
        _orbits.push_back(std::vector<nvis::vec2>());
        iterate<map_type>(_last, _orbits.back(), pmap, 200);
        
        for (int i=0 ; i<_orbits.back().size() ; ++i) {
            nvis::vec2& x = _orbits.back()[i];
            x = _plane_metric.modulo(x);
        }
        display();
    } else if (key == 'j') {
        check_jacobian(_last);
        display();
    } else if (key == 'c') {
        check_cubic(_last);
        display();
    } else if (key == 'q') {
        check_orbit(wc);
    }
    glutPostRedisplay();
}



// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    _last = nvis::vec2(-1,-1);
    
    spurt::record_newton_steps = true;
    spurt::record_search_steps = true;
    
    width = 800;
    height = 800;
    
    glutInit(&argc, argv);
    glutInitDisplayString("samples rgba double alpha");
    glutInitWindowSize(width, height);
    glutInitWindowPosition(20, 20);
    main_window = glutCreateWindow(argv[0]);
    
    // configure OpenGL for aesthetic results
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_POLYGON_SMOOTH);
    glHint(GL_LINE_SMOOTH, GL_NICEST);
    glHint(GL_POINT_SMOOTH, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH, GL_NICEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    initialize(argc, argv);
    init();
    
    int npoints = _plane_res[0] * _plane_res[1];
    
    map_type pmap(*_field);
    pmap.precision(eps);
    
    fixpoint fp;
    std::vector<nvis::vec2> dummy;
    nvis::ivec2 cell(cellid[0], cellid[1]);
    nvis::bbox2 bounds;
    bounds.min() = (*_plane)(cell);
    bounds.max() = (*_plane)(cell + nvis::ivec2(1,1));
    
    compute_cell_index(cell, per);
    nvis::vec2 x = (*_plane)(cell) + 0.5*_plane_spacing;
    
    std::cerr << "running meta_newton in cell " << cell << std::endl;
    newton_steps.clear();
    try {
        meta_newton(pmap, _plane_metric, bounds,
                    x, 5, per, fp, dummy, 1.0e-3, true, 50);
    } catch(std::runtime_error& e) {
        std::cerr << "caught: " << e.what() << std::endl;
    }
    
    display();
    
    // set GLUT callback functions
    glutDisplayFunc(display);
    glutReshapeFunc(GLUT_helper::glut_helper_reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(GLUT_helper::glut_helper_motion);
    glutKeyboardFunc(keyboard);
    
    // adjust camera to mesh
    GLUT_helper::resetCamera();
    
    // Enter GLUT event loop
    glutMainLoop();
    
    
    return 0;
}
