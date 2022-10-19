#include "period_analysis.hpp"
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <iomanip>

#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <boost/rational.hpp>
#include <kdtree++/kdtree.hpp>

#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>

#include <teem/hest_helper.hpp>

#include <data/grid.hpp>
#include <data/raster_data.hpp>
#include <data/kdtree.hpp>
#include <graphics/colors.hpp>
#include <image/nrrd_wrapper.hpp>
#include <maps_lib/display.hpp>
#include <maps_lib/period.hpp>
#include <math/fixed_vector.hpp>
#include "definitions.hpp"
#include "xmt_poincare_map.hpp"
#include "map_field_wrapper.hpp"
#include <math/rational.hpp>
#include <data/adaptive_quad_mesh.hpp>
#include <maps_lib/logical2physical.hpp>

// GLUT business
#include <graphics/GLUT_helper.hpp>
#include <graphics/GUI/GLUI_Wrapper.h>
#include <graphics/colors.hpp>


using namespace xavier;
using namespace map_analysis;
using namespace map_display;
// using namespace div_cleaning;

typedef grid<double, 3>                         grid_type;
typedef raster_data<nvis::vec3, double, 3>      trilinear_field_type;
typedef nvis::ivec3                             ivec_type;

xavier::map_metric  metric2d;

typedef xmt_poincare_map<xavier::map::wrapper<trilinear_field_type> >   trilinear_map_type;

trilinear_field_type*   _field;
std::vector<nvis::vec3> _vectors;
grid_type*              _domain;
nvis::ivec2             _resolution;
nvis::vec2              _last;
nvis::vec2              _step;

nvis::bbox2 _bounds;

char*    in, *phys;
int     niter, width, height, maxp;
double  eps;
float   pt_sz;
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
    hestOptAdd(&hopt, "n",      "# iter",               airTypeInt,     0,  1,  &niter,     "100",          "number of iterations");
    hestOptAdd(&hopt, "mp",     "max period",           airTypeInt,     0,  1,  &maxp,      "50",           "max considered period");
    hestOptAdd(&hopt, "e",      "eps",                  airTypeDouble,  0,  1,  &eps,       "1.0e-8",       "integration precision");
    hestOptAdd(&hopt, "ps",     "point size",           airTypeFloat,   0,  1,  &pt_sz,     "2",            "point size for display");
    hestOptAdd(&hopt, "phys",   "phys space",           airTypeString,  0,  1,  &phys,      "none",         "HDF5 reference file for mapping");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Safety factor probe",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}
// --------------------------------------------------------------------------------

typedef std::vector<nvis::vec2>                     orbit_type;
std::vector<std::pair<orbit_type, nvis::fvec3> >    _orbits;
std::vector<nvis::vec2>                             mesh_edges;

std::vector<std::pair<nvis::vec2, int> >    vf;
xavier::discrete_color_map<int>*             cmap;

inline double average_distance(const std::vector<nvis::vec2>& steps, int p)
{
    if (steps.size() < p) {
        return std::numeric_limits<double>::max();
    }
    double d=0;
    int count = 0;
    for (int i=0 ; i<steps.size() ; ++i) {
        if (i+p<steps.size()) {
            d += metric2d.distance(steps[i], steps[i+p]);
            ++count;
        }
        if (i>=p) {
            d += metric2d.distance(steps[i-p], steps[i]);
            ++count;
        }
    }
    return d/(double)count;
}

inline int best_period(const std::vector<nvis::vec2>& steps, int maxp)
{
    std::map<double, int> period_to_dist;
    for (int p=1 ; p<=maxp ; ++p) {
        period_to_dist[average_distance(steps, p)] = p;
    }
    
    return period_to_dist.begin()->second;
}

struct cell_data_type {
};

typedef std::pair<nvis::vec2, int>                          data_type;
typedef xavier::AQM_root<double, data_type, cell_data_type> quad_mesh_type;
typedef quad_mesh_type::node_type                           node_type;

xavier::logical2physical*   converter;
quad_mesh_type*             mesh;

struct edge {
    edge(const nvis::ivec3& i0, const nvis::ivec3& i1) {
        nvis::lexicographical_order less;
        if (less(i1, i0)) {
            id0 = i1;
            id1 = i0;
        } else {
            id0 = i0;
            id1 = i1;
        }
    }
    
    edge(const edge& e) {
        id0 = e.id0;
        id1 = e.id1;
    }
    
    bool operator<(const edge& e) const {
        nvis::lexicographical_order less;
        if (less(id0, e.id0)) {
            return true;
        } else if (less(e.id0, id0)) {
            return false;
        } else {
            return less(id1, e.id1);
        }
    }
    
    nvis::ivec3 id0, id1;
};

struct collect_leaves {
    collect_leaves() : __leaves() {}
    void reset() {
        __leaves.clear();
    }
    void operator()(node_type* node) {
        __leaves.push_back(node);
    }
    std::vector<node_type*> __leaves;
};

struct collect_edges {
    collect_edges() : __edges() {}
    void reset() {
        __edges.clear();
    }
    void operator()(const node_type* node) {
        nvis::ivec3 id = node->index();
        __edges.insert(edge(id, id+nvis::ivec3(1,0,0)));
        __edges.insert(edge(id+nvis::ivec3(1,0,0), id+nvis::ivec3(1,1,0)));
        __edges.insert(edge(id+nvis::ivec3(1,1,0), id+nvis::ivec3(0,1,0)));
        __edges.insert(edge(id+nvis::ivec3(0,1,0), id));
    }
    std::set<edge>  __edges;
};

struct collect_leaf_corners {
    collect_leaf_corners() : __vertices() {}
    void reset() {
        __vertices.clear();
    }
    void operator()(const node_type* node) {
        nvis::ivec3 id = node->index();
        __vertices.insert(id);
        __vertices.insert(id+nvis::ivec3(1,0,0));
        __vertices.insert(id+nvis::ivec3(1,1,0));
        __vertices.insert(id+nvis::ivec3(0,1,0));
    }
    std::set<nvis::ivec3, nvis::lexicographical_order>  __vertices;
};

inline nvis::vec2 to_pos(const nvis::ivec3& id, const nvis::bbox2& bounds,
                         const nvis::ivec2& res)
{
    int d = id[2];
    int f = (1 << d);
    nvis::ivec2 r = f*res;
    nvis::vec2 h = bounds.size() / nvis::vec2(r);
    return bounds.min() + nvis::vec2(id[0], id[1])*h;
}

enum {
    FIXED_POINT = 0,
    CONVERGED,
    DEGENERATE
};

int check_quad(node_type* node)
{
    node_type::quad_type q = node->quad();
    std::set<int> available_periods;
    for (int i=0 ; i<4 ; ++i) {
        available_periods.insert(q.point_val[i].second);
    }
    // loop over periods in ascending order
    for (std::set<int>::const_iterator it=available_periods.begin() ;
            it!=available_periods.end() ; ++it) {
            
    }
}

static void init()
{
    Nrrd* nin = nrrdNew();
    nin = xavier::readNrrd(in);
    
    // verify data type
    if(nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        exit(-1);
    }
    
    std::vector<double> _array;
    xavier::to_vector(_array, nin);
    ivec_type dims(nin->axis[1].size, nin->axis[2].size, nin->axis[3].size);
    nvis::vec3 spc(nin->axis[1].spacing, nin->axis[2].spacing, nin->axis[3].spacing);
    _domain = new grid_type(dims, spc, nvis::fixed_vector<bool, 3>(false, true, true));
    _vectors.resize(_array.size() / 3);
    for(int i = 0 ; i < _array.size() / 3 ; ++i) {
        _vectors[i][0] = _array[3*i  ];
        _vectors[i][1] = _array[3*i+1];
        _vectors[i][2] = _array[3*i+2];
    }
    _field = new trilinear_field_type(*_domain, _vectors);
    _field->verbose(false);
    
    double h = eps;
    
    grid_type::bounds_type bbox(_domain->bounds());
    _bounds = nvis::bounding_box<nvis::vec2> (nvis::vec2(bbox.min()[1], bbox.min()[0]),
              nvis::vec2(bbox.max()[1], bbox.max()[0]));
              
    metric2d.bounds() = _bounds;
    metric2d.periodic(0) = true;
    metric2d.periodic(1) = false;
    
    if (logical) {
        GLUT_helper::box = _bounds;
    } else {
        GLUT_helper::box = converter->bounds();
    }
    std::cerr << "box set to " << GLUT_helper::box << std::endl;
    
    mesh = new quad_mesh_type(dims[1]-1, dims[0]-1, _bounds);
    
    _step = _bounds.size() / nvis::vec2(dims[1]-1, dims[0]-1);
    _orbits.clear();
    trilinear_map_type base_map(*_field);
    base_map.precision(eps);
    _resolution = nvis::ivec2(dims[1], dims[0]);
    
    std::vector<data_type>  _buffer(dims[0]*dims[1]);
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n = 0 ; n < dims[0]*dims[1] ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            int j = n / dims[1];
            int i = n % dims[1];
            nvis::vec2 x0 = _bounds.min() + _step * nvis::vec2(i,j);
            
            trilinear_map_type* amap = base_map.clone();
            std::vector<nvis::vec2> tmp;
            try {
                amap->map(x0, tmp, niter);
            } catch(...) {
                _buffer[n] = data_type(nvis::vec2(0,0), 0);
                continue;
            }
            std::vector<nvis::vec2> steps;
            steps.push_back(x0);
            std::copy(tmp.begin(), tmp.end(), std::back_inserter(steps));
            
            int period = best_period(steps, maxp);
            nvis::vec2 v = metric2d.displacement(steps[0], steps[period]);
            if (nvis::norm(v) > 0.1*metric2d.bounds().size()[0]) {
                _buffer[n] = data_type(nvis::vec2(0, 0), 0);
            } else {
                _buffer[n] = data_type(v, period);
            }
            
            if (thread_id == 0) {
                std::ostringstream os;
                os << "\rcompleted " << n << " / " << dims[0]*dims[1] << " ("
                   << 100*n/dims[0]/dims[1] << "%), added vector " << v << " of period " << period << "      \r";
                std::cout << os.str();
            }
        }
    }
    std::cout << '\n';
    for (int n=0 ; n<_buffer.size() ; ++n) {
        int i = n%dims[1];
        int j = n/dims[1];
        mesh->point_data().add(nvis::ivec3(i,j,0), _buffer[n]);
    }
    
    // draw mesh
    std::cerr << "collecting edges... " << std::flush;
    collect_edges __ce;
    mesh->for_each_leaf_const(__ce);
    std::cerr << "done. " << __ce.__edges.size() << " unique edges\n";
    mesh_edges.clear();
    
    for (std::set<edge>::const_iterator it=__ce.__edges.begin() ;
            it!=__ce.__edges.end() ; ++it) {
        const edge& e = *it;
        mesh_edges.push_back(to_pos(e.id0, _bounds, _resolution));
        mesh_edges.push_back(to_pos(e.id1, _bounds, _resolution));
    }
    
    //
    
    std::vector<nvis::fvec3> colors;
    xavier::spiral_scale(colors, maxp, 1);
    std::vector<int> reference_values(maxp);
    for (int i=1 ; i<=maxp ; ++i) {
        reference_values[i] = i;
    }
    cmap = new xavier::discrete_color_map<int>(reference_values, colors);
}

// --------------------------------------------------------------------------------


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
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glColor3f(0.2, 0.2, 0.2);
    glLineWidth(1.0);
    for (int i=0 ; i<mesh_edges.size()/2 ; ++i) {
        glBegin(GL_LINES);
        nvis::vec2 x = mesh_edges[2*i  ];
        nvis::vec2 y = mesh_edges[2*i+1];
        if (!logical) {
            x = (*converter)(x);
            y = (*converter)(y);
        }
        glVertex2f(x[0], x[1]);
        glVertex2f(y[0], y[1]);
        glEnd();
    }
    
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(pt_sz);
    glBegin(GL_POINTS);
    for (int i=0 ; i<_orbits.size() ; ++i) {
        const nvis::fvec3& c = _orbits[i].second;
        glColor3f(c[0], c[1], c[2]);
        for (int j=0 ; j<_orbits[i].first.size() ; ++j) {
            nvis::vec2 x = metric2d.modulo(_orbits[i].first[j]);
            if (!logical) {
                x = (*converter)(x);
            }
            glVertex2f(x[0], x[1]);
        }
    }
    glEnd();
    
    collect_leaf_corners clc;
    mesh->for_each_leaf_const(clc);
    std::cerr << "done. " << clc.__vertices.size() << " unique vertices\n";
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(1.0);
    const double cos_alpha = cos(M_PI / 12.);
    const double sin_alpha = sin(M_PI / 12.);
    glBegin(GL_LINES);
    for (std::set<nvis::ivec3, nvis::lexicographical_order>::const_iterator it=clc.__vertices.begin() ;
            it!=clc.__vertices.end() ; ++it) {
        const data_type& pd = mesh->get_point_data(*it);
        int period = pd.second;
        nvis::fvec3 col = (*cmap)(period);
        glColor3f(col[0], col[1], col[2]);
        int i = (*it)[0];
        int j = (*it)[1];
        nvis::vec2 p0 = _bounds.min() + nvis::vec2(i,j)*_step;
        nvis::vec2 p1 = p0 + pd.first;
        if (!logical) {
            p0 = (*converter)(p0);
            p1 = (*converter)(p1);
        }
        glVertex2f(p0[0], p0[1]);
        glVertex2f(p1[0], p1[1]);
        nvis::vec2 e0 = 0.2 * (p0 - p1);
        nvis::vec2 e1(-e0[1], e0[0]);
        glVertex2f(p1[0], p1[1]);
        nvis::vec2 y = p1 + cos_alpha * e0 + sin_alpha * e1;
        glVertex2f(y[0], y[1]);
        glVertex2f(p1[0], p1[1]);
        y = p1 + cos_alpha * e0 - sin_alpha * e1;
        glVertex2f(y[0], y[1]);
    }
    glEnd();
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

std::pair<double, double>
period(const nvis::vec2& x0, std::vector<nvis::vec2>& steps,
       const xavier::map_metric& metric, int n)
{
    trilinear_map_type _map(*_field);
    _map.precision(eps);
    try {
        _map.map(x0, steps, n);
    } catch (...) {
        return std::pair<double, double>(-1, -1);
    }
    
    return map_analysis::period_x_periodic(steps, metric);
}

bool save_to_file;
void mouse(int button, int state, int x, int y)
{
    nvis::vec3 _wc = GLUT_helper::world_coordinates(x, y);
    wc = nvis::vec2(_wc[0], _wc[1]);
    
    if (nvis::all(wc == _last)) {
        return;
    } else {
        _last = wc;
    }
    
    _orbits.push_back(std::pair<orbit_type, nvis::fvec3>());
    std::pair<double, double> sf = period(wc, _orbits.back().first, metric2d, 10*niter);
    
    _orbits.back().second = (*cmap)(best_period(_orbits.back().first, maxp));
    _orbits.back().first.push_back(wc);
    display();
    
    std::cout << "safety factor at " << wc << " = "
              << sf.first << std::endl;
              
    std::ostringstream os;
    if (save_to_file) {
        std::fstream fs;
        os << wc[0] << "_" << wc[1] << "sf.txt";
        std::cerr << "exporting " << os.str() << '\n';
        
        fs.open(os.str().c_str(), std::ios::out);
        const std::vector<nvis::vec2>& steps = _orbits.back().first;
        for (int i=0 ; i<steps.size() ; ++i) {
            nvis::vec2 dist = steps[i] - wc;
            double q = (double)(i + 1) / dist[0] * metric2d.width();
            nvis::vec2 s = metric2d.displacement(wc, steps[i]);
            fs << q << "\t" << s[0] << "\t" << s[1] << "\t" << nvis::norm(s) << '\n';
        }
        fs.close();
        os.clear();
        os.str("");
    }
    
    os << "q(" << wc << ")=" << sf.first;
    
    save_to_file = false;
    
    glutSetWindowTitle(os.str().c_str());
    GLUT_helper::glut_helper_mouse(button, state, x, y);
}

void keyboard(unsigned char key, int x, int y)
{
    if(key == 'r') {
        GLUT_helper::resetCamera();
    } else if (key == 's') {
        save_to_file = true;
    } else if(key == 'x') {
        // output position coordinates
    }
    glutPostRedisplay();
}

// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    _last = nvis::vec2(-1,-1);
    
    initialize(argc, argv);
    save_to_file = false;
    
    if(strcmp(phys, "none")) {
        logical = false;
        converter = new xavier::logical2physical(phys);
    } else {
        logical = true;
    }
    
    width = 1200;
    height = 800;
    
    // initialize GLUT
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
    
    init();
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
