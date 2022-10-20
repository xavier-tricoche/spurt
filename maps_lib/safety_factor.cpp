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

#include <teem/nrrd.h>

#include <data/grid.hpp>
#include <data/raster.hpp>
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

// GLUT business
#include <graphics/GLUT_helper.hpp>
#include <graphics/GUI/GLUI_Wrapper.h>
#include <graphics/colors.hpp>


using namespace spurt;
using namespace map_analysis;
using namespace map_display;
// using namespace div_cleaning;

template<typename Val_> 
using image_type = image<Val_, 3, double, size_t>;

template<typename Val_>
class legacy_wrapper : public image_type<Val_> {
public:
    typedef image_type<Val_> base_type;
    typedef Val_             value_type;
    
    legacy_wrapper(base_type& image, const map_metric& metric) 
        : base_type(image), _metric(metric), _verbose(false) {}
        
    legacy_wrapper(const grid_type& grid, const std::vector<value_type>& val) 
        : base_type(grid, val), _verbose(false) {}
        
    grid_type get_grid() const {
        return grid_type(this->grid());
    }
    
    value_type interpolate(const typename base_type::point_type& p) const {
        return this->value(p);
    }
    
    void verbose(bool v) const {
        _verbose = v;
    }
    
private:
    spurt::map_metric    _metric;
    mutable bool          _verbose;
};


typedef grid<double, 3>                         grid_type;
typedef legacy_wrapper<nvis::vec3>              trilinear_field_type;
typedef nvis::ivec3                             ivec_type;

spurt::map_metric  metric2d;

typedef xmt_poincare_map<spurt::map::wrapper<trilinear_field_type> >   trilinear_map_type;

trilinear_field_type*   _field;
std::vector<nvis::vec3> _vectors;
grid_type*              _domain;

nvis::bbox2 _bounds;

char*    in;
int     niter, width, height;
double  eps;
float   pt_sz;

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
    hestOptAdd(&hopt, "i",  "input",                airTypeString,  1,  1,  &in,        NULL,           "input file name (NRRD)");
    hestOptAdd(&hopt, "m",  "# iter",               airTypeInt,     0,  1,  &niter,     "100",          "number of iterations");
    hestOptAdd(&hopt, "e",  "eps",                  airTypeDouble,  0,  1,  &eps,       "1.0e-8",       "integration precision");
    hestOptAdd(&hopt, "ps", "point size",           airTypeFloat,   0,  1,  &pt_sz,     "2",            "point size for display");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Safety factor probe",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}
// --------------------------------------------------------------------------------

typedef std::vector<nvis::vec2>     orbit_type;
std::vector<std::pair<orbit_type, nvis::fvec3> >    _orbits;
std::vector<nvis::vec2> mesh_edges;

static void init()
{
    Nrrd* nin = nrrdNew();
    nin = spurt::readNrrd(in);
    
    // verify data type
    if(nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        exit(-1);
    }
    
    std::vector<double> _array;
    spurt::to_vector(_array, nin);
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
    
    GLUT_helper::box = _bounds;
    std::cerr << "box set to " << GLUT_helper::box << std::endl;
    
    nvis::vec2 step = _bounds.size() / nvis::vec2(dims[1]-1, dims[0]-1);
    _orbits.clear();
    
    // draw mesh
    std::cerr << "bounds = " << _bounds << std::endl;
    std::cerr << "step = " << step << std::endl;
    mesh_edges.clear();
    double x = _bounds.min()[0];
    for (int i=0 ; i<dims[1] ; ++i, x+=step[0]) {
        mesh_edges.push_back(nvis::vec2(x, _bounds.min()[1]));
        mesh_edges.push_back(nvis::vec2(x, _bounds.max()[1]));
    }
    double y = _bounds.min()[1];
    for (int j=0 ; j<dims[0] ; ++j, y+=step[1]) {
        mesh_edges.push_back(nvis::vec2(_bounds.min()[0], y));
        mesh_edges.push_back(nvis::vec2(_bounds.max()[0], y));
    }
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
        glVertex2f(mesh_edges[2*i  ][0], mesh_edges[2*i  ][1]);
        glVertex2f(mesh_edges[2*i+1][0], mesh_edges[2*i+1][1]);
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
            glVertex2f(x[0], x[1]);
        }
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
       const spurt::map_metric& metric, int n)
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
    
    _orbits.push_back(std::pair<orbit_type, nvis::fvec3>());
    _orbits.back().second = nvis::fvec3(drand48(), drand48(), drand48());
    
    std::pair<double, double> sf = period(wc, _orbits.back().first, metric2d, niter);
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
    initialize(argc, argv);
    save_to_file = false;
    
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
