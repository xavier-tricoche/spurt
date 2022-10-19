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
typedef divfree_field<trilinear_field_type>     divfree_field_type;
typedef nvis::ivec3                             ivec_type;

xavier::map_metric  metric2d;

typedef xmt_poincare_map<xavier::map::wrapper<trilinear_field_type> >   trilinear_map_type;
typedef xmt_poincare_map<xavier::map::wrapper<divfree_field_type> >     divfree_map_type;

nvis::bbox2 _bounds;

char*    in, *seed_l, *seed_s, *out;
int     niter, rx, dolink, width, height, maxdepth, minp, maxp;
double  eps, dq, xf;
bool    div_free;
float   pt_sz, ln_w;

std::vector<std::vector<nvis::vec2> > xavier::broken_manifolds;

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
    hestOptAdd(&hopt, "o",  "output",               airTypeString,  0,  1,  &out,       "none",         "output file name (NRRD)");
    hestOptAdd(&hopt, "m",  "# iter",               airTypeInt,     0,  1,  &niter,     "100",          "number of iterations");
    hestOptAdd(&hopt, "d",  "max depth",            airTypeInt,     0,  1,  &maxdepth,  "4",            "max refinement depth");
    hestOptAdd(&hopt, "x",  "magnification",        airTypeDouble,  0,  1,  &xf,        "1",            "magnification coefficient");
    hestOptAdd(&hopt, "dq", "delta q",              airTypeDouble,  0,  1,  &dq,        "0.05",         "max q discontinuity");
    hestOptAdd(&hopt, "min","min period",           airTypeInt,     0,  1,  &minp,      "1",            "minimum period to consider");
    hestOptAdd(&hopt, "max","max period",           airTypeInt,     0,  1,  &maxp,      "15",           "maximum period to consider");
    hestOptAdd(&hopt, "e",  "eps",                  airTypeDouble,  0,  1,  &eps,       "1.0e-8",       "integration precision");
    hestOptAdd(&hopt, "ps", "point size",           airTypeFloat,   0,  1,  &pt_sz,     "2",            "point size for display");
    hestOptAdd(&hopt, "lw", "line width",           airTypeFloat,   0,  1,  &ln_w,      "1",            "line width for display");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Poincare map visualization",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}
// --------------------------------------------------------------------------------

std::vector<nvis::vec2> edges, interesting_edges;
std::vector<std::pair<double, std::vector<nvis::vec2> > > orbits;
std::vector<std::pair<nvis::ivec3, nvis::vec3> > quads;
xavier::band_color_map<double>* cmap;

static void init()
{
    Nrrd* nin = nrrdNew();
    nin = xavier::readNrrd(in);
    
    // verify data type
    if(nin->dim != 4 || nin->axis[0].size != 3) {
        std::cerr << "invalid input NRRD file.\n";
        exit(-1);
    }
    
    std::vector<double> __array;
    xavier::to_vector(__array, nin);
    ivec_type dims(nin->axis[1].size, nin->axis[2].size, nin->axis[3].size);
    nvis::vec3 spc(nin->axis[1].spacing, nin->axis[2].spacing, nin->axis[3].spacing);
    grid_type domain(dims, spc, nvis::fixed_vector<bool, 3>(false, true, true));
    std::vector<nvis::vec3> vectors(__array.size() / 3);
    for(int i = 0 ; i < __array.size() / 3 ; ++i) {
        vectors[i][0] = __array[3*i  ];
        vectors[i][1] = __array[3*i+1];
        vectors[i][2] = __array[3*i+2];
    }
    trilinear_field_type basic_field(domain, vectors);
    basic_field.verbose(false);
    
    double h = eps;
    
    trilinear_map_type basic_map(basic_field);
    basic_map.precision(eps);
    
    grid_type::bounds_type bbox(domain.bounds());
    _bounds = nvis::bounding_box<nvis::vec2> (nvis::vec2(bbox.min()[1], bbox.min()[0]),
              nvis::vec2(bbox.max()[1], bbox.max()[0]));
              
    metric2d.bounds() = _bounds;
    metric2d.periodic(0) = true;
    metric2d.periodic(1) = false;
    
    GLUT_helper::box = _bounds;
    std::cerr << "box set to " << GLUT_helper::box << std::endl;
    
    nvis::ivec2 resolution(floor(dims[1]*xf), floor(dims[0]*xf));
    
    std::vector<double> qs(orbits.size());
    for (int i=0 ; i<orbits.size() ; ++i) {
        qs[i] = orbits[i].first;
    }
    
    interesting_edges.clear();
    std::set<double> vals;
    typedef boost::rational<int> rational;
    for (int i=minp ; i<=maxp ; ++i) {
        for (int j=1 ; j<=std::min(i,3) ; ++j) {
            rational q(i,j);
            vals.insert(xavier::value<int, double>(q));
        }
    }
    std::cerr << "there are " << vals.size() << " relevant rational periods for considered periods\n";
    std::copy(vals.begin(), vals.end(), std::ostream_iterator<double>(std::cerr, " "));
    
    std::vector<double> copy(vals.begin(), vals.end());
    xavier::period_convergence_predicate predicate(3, copy, dq, 0.1*dq);
    
    xavier::period_analysis<trilinear_map_type, xavier::period_convergence_predicate>
    (basic_map, metric2d, resolution, predicate,
     niter, maxdepth, edges, orbits, quads);
     
     
    std::vector<nvis::fvec3> scale(vals.size());
    std::vector<double> cp(vals.size()-2);
    for (int i=0 ; i<copy.size()-1 ; ++i) {
        cp[i] = 0.5*(copy[i] + copy[i+1]);
    }
    std::vector<nvis::fvec3> colors;
    xavier::spiral_scale(colors, vals.size(), 0.5);
    cmap = new xavier::band_color_map<double>(cp, colors);
    
    for (int i=0 ; i<quads.size() ; ++i) {
        double min = (quads[i].second)[0];
        double max = (quads[i].second)[2];
        bool selected = false;
        for (std::set<double>::const_iterator it=vals.begin() ; it!=vals.end() && !selected; ++it) {
            double q = *it;
            if (q > min-dq || q < max+dq) {
                selected = true;
            }
        }
        if (selected) {
            nvis::ivec3 idx = quads[i].first;
            nvis::ivec2 r = resolution * (1 << idx[2]);
            nvis::vec2 h = _bounds.size() / nvis::vec2(r[0], r[1]);
            nvis::vec2 x(h[0], 0);
            nvis::vec2 y(0, h[1]);
            nvis::vec2 p0 = _bounds.min() + h*nvis::vec2(idx[0], idx[1]);
            interesting_edges.push_back(p0);
            interesting_edges.push_back(p0+x);
            interesting_edges.push_back(p0+x+y);
            interesting_edges.push_back(p0+y);
        }
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
    glClearColor(1, 1, 1, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glColor3f(0., 0., 0.);
    glLineWidth(ln_w);
    for (int i=0 ; i<edges.size()/4 ; ++i) {
        glBegin(GL_LINE_STRIP);
        glVertex2f(edges[4*i][0], edges[4*i][1]);
        glVertex2f(edges[4*i+1][0], edges[4*i+1][1]);
        glVertex2f(edges[4*i+2][0], edges[4*i+2][1]);
        glVertex2f(edges[4*i+3][0], edges[4*i+3][1]);
        glVertex2f(edges[4*i][0], edges[4*i][1]);
        glEnd();
    }
    
    glColor3f(1., 0., 0.);
    glLineWidth(3*ln_w);
    for (int i=0 ; i<interesting_edges.size()/4 ; ++i) {
        glBegin(GL_LINE_STRIP);
        glVertex3f(interesting_edges[4*i][0], interesting_edges[4*i][1], 1);
        glVertex3f(interesting_edges[4*i+1][0], interesting_edges[4*i+1][1], 1);
        glVertex3f(interesting_edges[4*i+2][0], interesting_edges[4*i+2][1], 1);
        glVertex3f(interesting_edges[4*i+3][0], interesting_edges[4*i+3][1], 1);
        glVertex3f(interesting_edges[4*i][0], interesting_edges[4*i][1], 1);
        glEnd();
    }
    
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(pt_sz);
    glBegin(GL_POINTS);
    for (int i=0 ; i<orbits.size() ; ++i) {
        nvis::fvec3 c = (*cmap)(orbits[i].first);
        glColor3f(c[0], c[1], c[2]);
        for (int j=0 ; j<orbits[i].second.size() ; ++j) {
            const nvis::vec2 x = orbits[i].second[j];
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

void mouse(int button, int state, int x, int y)
{
    nvis::vec3 _wc = GLUT_helper::world_coordinates(x, y);
    wc = nvis::vec2(_wc[0], _wc[1]);
    std::ostringstream os;
    os << "Test - (" << wc[0] << ", " << GLUT_helper::box.max()[1] - wc[1] << ")";
    glutSetWindowTitle(os.str().c_str());
    GLUT_helper::glut_helper_mouse(button, state, x, y);
}

void keyboard(unsigned char key, int x, int y)
{
    if(key == 'r') {
        GLUT_helper::resetCamera();
    } else if(key == 'x') {
        // output position coordinates
    }
    glutPostRedisplay();
}

// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    initialize(argc, argv);
    
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
