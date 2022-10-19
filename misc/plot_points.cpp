#include <iostream>
#include <sstream>
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

#include <graphics/colors.hpp>
#include <image/nrrd_wrapper.hpp>

// GLUT business
#include <graphics/GLUT_helper.hpp>
#include <graphics/GUI/GLUI_Wrapper.h>

#if _OPENMP
#include <omp.h>
#endif

// linear algebra
#define ARMA_USE_LAPACK
#include "armadillo"

using namespace xavier;
using namespace GLUT_helper;
int main_window;
nvis::vec2 wc;

char*    in;
int     npts, xaxis, yaxis;
float   ps;

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
    hestOptAdd(&hopt, "i",  "input",        airTypeString,  1,  1,  &in,        NULL,           "input file name (NRRD)");
    hestOptAdd(&hopt, "ps", "point size",   airTypeFloat,   0,  1,  &ps,        "1",            "point size for display");
    hestOptAdd(&hopt, "x",  "x axis",       airTypeInt,     0,  1,  &xaxis,     "0",            "column for x axis");
    hestOptAdd(&hopt, "y",  "y axis",       airTypeInt,     0,  1,  &yaxis,     "1",            "column for y axis");
    
    _hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                    me, "Plot 2D points",
                    AIR_TRUE, AIR_TRUE, AIR_TRUE);
}
void guiCallBack(int control);

// --------------------------------------------------------------------------------

nvis::bbox2 bbox;
std::vector<nvis::vec2> pts;
static void init()
{
    Nrrd* nin = nrrdNew();
    nin = xavier::nrrd_utils::readNrrd(in);
    
    assert(nin->dim == 2 &&
           xaxis >= 0 && yaxis >= 0 &&
           nin->axis[0].size == 6 &&
           std::max(xaxis, yaxis) < 6);
           
    std::vector<double> array;
    xavier::to_vector(array, nin);
    
    int npts = nin->axis[1].size;
    bbox.reset();
    pts.resize(npts);
    
    typedef nvis::fixed_vector<double, 6> vec6;
    
    vec6 mean(0);
    for (int i = 0 ; i < npts ; ++i) {
        vec6 pt;
        for (int n = 0 ; n < 6 ; ++n) {
            pt[n] = array[6*i+n];
        }
        mean += pt;
    }
    mean /= (double)npts;
    
    // set up covariance matrix
    arma::mat covar(6, 6);
    for (int i = 0 ; i < npts ; ++i) {
        vec6 pt;
        for (int n = 0 ; n < 6 ; ++n) {
            pt[n] = array[6*i+n];
        }
        pt -= mean;
        for (int r = 0 ; r < 6 ; ++r) {
            for (int c = r ; c < 6 ; ++c) {
                covar(r, c) += pt[r] * pt[c];
            }
        }
    }
    for (int r = 1 ; r < 6 ; ++r) {
        for (int c = 0 ; c < r ; ++c) {
            covar(r, c) = covar(c, r);
        }
    }
    covar *= 1. / (double)npts;
    arma::vec eigval;
    arma::mat eigvec;
    
    eig_sym(eigval, eigvec, covar);
    std::map<double, int> evalues;
    for (int i = 0 ; i < 6 ; ++i) {
        evalues.insert(std::pair<double, int>(eigval[i], i));
    }
    std::cerr << "sorted eigenvectors are:\n";
    vec6 ex, ey;
    int count = 0;
    for (std::map<double, int>::const_iterator it = evalues.begin() ; it != evalues.end() ; ++it, ++count) {
        int id = it->second;
        vec6 evec;
        for (int i = 0 ; i < 6 ; ++i) {
            evec[i] = eigvec(i, id);
        }
        std::cerr << "#" << count << ": eval = " << it->first << ", evec = " << evec << '\n';
        if (count == xaxis) {
            ex = evec;
        } else if (count == yaxis) {
            ey = evec;
        }
    }
    
    for (int i = 0 ; i < npts ; ++i) {
        vec6 pt;
        for (int n = 0 ; n < 6 ; ++n) {
            pt[n] = array[6*i+n];
        }
        pt -= mean;
        pts[i][0] = nvis::inner(pt, ex);
        pts[i][1] = nvis::inner(pt, ey);
        bbox.add(pts[i]);
    }
    
    GLUT_helper::box = bbox;
    std::cerr << "bounds are " << bbox << std::endl;
}

// --------------------------------------------------------------------------------

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
    
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(ps);
    glColor3f(1, 1, 1);
    
    glBegin(GL_POINTS);
    for (int i = 0 ; i < pts.size() ; ++i) {
        glVertex2f(pts[i][0], pts[i][1]);
    }
    glEnd();
}

void display(void)
{
    setup_display(draw);
    glutSwapBuffers();
}

void guiCallback(int)
{
    display();
}

void mouse(int button, int state, int x, int y)
{
    nvis::vec3 _wc = world_coordinates(x, y);
    wc = nvis::vec2(_wc[0], _wc[1]);
    std::ostringstream os;
    os << "Poincare plot - (" << wc[0] << ", " << wc[1] << ")";
    glutSetWindowTitle(os.str().c_str());
    glut_helper_mouse(button, state, x, y);
}

void keyboard(unsigned char key, int x, int y)
{
    if (key == 'r') {
        resetCamera();
    }
    // else if (key == 'x') {
    //  // output position coordinates
    // }
    glutPostRedisplay();
}

// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    initialize(argc, argv);
    
    width = 1200;
    height = 1200;
    
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
    glutReshapeFunc(glut_helper_reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(glut_helper_motion);
    glutKeyboardFunc(keyboard);
    
    // adjust camera to mesh
    resetCamera();
    resetCamera();
    
    // Enter GLUT event loop
    glutMainLoop();
    
    return 0;
}





