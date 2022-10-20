#include <iostream>
#include <list>

#include <math/fixed_vector.hpp>
#include <vector>
#include <maps_lib/display.hpp>
#include <maps_lib/definitions.hpp>
#include <data/grid.hpp>
#include <data/raster.hpp>
#include <math/divergence_cleaning.hpp>

#include <util/wall_timer.hpp>

#include <graphics/GLUT_helper.hpp>

#include <teem/hest.h>
#include <image/nrrd_wrapper.hpp>

#include <format/read_nathan_poincare.hpp>

#include <boost/random.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/prim_minimum_spanning_tree.hpp>

#if _OPENMP
#include <omp.h>
#endif

using namespace spurt;
using namespace map_analysis;
using namespace map_display;

char*    in;
int     n, rx;
int     _link;

nvis::bbox2     _bounds;


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
    hestOptAdd(&hopt, "i",  "input",                airTypeString,  1,  1,  &in,        NULL,           "input file name");
    hestOptAdd(&hopt, "n",  "number displayed",     airTypeInt,     0,  1,  &n,         "100",          "number of displayed particle paths");
    hestOptAdd(&hopt, "l",  "link strategy",        airTypeInt,     0,  1,  &_link,     "0",            "0: none, 1: MST, 2: best period");
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Poincare map visualization of particle trajectories",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

typedef std::vector<nvis::vec2>         orbit_type;
std::vector<orbit_type>                 pplot;      // poincare plot



void display();

// minimum spanning tree
void compute_mst(const std::vector<nvis::vec2>& pts, std::vector<std::pair<int, int> >& edges)
{
    using namespace boost;
    
    typedef adjacency_list<vecS, vecS, undirectedS, no_property, property<edge_weight_t, double> >  Graph;
    typedef graph_traits<Graph>::edge_descriptor                                                    Edge;
    
    Graph g(pts.size());
    
    for (unsigned int i = 0; i < pts.size(); ++i)
        for (unsigned int j = i + 1; j < pts.size(); ++j) {
            double l = nvis::norm(pts[i] - pts[j]);
            
            // if (l < 0.01)
            add_edge(i, j, l, g);
        }
        
    std::vector<Edge> spanning_tree;
    std::vector<graph_traits<Graph>::vertex_descriptor > p(num_vertices(g));
    prim_minimum_spanning_tree(g, &p[0]);
    
    std::vector<int> parents = std::vector<int>(p.begin(), p.end());
    edges.clear();
    for (int i = 0 ; i < parents.size() ; ++i) {
        if (parents[i] != i) {
            edges.push_back(std::pair<int, int>(parents[i], i));
        }
    }
}

void compute_best_period(const std::vector<nvis::vec2>& pts, std::vector<std::pair<int, int> >& edges)
{
    std::map<double, int> lengths;
    for (int i = 1 ; i < pts.size() / 2 ; ++i) {
        double l = 0;
        int count = 0;
        for (int n = 0 ; n < pts.size() ; ++n) {
            if (n + i < pts.size()) {
                l += nvis::norm(pts[n] - pts[n+i]);
                ++count;
            }
        }
        
        l /= (double)count;
        lengths[l] = i;
    }
    
    int period = lengths.begin()->second;
    edges.clear();
    for (int i = 0 ; i < pts.size() / 2 ; ++i) {
        if (i + period < pts.size()) {
            edges.push_back(std::pair<int, int>(i, i + period));
        }
    }
}


static void init()
{
    std::vector<std::vector<nvis::vec3> > __pplot;
    nvis::bbox3 bounds = spurt::read_nathan_poincare(__pplot, in);
    _bounds = nvis::bbox2(nvis::vec2(bounds.min()[0], bounds.min()[1]),
                          nvis::vec2(bounds.max()[0], bounds.max()[1]));
    nvis::vec2 __size = _bounds.size();
    _bounds.min()[0] += 0.65 * __size[0];
    _bounds.max()[1] -= 0.55 * __size[1];
    pplot.resize(__pplot.size());
    
    if (_link > 0) {
        for (int i = 0 ; i < __pplot.size() ; ++i) {
            std::vector<nvis::vec2> unsorted(__pplot[i].size());
            for (int j = 0 ; j < __pplot[i].size() ; ++j) {
                const nvis::vec3& x = __pplot[i][j];
                unsorted[j] = nvis::vec2(x[0], x[1]);
            }
            std::vector<std::pair<int, int> > edges;
            // std::cerr << "computing MST\n";
            if (_link == 1) {
                compute_mst(unsorted, edges);
            } else {
                compute_best_period(unsorted, edges);
            }
            
            pplot[i].clear();
            for (int j = 0 ; j < edges.size() ; ++j) {
                pplot[i].push_back(unsorted[edges[j].first]);
                pplot[i].push_back(unsorted[edges[j].second]);
            }
        }
    } else {
        for (int i = 0 ; i < pplot.size() ; ++i) {
            pplot[i].resize(__pplot[i].size());
            for (int j = 0 ; j < __pplot[i].size() ; ++j) {
                const nvis::vec3& x = __pplot[i][j];
                pplot[i][j] = nvis::vec2(x[0], x[1]);
            }
        }
    }
    
    std::random_shuffle(pplot.begin(), pplot.end());
    
//
//
// // initialize OpenGL
//
// if (strcmp(out, "none")) {
//  int nx = _bounds.size()[0] * rx;
//  int ny = _bounds.size()[1] * rx;
//  unsigned char *raster = (unsigned char*)calloc(3 * nx * ny, sizeof(unsigned char));
//  std::vector<orbit_type>::iterator it;
//  for (it = pplot.begin() ; it != pplot.end() ; ++it) {
//      nvis::fvec3 c(drand48(), drand48(), drand48());
//      c *= 256;
//      unsigned char col[3] = { floor(c[0]), floor(c[1]), floor(c[2]) };
//      orbit_type::const_iterator it2;
//      for (it2 = it->begin() ; it2 != it->end() ; ++it2) {
//          nvis::vec2 x =  rx * *it2;
//          int id = floor(x[0]) + nx * floor(x[1]);
//
//          raster[3*id] = col[0];
//          raster[3*id+1] = col[1];
//          raster[3*id+2] = col[2];
//      }
//  }
//
//  Nrrd *nout = nrrdNew();
//  size_t __size[] = {3, nx, ny};
//  nrrdWrap_nva(nout, raster, nrrdTypeUChar, 3, __size);
//  nrrdSave(out, nout, NULL);
//  std::cerr << "poincare plot exported\n";
// }


    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    
    glShadeModel(GL_FLAT);
    glClearColor(1.0, 1.0, 1.0, 0.0);
}

// --------------------------------------------------------------------------------

void display(void)
{
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    glEnableClientState(GL_VERTEX_ARRAY);
    
#if 0
    // draw mesh
    glEnable(GL_BLEND);
    glLineWidth(1.0);
    glColor3f(.5, .5, .5);
    glEnable(GL_LINE_SMOOTH);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glVertexPointer(2, GL_DOUBLE, 0, &grid.front());
    glDrawArrays(GL_QUADS, 0, grid.size());
    
    // draw cursor
    glEnable(GL_BLEND);
    glPointSize(10.0);
    glColor3f(1.0, 0.7, 0);
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
    glVertex2dv((GLdouble*)&cursor);
    glEnd();
    
    // draw seed candidates
    glEnable(GL_BLEND);
    glPointSize(3.0);
    glColor3f(0.0, 0.0, 0.0);
    glEnable(GL_POINT_SMOOTH);
    glBegin(GL_POINTS);
    for (std::vector<seedpoint>::iterator si = spts.begin(); si != spts.end(); ++si) {
        glVertex2dv((GLdouble*)&si->pos);
    }
    glEnd();
#endif
    
    // draw poincare plot - each orbit in a different color
    
    unsigned int count = 0;
    glDisable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glPointSize(2);
    glLineWidth(1);
    for (int i = 0 ; i < n ; ++i) {
        glEnableClientState(GL_VERTEX_ARRAY);
        glDisable(GL_BLEND);
        nvis::fvec3 c(drand48(), drand48(), drand48());
        glColor3f(c[0], c[1], c[2]);
        glVertexPointer(2, GL_DOUBLE, 0, &pplot[i].front());
        glDrawArrays(GL_POINTS, 0, pplot[i].size());
        if (_link > 0) {
            glVertexPointer(2, GL_DOUBLE, 0, &pplot[i].front());
            glDrawArrays(GL_LINES, 0, pplot[i].size());
            glDisableClientState(GL_VERTEX_ARRAY);
        }
    }
    
#if 0
    glDisable(GL_BLEND);
    glPointSize(1.0);
    glColor3f(0.0, 0.0, 0.5);
    glDisable(GL_POINT_SMOOTH);
    glVertexPointer(2, GL_DOUBLE, 0, &dpts.front());
    glDrawArrays(GL_POINTS, 0, dpts.size());
    
    glDisable(GL_BLEND);
    glLineWidth(1.0);
    glColor3f(0.4, 0.2, 0.2);
    glDisable(GL_LINE_SMOOTH);
    glVertexPointer(2, GL_DOUBLE, 0, &lpts.front());
    glDrawArrays(GL_LINES, 0, lpts.size());
    
    // draw newton iteration points
    glEnable(GL_BLEND);
    glLineWidth(2.0);
    glColor3f(0.5, 0.5, 0.5);
    glEnable(GL_LINE_SMOOTH);
    glVertexPointer(2, GL_DOUBLE, 0, &npts.front());
    glDrawArrays(GL_LINE_STRIP, 0, npts.size());
    
    glDisableClientState(GL_VERTEX_ARRAY);
    
    // draw fixpoints
    glEnable(GL_BLEND);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_LINE_SMOOTH);
    glPointSize(8.0);
    glLineWidth(1.0);
    for (std::vector<fixpoint>::iterator fi = fpts.begin(); fi != fpts.end(); ++fi) {
        if (fi->saddle) {
            glColor3f(1, 0, 0);
            
            vec2 p00 = fi->pos + 0.01 * fi->evec[0] / norm(fi->evec[0]);
            vec2 p01 = fi->pos - 0.01 * fi->evec[0] / norm(fi->evec[0]);
            vec2 p10 = fi->pos + 0.01 * fi->evec[1] / norm(fi->evec[1]);
            vec2 p11 = fi->pos - 0.01 * fi->evec[1] / norm(fi->evec[1]);
            
            glBegin(GL_LINES);
            glVertex2dv((GLdouble*)&p00);
            glVertex2dv((GLdouble*)&p01);
            glVertex2dv((GLdouble*)&p10);
            glVertex2dv((GLdouble*)&p11);
            glEnd();
        } else {
            glColor3f(0, 0, 1);
        }
        
        glBegin(GL_POINTS);
        glVertex2dv((GLdouble*)&fi->pos);
        glEnd();
    }
#endif
    
    glutSwapBuffers();
}

// --------------------------------------------------------------------------------

void reshape(int w, int h)
{
    glViewport(0, 0, (GLsizei)w, (GLsizei)h);
    
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(_bounds.min()[0], _bounds.max()[0],
            _bounds.min()[1], _bounds.max()[1],
            0, 1);
            
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}

// --------------------------------------------------------------------------------

void motion(int button, int state, int x, int y)
{
    // if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
    //  GLint vp[4];
    //  glGetIntegerv(GL_VIEWPORT, vp);
    //
    //  GLdouble mv[16], pr[16];
    //
    //  glGetDoublev(GL_MODELVIEW_MATRIX, mv);
    //  glGetDoublev(GL_PROJECTION_MATRIX, pr);
    //
    //  double tmp;
    //  gluUnProject(x, vp[3] - (y - vp[1]), 0.0, mv, pr, vp, &cursor[0], &cursor[1], &tmp);
    //
    //  glutPostRedisplay();
    // }
}

// --------------------------------------------------------------------------------

void keyboard(unsigned char key, int, int)
{
    std::cout << key << '\n';
    
#if 0
    switch (key) {
        case 'p':
            poincare_plot(cursor);
            break;
        case 'f':
            find_fixpoint(cursor);
            break;
        case 's':
            find_seedpoints();
            break;
        case 'i':
            for (int i = 0; i < spts.size(); ++i) {
                cursor = spts[i].pos;
                find_fixpoint(cursor);
            }
            break;
        case 'C':
            lpts.clear();
            fpts.clear();
        case 'c':
            ppts.clear();
            npts.clear();
            spts.clear();
            break;
        case 'k':
            ++K;
            std::cout << "K=" << K << std::endl;
            break;
        case 'j':
            if (K > 1) {
                --K;
                std::cout << "K=" << K << std::endl;
            }
            break;
    }
#endif
    
    glutPostRedisplay();
}

// --------------------------------------------------------------------------------

int main(int argc, char** argv)
{
    initialize(argc, argv);
    
    glutInit(&argc, argv);
    glutInitDisplayString("samples rgba double alpha");
    glutInitWindowSize(768, 768);
    glutInitWindowPosition(1400 - 768, 20);
    glutCreateWindow(argv[0]);
    
    init();
    
    glutReshapeFunc(reshape);
    glutDisplayFunc(display);
    glutKeyboardFunc(keyboard);
    glutMouseFunc(motion);
    glutMainLoop();
    return 0;
}












































































































