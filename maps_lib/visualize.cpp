#include "adaptive_map_sampling.hpp"

#include <iostream>
#include <list>

#include <glut.h>
#include <math/fixed_vector.hpp>
#include <vector>
#include <tokamak/map2d.hpp>
#include <tokamak/poincare_map.hpp>
#include "display.hpp"
#include <tokamak/tokamak_nimrod_parametric.hpp>
#include "definitions.hpp"

#include <util/wall_timer.hpp>

#if _OPENMP
#include <omp.h>
#endif

size_t m = 0, n = 0;

using namespace map_analysis;
using namespace map_display;

nvis::bbox2 _bounds;
metric_type metric;

typedef std::list< nvis::vec2 >     orbit_type;

// graphics stuff
std::vector<nvis::vec3>     colors;     // random colors
std::list<orbit_type>       pplot;      // poincare plot
triangle_list               mesh;       // sampling triangulation
std::vector<triangle_list>  p_mesh;     // period-specific triangulation
std::list<triangle_list>    singular;   // singular triangles
std::vector<line_list>      p_vectors;  // period-specific vector field

// --------------------------------------------------------------------------

void display();

template<typename M>
void poincare_plot(std::list<orbit_type>& points, const M& pmap,
                   size_t nseeds, size_t niter)
{
    typedef M       map_type;
    
    size_t nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    std::vector<std::list< orbit_type> > orbits(nbthreads);
    
    srand48(time(0));
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for (int i = 0 ; i < nseeds ; ++i) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            if (!thread_id) {
                std::cerr << i << "/" << nseeds << "          \r";
            }
            
            nvis::vec2 loc(drand48(), drand48());
            nvis::vec2 x0 = _bounds.min() + loc * _bounds.size();
            
            std::vector<nvis::vec2> steps;
            map_type* my_map = pmap.clone();
            try {
                my_map->map(x0, steps, niter);
            } catch (...) {
                continue;
            }
            
            orbits[thread_id].push_back(orbit_type());
            orbit_type& orb = orbits[thread_id].back();
            for (unsigned int n = 0 ; n < steps.size() ; ++n) {
                orb.push_back(metric.modulo(steps[n]));
            }
        }
        std::cerr << '\n';
        
        for (int i = 0 ; i < orbits.size() ; ++i) {
            std::copy(orbits[i].begin(), orbits[i].end(), std::back_inserter(points));
        }
    }
}

// --------------------------------------------------------------------------------

static void init()
{
    adaptive_map_sampling_params params;
    tokamak_nimrod_parametric* field;
    field = new tokamak_nimrod_parametric(std::string(params.hdf_name), params.time_step);
    field->periodic_coordinates(false);
    bool per[2] = {true, false};
    metric = metric_type(field->bounds(), per);
    
    static_data::metric = metric;
    
    double h = params.eps;
    
    poincare_map pmap(field);
    pmap.precision(params.eps);
    
    _bounds = field->bounds();
    if (!n) {
        n = params.n;
    }
    if (!m) {
        m = params.m;
    }
    poincare_plot(pplot, pmap, n, m);
    
    // initialize OpenGL
    
    colors.resize(256);
    for (int i = 0 ; i < colors.size() ; ++i) {
        colors[i] = nvis::vec3(drand48(), drand48(), drand48());
    }
    
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
    std::list<orbit_type>::iterator it;
    unsigned int count = 0;
    glDisable(GL_POINT_SMOOTH);
    for (it = pplot.begin() ; it != pplot.end() ; ++it, ++count) {
        glDisable(GL_BLEND);
        glPointSize(0.5);
        nvis::vec3 c = colors[count%256];
        glColor3f(c[0], c[1], c[2]);
        glVertexPointer(2, GL_DOUBLE, 0, &it->front());
        glDrawArrays(GL_POINTS, 0, it->size());
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
    if (argc >= 2) {
        n = atoi(argv[1]);
        if (argc == 3) {
            m = atoi(argv[2]);
        }
    }
    
    glutInit(&argc, argv);
    glutInitDisplayString("samples rgba double alpha");
    glutInitWindowSize(1200, 800);
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























































