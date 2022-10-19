#include <iostream>
#include <map>
#include <list>
#include <vector>

#include <math/fixed_vector.hpp>

#include "maps_lib/adaptive_map_sampling.hpp"
#include "maps_lib/definitions.hpp"
#include "maps_lib/orbits.hpp"

#include <teem/hest_helper.hpp>

// GLUT business
#include "opengl/Camera/CameraWrapper_Misc.h"

using namespace map_analysis;

typedef nvis::fvec3     color_type;

const color_type rainbow[] = {
    color_type(0, 0, 0), color_type(0, 0, 0.5), color_type(0, 0, 1), color_type(0, 0.5, 1), color_type(0, 1, 1),
    color_type(0, 1, 0.5), color_type(0, 1, 0), color_type(0.5, 1, 0), color_type(1, 1, 0), color_type(1, 0.5, 0),
    color_type(1, 0, 0), color_type(1, 0, 0.5), color_type(1, 0, 1), color_type(1, 0.5, 1), color_type(1, 1, 1)
};

template<typename T>
struct adaptive_color_map {
    adaptive_color_map(const std::vector<T>& vals, const std::vector<color_type>& scale)
        : t(scale.size()), colors(scale) {
        std::vector<T> tmp(vals);
        std::sort(tmp.begin(), tmp.end());
        unsigned int step = tmp.size() / (scale.size() - 1);
        for (int i = 0 ; i < t.size() ; ++i) {
            t[i] = tmp[i*step];
        }
        t.back() = tmp.back();
    }
    
    color_type operator()(const T& val) const {
        unsigned int bin = std::distance(t.begin(),
                                         std::lower_bound(t.begin(), t.end(), val));
        if (bin > t.size() - 2) {
            bin = t.size() - 2;
        }
        T u = (val - t[bin]) / (t[bin+1] - t[bin]);
        return (1. - u)*colors[bin] + u*colors[bin+1];
    }
    
    std::vector<T>          t;
    std::vector<color_type> colors;
};

// computation parameters
int n_samples, n_iter, max_tri, max_per;
char* fout, *fin, *ts;
double eps, dq, min_area, max_area, max_ratio, max_err;
void init(int argc, char* argv[])
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
    hestOptAdd(&hopt, "f",      "file",                 airTypeString,  1, 1, &fin,         NULL,       "input hdf5 file name");
    hestOptAdd(&hopt, "t",      "time",                 airTypeString,  1, 1, &ts,          NULL,       "time step string");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &fout,        NULL,       "output nrrd file name");
    hestOptAdd(&hopt, "eps",    "eps",                  airTypeDouble,  0, 1, &eps,         "1.0e-6",   "integration precision");
    hestOptAdd(&hopt, "n",      "# samples",            airTypeInt,     0, 1, &n_samples,   "5000",     "number of 1D samples");
    hestOptAdd(&hopt, "m",      "# rotations",          airTypeInt,     0, 1, &n_iter,      "50",       "number of rotations");
    hestOptAdd(&hopt, "dq",     "safety factor step",   airTypeDouble,  0, 1, &dq,          "0.05",     "maximum discrepancy between consecutive safety factors");
    hestOptAdd(&hopt, "mr",     "max aspect ratio",     airTypeDouble,  0, 1, &max_ratio,   "1.",       "maximum triangle aspect ratio");
    hestOptAdd(&hopt, "mt",     "max # triangles",      airTypeInt,     0, 1, &max_tri,     "150000",   "max number of triangles in adaptive sampling");
    hestOptAdd(&hopt, "ma",     "min triangle area",    airTypeDouble,  0, 1, &min_area,    "0.01",     "min triangle area in adaptive sampling");
    hestOptAdd(&hopt, "Ma",     "max triangle area",    airTypeDouble,  0, 1, &max_area,    "1.",       "max triangle area in adaptive sampling");
    hestOptAdd(&hopt, "mp",     "max period",           airTypeInt,     0, 1, &max_per,     "15",       "max considered period in fixed point search");
    hestOptAdd(&hopt, "err",    "max tolerance",        airTypeDouble,  0, 1, &max_err,     "0.1",      "max tolerance of approximation quality criterion");
    
    __hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                     me, "Adaptively sample phase portrait along orbits to achieve accurate piecewise linear interpolation on a per-period basis",
                     AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

template<typename T>
GLuint draw_poincare_plot(const T& color_map)
{
    GLuint plot_id = glGenLists(1);
    glNewList(plot_id, GL_COMPILE);
    glEnable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(1);
    glBegin(GL_POINTS);
    for (int i = 0 ; i < static_data::central_map_orbits.size() ; ++i) {
        const orbit& orb = static_data::central_map_orbits[i];
        color_type col = color_map(orb.period());
        glColor3f(col[0], col[1], col[2]);
        for (int j = 0 ; j < orb.size() ; ++j) {
            nvis::vec2 x = static_data::metric.modulo(orb[j]);
            glVertex2f(x[0], x[1]);
        }
    }
    glEnd();
    glEndList();
    return plot_id;
}

struct Edge {
    Edge() : i0(0), i1(0) {}
    Edge(size_t _i0, size_t _i1) : i0(_i0), i1(_i1) {
        if (i0 > i1) {
            size_t tmp = i0;
            i0 = i1;
            i1 = tmp;
        }
    }
    Edge(const Edge& e) : i0(e.i0), i1(e.i1) {}
    
    bool operator<(const Edge& e) const {
        if (i0 < e.i0) {
            return true;
        } else if (i0 > e.i0) {
            return false;
        }
        return i1 < e.i1;
    }
    
    size_t i0, i1;
};

template<typename T, typename I>
GLuint draw_mesh(const T& mesh, const color_type& color, const I& included)
{
    typedef typename T::triangle_type   triangle_type;
    
    // uniquify edges
    std::set<Edge> edges;
    for (int i = 0 ; i < mesh.get_nb_triangles() ; ++i) {
        if (!included(mesh, i)) {
            continue;
        }
        const triangle_type& tri = mesh.get_triangle_vertices(i);
        edges.insert(Edge(tri[0], tri[1]));
        edges.insert(Edge(tri[1], tri[2]));
        edges.insert(Edge(tri[2], tri[0]));
    }
    
    GLuint mesh_id = glGenLists(1);
    glNewList(mesh_id, GL_COMPILE);
    glBegin(GL_LINES);
    glColor3f(color[0], color[1], color[2]);
    for (std::set<Edge>::const_iterator it = edges.begin() ; it != edges.end() ; ++it) {
        const nvis::vec2& p0 = mesh.get_vertex(it->i0);
        const nvis::vec2& p1 = mesh.get_vertex(it->i1);
        glVertex2f(p0[0], p0[1]);
        glVertex2f(p1[0], p1[1]);
    }
    glEnd();
    glEndList();
    return mesh_id;
}

template<typename T>
struct default_mesh_functor {
    bool operator()(const T&, unsigned int) const {
        return true;
    }
};

template<typename T>
struct selection_mesh_functor {
    selection_mesh_functor(const std::list<unsigned int>& included)
        : _included(included.begin(), included.end()) {}
        
    bool operator()(const T&, unsigned int tri) const {
        return _included.find(tri) != _included.end();
    }
    
    std::set<unsigned int> _included;
};

template<typename T>
GLuint draw_mesh(const T& mesh, const color_type& color)
{
    return draw_mesh(mesh, color, default_mesh_functor<T>());
}

inline void draw_arrow(const nvis::vec2& p0, const nvis::vec2& p1)
{
    const double cos_alpha = cos(M_PI / 6.);
    const double sin_alpha = sin(M_PI / 6.);
    
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

template<typename T, typename F>
GLuint draw_hedgehog(const T& mesh, const color_type& color, const F& functor, double scaling)
{
    typedef typename T::data_type   data_type;
    
    size_t nb_pts = mesh.get_nb_vertices();
    
    GLuint hedgehog_id = glGenLists(1);
    glNewList(hedgehog_id, GL_COMPILE);
    glColor3f(color[0], color[1], color[2]);
    glBegin(GL_LINES);
    for (int i = 0 ; i < nb_pts ; ++i) {
        const data_type& d = mesh.get_data(i);
        nvis::vec2 v = functor(d);
        if (nvis::norm(v) == 0) {
            continue;
        }
        v *= scaling;
        const nvis::vec2& p0 = mesh.get_vertex(i);
        nvis::vec2 p1 = p0 + v;
        draw_arrow(p0, p1);
    }
    glEnd();
    glEndList();
    return hedgehog_id;
}

template<typename T>
struct p_map_functor {
    typedef typename T::data_type data_type;
    
    p_map_functor(unsigned int period, double dq)
        : _p(period), _dq(dq) {
        for (unsigned int d = 1 ; d <= _p ; ++d) {
            rational_type q(_p, d);
            if (q.numerator() != _p) {
                continue;
            } else {
                _qs.push_back(xavier::value(q));
            }
        }
    }
    
    nvis::vec2 operator()(const data_type& d) const {
        double q = d.period();
        std::vector<double> dist(_qs.size());
        for (int i = 0 ; i < _qs.size() ; ++i) {
            dist[i] = fabs(_qs[i] - q);
        }
        if (*std::min_element(dist.begin(), dist.end()) > _dq) {
            return nvis::vec2(0, 0);
        } else {
            return map_analysis::vector_value(d, _p, static_data::metric);
        }
    }
    
    std::vector<double> _qs;
    double              _dq;
    unsigned int        _p;
};

// GLUT code

// static variables
CameraWrapper       camera;
float               screen_ratio;
int                 width, height;
nvis::bbox2         box;
std::vector<GLuint> list_ids;

// callback functions
void display()
{
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);
    
    camera.setProjectionMatrix();
    camera.setModelviewMatrix();
    
    for (int i = 0 ; i < list_ids.size() ; ++i) {
        glCallList(list_ids[i]);
    }
    
    camera.markUnchanged();
    
    glutSwapBuffers();
}

void resetCamera()
{
    nvis::vec2 dim2D = box.size();
    nvis::vec2 center2D = box.center();
    
    Vector3d dim(dim2D[0], dim2D[1], min(dim2D[0], dim2D[1]));
    Vector3d center(center2D[0], center2D[1], 0);
    
    camera.scaleCameraToObject(center, dim, Vector3f(0, 1, 0), 90, 1);
    
    //camera.setLeftMouse_RotateAboutLookAtCenter();
    camera.setMiddleMouse_Pan();
    camera.setRightMouse_Dolly();
    
    // I should have a function that turns off the rotate complete, but this little hack does the trick
    camera.setRotateAboutLookAtCenterSensitivityX(0);
    camera.setRotateAboutLookAtCenterSensitivityY(0);
    
    // use functions like below to control sensitivity of the other movements
    //  though camera.scaleCameraToObject(...) does a pretty good job of figuring out a sensitivity for you
    //camera.setPanSensitivityX(10000);
    
    camera.setViewport(0, 0, width, height);
    
    glutPostRedisplay();
}

void reshape(int w, int h)
{
    width = w;
    height = h;
    
    if (h == 0) {
        h = 1;
    }
    
    if (screen_ratio > ((float)w) / h) {
        glViewport(0, 0, w, w / screen_ratio);
    } else if (screen_ratio < ((float)w) / h) {
        glViewport(0, 0, h*screen_ratio, h);
    } else {
        glViewport(0, 0, w, h);
    }
    
    CameraWrapper_GUI::glutResizeCallback(&camera, w, h);
    glutPostRedisplay();
}

void mouse(int button, int state, int x, int y)
{
    CameraWrapper_GUI::glutMouseCallback(&camera, button, state, x, y);
    
    // example of getting actual coordinates of pixel
    //   when third parameter is false, the function uses gluProject to determine the z "pixel" value
    //   when it is true, it reads the depth buffer to get the z value
    Vector3f pos = camera.unProject(x, y, false);
    printf("World Space x/y = %f %f\n", pos.x, pos.y);
    
    // camera.project(...): projects a world position onto the image plane
    
    glutPostRedisplay();
}

void motion(int x, int y)
{
    CameraWrapper_GUI::glutMouseMotionCallback(&camera, x, y);
    
    glutPostRedisplay();
}

void keyboard(unsigned char key, int x, int y)
{
    if (key == 'r') {
        resetCamera();
    }
    
    glutPostRedisplay();
}

int main(int argc, char* argv[])
{
    typedef mesh_type::data_type    data_type;
    
    init(argc, argv);
    
    adaptive_map_sampling_params params;
    params.hdf_name     = fin;
    params.time_step    = ts;
    params.out_name     = fout;
    params.eps          = eps;
    params.dq           = dq;
    params.mr           = max_ratio;
    params.ma           = min_area;
    params.Ma           = max_area;
    params.err          = max_err;
    params.n            = n_samples;
    params.m            = n_iter;
    params.mt           = max_tri;
    params.mp           = max_per;
    
    adaptive_map_sampling_output output;
    map_analysis::adaptive_map_sampling(output, params);
    
    box = output.base_mesh.get_bounds();
    nvis::vec2 total_diag = box.size();
    
    // loop over orbits to determine q-range
    std::vector<double> periods;
    for (int i = 0 ; i < static_data::central_map_orbits.size() ; ++i) {
        const orbit& orb = static_data::central_map_orbits[i];
        periods.push_back(orb.period());
    }
    // create color map
    std::vector<color_type> scale(rainbow, &rainbow[15]);
    adaptive_color_map<double> col_map(periods, scale);
    
    width = 1200;
    screen_ratio = total_diag[0] / total_diag[1];
    height = (int)((float)width / screen_ratio);
    std::cerr << "height = " << height << '\n';
    
    // initialize GLUT
    glutInit(&argc, argv);
    glutInitDisplayString("samples rgba double alpha");
    // glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(width, height);
    glutInitWindowPosition(20, 20);
    glutCreateWindow(argv[0]);
    
    // draw poincare plot
    GLuint plot_id = draw_poincare_plot(col_map);
    list_ids.push_back(plot_id);
    // draw base mesh
    list_ids.push_back(draw_mesh(output.base_mesh, color_type(0.1, 0.1, 0.1)));
    // draw last p-mesh
    // list_ids.push_back(draw_mesh(output.p_meshes.back(), color_type(0.6, 0.6, 0)));
    // draw candidate triangles
    for (int i = 0 ; i < output.p_sing_tris.size() ; ++i) {
        color_type c(0, 0, 1);
        list_ids.push_back(draw_mesh(output.p_meshes[i], c,
                                     selection_mesh_functor<mesh_type>(output.p_cand_tris[i])));
    }
    // draw singular triangles for each period
    for (int i = 0 ; i < output.p_sing_tris.size() ; ++i) {
        color_type c(1, 0, 0);
        list_ids.push_back(draw_mesh(output.p_meshes[i], c,
                                     selection_mesh_functor<mesh_type>(output.p_sing_tris[i])));
    }
    // draw hedgehogs
    for (int i = 0 ; i < output.p_meshes.size() ; ++i) {
        color_type c(drand48(), drand48(), drand48());
        p_map_functor<mesh_type> pmap(i + 1, params.dq);
        list_ids.push_back(draw_hedgehog(output.p_meshes[i], c, pmap, 0.05));
    }
    
    std::cerr << "list ids are:\n";
    for (int i = 0 ; i < list_ids.size() ; ++i) {
        std::cerr << list_ids[i] << " ";
    }
    
    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutMouseFunc(mouse);
    glutMotionFunc(motion);
    
    resetCamera();
    
    /* Enter GLUT event loop */
    glutMainLoop();
    
    return 0;
}




































































































