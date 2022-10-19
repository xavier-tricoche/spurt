#ifndef __GLUT_HELPER_HPP__
#define __GLUT_HELPER_HPP__

#include <iostream>
#include <map>
#include <list>
#include <vector>

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <GL/glut.h>
#include "Camera/CameraWrapper_Misc.h"

namespace GLUT_helper {

typedef spurt::fvec3             color_type;
typedef spurt::vec2              point2d;
typedef spurt::vec3              point3d;

static float                screen_ratio;
static int                  width, height;
static nvis::bbox2          box;
static CameraWrapper        camera;

// API to set user defined callback functions
extern void set_my_display(void (*f)(void));      /* drawing instructions*/
extern void set_my_idle(void (*f)(void));         /* idle instructions */

inline double max_dim(const nvis::bbox2& abox)
{
    spurt::vec2 span = abox.size();
    return std::max(span[0], span[1]);
}

inline nvis::bbox2 current_bounds()
{
    nvis::bbox2 cur_box;
    cur_box.min()[0] = camera.getLeftClipPlane();
    cur_box.max()[0] = camera.getRightClipPlane();
    cur_box.min()[1] = camera.getBottomClipPlane();
    cur_box.max()[1] = camera.getTopClipPlane();
}

inline double current_size()
{
    return max_dim(current_bounds());
}

inline void update_panning_sensitivity(float s = 2.)
{
    nvis::bbox2 cur_box = current_bounds();
    
    float sensitivityScale = max_dim(cur_box) * s;
    
    if (sensitivityScale == 0) {
        // std::cerr << "invalid sensitivity value. unchanged\n";
        return;
    }
    
//  std::cerr << "sensitivity updated to " << sensitivityScale << std::endl;
    camera.setPanSensitivityX(sensitivityScale);
    camera.setPanSensitivityY(sensitivityScale);
}

inline spurt::vec3 world_coordinates(int x, int y)
{
    camera.getActualPixelCoordinates(x,y);
    Vector3f v = camera.unProject(x,y);
    return spurt::vec3(v(0), v(1), v(2));
}

inline void resetCamera()
{
#if 0
    point2d dim2D = box.size();
    point2d center2D = box.center();
    
    Vector3d dim(dim2D[0], dim2D[1], min(dim2D[0], dim2D[1]));
    Vector3d center(center2D[0], center2D[1], 0);
    
    camera.scaleCameraToObject(center, dim, Vector3f(0, 1, 0), 90, dim2D[0] / dim2D[1]);
    
    float maxDim = max(dim2D[0], dim2D[1]);
    // camera.setOrthographic(-maxDim, maxDim, -maxDim, maxDim, camera.getNearZ(),
    //                        camera.getFarZ());
#endif
    
    // std::cerr << "znear = " << camera.getNearZ() << ", zfar = " << camera.getFarZ()
    //  << std::endl;
    
    camera.setOrthographic(box.min()[0], box.max()[0], box.min()[1], box.max()[1],
                           camera.getNearZ(), camera.getFarZ());
                           
                           
    // I should have a function that turns off the rotate complete, but this little hack does the trick
    camera.setRotateAboutLookAtCenterSensitivityX(0);
    camera.setRotateAboutLookAtCenterSensitivityY(0);
    
    // use functions like below to control sensitivity of the other movements
    //  though camera.scaleCameraToObject(...) does a pretty good job of figuring out a sensitivity for you
    //camera.setPanSensitivityX(10000);
    
    // float sensitivityScale = max_dim(box) * 2.0;
    //
    // camera.setDollySensitivity(sensitivityScale);
    // camera.setPanSensitivityX(sensitivityScale);
    // camera.setPanSensitivityY(sensitivityScale);
    
    camera.setViewport(0, 0, width, height);
    
    glutPostRedisplay();
}

inline void setSensitivity(double s)
{
    camera.setDollySensitivity(s);
    camera.setPanSensitivityX(s);
    camera.setPanSensitivityY(s);
    camera.setMiddleMouse_Pan();
}

// setup_display takes as input pointer to actual drawing instructions function
inline void setup_display(void (*f)(void))
{
    camera.setProjectionMatrix();
    camera.setModelviewMatrix();
    
    f();
    
    camera.markUnchanged();
}

// following functions must be static since their address will be passed to GLUT
static void glut_helper_reshape(int w, int h)
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
    
    // CameraWrapper_GUI::glutResizeCallback(&camera, w, h);
    
    update_panning_sensitivity();
    glutPostRedisplay();
}

static void glut_helper_mouse(int button, int state, int x, int y)
{
    // std::cerr << "mouse button #" << button << " is " << (state==0 ? "DOWN" : "UP") << std::endl;
    
    CameraWrapper_GUI::glutMouseCallback(&camera, button, state, x, y);
    
    // example of getting actual coordinates of pixel
    //   when third parameter is false, the function uses gluProject to determine the z "pixel" value
    //   when it is true, it reads the depth buffer to get the z value
    Vector3f pos = camera.unProject(x, y, false);
    
    // camera.project(...): projects a world position onto the image plane
    
    glutPostRedisplay();
}

static void glut_helper_motion(int x, int y)
{
    CameraWrapper_GUI::glutMouseMotionCallback(&camera, x, y);
    glutPostRedisplay();
}

static void glut_helper_keyboard(unsigned char key, int x, int y)
{
    if (key == 'r') {
        resetCamera();
    }
    
    glutPostRedisplay();
}

inline void __draw_vector(const spurt::vec2& x, const spurt::vec2& v)
{
    const double cos_alpha = cos(M_PI/12.);
    const double sin_alpha = sin(M_PI/12.);
    
    if (nvis::norm(v) == 0) {
        return;
    }
    spurt::vec2 y = x + v;
    glVertex2f(x[0], x[1]);
    glVertex2f(y[0], y[1]);
    spurt::vec2 e0 = -0.2 * v;
    spurt::vec2 e1(-e0[1], e0[0]);
    glVertex2f(y[0], y[1]);
    spurt::vec2 z = y + cos_alpha * e0 + sin_alpha * e1;
    glVertex2f(z[0], z[1]);
    glVertex2f(y[0], y[1]);
    z = y + cos_alpha * e0 - sin_alpha * e1;
    glVertex2f(z[0], z[1]);
}

inline void draw_vector(const spurt::vec2& x, const spurt::vec2& v, const spurt::fvec3& col, float width=1)
{
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(width);
    glColor3f(col[0], col[1], col[2]);
    glBegin(GL_LINES);
    __draw_vector(x, v);
    glEnd();
}

inline void draw_vectors(const std::vector<spurt::vec2>& xs,
                         const std::vector<spurt::vec2>& vs, const spurt::fvec3& col, float width=1)
{

    assert(xs.size() == vs.size());
    
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(width);
    glColor3f(col[0], col[1], col[2]);
    glBegin(GL_LINES);
    for (int i=0 ; i<xs.size() ; ++i) {
        const spurt::vec2& x = xs[i];
        const spurt::vec2& v = vs[i];
        __draw_vector(x, v);
    }
    glEnd();
}

inline void draw_vectors(const std::vector<std::pair<spurt::vec2, spurt::vec2> >& ps, const spurt::fvec3& col, float width=1)
{
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glLineWidth(width);
    glColor3f(col[0], col[1], col[2]);
    glBegin(GL_LINES);
    for (int i=0 ; i<ps.size() ; ++i) {
        const spurt::vec2& x = ps[i].first;
        const spurt::vec2& v = ps[i].second;
        __draw_vector(x, v);
    }
    glEnd();
}

inline void draw_curve(const std::vector<spurt::vec2>& xs, const spurt::fvec3& col, float width=1)
{
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glColor3f(col[0], col[1], col[2]);
    glLineWidth(width);
    glBegin(GL_LINE_STRIP);
    for (int i=0 ; i<xs.size() ; ++i) {
        glVertex2f(xs[i][0], xs[i][1]);
    }
    glEnd();
}

inline void draw_quad(const nvis::bbox2& box, const spurt::fvec3& col, float width=1)
{
    glEnable(GL_BLEND);
    glEnable(GL_LINE_SMOOTH);
    glColor3f(col[0], col[1], col[2]);
    glLineWidth(width);
    const spurt::vec2& minp = box.min();
    const spurt::vec2& maxp = box.max();
    glBegin(GL_LINE_STRIP);
    glVertex2f(minp[0], minp[1]);
    glVertex2f(maxp[0], minp[1]);
    glVertex2f(maxp[0], maxp[1]);
    glVertex2f(minp[0], maxp[1]);
    glVertex2f(minp[0], minp[1]);
    glEnd();
}

inline void draw_dots(const std::vector<spurt::vec2>& xs, const spurt::fvec3& col, float sz=1)
{
    glDisable(GL_POINT_SMOOTH);
    glDisable(GL_BLEND);
    glPointSize(sz);
    glColor3f(col[0], col[1], col[2]);
    glBegin(GL_POINTS);
    for (int i=0 ; i<xs.size() ; ++i) {
        spurt::vec2 x = xs[i];
        glVertex2f(x[0], x[1]);
    }
    glEnd();
}

} // GLUT_helper

#endif



























