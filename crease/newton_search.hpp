#ifndef __NEWTON_SEARCH_HPP__
#define __NEWTON_SEARCH_HPP__

#include <math/types.hpp>
#include <face.hpp>

namespace spurt::crease {
inline vec2 newton_step_2d(const vec2& g, const vec3& h)
{
    // solve 2x2 system: Hdx = -g
    double denom = h[0] * h[2] - h[1] * h[1];
    double du = (-g[0] * h[2] + g[1] * h[1]) / denom;
    double dv = (-h[0] * g[1] + h[1] * g[0]) / denom;
    return vec2(du, dv);
}

bool newton(const face_type& face, vec3& guess, double eps, 
            unsigned int maxit = 10)
{
    face.set_basis();
    const vec3& e0 = face.e0;
    const vec3& e1 = face.e1;
    
    vec2 x = local_coordinates(face, guess);
    vec2 g = project_on_face(measure.gradient(guess), face);
    vec3 h = project_on_face(measure.hessian(guess), face);
    
    double gm = norm(g);
    std::cout << "newton init: gm=" << gm << " at " << x << std::endl;
    
    for (unsigned int i = 0 ; i < maxit ; i++) {
        vec2 dx = newton_step_2d(g, h);
        x += dx;
        g = project_on_face(spurt::crease::the_wrapper->gradient(global_coord(face, x)), face);
        h = project_on_face(spurt::crease::the_wrapper->hessian(global_coord(face, x)), face);
        gm = norm(g);
        std::cout << "newton #" << i << ": gm=" << gm << " at " << x << std::endl;
        
        if (gm < eps) {
            guess = global_coord(face, x);
            return true;
        }
    }
    
    return false;
}

struct search_face_newton {
    int operator()(std::vector< vec3 >& xing,
                   const vec3& p0, const vec3& p1,
                   const vec3& p2, const vec3& p3,
                   unsigned int maxdepth, bool& something_found) const {
        face_type face(p0, p1, p2, p3);
        vec3 guess = 0.5 * (p0 + p3);
        if (newton(face, guess, spurt::crease::gradient_eps, 20)) {
            xing.push_back(guess);
            return 1;
        }
        
        return 0;
        
    }
};

}


#endif









