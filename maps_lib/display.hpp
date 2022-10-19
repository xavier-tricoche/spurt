#ifndef __DISPLAY_HPP__
#define __DISPLAY_HPP__

#include <vector>
#include <map>
#include <list>
#include <math/fixed_vector.hpp>

namespace map_display {
typedef nvis::vec2          point2d;
typedef nvis::fvec3         color;

struct vertex {
    vertex() : x(0, 0), col(0, 0, 0) {}
    vertex(const nvis::vec2& _x) : x(_x), col(0, 0, 0) {}
    vertex(const nvis::vec2& _x, const nvis::vec3& _col) : x(_x), col(_col) {}
    
    point2d         x;
    color           col;
};

typedef std::list<vertex>           vertex_list;
typedef std::pair<vertex, vertex>   line;
typedef std::list<line>             line_list;
typedef std::list<vertex>           polyline;
typedef std::list<polyline>         polyline_list;

struct triangle {
    triangle() {}
    triangle(const nvis::vec2& x0, const nvis::vec2& x1, const nvis::vec2& x2) {
        v[0] = vertex(x0);
        v[1] = vertex(x1);
        v[2] = vertex(x2);
    }
    triangle(const vertex& v0, const vertex& v1, const vertex& v2) {
        v[0] = v0;
        v[1] = v1;
        v[2] = v2;
    }
    
    vertex v[3];
};

typedef std::list<triangle>         triangle_list;
typedef std::list<vertex>           triangle_strip;

inline line_list arrow(const line& l, double a = 0.25)
{
    const double cos_alpha = cos(M_PI / 6.);
    const double sin_alpha = sin(M_PI / 6.);
    
    line_list _lines;
    _lines.push_back(l);
    color c = l.first.col;
    
    nvis::vec2 e0 = a * (l.second.x - l.first.x);
    nvis::vec2 e1(-e0[1], e0[0]);
    point2d x = l.second.x + cos_alpha * e0 + sin_alpha * e1;
    _lines.push_back(line(l.second, vertex(x, c)));
    x = l.second.x + cos_alpha * e0 - sin_alpha * e1;
    _lines.push_back(line(l.second, vertex(x, c)));
    
    return _lines;
}




}


#endif










