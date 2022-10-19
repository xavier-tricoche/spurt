#ifndef __XAVIER_FACE_HPP__
#define __XAVIER_FACE_HPP__

#include <math/fixed_vector.hpp>
#include <vector>
#include <iostream>
#include "crease.hpp"


namespace xavier {
namespace crease {

// sample points in local coordinates
// 3 --- 6 --- 2
// |     |     |
// 7 --- 8 --- 5
// |     |     |
// 0 --- 4 --- 1
const double coef[][2] = {
    { 0, 0 },           // 0
    { 1, 0 },           // 1
    { 1, 1 },           // 2
    { 0, 1 },           // 3
    { 0.5, 0 },         // 4
    { 1, 0.5 },         // 5
    { 0.5, 1 },         // 6
    { 0, 0.5 },         // 7
    { 0.5, 0.5 }        // 8
};

// subdivision scheme
const unsigned int indices[][4] = {
    { 0, 4, 8, 7 },
    { 4, 1, 5, 8 },
    { 8, 5, 2, 6 },
    { 7, 8, 6, 3 }
};

nvis::vec3 average(nvis::vec3 dirs[4]);
nvis::vec3 average(const std::vector< nvis::vec3 >& dirs);

struct face_information {
    face_information() : e0(0, 0, 0), e1(0, 0, 0), x0(0), y1(0) {}
    
    void set_info(const nvis::vec3& p0, const nvis::vec3& p1,
                  const nvis::vec3& p2, const nvis::vec3& p3,
                  double eps = 1.0e-3) {
        omega = p0;
        e0 = p1 - p0;
        e1 = p3 - p0;
        x0 = nvis::norm(e0);
        y1 = nvis::norm(e1);
        e0 /= x0;
        e1 /= y1;
        tol = fabs(eps);
        cps.clear();
    }
    
    bool add_crease_point(const nvis::vec3& p) {
        if (inside(p)) {
            cps.push_back(p);
            valid = true;
            return true;
        }
        return false;
    }
    
    bool contains_point() const {
        return cps.size();
    }
    
    bool inside(const nvis::vec3& p) {
        nvis::vec3 q = p - omega;
        double x = nvis::inner(q, e0) / x0;
        double y = nvis::inner(q, e1) / y1;
        return (x > -tol && x < 1 + tol && y > -tol && y < 1 + tol);
    }
    
    nvis::vec3 omega, e0, e1;
    double x0, y1;
    double tol;
    bool valid;
    std::vector< nvis::vec3 > cps;
};

struct face_type {
    face_type() : basis_set(false) {}
    
    face_type(const nvis::vec3& p0, const nvis::vec3& p1,
              const nvis::vec3& p2, const nvis::vec3& p3) {
        p[0] = p0;
        p[1] = p1;
        p[2] = p2;
        p[3] = p3;
        
        e0 = p[1] - p[0];
        e1 = p[3] - p[0];
        l0 = nvis::norm(e0);
        l1 = nvis::norm(e1);
        
        for (unsigned int i = 0 ; i < 4 ; i++) {
            g[i] = xavier::crease::the_wrapper->gradient(p[i]);
            Hg[i] = xavier::crease::the_wrapper->Hgradient(p[i]);
        }
        
        basis_set = true;
        depth = 0;
        reference.clear();
    }
    
    void set_basis() const {
        if (!basis_set) {
            e0 = p[1] - p[0];
            e1 = p[3] - p[0];
            basis_set = true;
        }
    }
    
    nvis::vec3 operator()(double u, double v) const {
        return p[0] + u*e0 + v*e1;
    }
    
    nvis::vec3 operator()(double u, double v, const nvis::vec3* vals) const {
        return (1 - u)*(1 - v)*vals[0] + u*(1 - v)*vals[1] + u*v*vals[2] + (1 - u)*v*vals[3];
    }
    
    double average_norm(const nvis::vec3* vals) const {
        double n=0;
        for (int i=0 ; i<4 ; ++i) {
            n += nvis::norm(vals[i]);
        }
        return 0.25*n;
    }
    
    nvis::vec3 p[4];
    nvis::vec3 g[4];
    nvis::vec3 Hg[4];
    mutable nvis::vec3 e0, e1;
    mutable double l0, l1;
    mutable bool basis_set;
    std::vector< int > reference;
    unsigned int depth;
};

void refine_face(face_type& out, const face_type& face,
                 const nvis::vec3& q, double h = 0.25);
                 
inline nvis::vec3 project_on_face(const nvis::vec6& m3d, const face_type& face)
{
    nvis::vec3 m2d;
    nvis::vec3 m3d_dot_e0 = prod(m3d, face.e0)/face.l0;
    nvis::vec3 m3d_dot_e1 = prod(m3d, face.e1)/face.l1;
    m2d[0] = nvis::inner(face.e0, m3d_dot_e0) / face.l0;
    m2d[1] = nvis::inner(face.e1, m3d_dot_e0) / face.l1;
    m2d[1] = nvis::inner(face.e1, m3d_dot_e1) / face.l1;
    return m2d;
}

inline nvis::vec2 project_on_face(const nvis::vec3& v3d, const face_type& face)
{
    return nvis::vec2(nvis::inner(v3d, face.e0) / (face.l0*face.l0),
                      nvis::inner(v3d, face.e1) / (face.l1*face.l1));
}

std::ostream& operator<<(std::ostream& out, const face_type& face);

inline nvis::vec2 local_coord(const face_type& face, const nvis::vec3& p)
{
    nvis::vec3 q = p - face.p[0];
    double u = nvis::inner(q, face.e0) / face.l0;
    double v = nvis::inner(q, face.e1) / face.l1;
    return nvis::vec2(u, v);
}

inline nvis::vec3 global_coord(const face_type& face, const nvis::vec2& x)
{
    return face.p[0] + x[0]*face.e0 + x[1]*face.e1;
}

inline bool inside(const face_type& face, const nvis::vec3& p)
{
    nvis::vec2 x = local_coord(face, p);
    return (x[0] >= 0 && x[0] <= 1 && x[1] >= 0 && x[1] <= 1);
}

inline nvis::vec3 position(unsigned int pid, const face_type& face)
{
    if (pid < 4) {
        return face.p[pid];
    } else {
        return face(coef[pid][0], coef[pid][1], face.p);
    }
}

inline nvis::vec3 gradient(unsigned int pid, const face_type& face, bool interpolate)
{
    nvis::vec3 g;
    if (pid < 4) {
        return face.g[pid];
    } else if (interpolate) {
        face(coef[pid][0], coef[pid][1], face.g);
    } else {
        nvis::vec3 p = position(pid, face);
        g = xavier::crease::the_wrapper->gradient(p);
    }
    
    return g;
}

inline nvis::vec3 Hgradient(unsigned int pid, const face_type& face, bool interpolate)
{
    nvis::vec3 Hg;
    if (pid < 4) {
        return face.Hg[pid];
    } else if (interpolate) {
        face(coef[pid][0], coef[pid][1], face.Hg);
    } else {
        nvis::vec3 p = position(pid, face);
        Hg = xavier::crease::the_wrapper->Hgradient(p);
    }
    
    return Hg;
}
} // namespace xavier
} // namespace crease

#endif













