#include "face.hpp"


namespace spurt {
namespace crease {

vec3 average(vec3 dirs[4]);
vec3 average(const std::vector< vec3 >& dirs);


face_type::facetype(const vec3& p0, const vec3& p1,
                    const vec3& p2, const vec3& p3) {
    p[0] = p0;
    p[1] = p1;
    p[2] = p2;
    p[3] = p3;
    
    e0 = p[1] - p[0];
    e1 = p[3] - p[0];
    l0 = norm(e0);
    l1 = norm(e1);
    
    for (unsigned int i = 0 ; i < 4 ; i++) {
        g[i] = spurt::crease::the_wrapper->gradient(p[i]);
        Hg[i] = spurt::crease::the_wrapper->Hgradient(p[i]);
    }
    
    basis_set = true;
    depth = 0;
    reference.clear();
}

vec3 gradient(unsigned int pid, const face_type& face, bool interpolate)
{
    vec3 g;
    if (pid < 4) {
        return face.g[pid];
    } else if (interpolate) {
        face(coef[pid][0], coef[pid][1], face.g);
    } else {
        vec3 p = position(pid, face);
        g = spurt::crease::the_wrapper->gradient(p);
    }
    
    return g;
}

vec3 Hgradient(unsigned int pid, const face_type& face, bool interpolate)
{
    vec3 Hg;
    if (pid < 4) {
        return face.Hg[pid];
    } else if (interpolate) {
        face(coef[pid][0], coef[pid][1], face.Hg);
    } else {
        vec3 p = position(pid, face);
        Hg = spurt::crease::the_wrapper->Hgradient(p);
    }
    
    return Hg;
}

void refine_face(face_type& out, const face_type& face,
                 const vec3& q, double h)
{
    // compute local coordinates
    vec2 x = local_coord(face, q);
    double xmin = std::max(0., x[0] - h);
    double xmax = std::min(1., x[0] + h);
    double ymin = std::max(0., x[1] - h);
    double ymax = std::min(1., x[1] + h);
    out = face_type(face(xmin, ymin), face(xmax, ymin),
                    face(xmax, ymax), face(xmin, ymax));
    std::copy(face.reference.begin(), face.reference.end(), std::back_inserter(out.reference));
    out.reference.push_back(8);
    
    out.depth = face.depth;
    // we do not increment the depth when we refine around a PVO solution
}


} // namespace spurt
} // namespace crease

#endif













