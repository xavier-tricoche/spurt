#ifndef __CREASE_HPP__
#define __CREASE_HPP__

#include <vector>
#include <math/fixed_vector.hpp>

namespace spurt {
typedef fixed_vector<double, 6>   vec6;
typedef fixed_vector<double, 4>   vec4;
}

namespace spurt {
class MeasureWrapper;

namespace crease {

// visual debugging
extern std::vector< std::vector< vec3 > > vertices;
extern std::vector< vec3 > current_vertices;
extern std::vector< vec3 > problematic_voxels;
extern std::vector< vec3 > fixed_voxels;
extern std::vector< vec3 > ok_faces;
extern std::vector< unsigned int > added_vertices;
extern std::vector< vec3 > crossing_faces;
extern std::vector< vec3 > pvo_faces;
extern std::vector< vec3 > show_all;
extern std::vector< vec3 > pvo_triangles;
extern std::vector< vec3 > round1, round2, round12;

typedef std::vector< vec3 > path;
extern std::vector< path > paths;
extern std::vector< vec3 > intermediate_steps;
extern std::string flag_file_name;

extern bool bold_move;

extern double value_threshold;
extern double value_threshold_select;
extern double strength_threshold;
extern double strength_threshold_select;
extern double confidence_threshold;
extern double gradient_eps;
extern double gradient_eps_rel;
extern unsigned int max_depth;
extern unsigned int max_depth_fix;
extern double max_int_error;
extern double max_align_error;
extern bool is_ridge;
extern int crease_kind;
extern unsigned int upsample;
extern unsigned int failed_conv;
extern bool apply_filter;
extern unsigned int nb_pvo;
extern bool speedup;
extern bool fixing_voxel;
extern bool display_debug_info;
extern bool read_info;

extern unsigned int nb_crossings;

extern spurt::MeasureWrapper* the_wrapper;

inline bool is_ok(double val, double str)
{
    if (is_ridge) {
        return (val > value_threshold_select &&
                str < strength_threshold_select);
    } else {
        return (val < value_threshold_select &&
                str > strength_threshold_select);
    }
}

// matrix-vector product between 3D symmetric matrix represented as
// 6D vector and 3D vector
inline vec3 prod(const vec6& H, const vec3& g)
{
    vec3 Hg;
    Hg[0] = H[0] * g[0] + H[1] * g[1] + H[2] * g[2];
    Hg[1] = H[1] * g[0] + H[3] * g[1] + H[4] * g[2];
    Hg[2] = H[2] * g[0] + H[4] * g[1] + H[5] * g[2];
    return Hg;
}

// matrix-vector product between 2D matrix represented as
// 4D vector and 2D vector
inline vec2 prod(const vec4& H, const vec2& g)
{
    vec2 Hg;
    Hg[0] = H[0] * g[0] + H[1] * g[1];
    Hg[1] = H[2] * g[0] + H[3] * g[1];
    return Hg;
}

}
}


#endif

