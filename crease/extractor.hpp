#pragma once

#include <vector>
#include <map>

#include "measure_wrapper.hpp"
#include "face.hpp"
#include "grid.hpp"
#include "crease.hpp"
#include <math/types.hpp>

namespace spurt {
namespace crease {

// some helper functions
static void PCA(std::vector< vec3 >& evec, const std::vector< vec3 >& dirs);
static void Furst_orientation(vec3& ev1, const vec3& ev0);
static void Naive_orientation(vec3& ev1, const vec3& ev0);


struct crease_point {
    crease_point() : p(), dot(0), ngm(1) {}

    vec3 p;
    double dot;
    double ngm;
};

struct point_on_face {
    point_on_face() {}

    vec3 p[4];
    vec3 q;
    vec3 e, g, Hg;
};

enum extraction_method { PVO, MC, ZG };

struct extractor 
{   
    // for debugging display
    std::vector< std::vector< vec3 > > vertices;
    std::vector< vec3 > problematic_voxels;
    std::vector< vec3 > fixed_voxels;
    std::vector< vec3 > current_vertices;
    std::vector< vec3 > ok_faces;
    std::vector< unsigned int > added_vertices;
    std::vector< vec3 > crossing_faces;
    std::vector< path > paths;
    std::vector< vec3 > pvo_faces;
    std::vector< vec3 > show_all;
    std::vector< vec3 > intermediate_steps;
    std::vector< vec3 > round1, round2, round12;
    std::string flag_file_name;
    bool bold_move;

    // spatial accuracy
    unsigned int max_depth = 4;
    unsigned int max_depth_fix = 6;
    unsigned int upsample;
    double max_int_error = 0.05;
    double max_align_error = 0.01;

    unsigned int nb_crossings;

    // crease type
    bool is_ridge;
    int crease_kind;
    unsigned int failed_conv;

    // interface to gage
    MeasureWrapper* the_wrapper;
    
    std::map< FaceId, vec3 > reference_points;
    FaceId current_face_id;
    face_information current_face_info;
    std::vector< point_on_face > all_points_on_face;
    std::map< FaceId, unsigned int > face_found;
    std::vector< vec3 > all_face_points;
    std::vector< std::pair< unsigned int, unsigned int > > all_edges;
    unsigned int nb_segments;

    // extraction parameters
    
    extraction_method ext_meth;

    // crease points found on voxel faces
    std::vector< vec3 > all_face_points;

    int filter(const std::vector< double >& vals, const std::vector< double >& strs);

    int filter(const std::vector< double >& vals,
               const std::vector< double >& strs,
               const std::vector< double >& cfds);

    typedef std::vector< unsigned int > line;

    void extract_lines(std::vector< line >& creases, const Nrrd* nrrd);

    void connect_segments(std::vector< line >& creases);

    bool good_value(double val);
    bool good_strength(double str);
    void initialize();
};

}
}

// ---------------------------------------------------------------------------
// filtering
inline bool spurt::crease::extractor::good_value(double val)
{
    return (is_ridge ?
            (val > crease::value_threshold) :
            (val < crease::value_threshold));
}

inline bool spurt::crease::extractor::good_strength(double str)
{
    return(is_ridge ?
           (str < crease::strength_threshold) :
           (str > crease::strength_threshold));
}
