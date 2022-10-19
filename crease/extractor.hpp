#ifndef __XAVIER_EXTRACTOR_HPP__
#define __XAVIER_EXTRACTOR_HPP__

#include <vector>
#include <map>

#include "measure_wrapper.hpp"
#include "math/fixed_vector.hpp"
#include "face.hpp"
#include "grid.hpp"
#include "crease.hpp"

namespace spurt {
namespace crease {

extern std::map< spurt::FaceId, vec3 > reference_points;
extern spurt::FaceId current_face_id;
extern face_information current_face_info;

struct point_on_face {
    point_on_face() {}
    
    vec3 p[4];
    vec3 q;
    vec3 e, g, Hg;
};
extern std::vector< point_on_face > all_points_on_face;

// extraction parameters
struct crease_point {
    crease_point() : p(), dot(0), ngm(1) {}
    
    vec3 p;
    double dot;
    double ngm;
};

// some helper functions
void PCA(std::vector< vec3 >& evec, const std::vector< vec3 >& dirs);
void Furst_orientation(vec3& ev1, const vec3& ev0);
void Naive_orientation(vec3& ev1, const vec3& ev0);

enum extraction_method { PVO, MC, ZG };
extern extraction_method ext_meth;

// crease points found on voxel faces
extern std::vector< vec3 > all_face_points;

int filter(const std::vector< double >& vals,
           const std::vector< double >& strs);
           
int filter(const std::vector< double >& vals,
           const std::vector< double >& strs,
           const std::vector< double >& cfds);
           
typedef std::vector< unsigned int > line;

void extract_lines(std::vector< line >& creases, const Nrrd* nrrd);

void connect_segments(std::vector< line >& creases);

bool good_value(double val);
bool good_strength(double str);
}
}

// ---------------------------------------------------------------------------
// filtering
inline bool spurt::crease::good_value(double val)
{
    return (crease::is_ridge ?
            (val > crease::value_threshold) :
            (val < crease::value_threshold));
}

inline bool spurt::crease::good_strength(double str)
{
    return(crease::is_ridge ?
           (str < crease::strength_threshold) :
           (str > crease::strength_threshold));
}

#endif



















































