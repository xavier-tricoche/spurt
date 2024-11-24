#ifndef __CREASEEXTRACTION_HPP__
#define __CREASEEXTRACTION_HPP__

#include <math/types.hpp>
#include <data/image.hpp>
#include <data/edge.hpp>
#include <image/convert.hpp>
#include <image/crease.hpp>
#include <image/probe.hpp>
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <teem/nrrd.h>

namespace spurt {
namespace crease {

extern std::vector< double > crease_strength;
extern std::vector< unsigned int > isolated;
extern std::vector< bool > skipped;
extern std::vector< bool > singular;
extern float* reconstructed_image;
extern std::vector< double > crease_strength;
typedef raster_grid<long, double, 2> grid_type;


void extract(const Nrrd* nrrd,
             const grid_type& grid,
             double threshold,
             bool ridge,
             std::vector< vec2 >& intersections,
             std::vector< std::list< unsigned int > >& creases);
             
// helper functions to navigate the cell structure
int next(const std::pair< int, int >& choices, unsigned int cur);
int right(unsigned int i, unsigned int j, const grid_type& grid);
int top(unsigned int i, unsigned int j, const grid_type& grid);
int left(unsigned int i, unsigned int j, const grid_type& grid);
int down(unsigned int i, unsigned int j, const grid_type& grid);
unsigned int cellid(unsigned int i, unsigned int j,
                    const grid_type& grid);
};
};

inline int spurt::crease::
next(const std::pair< int, int >& choices, unsigned int cur)
{
    // choices.first < 0 => choices.second < 0
    if (choices.first != (int)cur) return choices.first;
    else return choices.second;
}

inline int spurt::crease::
right(unsigned int i, unsigned int j, const grid_type& grid)
{
    if (i < grid.resolution()[0] - 2) return i + 1 + j*(grid.resolution()[0] - 1);
    else return -1;
}

inline int spurt::crease::
top(unsigned int i, unsigned int j, const grid_type& grid)
{
    if (j < grid.resolution()[1] - 2) return i + (j + 1)*(grid.resolution()[0] - 1);
    else return -1;
}

inline int spurt::crease::
left(unsigned int i, unsigned int, const grid_type&)
{
    if (i > 0) return i -1;
    else return -1;
}

inline int spurt::crease::
down(unsigned int i, unsigned int j, const grid_type& grid)
{
    if (j > 0) return i + (j - 1)*(grid.resolution()[0] - 1);
    else return -1;
}

inline unsigned int spurt::crease::
cellid(unsigned int i, unsigned int j, const grid_type& grid)
{
    return i + j*(grid.resolution()[0] - 1);
}

#endif

