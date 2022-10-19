#ifndef __CREASEEXTRACTION_HPP__
#define __CREASEEXTRACTION_HPP__

#include <math/fixed_vector.hpp>
#include <data/raster.hpp>
#include <data/edge.hpp>
#include <image/convert.hpp>
#include <image/crease.hpp>
#include <image/probe.hpp>
#include <limits>
#include <vector>
#include <map>
#include <set>
#include <teem/nrrd.h>

namespace xavier {
namespace crease {

extern std::vector< double > crease_strength;
extern std::vector< unsigned int > isolated;
extern std::vector< bool > skipped;
extern std::vector< bool > singular;
extern float* reconstructed_image;
extern std::vector< double > crease_strength;


void extract(const Nrrd* nrrd,
             const raster_grid<2>& grid,
             double threshold,
             bool ridge,
             std::vector< nvis::vec2 >& intersections,
             std::vector< std::list< unsigned int > >& creases);
             
// helper functions to navigate the cell structure
int next(const std::pair< int, int >& choices, unsigned int cur);
int right(unsigned int i, unsigned int j, const raster_grid<2>& grid);
int top(unsigned int i, unsigned int j, const raster_grid<2>& grid);
int left(unsigned int i, unsigned int j, const raster_grid<2>& grid);
int down(unsigned int i, unsigned int j, const raster_grid<2>& grid);
unsigned int cellid(unsigned int i, unsigned int j,
                    const raster_grid<2>& grid);
};
};

inline int xavier::crease::
next(const std::pair< int, int >& choices, unsigned int cur)
{
    // choices.first < 0 => choices.second < 0
    if (choices.first != (int)cur) return choices.first;
    else return choices.second;
}

inline int xavier::crease::
right(unsigned int i, unsigned int j, const xavier::raster_grid<2>& grid)
{
    if (i < grid.resolution()[0] - 2) return i + 1 + j*(grid.resolution()[0] - 1);
    else return -1;
}

inline int xavier::crease::
top(unsigned int i, unsigned int j, const xavier::raster_grid<2>& grid)
{
    if (j < grid.resolution()[1] - 2) return i + (j + 1)*(grid.resolution()[0] - 1);
    else return -1;
}

inline int xavier::crease::
left(unsigned int i, unsigned int, const xavier::raster_grid<2>&)
{
    if (i > 0) return i -1;
    else return -1;
}

inline int xavier::crease::
down(unsigned int i, unsigned int j, const xavier::raster_grid<2>& grid)
{
    if (j > 0) return i + (j - 1)*(grid.resolution()[0] - 1);
    else return -1;
}

inline unsigned int xavier::crease::
cellid(unsigned int i, unsigned int j, const xavier::raster_grid<2>& grid)
{
    return i + j*(grid.resolution()[0] - 1);
}

#endif

