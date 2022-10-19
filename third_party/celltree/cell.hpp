/*************************************************************************
* celltree: boundary interval hierarchy tree-based solution for point 
*           location and interpolation in unstructured meshes. This
*           technique was first proposed in:
*           "Fast, Memory-Efficient Cell location in Unstructured Grids
*            for Visualization"
*                by Christop Garth and Kenneth I Joy, VisWeek 2011.
* Author: Christoph Garth
* Copyright (C) 2011 Christoph Garth
* All rights reserved
**************************************************************************/

#ifndef __cell_hpp
#define __cell_hpp

namespace CellTree {

enum cell_kind
{
    TETRAHEDRON,
    HEXAHEDRON,
    PYRAMID,
    PRISM,
};

// -------------------------------------------------------------------------

inline unsigned int cell_size( cell_kind kind )
{
    switch( kind )
    {
    case TETRAHEDRON: return 4;
    case HEXAHEDRON:  return 8;
    case PRISM:       return 6;
    case PYRAMID:     return 5;
    }

    return 0;
}

// -------------------------------------------------------------------------

struct cell
{
    cell_kind    kind;
    unsigned int start;
    
    cell() {}
    cell( cell_kind k, unsigned int s ) : kind(k), start(s) {}
};

// -------------------------------------------------------------------------

bool invert_cell( double* c, const double* pts, const double* pos, cell_kind kind );

void interpolate_cell( double* result, const double* values, const double* c,
                       const unsigned int dim, cell_kind kind );

unsigned int intersect_cell( double* t, const double* pts, const double* origin, const double* direction,
                     cell_kind kind );
                     
} // namespace CellTree
                     
#include "detail/cell.hpp"

#endif // __cell_hpp
