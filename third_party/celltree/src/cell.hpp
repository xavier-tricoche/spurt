#ifndef __cell_hpp
#define __cell_hpp

namespace celltree {

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

bool invert_cell( float* c, const float* pts, const float* pos, cell_kind kind );

void interpolate_cell( float* result, const float* values, const float* c,
                       const unsigned int dim, cell_kind kind );

unsigned int intersect_cell( float* t, const float* pts, const float* origin, const float* direction,
                     cell_kind kind );

} // namespace celltree

#endif // __cell_hpp
