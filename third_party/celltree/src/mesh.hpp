#ifndef __mesh_hpp
#define __mesh_hpp

#include <vector>
#include <limits>
#include <cstdlib>
#include <stdint.h>

#include "cell.hpp"
#include "mesh_traits.hpp"

namespace celltree {

struct point
{
    float coord[3];
};

// -------------------------------------------------------------------------

struct mesh
{
    struct cell
    {
        cell_kind    kind;
        unsigned int start;
    
        cell() {}
        cell( cell_kind k, unsigned int s ) : kind(k), start(s) {}
    };

    float*         points;
    cell*          cells;
    unsigned int*  indices;
    
    unsigned int   npoints;
    unsigned int   ncells;
    unsigned int   nindices;
    
    mesh( unsigned int np, unsigned int nc, unsigned int ni,
          float* p, cell* c, unsigned int* i ) : 
        npoints(np), ncells(nc), nindices(ni), points(p), cells(c), indices(i)
    {
    }
    
    ~mesh()
    {
        delete[] points;
        delete[] cells;
        delete[] indices;
    }
    
private:
    
    mesh( const mesh& );
};

// -------------------------------------------------------------------------

template<>
struct mesh_traits<mesh>
{
    typedef unsigned int size_type;
    typedef float        coord_type;
   
    static size_t memsize( const mesh& m ) 
    {
        return m.npoints * sizeof(float)*3 + 
               m.nindices * sizeof(uint32_t) + 
               m.ncells * sizeof(mesh::cell);
    }
   
    static void extents( const mesh& m, coord_type* min, coord_type* max )
    {
        const float* pptr = m.points;
        
        for( unsigned int d=0; d<3; ++d )
            min[d] = max[d] = *(pptr)++;

        for( size_type i=1; i<m.npoints; ++i )
        {
            for( unsigned int d=0; d<3; ++d, ++pptr )
            {
                if( *pptr < min[d] )    min[d] = *pptr;
                if( *pptr > max[d] )    max[d] = *pptr;
            }
        }
    }
   
    static void minmax( const mesh& m, size_type index, coord_type* min, coord_type* max )
    {
        const mesh::cell& c = m.cells[index];
        const size_type* ii = m.indices + c.start;

        for( unsigned int d=0; d<3; ++d )
        {
            min[d] =  std::numeric_limits<coord_type>::max();
            max[d] = -std::numeric_limits<coord_type>::max();
        }

        for( size_type i=0; i<cell_size( c.kind ); ++i, ++ii )
        {
            const float* p = m.points + 3*(*ii);

            for( unsigned int d=0; d<3; ++d )
            {
                if( p[d] < min[d] )   min[d] = p[d];
                if( p[d] > max[d] )   max[d] = p[d];
            }
        }
    }

    static void center( const mesh& m, size_type index, coord_type* center )
    {
        const mesh::cell& c = m.cells[index];
        const size_type* ii = m.indices + c.start;

        for( unsigned int d=0; d<3; ++d )
            center[d] = 0;
        
        for( size_type i=0; i<cell_size( c.kind ); ++i, ++ii )
        {
            const float* p = m.points + 3*(*ii);

            for( unsigned int d=0; d<3; ++d )
                center[d] += p[d];
        }

        for( unsigned int d=0; d<3; ++d )
            center[d] /= cell_size( c.kind );
    }

    static unsigned int cells_size( const mesh& m ) 
    {
        return m.ncells;
    }

    template<typename OIter>
    static size_type get_indices( const mesh& m, size_type index, OIter out )
    {
        const mesh::cell& c = m.cells[index];
        const size_type* ii = m.indices + c.start;

        for( unsigned int i=0; i<cell_size(c.kind); ++i )
            *(out++) = *(ii++);
        
        return cell_size(c.kind);
    }

    static unsigned int 
    intersect( const mesh& m, size_type index, float* t, const coord_type* orig, const coord_type* dir )
    {
        const mesh::cell& c = m.cells[index];
        const size_type* ii = m.indices + c.start;

        float tmp[24];
        
        for( unsigned int i=0; i<cell_size(c.kind); ++i, ++ii )
        {
            const float* p = m.points + 3*(*ii);
            
            for( unsigned int d=0; d<3; ++d )
                tmp[3*i+d] = p[d];
        }
        
        return intersect_cell( t, tmp, orig, dir, c.kind );
    }
    
    static unsigned int 
    invert( const mesh& m, size_type index, const coord_type* pos, coord_type* param )
    {
        const mesh::cell& c = m.cells[index];
        const size_type* ii = m.indices + c.start;

        float tmp[24];

        for( unsigned int i=0; i<cell_size(c.kind); ++i, ++ii )
        {
            const float* p = m.points + 3*(*ii);
            
            for( unsigned int d=0; d<3; ++d )
                tmp[3*i+d] = p[d];
        }
        
        if( invert_cell( param, tmp, pos, m.cells[index].kind ) )
            return cell_size(c.kind);
        
        return 0;
    }

    static void interpolate( const mesh& m, size_type index, const coord_type* param, const coord_type* values, 
                             unsigned int dim, coord_type* result )
    {
        const mesh::cell& c = m.cells[index];
        interpolate_cell( result, values, param, dim, c.kind );
    }

    static void random( const mesh& m, size_type index, coord_type* result )
    {
        const mesh::cell& c = m.cells[index];
        const size_type* ii = m.indices + c.start;
        
        coord_type param[3] = { (float)drand48(), (float)drand48(), (float)drand48() };
        
        float tmp[24];

        for( unsigned int i=0; i<cell_size(c.kind); ++i, ++ii )
        {
            const float* p = m.points + 3*(*ii);
            
            for( unsigned int d=0; d<3; ++d )
                tmp[3*i+d] = p[d];
        }

        interpolate_cell( result, tmp, param, 3, c.kind );
    }
};

} // namespace celltree

// -------------------------------------------------------------------------

#endif // __mesh_hpp
