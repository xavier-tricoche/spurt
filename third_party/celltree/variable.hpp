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

#ifndef __variable_hpp
#define __variable_hpp

#include <vector>
#include "mesh_traits.hpp"

namespace CellTree {

struct variable
{
    double*      data;
    unsigned int dim;
    unsigned int size;
    
    variable( unsigned int _dim, unsigned int _size, double* _data ) :
        dim(_dim), size(_size), data(_data)
    {
    }
    
    ~variable()
    {
        delete[] data;
    }
    
private:
    
    variable( const variable& );
};

// -------------------------------------------------------------------------

template<> 
struct variable_traits<variable>
{
    typedef double value_type;
    
    static size_t memsize( const variable& v )
    {
        return v.dim * v.size * sizeof(double);
    }
    
    template<typename IIter, typename OIter>
    static OIter copy_values( const variable& v, IIter begin, IIter end, OIter out )
    {
        for( ; begin != end; ++begin )
        {
            const double* vi = v.data + v.dim * (*begin);
            
            for( size_t d=0; d<v.dim; ++d )
                *(out++) = *(vi++);
        }
        
        return out;
    }

    static size_t size( const variable& v )
    {
        return v.size;
    }

    static unsigned int dim( const variable& v )
    {
        return v.dim;
    }
};

} // namespace CellTree

#endif // __variable_hpp
