#ifndef __variable_hpp
#define __variable_hpp

#include <vector>
#include "mesh_traits.hpp"

namespace celltree {

struct variable
{
    float*       data;
    unsigned int dim;
    unsigned int size;
    
    variable( unsigned int _dim, unsigned int _size, float* _data ) :
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
    typedef float value_type;
    
    static size_t memsize( const variable& v )
    {
        return v.dim * v.size * sizeof(float);
    }
    
    template<typename IIter, typename OIter>
    static OIter copy_values( const variable& v, IIter begin, IIter end, OIter out )
    {
        for( ; begin != end; ++begin )
        {
            const float* vi = v.data + v.dim * (*begin);
            
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

} // namespace celltree

#endif // __variable_hpp