#ifndef __interpolator_hpp
#define __interpolator_hpp

#include <locator.hpp>
#include <mesh_traits.hpp>
#include "variable.hpp"

namespace celltree {

// -------------------------------------------------------------------------

class interpolator: public locator
{
    typedef mesh_traits<mesh>          mtraits;
    typedef variable_traits<variable>  vtraits;

public:
    
    typedef float coord_type;
    typedef float value_type;
 
    interpolator( const mesh* m, const variable* v, const celltree& ct ) :
        locator( m, ct ), m_var(v)
    {
    }

    unsigned int dim() const 
    {
        return vtraits::dim( *m_var );
    }

    bool operator()( const float time, const float* pos, float* result ) const
    {
        mtraits::size_type  cell;
        mtraits::coord_type coord[3];
        
        if( find_cell( pos, cell, coord ) )
        {
            unsigned int ind[8], nind;
            nind = mtraits::get_indices( *(this->m_mesh), cell, ind );
            
            float tmp[24];
            vtraits::copy_values( *m_var, ind, ind+nind, tmp );
            mtraits::interpolate( *(this->m_mesh), cell, coord, tmp, m_var->dim, result );

            return true;
        }
        
        return false;
    }

    const variable* m_var;
};

}

#endif // __interpolator_hpp
