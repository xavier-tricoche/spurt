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

#ifndef __interpolator_hpp
#define __interpolator_hpp

#include "locator.hpp"
#include "mesh_traits.hpp"

namespace CellTree {

class interpolator: public locator
{
    typedef mesh_traits<mesh>          mtraits;
    typedef variable_traits<variable>  vtraits;

public:

    typedef double coord_type;
    typedef double value_type;

    interpolator( const mesh* m, const variable* v, const celltree& ct ) :
        locator( m, ct ), m_var(v)
    {
    }

    ~interpolator()
    {
        SafeDelete(m_var);
    }

    unsigned int dim() const
    {
        return vtraits::dim( *m_var );
    }

    bool operator()( const double time, const double* pos, double* result ) const
    {
        mtraits::size_type  cell;
        mtraits::coord_type coord[3];

        if( find_cell( pos, cell, coord ) )
        {

            // std::cout << "functor: interpolation successful at (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
            unsigned int ind[8], nind;
            nind = mtraits::get_indices( *(this->m_mesh), cell, ind );

            double tmp[24];
            vtraits::copy_values( *m_var, ind, ind+nind, tmp );
            mtraits::interpolate( *(this->m_mesh), cell, coord, tmp, m_var->dim, result );

            return true;
        }

        return false;
    }

    void evaluate(const mtraits::size_type cell, mtraits::coord_type *local, double* result) const {
        unsigned int ind[8], nind;
        nind = mtraits::get_indices( *(this->m_mesh), cell, ind );

        // std::cout << "evaluate: interpolation successful at (" << local[0] << ", " << local[1] << ", " << local[2] << ")\n";


        double tmp[24];
        vtraits::copy_values( *m_var, ind, ind+nind, tmp );
        mtraits::interpolate( *(this->m_mesh), cell, local, tmp, m_var->dim, result );
    }

    bool find(const double* pos, mtraits::size_type& cell, mtraits::coord_type* local) const {
        // std::cout << "calling findcell at (" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")\n";
        return find_cell( pos, cell, local );
    }

    const variable* m_var;
};

} // namespace CellTree


#endif // __interpolator_hpp
