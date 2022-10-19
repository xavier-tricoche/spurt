#ifndef __vtk_mesh_traits_h

#include <vtkDataSet.h>
#include <vtkGenericCell.h>

#include "mesh_traits.hpp"

template<>
struct mesh_traits<vtkDataSet*>
{
    typedef vtkIdType size_type;
    typedef double    coord_type;
   
    static size_t memsize( vtkDataSet* ds ) 
    {
        return ds->GetActualMemorySize() * 1024;
    }
   
    template<typename V>
    static void extents( vtkDataSet* ds, V* min, V* max )
    {
        double bounds[6];
        ds->GetBounds( bounds );
        
        min[0] = bounds[0]; max[0] = bounds[1];
        min[1] = bounds[2]; max[1] = bounds[3];
        min[2] = bounds[4]; max[2] = bounds[5];
    }
   
    template<typename V>
    static void minmax( vtkDataSet* ds, size_type index, V* min, V* max )
    {
        double bounds[6];
        ds->GetCellBounds( index, bounds );

        min[0] = bounds[0]; max[0] = bounds[1];
        min[1] = bounds[2]; max[1] = bounds[3];
        min[2] = bounds[4]; max[2] = bounds[5];
    }

    template<typename V>
    static void center( vtkDataSet* ds, size_type index, V* center )
    {
        vtkGenericCell* cell = vtkGenericCell::New();
        ds->GetCell( index, cell );
        
        double pcenter[3];
        int subid = cell->GetParametricCenter( pcenter );
        
        double wcenter[3], w[8];
        cell->EvaluateLocation( subid, pcenter, wcenter, w );
        
        center[0] = wcenter[0];
        center[1] = wcenter[1];
        center[2] = wcenter[2];
        
        cell->Delete();
    }

    template<typename V>
    static void random( vtkDataSet* ds, size_type index, V* point )
    {
        vtkGenericCell* cell = vtkGenericCell::New();
        ds->GetCell( index, cell );
        
        double pcoord[3];
        
        pcoord[0] = drand48();
        pcoord[1] = drand48();
        pcoord[2] = drand48();
        
        double wpoint[3], w[8];
        
        int subid = 0;
        cell->EvaluateLocation( subid, pcoord, wpoint, w );
        
        point[0] = wpoint[0];
        point[1] = wpoint[1];
        point[2] = wpoint[2];
        
        cell->Delete();
    }

    static std::size_t points_size( vtkDataSet* ds ) 
    {
        return ds->GetNumberOfPoints();
    }
    
    static std::size_t cells_size( vtkDataSet* ds ) 
    {
        return ds->GetNumberOfCells();
    }

    // template<typename OIter>
    // static OIter copy_indices( const mesh& m, size_type index, OIter out )
    // {
    //     const mesh::cell& c = m.cells[index];
    //     const size_type* ii = m.indices + c.start;
    // 
    //     for( unsigned int i=0; i<cell_size(c.kind); ++i )
    //         *(out++) = *(ii++);
    //     
    //     return out;
    // }
    // 
    // static unsigned int 
    // intersect( const mesh& m, size_type index, float* t, const coord_type* orig, const coord_type* dir )
    // {
    //     const mesh::cell& c = m.cells[index];
    //     const size_type* ii = m.indices + c.start;
    // 
    //     float tmp[24];
    //     
    //     for( unsigned int i=0; i<cell_size(c.kind); ++i, ++ii )
    //     {
    //         const float* p = m.points + 3*(*ii);
    //         
    //         for( unsigned int d=0; d<3; ++d )
    //             tmp[3*i+d] = p[d];
    //     }
    //     
    //     return intersect_cell( t, tmp, orig, dir, c.kind );
    // }
};

#endif // __vtk_mesh_traits_h