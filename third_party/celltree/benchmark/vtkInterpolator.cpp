#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkAbstractCellLocator.h>
#include <vtkGenericCell.h>

#include "vtkInterpolator.hpp"
#include "vtkDLRReader.hpp"

// -------------------------------------------------------------------------

vtkInterpolator::vtkInterpolator( vtkDataSet* ds, vtkAbstractCellLocator* loc ) :
    m_ds(ds), m_loc(loc), m_cell(vtkGenericCell::New())
{
    m_var  = m_ds->GetPointData()->GetVectors();
}

// -------------------------------------------------------------------------

bool vtkInterpolator::operator()( double time, const double* pos, double* res ) const
{
    double x[3] = { pos[0], pos[1], pos[2] }, tol2, pcoords[3], weights[8];
    
    vtkIdType cellid = m_loc->FindCell( x, tol2, m_cell, pcoords, weights );
  
    if( cellid != -1 )
    {
        double vector[3];

        if( m_cell->GetNumberOfPoints() == 0 )
        {
            fprintf( stderr, "warning: cell %ld has zero points\n", (long int)cellid );
            
            m_cell->PrintSelf( std::cout, vtkIndent() );
            abort();
        }

        
        vtkIdType id = m_cell->GetPointId(0);
        m_var->GetTuple( id, vector );

        for( int j=0; j<m_var->GetNumberOfComponents(); ++j )
            res[j] = vector[j] * weights[0];

        for( int i=1; i<m_cell->GetNumberOfPoints(); ++i )
        {
            id = m_cell->GetPointId(i);
            m_var->GetTuple( id, vector );
            
            for( int j=0; j<m_var->GetNumberOfComponents(); ++j )
                res[j] += vector[j] * weights[i];
        }
        
        return true;
    }
    
    return false;
}
