#ifndef __vtkinterpolator_h
#define __vtkinterpolator_h

class vtkAbstractCellLocator;
class vtkDataArray;
class vtkDataSet;
class vtkGenericCell;

class vtkInterpolator
{
public:
    
    typedef double coord_type;
    typedef double value_type;

    typedef vtkDataSet* mesh_type;
    
    vtkInterpolator( vtkDataSet* ds, vtkAbstractCellLocator* loc );

    bool operator()( double time, const double* pos, double* result ) const;

protected:    

    vtkDataSet*             m_ds;
    vtkAbstractCellLocator* m_loc;
    vtkDataArray*           m_var;
    vtkGenericCell*         m_cell;
};

#endif // __vtkinterpolator_h
