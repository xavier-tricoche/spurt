#include <vtkPolyData.h>
#include <vtkObjectFactory.h>

#include "vtkCellTree.hpp"

#include "vtk_mesh_traits.hpp"
#include "celltree.hpp"
#include "celltree_builder.hpp"

//----------------------------------------------------------------------------

vtkCxxRevisionMacro(vtkCellTree, "$Revision: 1.0 $");
vtkStandardNewMacro(vtkCellTree);

//----------------------------------------------------------------------------

//////////////////////////////////////////////////////////////////////////////
// Main management and support for tree
//////////////////////////////////////////////////////////////////////////////

vtkCellTree::vtkCellTree(void) 
{
    this->ct             = NULL;
    this->LazyEvaluation = 1;
}

//---------------------------------------------------------------------------

vtkCellTree::~vtkCellTree(void) 
{
    this->FreeSearchStructure();
    this->FreeCellBounds();
}

//---------------------------------------------------------------------------

void vtkCellTree::FreeSearchStructure(void) 
{
    if( this->ct ) 
        delete this->ct;
}

//---------------------------------------------------------------------------

void vtkCellTree::BuildLocator()
{
    if( this->LazyEvaluation ) 
        return;
    
    this->ForceBuildLocator();
}

//---------------------------------------------------------------------------

void vtkCellTree::BuildLocatorIfNeeded()
{
    if( this->LazyEvaluation ) 
    {
        if( !ct || (ct && (this->MTime>this->BuildTime)) ) 
        {
            this->Modified();
            vtkDebugMacro(<< "Forcing BuildLocator");
            this->ForceBuildLocator();
        }
    }
}

//---------------------------------------------------------------------------

void vtkCellTree::ForceBuildLocator()
{
  //
  // don't rebuild if build time is newer than modified and dataset modified time
  if( ct && (this->BuildTime>this->MTime) && (this->BuildTime>DataSet->GetMTime()) ) 
      return;

  // don't rebuild if UseExistingSearchStructure is ON and a tree structure already exists
  if( (ct) && this->UseExistingSearchStructure ) 
  {
      this->BuildTime.Modified();
      vtkDebugMacro(<< "BuildLocator exited - UseExistingSearchStructure");
      return;
  }

  this->BuildLocatorInternal();
}

//---------------------------------------------------------------------------

void vtkCellTree::BuildLocatorInternal()
{
    //
    vtkIdType numCells;
    if( !this->DataSet || (numCells = this->DataSet->GetNumberOfCells()) < 1 ) 
    {
        vtkDebugMacro( << "No Cells to divide");
        numCells = 0;
    }

    vtkDebugMacro( << "Creating BSPTree for " << numCells << " cells");

    this->FreeSearchStructure();
    this->FreeCellBounds();

    ct = new celltree;

    celltree_builder builder;

    builder.m_leafsize = 8;
    builder.m_buckets  = 13;

    builder.build( *ct, this->DataSet );

    this->BuildTime.Modified();
}

void vtkCellTree::GenerateRepresentation(int level, vtkPolyData *vtkNotUsed(pd)) 
{
}

void vtkCellTree::GenerateRepresentationLeafs(vtkPolyData *pd) 
{
    GenerateRepresentation(-1,pd);
}

//---------------------------------------------------------------------------

int vtkCellTree::IntersectWithLine( double p1[3], double p2[3], double tol,
                                    double &t, double x[3], double pcoords[3], 
                                    int &subId, vtkIdType &cellId, vtkGenericCell *cell )
{
    int hit = this->IntersectWithLine(p1, p2, tol, t, x, pcoords, subId, cellId);

    if(hit)
        this->DataSet->GetCell( cellId, cell );

    return hit;
}

//---------------------------------------------------------------------------

int vtkCellTree::IntersectWithLine( double p1[3], double p2[3], double tol,
                                    double &t, double x[3], double pcoords[3], 
                                    int &subId, vtkIdType &cellId )
{
}

//---------------------------------------------------------------------------

// int vtkCellTree::IntersectCellInternal( vtkIdType cell_ID, double p1[3], double p2[3], 
//                                         double tol, double &t, double ipt[3], 
//                                         double pcoords[3], int &subId )
// {
//     this->DataSet->GetCell( cell_ID, this->GenericCell );
//     return this->GenericCell->IntersectWithLine( p1, p2, tol, t, ipt, pcoords, subId );
// }

//---------------------------------------------------------------------------

struct leaf_functor
{
    leaf_functor( vtkDataSet* ds, const double pos[3] ) : 
        m_ds(ds), m_id(-1)
    {
        m_pos[0] = pos[0];
        m_pos[1] = pos[1];
        m_pos[2] = pos[2];
        
        m_cell = vtkGenericCell::New();
    }
    
    ~leaf_functor()
    {
        m_cell->Delete();
    }
    
    bool operator()( const uint32_t* begin, const uint32_t* end )
    {
        double pcoord[3], w[8], dist2;
        int    subid;
        
        for( ; begin != end ; ++begin )
        {
            // double bnd[6];
            // m_ds->GetCellBounds( *begin, bnd );
            // 
            // if( m_pos[0] < bnd[0] || m_pos[0] > bnd[1] ||
            //     m_pos[1] < bnd[2] || m_pos[1] > bnd[3] ||
            //     m_pos[2] < bnd[4] || m_pos[2] > bnd[5] )
            //     continue;
                
            m_ds->GetCell( *begin, m_cell );
            
            
            if( m_cell->EvaluatePosition( m_pos, 0, subid, pcoord, dist2, w ) )
            {
                m_id = *begin;
                return true;
            }
        }

        return false;
    }
    
    vtkDataSet*     m_ds;
    double          m_pos[3];
    vtkIdType       m_id;
    vtkGenericCell* m_cell;
};

vtkIdType vtkCellTree::FindCell( double x[3], double, vtkGenericCell *cell, 
                                 double pcoords[3], double *weights )
{
    this->BuildLocatorIfNeeded();

    double pcoord[3], w[8], bnd[6], dist2;
    int    subid;
    unsigned int nodes;

    const float _x[3] = { x[0], x[1], x[2] };
    celltree::point_traversal pt( *ct, _x, &nodes );

    while( const celltree::node* n = pt.next() )
    {
        const unsigned int* begin = &ct->leaves[n->start()];
        const unsigned int* end   = begin + n->size();

        for( ; begin!=end; ++begin )
        {
            this->DataSet->GetCellBounds( *begin, bnd );
            
            if( x[0] < bnd[0] || x[0] > bnd[1] ||
                x[1] < bnd[2] || x[1] > bnd[3] ||
                x[2] < bnd[4] || x[2] > bnd[5] )
                continue;
                
            this->DataSet->GetCell( *begin, cell );
                        
            if( cell->EvaluatePosition( x, 0, subid, pcoord, dist2, w ) )
                return *begin;
        }    
    }
    
    return -1;
}

//---------------------------------------------------------------------------

void vtkCellTree::PrintSelf( ostream& os, vtkIndent indent )
{
  this->Superclass::PrintSelf( os,indent );
}

//---------------------------------------------------------------------------

unsigned long vtkCellTree::GetMemorySize()
{
    this->BuildLocatorIfNeeded();
    return ct->memsize();
}

