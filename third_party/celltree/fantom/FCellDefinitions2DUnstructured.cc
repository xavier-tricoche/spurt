//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions2DUnstructured.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:02 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#include "FCellDefinitions2DUnstructured.hh"
#include "FNeighborhoodDataUnstructured.hh"

#include "FQuadrilateralCell2D.hh"
#include "FAxisParallelQuadCell2D.hh"
#include "FAxisParallelTriangleCell2D.hh"
#include "FTriangleCell2D.hh"

#include "FException.hh"

#include <cassert>

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions2DUnstructured::FCellDefinitions2DUnstructured( positive nb_pos, 
								const std::string& name,
								vector< pair<FCell::CellType,unsigned int> >& types,
								vector<FIndex> & vertices )
    : FCellDefinitions( name ), nbQuadrilateralCells( 0 ), nbTriangleCells( 0 )
{
    assert( types.size() );
    assert( vertices.size() );
    
    nbPos = nb_pos;
    cell_types.swap( types );
    cell_vertices.swap( vertices );

    // cell_types contains a dummy last cell   
    nbCells = cell_types.size()-1;

    neighborData = new FNeighborhoodDataUnstructured( this, 
 						      cell_types,
 						      cell_vertices );
}

FCellDefinitions2DUnstructured::FCellDefinitions2DUnstructured( positive nb_pos, 
				      const std::string& name)
: FCellDefinitions( name )
{
  nbPos = nb_pos;
}

//--------------------------------------------------------------------------- 

FCellDefinitions2DUnstructured::~FCellDefinitions2DUnstructured()
{
    // base class will free neighborData
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions2DUnstructured::getClassName() const
{
  static FString name("FCellDefinitions2DUnstructured");
  return name;
}


//--------------------------------------------------------------------------- 
#ifdef OLD_CODE
void FCellDefinitions2DUnstructured::addNewCell( FCell::CellType myCellType, 
						 vector< FIndex > vertices )
{
    cell_types.push_back( 
	pair<FCell::CellType,positive>( myCellType, cell_vertices.size() ) );

    for (unsigned int i=0 ; i<vertices.size() ; i++) 
	cell_vertices.push_back( vertices[i] );

    if (myCellType == FCell::QUADRILATERAL_2D)
	nbQuadrilateralCells++;
    else
	nbTriangleCells++;
}
#endif
//--------------------------------------------------------------------------- 

void FCellDefinitions2DUnstructured::getCellVerticesIndices( const FIndex& cellId, 
							     vector<FIndex>& vertices ) const
{
    assert( cellId < nbCells );
    vertices.clear();

    // cellIndex + 1 is always valid thanks an additional last element 
    // in cell_types
    for( positive i = cell_types[cellId].second;
	 i < cell_types[(positive)cellId+1].second; i++)
	vertices.push_back( cell_vertices[i] );
}

//--------------------------------------------------------------------------- 

void FCellDefinitions2DUnstructured::getCellType( const FIndex& cellId, 
						  FCell::CellType& cell_type ) const
{
    assert( cellId < nbCells );
    cell_type = cell_types[ cellId ].first;
}

//--------------------------------------------------------------------------- 
positive FCellDefinitions2DUnstructured::memSize() const
{        
    return
      cell_types.capacity()*sizeof(cell_types[0])
      +
      cell_vertices.capacity()*sizeof(cell_vertices[0])
      +
      neighborData->memSize()
      +
      sizeof(*this);
}
 
//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions2DUnstructured::getCellTorso( const FIndex& cellId ) const 
{
    assert( cellId < nbCells );
    
    FCell *cell;

    if( cell_types[cellId].first == FCell::TRIANGLE_2D )
	cell = new FTriangleCell2D(
	    vector<FIndex>( cell_vertices.begin() + cell_types[(positive)cellId  ].second,
			    cell_vertices.begin() + cell_types[(positive)cellId+1].second ) );
    else if( cell_types[cellId].first == FCell::QUADRILATERAL_2D )
	cell = new FQuadrilateralCell2D(     
	    vector<FIndex>( cell_vertices.begin() + cell_types[(positive)cellId  ].second,
			    cell_vertices.begin() + cell_types[(positive)cellId+1].second ) );
    else if( cell_types[cellId].first == FCell::AXIS_PARALLEL_QUAD_2D )
	cell = new FAxisParallelQuadCell2D(     
	    vector<FIndex>( cell_vertices.begin() + cell_types[(positive)cellId  ].second,
			    cell_vertices.begin() + cell_types[(positive)cellId+1].second ) );
    
    else if( cell_types[cellId].first == FCell::AXIS_PARALLEL_TRI_2D )
	cell = new FAxisParallelTriangleCell2D(     
	    vector<FIndex>( cell_vertices.begin() + cell_types[(positive)cellId  ].second,
			    cell_vertices.begin() + cell_types[(positive)cellId+1].second ) );
    else
    {
      cout<<cell_types[cellId].first<<endl;
      FException e("Until now FCellDefinitions2DUnstructured::getCellTorso only supports:\
                    \nFTriangleCell2D\nFQuadrilateralCell2D\nFAxisParallelQuadCell2D\nFAxisParallelTriangleCell2D");
      throw e;
    }

    return shared_ptr<FCell>( cell );
}

