//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions2Din3DUnstructured.cc,v $
// Language:  C++
// Date:      $Date: 2003/09/25 13:16:04 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FCellDefinitions2Din3DUnstructured.hh"
#include "FNeighborhoodDataUnstructured.hh"

#include "FTriangleCell3D.hh"
#include "FQuadrilateralCell3D.hh"

#include <cassert>

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions2Din3DUnstructured::FCellDefinitions2Din3DUnstructured( positive nb_pos, 
									const std::string& newname,
									std::vector< pair<FCell::CellType, unsigned int> >& types,
									std::vector<FIndex> & vertices ,
									bool checkConsistency)
    : FCellDefinitions( newname ),
      nbQuadrilateralCells( 0 ), nbTriangleCells( 0 )
{


  if(checkConsistency){

    vector< pair<FCell::CellType, unsigned int> >::iterator typeIt;
  
    vector<FIndex>::iterator vertIt;

    for(typeIt=types.begin();typeIt+1<types.end();typeIt++){

      assert( typeIt->second < vertices.size() );

      if( typeIt->first == FCell::QUADRILATERAL_3D )
	assert( (typeIt+1)->second - typeIt->second  == 4 );
      else if( typeIt->first == FCell::TRIANGLE_3D )
	assert( (typeIt+1)->second - typeIt->second  == 3 );
      else
	assert(false);
    }
    assert( typeIt->second == vertices.size() );
    
    for(vertIt=vertices.begin();vertIt!=vertices.end();vertIt++)
      assert( *vertIt<nb_pos );
      
  }
    

  assert( types.size() );
  assert( vertices.size() );
  assert( nb_pos );
  
  cell_types.swap( types );
  cell_vertices.swap( vertices );
  
  nbCells = cell_types.size()-1;
  nbPos = nb_pos;

  neighborData 
    = new FNeighborhoodDataUnstructured( this, cell_types, cell_vertices );
}

//--------------------------------------------------------------------------- 

FCellDefinitions2Din3DUnstructured::FCellDefinitions2Din3DUnstructured( positive nb_pos, 
									const std::string& newname)
    : FCellDefinitions( newname ),
      nbQuadrilateralCells( 0 ), nbTriangleCells( 0 )
{
  nbPos = nb_pos;
}

//--------------------------------------------------------------------------- 

FCellDefinitions2Din3DUnstructured::~FCellDefinitions2Din3DUnstructured()
{
    // neighborData is deleted by base class
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions2Din3DUnstructured::getClassName() const
{
  static FString name("FCellDefinitions2Din3DUnstructured");
  return name;
}


//--------------------------------------------------------------------------- 

void FCellDefinitions2Din3DUnstructured::getCellVerticesIndices( const FIndex& cellId, 
								 std::vector< FIndex >& vertices ) const
{
    assert( cellId < nbCells );

    vertices.clear();

    for( positive i = cell_types[cellId].second;
	 i < cell_types[(positive)cellId+1].second; 
	 // cellIndex + 1 is always valid thanks an additional last element 
	 // in cell_types
	 i++)
	
	vertices.push_back( cell_vertices[i] );
}

//--------------------------------------------------------------------------- 

void FCellDefinitions2Din3DUnstructured::getCellType( const FIndex& cellId, 
						      FCell::CellType& cell_type ) const
{
    assert( cellId < nbCells );

    cell_type = cell_types[cellId].first;
}

//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions2Din3DUnstructured::getCellTorso(const FIndex& cellId) const 
{
    assert( cellId < nbCells );

    switch( cell_types[cellId].first )
    {
    case FCell::TRIANGLE_3D:
	return shared_ptr<FCell>
	    ( new FTriangleCell3D 
	      ( vector< FIndex >
		(cell_vertices.begin()+cell_types[(positive)cellId   ].second,
		 cell_vertices.begin()+cell_types[(positive)cellId+1 ].second ) 
		  ) 
		);
	break;
    case FCell::QUADRILATERAL_3D:
	return shared_ptr<FCell>
	    ( new FQuadrilateralCell3D 
	      ( vector< FIndex >
		(cell_vertices.begin()+cell_types[(positive)cellId   ].second,
		 cell_vertices.begin()+cell_types[(positive)cellId+1 ].second ) 
		  ) 
		);
	break;
    default:
      FException e("FCellDefinitions2Din3DUnstructured::getCellTorso only suitable for  FCell::TRIANGLE_3D  FCell::QUADRILATERAL_3D");
      throw e;
    }
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions2Din3DUnstructured::memSize() const
{
    return 
      cell_types.size()*sizeof(cell_types[0])
      +
      cell_vertices.size()*sizeof(cell_vertices[0])
      +
      neighborData->memSize()
      +
      sizeof(*this);
}
