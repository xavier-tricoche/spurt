//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile$
// Language:  C++
// Date:      $Date$
// Author:    $Author$
// Version:   $Revision$
//
//--------------------------------------------------------------------------- 

#include "FCellDefinitions2DPoints.hh"
#include "FNeighborhoodDataUnstructured.hh"

#include "FPointCell2D.hh"

#include <cassert>

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions2DPoints::
FCellDefinitions2DPoints( positive nb_pos, 
				 const std::string& newname )
    : FCellDefinitions( newname )
{
    assert( nb_pos );
  
    nbCells = nb_pos;
    nbPos = nb_pos;
    neighborData = new FNeighborhoodDataUnstructured(this );
}

//--------------------------------------------------------------------------- 

FCellDefinitions2DPoints::~FCellDefinitions2DPoints()
{
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions2DPoints::getClassName() const
{
  static FString name("FCellDefinitions2DPoints");
  return name;
}


//--------------------------------------------------------------------------- 

void FCellDefinitions2DPoints::getCellVerticesIndices( const FIndex& cellId, 
							      std::vector< FIndex >& vertices ) const
{
    assert( cellId < nbCells );

    vertices.clear();

	vertices.push_back( cellId.getIndex() );
}

//--------------------------------------------------------------------------- 

void FCellDefinitions2DPoints::getCellType( const FIndex& cellId, 
						   FCell::CellType& cell_type ) const
{
    assert( cellId < nbCells );

    cell_type = FCell::POINT_2D;
}

//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions2DPoints::getCellTorso(const FIndex& cellId) const 
{
    assert( cellId < nbCells );

    return shared_ptr<FCell> ( new FPointCell2D ( cellId ) );
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions2DPoints::memSize() const
{
  if ( neighborData )
  {
    return 
      neighborData->memSize()
      +
      sizeof(*this);
  }
  else
  {
    return sizeof( *this );
  }
}
