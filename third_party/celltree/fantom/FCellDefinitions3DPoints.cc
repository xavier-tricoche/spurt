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

#include "FCellDefinitions3DPoints.hh"
#include "FNeighborhoodDataUnstructured.hh"

#include "FPointCell3D.hh"

#include <cassert>

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions3DPoints::
FCellDefinitions3DPoints( positive nb_pos, 
				 const std::string& newname )
    : FCellDefinitions( newname )
{

    assert( nb_pos );
  
    nbCells = nb_pos;
    nbPos = nb_pos;

    neighborData = new FNeighborhoodDataUnstructured(this );
}

//--------------------------------------------------------------------------- 

FCellDefinitions3DPoints::~FCellDefinitions3DPoints()
{
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions3DPoints::getClassName() const
{
  static FString name("FCellDefinitions3DPoints");
  return name;
}


//--------------------------------------------------------------------------- 

void FCellDefinitions3DPoints::getCellVerticesIndices( const FIndex& cellId, 
							      std::vector< FIndex >& vertices ) const
{
    assert( cellId < nbCells );

    vertices.clear();

	vertices.push_back( cellId.getIndex() );
}

//--------------------------------------------------------------------------- 

void FCellDefinitions3DPoints::getCellType( const FIndex& cellId, 
						   FCell::CellType& cell_type ) const
{
    assert( cellId < nbCells );

    cell_type = FCell::POINT_3D;
}

//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions3DPoints::getCellTorso(const FIndex& cellId) const 
{
    assert( cellId < nbCells );

    return shared_ptr<FCell>
	( new FPointCell3D ( cellId ) );
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions3DPoints::memSize() const
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
