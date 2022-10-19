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

#include "FCellDefinitions2DTriangulation.hh"
#include "FNeighborhoodDataTriangulation.hh"

#include "FTriangleCell2D.hh"

#include <cassert>

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions2DTriangulation::
FCellDefinitions2DTriangulation( positive nb_pos, 
				 const std::string& newname,
				 std::vector<FIndex> & vertices,
				 bool buildneighborhooddata )
    : FCellDefinitions2DUnstructured( nb_pos ,newname )
{
#ifndef NODEBUG
    cout << "nb of vertices in cell def= " << vertices.size() << endl;
#endif

    assert( vertices.size() );
    assert( !( vertices.size()%3 ) ); 
    assert( nb_pos );
  
    nbCells = vertices.size()/3;
    nbPos = nb_pos;
    
    cell_vertices.swap( vertices );

    if ( buildneighborhooddata )
	neighborData 
	    = new FNeighborhoodDataTriangulation( this, cell_vertices );
}

//--------------------------------------------------------------------------- 

FCellDefinitions2DTriangulation::~FCellDefinitions2DTriangulation()
{
    // neighborData is deleted by base class
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions2DTriangulation::getClassName() const
{
  static FString name("FCellDefinitions2DTriangulation");
  return name;
}


//--------------------------------------------------------------------------- 

void FCellDefinitions2DTriangulation::getCellVerticesIndices( const FIndex& cellId, 
							      std::vector< FIndex >& vertices ) const
{
    assert( cellId < nbCells );

    vertices.clear();

    for( positive i=cellId.getIndex()*3 ; i<(cellId.getIndex()+1)*3 ; ++i )
	vertices.push_back( cell_vertices[i] );
}

//--------------------------------------------------------------------------- 

void FCellDefinitions2DTriangulation::getCellType( const FIndex& cellId, 
						   FCell::CellType& cell_type ) const
{
    assert( cellId < nbCells );

    cell_type = FCell::TRIANGLE_2D;
}

//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions2DTriangulation::getCellTorso(const FIndex& cellId) const 
{
    assert( cellId < nbCells );

    return shared_ptr<FCell>
	( new FTriangleCell2D 
	  ( vector< FIndex >
	    (cell_vertices.begin()+3*cellId.getIndex(),
	     cell_vertices.begin()+3*( cellId.getIndex()+1 ) ) 
	      ) 
	    );
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions2DTriangulation::memSize() const
{
    return 
      cell_vertices.size()*sizeof(cell_vertices[0])
      +
      neighborData->memSize()
      +
      sizeof(*this);
}
