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

#include "FCellDefinitions3DTriangulation.hh"
#include "FNeighborhoodDataTriangulation.hh"

#include "FTriangleCell3D.hh"

#include <cassert>

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions3DTriangulation::
FCellDefinitions3DTriangulation( positive nb_pos, 
				 const std::string& newname,
				 std::vector<FIndex> & vertices,
				 bool buildneighborhooddata )
    : FCellDefinitions2Din3DUnstructured( nb_pos ,newname )
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

FCellDefinitions3DTriangulation::~FCellDefinitions3DTriangulation()
{
    // neighborData is deleted by base class
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions3DTriangulation::getClassName() const
{
  static FString name("FCellDefinitions3DTriangulation");
  return name;
}

//--------------------------------------------------------------------------- 

void FCellDefinitions3DTriangulation::getCellVerticesIndices( const FIndex& cellId, 
							      std::vector< FIndex >& vertices ) const
{
    assert( cellId < nbCells );

    vertices.clear();

    for( positive i=cellId.getIndex()*3 ; i<(cellId.getIndex()+1)*3 ; ++i )
	vertices.push_back( cell_vertices[i] );
}

//--------------------------------------------------------------------------- 

void FCellDefinitions3DTriangulation::getCellType( const FIndex& cellId, 
						   FCell::CellType& cell_type ) const
{
    assert( cellId < nbCells );

    cell_type = FCell::TRIANGLE_3D;
}

//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions3DTriangulation::getCellTorso(const FIndex& cellId) const 
{
    assert( cellId < nbCells );

    return shared_ptr<FCell>
	( new FTriangleCell3D 
	  ( vector< FIndex >
	    (cell_vertices.begin()+3*cellId.getIndex(),
	     cell_vertices.begin()+3*( cellId.getIndex()+1 ) ) 
	      ) 
	    );
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions3DTriangulation::memSize() const
{
    return 
      cell_vertices.size()*sizeof(cell_vertices[0])
      +
      neighborData->memSize()
      +
      sizeof(*this);
}
