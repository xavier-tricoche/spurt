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

#include "FCellDefinitions2DLines.hh"
#include "FNeighborhoodDataUnstructured.hh"

#include "FLineCell3D.hh"

#include "eassert.hh"

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions2DLines::
FCellDefinitions2DLines( positive nb_pos, 
				 const std::string& newname,
				 std::vector<FIndex> & vertices,
				 bool /*buildneighborhooddata*/ )
    : FCellDefinitions( newname )
{
#ifndef NODEBUG
    cout << "nb of vertices in cell def= " << vertices.size() << endl;
#endif

    eassert( vertices.size() );
    eassert( !( vertices.size()%2) ); 
    eassert( nb_pos );
  
    nbCells = vertices.size()/2;
    nbPos = nb_pos;
    
    cell_vertices.swap( vertices );

//    if ( buildneighborhooddata )
//	neighborData 
//	    = new FNeighborhoodDataLines( this, cell_vertices );
    neighborData = new FNeighborhoodDataUnstructured(this );
}

//--------------------------------------------------------------------------- 

FCellDefinitions2DLines::~FCellDefinitions2DLines()
{
    // neighborData is deleted by base class
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions2DLines::getClassName() const
{
  static FString name("FCellDefinitions2DLines");
  return name;
}


//--------------------------------------------------------------------------- 

void FCellDefinitions2DLines::getCellVerticesIndices( const FIndex& cellId, 
							      std::vector< FIndex >& vertices ) const
{
    eassert( cellId < nbCells );

    vertices.clear();

    for( positive i=cellId.getIndex()*2 ; i<(cellId.getIndex()+1)*2 ; ++i )
	vertices.push_back( cell_vertices[i] );
}

//--------------------------------------------------------------------------- 

void FCellDefinitions2DLines::getCellType( const FIndex& cellId, 
						   FCell::CellType& cell_type ) const
{
    eassert( cellId < nbCells );

    cell_type = FCell::LINE_2D;
}

//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions2DLines::getCellTorso(const FIndex& cellId) const 
{
    eassert( cellId < nbCells );

    return shared_ptr<FCell>
	( new FLineCell3D 
	  ( vector< FIndex >
	    (cell_vertices.begin()+2*cellId.getIndex(),
	     cell_vertices.begin()+2*( cellId.getIndex()+1 ) ) 
	      ) 
	    );
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions2DLines::memSize() const
{
    return 
      cell_vertices.size()*sizeof(cell_vertices[0])
      +
      neighborData->memSize()
      +
      sizeof(*this);
}
