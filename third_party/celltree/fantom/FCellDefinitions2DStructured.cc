//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions2DStructured.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:01 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#include "FCellDefinitions2DStructured.hh"
#include "FNeighborhoodDataStructured2D.hh"

#include "FTriangleCell2D.hh"
#include "FQuadrilateralCell2D.hh"

#include <cassert>

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions2DStructured::FCellDefinitions2DStructured( positive sizeX, 
							    positive sizeY,
							    const std::string& newname, 
							    bool triang )
    : FCellDefinitions( newname ), triangulated( triang ),
      nbX ( sizeX ), nbY ( sizeY )
{
    nbCells = nbX * nbY;
    
    if( triangulated )
	nbCells*=2;

    nbPos = ( nbX + 1 ) * ( nbY + 1 );

    neighborData = new FNeighborhoodDataStructured2D( sizeX + 1, sizeY + 1,
						      triangulated );
}

//--------------------------------------------------------------------------- 

FCellDefinitions2DStructured::~FCellDefinitions2DStructured()
{
    // base class deletes neighborData
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions2DStructured::getClassName() const
{
  static FString name("FCellDefinitions2DStructured");
  return name;
}


//--------------------------------------------------------------------------- 

void FCellDefinitions2DStructured::getCellVerticesIndices( const FIndex& cellId, 
							   std::vector< FIndex >& vertices ) const
{
    assert( cellId < nbCells );

    vertices.clear();

    if( !triangulated ) 
    {
	positive vertexIndex;
      
	// left lower edge 
	vertexIndex = ( cellId % nbX ) 
	    + ( cellId / nbX ) * ( nbX + 1 );
	vertices.push_back( vertexIndex );
      
	// right lower edge
	vertexIndex++;
	vertices.push_back( vertexIndex );
      
	// right upper edge
	vertexIndex += nbX + 1;
	vertices.push_back( vertexIndex );
      
	// left upper edge
	vertexIndex--;
	vertices.push_back( vertexIndex );
    }
    else 
    {
	if( ( cellId % 2 ) == 0 ) 
	{ 
            // down 
	    FIndex tmpCellId = cellId / 2;
	    positive vertexIndex;
	
	    // left lower edge 
	    vertexIndex = ( tmpCellId % nbX ) 
		+ ( tmpCellId / nbX ) * ( nbX + 1 );    
	    vertices.push_back( vertexIndex );
	
	    // right lower edge
	    vertexIndex++;
	    vertices.push_back( vertexIndex );

	    // left upper edge
	    vertexIndex += nbX;
	    vertices.push_back( vertexIndex );
	}
	else 
	{ 
	    // up
	    FIndex tmpCellId = cellId / 2;
	    positive vertexIndex;
	    
	    // right lower edge 
	    vertexIndex = ( tmpCellId % nbX ) 
		+ ( tmpCellId / nbX ) * ( nbX + 1 ) + 1;
	    vertices.push_back( vertexIndex );
	    
	    // right upper edge
	    vertexIndex += nbX + 1;
	    vertices.push_back( vertexIndex );
	    
	    // left upper edge
	    vertexIndex--;
	    vertices.push_back( vertexIndex );
	}
    }
}

//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions2DStructured::getCellTorso(const FIndex& cellId) const 
{
    vector<FIndex> vertices;
    getCellVerticesIndices(cellId,vertices);
    
    FCell *cell;

    if( triangulated )
	cell = new FTriangleCell2D(vertices);
    else
	cell = new FQuadrilateralCell2D(vertices);    

    return shared_ptr<FCell>( cell );
}

//--------------------------------------------------------------------------- 

bool FCellDefinitions2DStructured::getTriangulated() const 
{
    return triangulated;
}

//--------------------------------------------------------------------------- 

void FCellDefinitions2DStructured::getCellType( const FIndex& /*cellId*/, 
						FCell::CellType& cell_type ) const
{
    if (!triangulated)
	cell_type = FCell::QUADRILATERAL_2D;
    else
	cell_type = FCell::TRIANGLE_2D;
} 

//--------------------------------------------------------------------------- 

positive FCellDefinitions2DStructured::memSize() const
{
    return sizeof(*this);
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions2DStructured::getDimensionX() const
{
  return nbX;
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions2DStructured::getDimensionY() const
{
  return nbY;
}

//--------------------------------------------------------------------------- 
