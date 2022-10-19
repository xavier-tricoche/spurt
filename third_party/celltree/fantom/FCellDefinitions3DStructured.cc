//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions3DStructured.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:02 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#include "FCellDefinitions3DStructured.hh"
#include "FNeighborhoodDataStructured3D.hh"

#include "FArbitraryHexahedronCell.hh"
#include "FTetrahedronCell.hh"

#include <cassert>

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions3DStructured::FCellDefinitions3DStructured( positive sizeX, 
							    positive sizeY, 
							    positive sizeZ, 
							    const std::string& newname,
							    bool triang )
    : FCellDefinitions( newname ), triangulated( triang ), 
      nbX ( sizeX ), nbY ( sizeY ), nbZ( sizeZ )

{
    nbCells = nbX * nbY * nbZ;

    if( triangulated ) 
	nbCells*=6;

    nbPos = ( nbX + 1 ) * ( nbY + 1 ) * ( nbZ + 1 );

    neighborData = new FNeighborhoodDataStructured3D( sizeX + 1, 
						      sizeY + 1,
						      sizeZ + 1,
						      triangulated );
}

//--------------------------------------------------------------------------- 

FCellDefinitions3DStructured::~FCellDefinitions3DStructured()
{
    // neighborData is deleted in base class
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions3DStructured::getClassName() const
{
  static FString name("FCellDefinitions3DStructured");
  return name;
}

//--------------------------------------------------------------------------- 

void FCellDefinitions3DStructured::getCellVerticesIndices( const FIndex& cellId, 
							   std::vector<FIndex>& vertices ) const
{
    assert( cellId < nbCells );

    if( triangulated ) 
    {
	// tetrahedrized case:
        // cube is counter-clockwise enumerated,i.e.
	// first cut goes through plane {0,2,6,4}
        
	static positive tetraVertices[6][4] = {
	    {0,1,2,6},{0,5,1,6},{0,2,3,6},
	    {0,3,7,6},{0,4,5,6},{0,7,4,6}
	};
    
	// stepsize of vertex index for
	// one step in Y, Z-direction
	// (in X dir it is one)
	positive stepY = nbX + 1, 
	         stepZ = stepY * (nbY + 1);
    
	// hexahedron index
	positive hexInd = cellId/6, baseaddr;
    	  
	// compute index of lower left edge of cell out of cell Index
	// xCoords-1 : because in one direction ,
	// there are xCoords pos, but only xCoords-1 cells
	baseaddr =        (hexInd%(nbX)); 
	hexInd/=(nbX);
	baseaddr += stepY*(hexInd%(nbY)); 
	hexInd/=(nbY);
	baseaddr += stepZ*hexInd;
    
	positive *ind=tetraVertices[cellId%6];
    
	// compute vertex indices from tetraeder-in-voxel-indices
	vertices.resize(4);
    
	for(int j=0;j<4;j++) 
	{
	    int i = ind[j];
	    int index = baseaddr;
      
	    // if least significant bit of i is one, add one to baseaddr
	    index += i&1; 
      
            // right shift bits of i by one
	    i>>=1;
      
	    // if actual LSB of i = 2nd LSB of ind[j]  
	    // is one, add stepY
	    index += stepY * (i&1);
      
            // right shift bits of i by one 
	    i>>=1;
      
	    // if actual LSB of i = 3rd LSB of ind[j] 
	    // is one, add stepZ
	    index += stepZ * (i&1);
      
	    vertices[j] = index;
	}
    }
    else 
    {
	// non-tetrahedrized case
	positive baseaddr, hexInd = cellId;
	positive stepY = (nbX + 1), stepZ = stepY * (nbY + 1);
    
	//compute index of lower left edge of cell out of hexind
	baseaddr =          (hexInd%(nbX)); hexInd /= (nbX);
	baseaddr += stepY * (hexInd%(nbY)); hexInd /= (nbY);
	baseaddr += stepZ *  hexInd ;
    
	//compute VertexIndex from cellIndex
	vertices.resize(8);
	vector<FIndex>::iterator v = vertices.begin();
    
	v[0] = baseaddr;
	v[1] = baseaddr + 1;
	v[2] = baseaddr + stepY + 1;
	v[3] = baseaddr + stepY;
    
	v[4] = v[0] + stepZ;
	v[5] = v[1] + stepZ;
	v[6] = v[2] + stepZ;
	v[7] = v[3] + stepZ;
    }  
}

//--------------------------------------------------------------------------- 

bool FCellDefinitions3DStructured::getTriangulated() const 
{
  return triangulated;
}

//--------------------------------------------------------------------------- 

void FCellDefinitions3DStructured::getCellType( const FIndex& /*cellId*/, 
						FCell::CellType& cell_type ) const
{
  if( !triangulated )
      cell_type = FCell::ARBITRARY_HEX;
  else
      cell_type = FCell::TETRAHEDRON;
} 

//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions3DStructured::getCellTorso( const FIndex& cellId ) const
{
    vector<FIndex> vertices;
    getCellVerticesIndices(cellId,vertices);

    FCell *cell;

    if( triangulated )
	cell = new FTetrahedronCell(vertices);
    else
	cell = new FArbitraryHexahedronCell(vertices);
    
    return shared_ptr<FCell>( cell );
}

positive FCellDefinitions3DStructured::memSize() const
{
    return sizeof(*this);
}
  
positive FCellDefinitions3DStructured::getDimensionX() const
{
    return nbX;
}

positive FCellDefinitions3DStructured::getDimensionY() const
{
    return nbY;
}

positive FCellDefinitions3DStructured::getDimensionZ() const
{
    return nbZ; 
}
