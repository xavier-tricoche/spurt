//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions3DRectilinear.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:02 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#include "FCellDefinitions3DRectilinear.hh"

#include "FAxisParallelHexCell.hh"
#include "FTetrahedronCell.hh"

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions3DRectilinear::FCellDefinitions3DRectilinear( positive sizeX, 
							      positive sizeY, 
							      positive sizeZ, 
							      const std::string& newname,
							      bool triang )
    : FCellDefinitions3DStructured( sizeX,sizeY,sizeZ,newname,triang ) 
{
}

//--------------------------------------------------------------------------- 

FCellDefinitions3DRectilinear::~FCellDefinitions3DRectilinear()
{
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions3DRectilinear::getClassName() const
{
  static FString name("FCellDefinitions3DRectilinear");
  return name;
}


//--------------------------------------------------------------------------- 

void FCellDefinitions3DRectilinear::getCellType( const FIndex& /*cellId*/, 
						 FCell::CellType& cell_type ) const
{
    if( !triangulated )
	cell_type = FCell::AXIS_PARALLEL_HEX;
    else
	cell_type = FCell::TETRAHEDRON;
}

//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions3DRectilinear::getCellTorso( const FIndex& cellId ) const
{
    vector<FIndex> vertices;
    getCellVerticesIndices(cellId,vertices);

    FCell *cell;

    if( triangulated )
	cell = new FTetrahedronCell(vertices);
    else
	cell = new FAxisParallelHexCell(vertices);

    return shared_ptr<FCell>( cell );
}
