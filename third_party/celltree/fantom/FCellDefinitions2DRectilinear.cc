//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions2DRectilinear.cc,v $
// Language:  C++
// Date:      $Date: 2003/09/10 10:35:18 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FCellDefinitions2DRectilinear.hh"
#include "FAxisParallelTriangleCell2D.hh"
#include "FAxisParallelQuadCell2D.hh"

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions2DRectilinear::FCellDefinitions2DRectilinear( positive sizeX, 
							      positive sizeY,
							      const std::string& newname, 
							      bool triang )
    : FCellDefinitions2DStructured( sizeX, sizeY, newname, triang )
{
}

//--------------------------------------------------------------------------- 

FCellDefinitions2DRectilinear::~FCellDefinitions2DRectilinear()
{
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions2DRectilinear::getClassName() const
{
  static FString name("FCellDefinitions2DRectilinear");
  return name;
}


//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions2DRectilinear::getCellTorso( const FIndex& cellId ) const 
{
    // index check in base class

    vector<FIndex> vertices;
    getCellVerticesIndices(cellId,vertices);

    FCell *cell;

    if(triangulated)
	cell = new FAxisParallelTriangleCell2D(vertices);
    else
	cell = new FAxisParallelQuadCell2D(vertices);    

    return shared_ptr<FCell>( cell );
}

//--------------------------------------------------------------------------- 

void FCellDefinitions2DRectilinear::getCellType( const FIndex& /*cellId*/, 
						 FCell::CellType& cell_type ) const
{
    if( triangulated )
	cell_type = FCell::AXIS_PARALLEL_TRI_2D;
    else
	cell_type = FCell::AXIS_PARALLEL_QUAD_2D;
} 

//--------------------------------------------------------------------------- 
