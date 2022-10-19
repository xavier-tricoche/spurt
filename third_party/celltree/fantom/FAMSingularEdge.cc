//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularEdge.cc,v $
// Language:  C++
// Date:      $Date: 2003/03/27 08:21:25 $
// Author:    $Author: garth $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#include "FAMSingularEdge.hh"

//--------------------------------------------------------------------------- 

FAMSingularEdge::FAMSingularEdge() 
  : pos1(), pos2(), type( NONE ), cellId1(), cellId2()
{
}

//--------------------------------------------------------------------------- 

FAMSingularEdge::FAMSingularEdge( const FAMSingularEdge& edge ) 
  : FObject(), pos1(edge.pos1), pos2(edge.pos2), type(edge.type), 
    cellId1(edge.cellId1), cellId2(edge.cellId2)
{
}

//--------------------------------------------------------------------------- 

FAMSingularEdge::FAMSingularEdge(const FPosition& pos1, 
				 const FPosition& pos2,
				 const FIndex& cellId1, 
				 const FIndex& cellId2,
				 const edgeType& type)
{
  this->pos1 = pos1;
  this->pos2 = pos2;
  this->cellId1 = cellId1;
  this->cellId2 = cellId2;
  this->type = type;
}

//--------------------------------------------------------------------------- 

void FAMSingularEdge::getEdgeVertices( FPosition& pos1, FPosition& pos2 )
{
  pos1 = this->pos1;
  pos2 = this->pos2;
}

//--------------------------------------------------------------------------- 

void FAMSingularEdge::getEdgeType( edgeType& type )
{
  type = this->type;
}

//--------------------------------------------------------------------------- 

void FAMSingularEdge::getEdgeCells( FIndex& cell1, FIndex& cell2 )
{
  cell1 = cellId1;
  cell2 = cellId2;
}

//--------------------------------------------------------------------------- 

const FString& FAMSingularEdge::getClassName() const
{
  static FString className("FAMSingularEdge");

  return className;
}

//--------------------------------------------------------------------------- 
