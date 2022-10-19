//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPointCell2D.cc,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:02 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.31 $
//
//--------------------------------------------------------------------------- 

#include "FPointCell2D.hh"

#include "FException.hh"

#include "FArray.hh"
#include "FMatrix.hh"
#include "FPosition.hh"

#include <iostream>
#include <utility>
#include <vector>

//==========================================================================

#define FTRIANGLECELL3D_epsilon 1.0e-8

//==========================================================================

// first the cell's geometry description
const FCell::geoInfo FPointCell2D::myGeoDescription =
  {
    // dimension
    2,
    // # vertices
    1, 
    // # edges
    0, 
    // # faces
    0, 
    // edges
    {},
    // face sizes
    {},
    // faces
    {}
  };

//--------------------------------------------------------------------------- 
// now comes the fun stuff ...
//--------------------------------------------------------------------------- 


//--------------------------------------------------------------------------- 

FPointCell2D::FPointCell2D()
  : FCell(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
}

//--------------------------------------------------------------------------- 

FPointCell2D::FPointCell2D( const FIndex& vertId)
  : FCell(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
  vertexIndices[0] = vertId;
}

//--------------------------------------------------------------------------- 

FPointCell2D::FPointCell2D( const vector<FIndex>& vertIds )
  : FCell(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
  vertexIndices[0]=vertIds[0] ;
}

//--------------------------------------------------------------------------- 

FPointCell2D::FPointCell2D( const FIndex*vertIds )
  : FCell(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
    vertexIndices[0]=vertIds[0] ;
}

//--------------------------------------------------------------------------- 

FPointCell2D::~FPointCell2D()
{}

//--------------------------------------------------------------------------- 

FCell*  FPointCell2D::getClone() const
{
  try 
    {
      return new FPointCell2D( vertexIndices );
    }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

bool FPointCell2D::isInside( const FPosition& /*pos*/ ) const
{
  return false;
}

//--------------------------------------------------------------------------- 

void FPointCell2D::interpolate(FTensor& /*result*/,
				  const FPosition& /*position*/) const
{
  THROW_EXCEPTION( FException, "Cannot interpolate in 2D Point Cell" );
}

//--------------------------------------------------------------------------- 

void FPointCell2D::getZeros(list<FAMSingularPoint>& /*result*/) const
{
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//--------------------------------------------------------------------------- 

FCell::CellType FPointCell2D::getCellType(void) const
{
  return FCell::POINT_2D;
}



//--------------------------------------------------------------------------- 

void FPointCell2D::derivatives(FTensor&, const FPosition&) const
{
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//--------------------------------------------------------------------------- 

positive FPointCell2D::memSize() const
{
  /**
   * \todo FIXME
   */
  return
    (
     3//basis
     +3//locPos
    ) *3*sizeof(double) 
    +
    (tensorData?
     tensors[0].size()*sizeof(double)*
     sizeOfCellType() //tensorData
     + (sizeof(FTensor)+tensors[0].size()*sizeof(double))
     *(3//projTensors
       +3//locTensors
       )
     :0)
    +
    sizeof(*this);
}
//----------------------------------------------------------------------------
