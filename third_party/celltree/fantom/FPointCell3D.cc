//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPointCell3D.cc,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:02 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.31 $
//
//--------------------------------------------------------------------------- 

#include "FPointCell3D.hh"

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
const FCell::geoInfo FPointCell3D::myGeoDescription =
  {
    // dimension
    3,
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

FPointCell3D::FPointCell3D()
  : FCell(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
}

//--------------------------------------------------------------------------- 

FPointCell3D::FPointCell3D( const vector<FIndex>& vertIds )
  : FCell(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
  vertexIndices[0]=vertIds[0] ;
}

//--------------------------------------------------------------------------- 

FPointCell3D::FPointCell3D( const FIndex& vertId )
  : FCell(myGeoDescription, myIndices, myPositions, myPositionData, myTensors )
{
  vertexIndices[0] = vertId;
}

//--------------------------------------------------------------------------- 

FPointCell3D::FPointCell3D( const FIndex*vertIds )
  : FCell(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
    vertexIndices[0]=vertIds[0] ;
}

//--------------------------------------------------------------------------- 

FPointCell3D::~FPointCell3D()
{}

//--------------------------------------------------------------------------- 

FCell*  FPointCell3D::getClone() const
{
  try 
    {
      return new FPointCell3D( vertexIndices );
    }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

bool FPointCell3D::isInside( const FPosition& /*pos*/ ) const
{
  return false;
}

//--------------------------------------------------------------------------- 

void FPointCell3D::interpolate(FTensor& /*result*/,
				  const FPosition& /*position*/) const
{
  THROW_EXCEPTION( FException, "Cannot interpolate in 3D Point Cell" );
}

//--------------------------------------------------------------------------- 

void FPointCell3D::getZeros(list<FAMSingularPoint>& /*result*/) const
{
  THROW_EXCEPTION( FException, "Should return points that are zero" );
}

//--------------------------------------------------------------------------- 

FCell::CellType FPointCell3D::getCellType(void) const
{
  return FCell::POINT_3D;
}



//--------------------------------------------------------------------------- 

void FPointCell3D::derivatives(FTensor& /*result*/, const FPosition&) const
{
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//--------------------------------------------------------------------------- 


positive FPointCell3D::memSize() const
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
