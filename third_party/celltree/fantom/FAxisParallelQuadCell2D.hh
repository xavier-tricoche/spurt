//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAxisParallelQuadCell2D.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:29 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 


#ifndef __FAxisParallelQuadCell2D_hh
#define __FAxisParallelQuadCell2D_hh

#include "FCell2D.hh"
#include "FObject.hh"
#include "FString.hh"
#include "FPosition.hh"
#include "FAMSingularPoint.hh"
#include "FArray.hh"

#include <list>
#include <vector>
#include <complex>

//===========================================================================

/**
 *The FAxisParallelQuadCell2D class is a derived type of FCell2D that 
 *represents rectangle cells in 2D rectilinear grids.
 */
class FAxisParallelQuadCell2D : public FCell2D
{
public:

  /**
   *\par Description:
   *Constructor: gives an empty rectilinear.
   */
  FAxisParallelQuadCell2D();

  /**
   *\par Description:
   *Constructor: gives a 2D rectangle.
   *\pre
   *\exception
   *none
   *\param
   *vert1Id,\b vert2Id,\b vert3Id,\b vert4Id: indices of the rectangle 
   *vertices.
   */
  FAxisParallelQuadCell2D(const vector<FIndex>& vertIds);
  FAxisParallelQuadCell2D(const FIndex* vertIds);
  
  /**
   *\par Description:
   *Destructor.
   */
  ~FAxisParallelQuadCell2D();

  /** 
   *\par Description:
   *Returns a pointer on a copy of the rectangle.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  FCell* getClone() const;

  /** 
   *\par Description:
   *Returns the interpolated tensor value at \b position.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned interpolated tensor.
   *\param
   *position: position at which to compute the interpolant
   */
  void interpolate(FTensor& result, const FPosition& position) const;

  /** 
   *\par Description:
   *Returns the derivatives at \b position.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned derivatives.
   *\param
   *position: position at which to compute the derivatives
   */
  void derivatives(FTensor& result, const FPosition& position) const;

  /** 
   *\par Description:
   *Returns zeros lying in the cell.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned zeros.
   */
  void getZeros(list<FAMSingularPoint>& result) const;

  /** 
   *\par Description:
   *Indicates if the given \b position lies in the cell.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *position: position to test.
   */
  bool isInside(const FPosition& position) const;

  /**
   *\par Description:
   *same as function above... but static!
   */
  static bool isInside(const FPosition& pos, 
		       const vector< FIndex >& verticesId,
		       FPositionSet *posSet);

  /** 
   *\par Description:
   *Builds the cell bounding box.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  void buildBoundingBox(void) const;

  /** 
   *\par Description:
   *returns the type of the cell. (similar to getClassName, but more efficient
   */
  FCell::CellType getCellType(void) const;

  /** 
   *\par Description:
   *Returns the size of the FIndex array that is needed by the class.
   *\param
   *none
   *\return
   *number of needed FIndexes.
   */
  positive sizeOfCellType() const;

  static const geoInfo myGeoDescription;

  positive memSize() const;

private:
  /**
     my... : arrays for FCell superclass
  */
  FIndex myIndices[4];
  FRefArray myPositions[4];
  double myPositionData[8];
  FRefTensor myTensors[4];


};


#endif // __FAxisParallelQuadCell2D_hh
