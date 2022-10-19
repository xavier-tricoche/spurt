//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FLineCell2D.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:29 $
// Author:    $Author: garth $
// Version:   $Revision: 1.12 $
//
//--------------------------------------------------------------------------- 

#ifndef __FLineCell2D_hh
#define __FLineCell2D_hh

#include "FCell1Din2D.hh"
#include "FObject.hh"
#include "FString.hh"
#include "FPosition.hh"
#include "FAMSingularPoint.hh"
#include "FArray.hh"

//===========================================================================

/** 
 *The FLineCell2D class is an implementation of the geometric triangle
 *in 2D space derived from the FCell2D abstract class.
 */
class FLineCell2D : public FCell1Din2D
{
public:

  /** 
   *\par Description:
   * Constructor: returns an triangle lying in a 2D space with indices
   * set to invalid
  */
  FLineCell2D();

  /** 
   *\par Description:
   * Constructor: returns a triangle lying in a 2D space 
   * which vertices are defined by given indices 
   *\param
   * vert1: indices of the line vertices (must be 2)
   */
  FLineCell2D( const vector<FIndex>& vertIds );
  FLineCell2D( const FIndex* vertIds );

  /** 
   *\par Description:
   * Destructor
   */
  ~FLineCell2D();


  /** 
   *\par Description:
   * returns a pointer on a copy of the triangle.
   *\exception
   * FException if operation fails
   *\return
   * FCell*
   */
  FCell* getClone() const;

  /** 
   *\par Description:
   *Returns the interpolated tensor value at \b position.
   *\pre
   *The result tensor has been allocated with the right size 
   *corresponding to the dimension and order of thetensorset
   *(i.e. with the same dimension and the same order).
   *The position lies in the triangle.
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
   *The result tensor has been allocated with the right size 
   *corresponding to the dimension and order of the current tensorset
   *(i.e. with the same dimension and one order more)
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
   *The position has dimension 2.
   *\post
   *none
   *\exception
   *FInvalidDimensionException
   *\param
   *position: position to test.
   */
  bool isInside(const FPosition& position) const;

  /**
   *\par Description:
   *same as function above... but static!
   */
  static bool isInside(const FPosition& pos, 
		       const vector< FPosition >& vertices);

  /** 
   *\par Description:
   *returns the type of the cell.
   *(similar to getClassName, but more efficient)
   */
  FCell::CellType getCellType(void) const;


  /** 
   *\par Description:
   *Builds the cell bounding box.
   *\pre
   *the position cache cache must be valid, !AND! this function assumes, that
   *the bounding box of the position cache represents the boundingbox for the
   *whole cell, if this is not the case, the concrete cells have to overload
   *this function!
   *\post
   *none
   *\exception
   *none
   */
  void buildBoundingBox(void) const;

  static const geoInfo myGeoDescription;

  /**
   *\return
   * approximate size of cell in bytes
   */
  positive memSize() const;
private:
  /**
     my... : arrays for FCell superclass
  */
  FIndex myIndices[2];
  FRefArray myPositions[2];
  double myPositionData[4];
  FRefTensor myTensors[2];


#if 0 
  // computation of barycentric coordinates (interpolation, isInside,...)
  void barycentricCoordinates(double b[3], 
			      const FPosition& position) const;
  
  // ...and of their derivatives
  void derivBaryCoord(double dbdx[3],
		      double dbdy[3]) const;

#endif
  // enable reuse of barycentric coordinates
  mutable FPosition lastPos;
  mutable double lastBs[3];
};

#endif // __FLineCell2D_hh
