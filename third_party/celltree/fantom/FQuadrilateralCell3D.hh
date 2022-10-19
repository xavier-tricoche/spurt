//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FQuadrilateralCell3D.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:29 $
// Author:    $Author: garth $
// Version:   $Revision: 1.17 $
//
//--------------------------------------------------------------------------- 


#ifndef __FQuadrilateralCell3D_hh
#define __FQuadrilateralCell3D_hh

#include "FCell2Din3D.hh"

#include "FBilinearSurface.hh"

#include <list>
#include <vector>
#include <complex>

//===========================================================================

/**
 *The FQuadrilateralCell3D class is a derived type of FCell2D that represents
 *quadrilterals in a 2D space.
 */
class FQuadrilateralCell3D : public FCell2Din3D
{
public:

  /**
   *\par Description:
   *Constructor: gives an empty 2D quadrilateral.
   */
  FQuadrilateralCell3D();

  /**
   *\par Description:
   *Constructor: gives a 2D quadrilateral.
   *\pre
   *\exception
   *none
   *\param
   *vert1: array with indices of the cell vertices
   *vertices.
   */
  FQuadrilateralCell3D( const vector< FIndex >& vert1 );
  FQuadrilateralCell3D( const FIndex * vert1 );
  
  /**
   *\par Description:
   *Destructor.
   */
  ~FQuadrilateralCell3D();

   /** 
   *\par Description:
   *Returns a pointer on a copy of the quadrilateral.
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
   *the result tensor has been allocated with the right dimension and 
   *order, corresponding to the current tensorset (i.e. it must have 
   *the same dimension and order as the tensorset).
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
   *the result tensor has been allocated with the right dimension and 
   *order, corresponding to the current tensorset (i.e. it must have 
   *the same dimension and one order more than the tensorset).
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
		       const vector< FPosition >& vertices);

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

  /**
   *\return
   * approximate size of cell in bytes
   */
  positive memSize() const;

  /** 
   *\par Description:
   *Compute the bouding box 
   *\pre
   *none
   *\post
   *bBox has been set to its new value
   *\param
   *none
   */
  void buildBoundingBox(void) const;

  static const geoInfo myGeoDescription;

private:
  /**
     my... : arrays for FCell superclass
  */
  FIndex myIndices[4];
  FRefArray myPositions[4];
  double myPositionData[4][3];
  FRefTensor myTensors[4];


  mutable FBilinearSurface surf;
  mutable double normal[3];
   typedef double double3[3];

  mutable double3 localCoords;
  mutable double3 lastPosition;
   
  inline void computeLocalCoords(const double p[3]) const;    

  //precondition:local coords were computed
  //compute derivation of position after local coords
  inline void dPosdLocal(double3 m[]) const;
  
  inline void dTensdLocal(FTensor&result) const;

public:
  template<class T>
  static void fromToPixelEnum(T*v)
  {
    std::swap( v[2], v[3] );
  }

  template<class T>
  static void fromToPixelEnum(vector<T>&v)
  {
    std::swap( v[2], v[3] );
  }
  
};


#endif // __FQuadrilateralCell3D_hh
