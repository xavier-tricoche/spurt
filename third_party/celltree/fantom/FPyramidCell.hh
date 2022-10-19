//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPyramidCell.hh,v $
// Language:  C++
// Date:      $Date: 2004/03/11 13:53:14 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef __FPyramidCell_hh
#define __FPyramidCell_hh
#include "FIndex.hh"
#include "FCell.hh"
#include "stdAliases.hh"
#include "FTensor.hh"
#include "FAMSingularPoint.hh"
#include "FMatrix.hh"
#include "FBilinearSurface.hh"

#include <list>
#include <vector>
#include <utility>
#include <utility>

//===========================================================================

/** 
 *The FPyramidCell class is an implementation of the tetrahedron
 *in 3D space derived from the FCell abstract class.
 */
class FPyramidCell : public FCell
{
public:
  /** 
   *\par Description:
   * returns an tetrahedron with indices set to invalid
   */
  FPyramidCell();

  /** 
   *\par Description:
   * Constructor: returns a tetrahedron which vertices are defined
   * by given indices
   *\param
   * vert: vector indices of the tetrahedron vertices
   */
  FPyramidCell( const vector<FIndex>& vert);
  FPyramidCell( const FIndex* vert);

  /** 
   *\par Description
   * Destructor
   */
  ~FPyramidCell();


  /** 
   *\pre Description:
   * returns a pointer on a copy of the tetrahedron.
   *\exception
   * FException if operation fails
   *\return
   * FCell*
   */
  FCell* getClone() const;

  /** 
   *\par Description:
   *implementation of FCell::interpolate
   */
  void interpolate(FTensor& result, 
  		   const FPosition& position) const;
  
  /** 
   *\par Description:
   *implementation of FCell::derivatives
   */
  void derivatives(FTensor& result, 
  		   const FPosition& position) const;

  /** 
   *\par Description:
   *implementation of FCell::getZeros
   */
  void getZeros(list<FAMSingularPoint>& result) const;

  /** 
   *\par Description:
   *implementation of FCell::isInside
   */
  bool isInside(const FPosition& position) const;

  /**
   *\par Description:
   *same as function above... but static!
   */
  //static bool isInside(const FPosition& pos, 
  //		       const vector< FPosition >& vertices);

  /** 
   *\par Description:
   *test
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
//   static bool isInside( const FPosition& position,
// 			FIndex* firstIndex,
// 			FPositionSet* posPtr );

  /** 
   *\par Description:
   *implementation of FCell::sizeOfCelltype
   */
  positive sizeOfCellType() const;
  
  /** 
   *\par Description:
   *implementation of FCell::getCellType
   *\return FCell::TETRAHEDRON
   */
  FCell::CellType getCellType(void) const;

  static const geoInfo myGeoDescription;

  /**
   *\return
   * approximate size of cell in bytes
   */
  positive memSize() const;

  /**
   *\par Description:
   *\see FCell::neighborFaceForPos
   */
  virtual bool neighborFaceForPos(const FArray& pos,FIndex&faceId ) const;

protected:

  mutable FBilinearSurface bottom;

private:
  /**
     my... : arrays for FCell superclass
  */
  FIndex myIndices[5];
  FRefArray myPositions[5];
  double myPositionData[5][3];
  FRefTensor myTensors[5];


  typedef double double3[3];

  mutable double3 localCoords;
  mutable double3 lastPosition;
   
  inline void computeLocalCoords(const double p[3]) const;    

  //precondition:local coords were computed
  //compute derivation of position after local coords
  void dPosdLocal(double3 m[]) const;

  void dTensdLocal(FTensor&result) const;

};

//=========================================================================== 

#endif // __FTetrahedronCell_hh
