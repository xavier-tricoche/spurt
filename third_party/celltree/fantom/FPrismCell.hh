//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPrismCell.hh,v $
// Language:  C++
// Date:      $Date: 2003/09/30 14:47:49 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef __FPrismCell_hh
#define __FPrismCell_hh
#include "FIndex.hh"
#include "FCell.hh"
#include "stdAliases.hh"
#include "FTensor.hh"
#include "FAMSingularPoint.hh"
#include "FMatrix.hh"

#include <list>
#include <vector>
#include <utility>
#include <utility>

//===========================================================================

/** 
 *The FPrismCell class is an implementation of the tetrahedron
 *in 3D space derived from the FCell abstract class.
 */
class FPrismCell : public FCell
{
public:
  /** 
   *\par Description:
   * returns an tetrahedron with indices set to invalid
   */
  FPrismCell();

  /** 
   *\par Description:
   * Constructor: returns a tetrahedron which vertices are defined
   * by given indices
   *\param
   * vert: vector indices of the tetrahedron vertices
   */
  FPrismCell( const vector<FIndex>& vert);
  FPrismCell( const FIndex* vert);

  /** 
   *\par Description
   * Destructor
   */
  ~FPrismCell();


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



private:
  /**
     my... : arrays for FCell superclass
  */
  FIndex myIndices[6];
  FRefArray myPositions[6];
  double myPositionData[6][3];
  FRefTensor myTensors[6];

  //--------------------------------------------------------------------------

  typedef double double3[3];

  mutable double3 a,b,c,d,e; //buffer vars

  //--------------------------------------------------------------------------

  //derivations after the different local coords in position[0] 
  mutable double3 pe0,pe1,pe2;

  //derivations of the above after local coord 2
  mutable double3 pde0,pde1;

  //params for cubic equation
  mutable double3 pa,pb,pc;

  //scalar product with e2
  mutable double pe2a,pe2b,pe2c;

  mutable bool paramsComputed;

  inline void computeParams() const;

  //--------------------------------------------------------------------------

  mutable double3 localCoords;
  mutable double3 lastPosition;
  void computeLocalCoords(const double p[3]) const;    

  //position which comes from local coords
  inline void interpolatePositionFromLocals(double result[3]) const;


  //precondition:local coords were computed
  //compute derivation of position after local coords
  inline void dPosdLocal(double3 m[]) const;


  mutable FTensor te0,te1,te2,tde0,tde1;
  inline void computeTensorParams() const;
  //precondition: computeLocalCoords,computeTensorParams have been invoked
  inline void dTensdLocal(FTensor&result) const;

  //pre: locals computed
  inline void interpolateTensorsFromLocals(FTensor & result) const;


};

//=========================================================================== 

#endif // __FTetrahedronCell_hh
