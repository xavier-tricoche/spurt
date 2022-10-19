//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTetrahedronCell.hh,v $
// Language:  C++
// Date:      $Date: 2003/11/05 20:57:09 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.15 $
//
//--------------------------------------------------------------------------- 

#ifndef __FTetrahedronCell_hh
#define __FTetrahedronCell_hh
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
 *The FTetrahedronCell class is an implementation of the tetrahedron
 *in 3D space derived from the FCell abstract class.
 */
class FTetrahedronCell : public FCell
{
public:
  /** 
   *\par Description:
   * returns an tetrahedron with indices set to invalid
   */
  FTetrahedronCell();

  /** 
   *\par Description:
   * Constructor: returns a tetrahedron which vertices are defined
   * by given indices
   *\param
   * vert: vector indices of the tetrahedron vertices
   */
  FTetrahedronCell( const vector<FIndex>& vert);
  FTetrahedronCell( const FIndex* vert);

  /** 
   *\par Description
   * Destructor
   */
  ~FTetrahedronCell();


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
   *\par {\bf Description:}\\
   * Returns zeros of tensor field
   * even if they do not lie in the tetrahedron 
   * but any where in world space.
   *\pre
   * none
   *\post
   * none
   *\exception
   * none
   *\param
   * result: returned zeros.
   *\return
   * bool: is zero inside cell?
   */
  bool getZerosArbitrary(list<FAMSingularPoint>& result) const;

  /** 
   *\par Description:
   *implementation of FCell::isInside
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
   *test
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  static bool isInside( const FPosition& position,
			FIndex* firstIndex,
			FPositionSet* posPtr );

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

  ///variable for baryCoord computation:
  /// (baryz0,baryz1,baryz2)=normalVectors*(x-p[3])
  mutable double normalVectors[3][3]; 
  mutable double denom;

  ///is the matrix normalVectors valid ?
  mutable bool normalVectorsComputed;

  ///derivation after global coordinates
  mutable FTensor deriv;

  ///compute the deriv tensor(which is then also used for interpolation)
  void buildDeriv() const;

  /// initializes the matrix for fast computation of barycentric coords
  void initBaryCoords() const;

  /// computation of barycentric coordinates (interpolation, isInside,...)
  void barycentricCoordinates(double b[4], const FPosition& position) const;

  /// ...and their derivatives
  void derivBaryCoord(double dbdx[4], double dbdy[4], double dbdz[4]) const;


private:
  /**
     my... : arrays for FCell superclass
  */
  FIndex myIndices[4];
  FRefArray myPositions[4];
  double myPositionData[12];
  FRefTensor myTensors[4];

   
  // enable reuse of barycentric coordinates computation
  mutable FPosition lastPos;
  mutable double lastBs[4];
  mutable bool reuseInit;
};

//=========================================================================== 

#endif // __FTetrahedronCell_hh
