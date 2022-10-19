//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAxisParallelTetCell.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:29 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAxisParallelTetCell_hh
#define __FAxisParallelTetCell_hh

#include "stdAliases.hh"
#include "FTensor.hh"
#include "FMatrix.hh"
#include "FAMSingularPoint.hh"
#include "FCell.hh"

//===========================================================================

/** 
 *The FAxisParallelTetCell class is an implementation of the tetrahedron
 *in 3D space derived from the FCell abstract class.
 */

class FAxisParallelTetCell : public FCell
{
public:
  /** 
   *\par Description:
   * returns an tetrahedron with indices set to invalid
   */
  FAxisParallelTetCell();

  /** 
   *\par Description:
   * this cell describes a cell which comes from the
   * decomposition of an axis parallel hexahedron cell.
   * it is bounded by the conditions:
   * 0<x[s]<x[m]<x[b]<1, where x denotes the local coordinate in the
   * hexahedron
   *
   *\pre 
   *vert is in the correct order:
   * vert[0]= hvert[0], vert[1] =hvert[1<<b], 
   * vert[2] =hvert[(1<<b)|(1<<s)], vert[3]=hvert[7],
   * where hvert are the vertices of the surrounding hexahedron
   * in vtk_voxel enumeration
   *
   *\post none
   *\exception none
   *\param vert: indices of vertex coordinates 
   *\param s: index of smallest coordinate 
   *\param m: index of middle   coordinate
   *\param b: index of biggest  coordinate
   */
  FAxisParallelTetCell(const vector<FIndex>&vert, int s, int m, int b);
  FAxisParallelTetCell(const FIndex*vert, int s, int m, int b);

  /** 
   *\par Description
   * Destructor
   */
  ~FAxisParallelTetCell(); 


  /** 
   *\par Description:
   * returns a pointer on a copy of the tetrahedron.
   *\exception
   * FException if operation fails
   *\return
   * FCell*
   */
  FCell* getClone() const;


  /** 
   *\par Description:
   *implementation of FCell::interpolate(FTensor& result, 
   *		   const FPosition& position)
   */
  void interpolate(FTensor& result, 
		   const FPosition& position) const;

  /** 
   *\par Description:
   *implementation of FCell::derivatives(FTensor& result, 
   *		   const FPosition& position)
   */ 
  void derivatives(FTensor& result, 
		   const FPosition& position) const;

  /** 
   *\par Description:
   *implementation of FCell::getZeros(list<FAMSingularPoint>& result)
   */
  void getZeros(list<FAMSingularPoint>& result) const;

  /** 
   *\par Description:
   *impementation of FCell::isInside(const FPosition& position)
   */
  bool isInside(const FPosition& position) const;

  /**
   *\par Description:
   *same as function above... but static!
   */
  static bool isInside(const FPosition& pos, 
		       const vector< FPosition >& vertices );

  /** 
   *\par Description:
   *implementation of FCell::getCellType(void)
   */
  FCell::CellType getCellType(void) const;

  positive sizeOfCellType () const;

  static const geoInfo myGeoDescription;

  positive memSize() const;

protected:

  ///see FCell.hh
  void buildBoundingBox(void) const;
  

private:
  /**
     my... : arrays for FCell superclass
  */
  FIndex myIndices[4];
  FRefArray myPositions[4];
  double myPositionData[12];
  FRefTensor myTensors[4];
  
  
  mutable double ps,pm,pb;          
  // position of lower left corner
  mutable double invds,invdm,invdb; // 1.0/length of tetrahedron
  mutable int s,m,b;

  mutable FTensor deriv;


};

//=========================================================================== 

#endif // __FAxisParallelTetCell_hh
