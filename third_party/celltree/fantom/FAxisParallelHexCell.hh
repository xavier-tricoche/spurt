//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAxisParallelHexCell.hh,v $
// Language:  C++
// Date:      $Date: 2004/03/11 13:53:13 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAxisParallelHexCell_hh
#define __FAxisParallelHexCell_hh

#include "stdAliases.hh"
#include "FTensor.hh"
#include "FMatrix.hh"
#include "FAMSingularPoint.hh"
#include "FHexahedronCell.hh"
#include "FTrilinearInterpolation.hh"

//===========================================================================

/** 
 *The FAxisParallelHexCell class is an implementation of the tetrahedron
 *in 3D space derived from the FCell abstract class.
 */

class FAxisParallelHexCell : public FHexahedronCell
{
public:
  /** 
   *\par Description:
   * returns an hexahedron with indices set to invalid
   */
  FAxisParallelHexCell();

  /** 
   *\par Description:
   * see FHexahedronCell::FHexahedronCell(const vector<FIndex>&vert,bool useVTKHexahedronEnumeration=true)
   */
  FAxisParallelHexCell(const vector<FIndex>&vert,bool useVTKHexahedronEnumeration=true);
  FAxisParallelHexCell(const FIndex*);

  /** 
   *\par Description
   * Destructor
   */
  ~FAxisParallelHexCell(); 


  /** 
   *\par Description:
   * returns a pointer on a copy of the hexahedron.
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
   *impementation of FCell::isInside(const FPosition& position)
   */
  bool isInside(const FPosition& position) const;

  /**
   *\par Description:
   *same as function above... but static!
   */
  static bool isInside(const FPosition& pos, 
		       const vector< FPosition >& vertices,
		       bool useVTKHexahedronEnumeration=true );

  /** 
   *\par Description:
   *implementation of FCell::getCellType(void)
   */
  FCell::CellType getCellType(void) const;

  /**
   *\return
   * approximate size of cell in bytes
   */
  positive memSize() const;
protected:

  //from FHexahedronCell.hh
  void computeSing(float minb[3],float maxb[3],list<FAMSingularPoint>& result) const;

  ///see FCell.hh
  void buildBoundingBox(void) const;

  void computeTensorParams() const;

  mutable FTrilinearInterpolation tensorParams;

  //inverse size of cell in the 3 coordinates
  mutable double invLen[3];


};

//=========================================================================== 

#endif // __FAxisParallelHexCell_hh
