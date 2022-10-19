//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FArbitraryHexahedronCell.hh,v $
// Language:  C++
// Date:      $Date: 2004/03/11 13:53:13 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef __FArbitraryHexahedronCell_hh
#define __FArbitraryHexahedronCell_hh

#include "stdAliases.hh"
#include "FTensor.hh"
#include "FMatrix.hh"
#include "FAMSingularPoint.hh"
#include "FHexahedronCell.hh"
#include "FTrilinearInterpolation.hh"

//===========================================================================

/** 
 *The FArbitraryHexahedronCell class is an implementation of the tetrahedron
 *in 3D space derived from the FCell abstract class.
 */

class FArbitraryHexahedronCell : public FHexahedronCell
{
public:
  /** 
   *\par Description:
   * returns an hexahedron with indices set to invalid
   */
  FArbitraryHexahedronCell();

  /** 
   *\par Description:
   * Constructor: returns a hexahedron which vertices are defined
   * by given indices
   * see constructor of FHexcahedroncell
   */
  FArbitraryHexahedronCell(const vector<FIndex>&,bool useVTKHexahedronEnumeration=true);
  FArbitraryHexahedronCell(const FIndex*);

  /** 
   *\par Description
   * Destructor
   */
  ~FArbitraryHexahedronCell(); 


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
   *Returns the interpolated tensor value at position.
   *\pre
   *The result tensor has been allocated with the right dimension
   *and order with respect to the current tensorset (ie same 
   *dimension and order as it).
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned interpolated tensor.
   *\param
   *position: position at which to compute the interpolant
   */
  void interpolate(FTensor& result, 
		   const FPosition& position) const;

  /** 
   *{\bf Description:}\\
   *Returns the derivatives at {\bf position}.
   *\\{\bf Precondition:}\\
   *The result tensor has been allocated with the right dimension
   *and order with respect to the current tensorset (ie same dimension
   *as it and one order more).
   *\post
   *none
   *\exception
   *none 
   *\param
   *result: returned derivatives.
   */ 
  void derivatives(FTensor& result, 
		   const FPosition& position) const;


 
  /** 
   *\par Description:
   *abstract function from FCell
   */
  bool isInside(const FPosition& position) const;

  /**
   *\par Description:
   *same as function above... but static!
   *\param useVTKHexahedronEnumeration
   * should vtk hexahedron or vtk voxel vertex enumeration
   * be used ?
   */
  static bool isInside(const FPosition& pos, 
		       const vector< FPosition >& vertices,
		       bool useVTKHexahedronEnumeration=true);
  
  /** 
   *\par Description:
   *returns the type of the cell. (similar to getClassName, but more efficient
   */
  FCell::CellType getCellType(void) const;

  //for statistical issues
  static int NbcomputeTril;
  static int Nbchecknormvect;

  /**
   *\par Description:
   *\see FCell::neighborFaceForPos
   */
  virtual bool neighborFaceForPos(const FArray& pos,FIndex&faceId ) const;

protected:
  
  //compute limiting faces of hexahedron
  //precondition:positionstoupdate has been invoked
  void computeNormalVectors() const;
  mutable bool normalVectorsComputed;
  mutable FVector normals[6];
  ///values for inner face equation
  mutable double loh[6];
  ///values for outer face equation  
  mutable double hih[6];

  ///compute local coords ( =trilCoords, s.d.)
  void computeTrilinearCoords(const FPosition&p,int coordToIgnore=-1) const ;
  ///coordToIgnore:
  ///this local coord is constant in the  the first step
  ///of the local coord computation
  ///if it ist set >=0, also the local coords are not reset


  ///has there been an invocation of computeTrilinearCoords with the new vertex positions?
  mutable bool trilinearCoordsComputed;
  mutable FPosition lastPositionEvaluated;

  ///  the trilinear coordinates which  c..Coords 
  ///computes and which is widely used
  mutable FVector trilCoords;


  void interpolate(FTensor&t) const;

  //compute derivation of tensorfield after 
  //the trilinear coordinate with number ind
  //global parameters used:trilCoords, tensorParams
  void computeTrilDerivTensor(FTensor&ret,int ind) const;
  void computeTrilDerivTensor(FTensor&ret) const;

  //----------------------------------------------

  void computeTensorParams() const;
  mutable FTrilinearInterpolation tensorParams;
  
  //----------------------------------------------

  /*
   *\par Description
   *compute global out of local coordinates
   */
  void interpolate(FPosition&p) const;

  /**
   *\par Description
   *compute derivation of position after local coordinates
   *global params used:trilCoords, pointParams
   *the trilinear coordinate with number ind
   */
  void computeTrilDeriv(FVector&ret,int ind) const;
  /**
   *\par Description:
   *see above
   */
  void computeTrilDeriv(FMatrix&ret) const;


  //----------------------------------------------

  void computePointParams() const;
  /** 
   *\par Description
   * true if pointParams are valid for the current positions.
   */  
  mutable bool pointParamsComputed;
  /** 
   *\par Description
   * interpolation parameters for current positions
   */
  mutable FTrilinearInterpolation pointParams;

  //----------------------------------------------

  //tests geometryOk and resets flags if false,sets geometryOk to true
  void resetGeometryFlags() const;

  //testing functions, throw them away if fully tested
  //they give true if error occurred


  //overloaded from FHexahedrinCell for getzeros function
  void computeSing(float minb[3],float maxb[3],list<FAMSingularPoint>& result) const;


  /**
   *\return
   * approximate size of cell in bytes
   */
  positive memSize() const;
};

//=========================================================================== 

#endif // __FArbitraryHexahedronCell_hh
