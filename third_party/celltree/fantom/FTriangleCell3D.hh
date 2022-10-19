 //---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTriangleCell3D.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:13 $
// Author:    $Author: garth $
// Version:   $Revision: 1.17 $
//
//--------------------------------------------------------------------------- 

#ifndef __FTriangleCell3D_hh
#define __FTriangleCell3D_hh

#include "FCell2Din3D.hh"
#include "FObject.hh"
#include "FString.hh"
#include "FPosition.hh"
#include "FAMSingularPoint.hh"
#include "FArray.hh"

//===========================================================================

/** 
 * The FTriangleCell3D class is an implementation of the geometric triangle
 * derived from the FCell2Din3D abstract class, lying in a 3D space.
 * \NOTE: A Constructor which takes a vector of positions is implemented in
 * FCell.
 */

class FTriangleCell3D : public FCell2Din3D
{
public: 

  //========================================================================
  // Standard: 3*Constructor , Destructur
  //========================================================================
  
  /** 
   *{\bf Description:}\\
   *Constructor: returns an empty triangle lying in a 3D space. Indices are
   * are set to invalid.
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   */
  FTriangleCell3D();
  
  /** 
   *{\bf Description:}\\
   *Constructor: Returns an triangle lying in a 3D space which is defined
   * thru the indices of the vertices.
   *\param
   * vertIds : Vector should have the lenght of 3 !
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   */
  FTriangleCell3D( const vector<FIndex>& vertIds );
  
  /** 
   *{\bf Description:}\\
   *Constructor: Returns an triangle lying in a 3D space.
   *\param
   * vertIds : Pointer to the first index of the triangle. Only
   * 3 indices are necessary !
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   */
  FTriangleCell3D( const FIndex* vertIds );

  /** 
   *\par Description:
   * Destructor
   */
  ~FTriangleCell3D();
  

  //========================================================================
  // Functions for general use
  //========================================================================
  
  /** 
   *{\bf Description:}\\
   * Returns a pointer to a copy of the triangle.
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   * FException if operation fails !
   *\return
   * FCell*
   */
  FCell* getClone() const;

  /** 
   * {\bf Description:}\\
   *Indicates if the given \b position lies in the cell.
   *\pre
   *The position has dimension 3.
   *\post
   *none
   *\exception
   *FInvalidDimensionException
   *\param
   *position: Position to test.
   */
  bool isInside(const FPosition& position) const;

  /** 
   * {\bf Description:}\\
   * Computes the intersection of a segment with the edges of the triangle
   * and returns the corresponding edge index
   *\pre
   * The first point lies in the cell
   *\pre
   * The second point lies outside 
   *\pre
   * The segment lies in the plane of the cell
   *\post
   *none
   *\exception
   * FInvalidPositionException
   *\param
   *result: returned intersection point
   *\param
   *pos: start position
   *\param
   *vec: intersection direction
   */
  char intersectEdge( FPosition& result, const FPosition& pos,
		      const FArray& vec ) const;

//   /**
//    *\par Description:
//    * Same as function above... but static!
//    */
//   static bool isInside(const FPosition& pos, 
// 		       const vector< FPosition >& vertices);
  
  /** 
   *\par {\bfDescription:} \\
   * Returns the interpolated tensor value at \b position.
   *\pre
   * The position lies in the triangle.
   *\post
   * none
   *\exception
   * none
   *\param
   * result: returned interpolated tensor.
   *\param
   * position: position at which to compute the interpolant
   */
  void interpolate(FTensor& result, const FPosition& position) const;
  
  /** 
   *\par {\bfDescription:} \\
   * Returns the interpolated tensor value at \b position with respect
   * to the projected field
   *\pre
   * The position lies in the triangle.
   *\post
   * none
   *\exception
   * none
   *\param
   * result: returned interpolated tensor.
   *\param
   * position: position at which to compute the interpolant
   */
  void interpolateInPlane(FTensor& result, const FPosition& position) const;

  /** 
   *\par {\bf Description:}\\
   * Returns zeros of tensor field (pojected into triangle plane)
   * which lie within the cell.
   *\pre
   * none
   *\post
   * none
   *\exception
   * none
   *\param
   * result: returned zeros.
   */
  void getZeros(list<FAMSingularPoint>& result) const;
  
  /** 
   *\par {\bf Description:}\\
   * Returns zeros of tensor field (pojected into triangle plane)
   * even if they do not lie in the triangle 
   * but any where in world space.
   *\pre
   * none
   *\post
   * none
   *\exception
   * none
   *\param
   * result: returned zeros.
   */
  void getZerosArbitrary(list<FAMSingularPoint>& result) const;

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
   *Returns vectors representing the original vectors at the triangle's 
   *vertices projected into the triangles plane
   */
  vector<FTensor> getReducedVectors();

  /** 
   *\par Description:
   *Returns vectors representing the original vectors at the triangle's 
   *vertices projected into the triangles plane and in local coordinates.
   */
  vector<FTensor> getLocalVectors();


    /** 
   *\par {\bf Description:}\\
   *returns the type of the cell.
   *(similar to getClassName, but more efficient)
   */
  FCell::CellType getCellType(void) const;

  static const geoInfo myGeoDescription;

    /** 
   *\par  Description:
   * Compute and return the triangles normalized normal
   */
  const FVector& normal() const;

  /**
   *\return
   * approximate size of cell in bytes
   */
  positive memSize() const;

protected:
  /**
     my... : arrays for FCell superclass
  */
  FIndex myIndices[3];
  FRefArray myPositions[3];
  double myPositionData[9];
  FRefTensor myTensors[3];

  // holds lokal triangle basis e0,e1 and e2
  // wherby e2 ^= plane normal
  mutable vector<FArray> basis;
  // holds tensors which are projected onto the triangle plane.
  mutable vector<FTensor> projTensors;
  // same as "projTensors",but only in lokal triangle coordinates.
  mutable vector<FTensor> locTensors;
  // positions of triangle in lokal triangle coordinates
  mutable vector<FArray> locPos;

  mutable bool isBasisSet;
  mutable bool isProjected;
  
  /** 
   * {\bf Description:}\\
   * Calculates the lokal basis of the triangle.Fkt.gets called
   * in the Constructors.
   * \pre
   * The member variable "vector<FArray> myTriBasis" holds
   * the basis.
   */
  void calculateTriBasis() const;
  
  /** 
   * {\bf Description:}\\
   * Projects the vectors onto the triangle plane.  
   * \pre
   * The member variables projTensors and locTensors hold
   * the projeced vectors.
   */
  void projectTensors() const;


    /** 
   *\par Description:
   * Calculate the barycentric coordinates of a given position.
   *\exception
   * FExeption
   *\param
   * bary: Coordinates, pos: given position
   */
  void baryCoord( double* bary, const FPosition& pos ) const;

   /** 
   *\par Description:
   * Calculate the positions of the triangle in lokal triangle 
   * coordinates.
   *\exception
   * FExeption
   * \post
   * vector<FArray> locPos holds coordinates 
   */
  void calculateLocTriPos() const;
  
  /** 
   *\par Description:
   * To a given position calculate the lokal triangle coordinates
   * \param
   * world: world coordinates, loc: lokal coordinates 
   */
  void getLocTriPos(FArray &world, FArray &loc) const;

  /** 
   *\par Description:
   * Given a singular point calculate the linear nature and set
   * the max norm (vector case!).
   *\exception
   * FExeption
   * \param
   * result: setLinearNature will be set
   */
  void getSingularType(FAMSingularPoint &result) const;

};

#endif // __FTriangleCell3D_hh
