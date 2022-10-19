//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FBoundingBox.hh,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:01 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.15 $
//
//--------------------------------------------------------------------------- 

#ifndef __FBoundingBox_hh
#define __FBoundingBox_hh

#include "FPosition.hh"
#include "stdAliases.hh"  

//===========================================================================

/** The FBoundingBox is a bounding box, that is parallel to the
 *  coordinate axes.  
 */

class FBoundingBox
{
  //=== Constructors ========================================================

public:
  /**
   *{\bf Description:}\\
   *Default constructor,
   *All range values are initialized to zero
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *\param
   *dimension defines the dimension of the bounding box
   */
  FBoundingBox(positive dimension = 3);

  
  /**
   *{\bf Description:}\\
   *copy constructor\\
   *{\bf Precondition:}\\
   *none
   *{\bf Postcondition:}\\
   *none
   *\exception
   *none
   * \param
   * bbox FBoundingBox to be copied
   */ 
  FBoundingBox( const FBoundingBox& bbox );

  /**
   *Constructs 2-dimensional bounding box from arguments.
   *  \param
   *  minX minimun x coordinate
   *  \param
   *  maxX maximum x coordinate
   *  \param
   *  minY minimun y coordinate
   *  \param
   *  maxY maximum y coordinate
   */
  FBoundingBox(double minX, double minY, double maxX, double maxY); 

  /**
   *Constructs 3-dimensional bounding box from arguments.
   *  \param
   *  minX minimun x coordinate
   *  \param
   *  maxX maximum x coordinate
   *  \param
   *  minY minimun y coordinate
   *  \param
   *  maxY maximum y coordinate
   *  \param
   *  minZ minimun z coordinate
   *  \param
   *  maxZ maximum z coordinate
   */
  FBoundingBox(double minX, double minY, double minZ, 
	       double maxX, double maxY, double maxZ);

  //=== Destructor  ==========================================================
public:
  /**
   *{\bf Description:}\\
   *Destructor
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *range is empty
   *\exception
   *none
   *\param
   *none
   */
  ~FBoundingBox();

  //=== Member Functions ====================================================
public:
  /**
   *sets bounding box to given arguments.
   *  \param
   *  minX minimun x coordinate
   *  \param
   *  maxX maximum x coordinate
   *  \param
   *  minY minimun y coordinate
   *  \param
   *  maxY maximum y coordinate
   */
  void setBoundingBox(double minX, double minY, double maxX, double maxY); 

  /**
   *sets bounding box to given arguments.
   *  \param
   *  minX minimun x coordinate
   *  \param
   *  maxX maximum x coordinate
   *  \param
   *  minY minimun y coordinate
   *  \param
   *  maxY maximum y coordinate
   *  \param
   *  minZ minimun z coordinate
   *  \param
   *  maxZ maximum z coordinate
   */
  void setBoundingBox(double minX, double minY, double minZ, double maxX,
                      double maxY, double maxZ);
  
  /** assignment operator.
   *  \param
   *  bbox bounding box to be assigned
   */
  FBoundingBox& operator= (const FBoundingBox& bbox);


  /** add-assignment operator.
   *  \param
   *  bbox bounding box to be added
   */
  FBoundingBox& operator+= (const FBoundingBox& bbox);

  /** comparison of two BoundingBoxes
   * \param
   * otherBBox the other box to compare with
   */
  bool operator==(const FBoundingBox& otherBBox) const;

  /** comparison of two BoundingBoxes
   * \param
   * otherBBox the other box to compare with
   */
  bool operator< (const FBoundingBox& otherBBox) const;
  

  /** Computes the center of this bounding box.
   *  \param
   *  center computed center position
   */
  void computeCenter3D( FVector& center ) const;

  /** Returns the center of the bounding box
   */
  FPosition center() const;

  /**
   * Returns bounding box diagonal length (ie the diameter
   * of the bounding box).
   */
  double diagonal() const;


  /** Computes the radius of the smallest sphere containing this
   *  bounding box. 
   *  \return radius of the bounding sphere
   */
  double computeBoundingSphereRadius() const;


  /**
   *tests if the given position is inside this boundingbox.
   * \pre
   *The Dimension of position and this have to be equal
   *\param
   *position the position to test
   */
  bool isInside(const FPosition& position) const; 


  /**
   *tests if the given BoundingBox is inside this boundingbox
   *\pre
   *The Dimension of BBox and this have to be equal
   *\param
   *inBBox the bounding box to test
   */
  bool isInside(const FBoundingBox& inBBox) const;

  /**
   *{\bf Description:}\\
   *tests if the given position/BoundingBox is included in this
   *\\{\bf Precondition:}\\
   *The Dimension of position/BBox and this have to be equal
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *
   *\param
   *position the position to test
   *\param
   *inBBox the bounding box to test
   */
  bool includes(const FBoundingBox& inBBox) const;

  /**
   *{\bf Description:}\\
   *tests if the given BoundingBox intersects the other in any way
   *\\{\bf Precondition:}\\
   *The Dimension of BBox and this have to be equal
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *\param
   *inBBox the bounding box to test
   */
  bool intersects( const FBoundingBox& inBBox ) const;

  FBoundingBox intersection( const FBoundingBox& rhs ) const;
  
  /**
   *{\bf Description:}\\
   *gets the values of the bounding box
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *\param
   *result vector in which the result is stored (either of size 4 or 8)
   */
  void getRange(std::vector<double>& result) const;            

  /**
   *{\bf Description:}\\
   *gets the values of the bounding box
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *\param
   *result variables in which the result is stored
   */
  void getRange( double& xmin, double& xmax, 
		 double& ymin, double& ymax, 
		 double& zmin, double& zmax ) const;

  /**
   *{\bf Description:}\\
   *gets the dimension of the bounding box
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *\param
   *none
   */
  positive getDimension() const;

  /**
   *{\bf Description:}\\
   *checks if the given Position is inside and enlarges the bounding
   *box if necessary
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *The bounding box is larger than before
   *\exception
   *none
   *\param
   *inPosition the position to check
   */
  
  void resize(const FPosition&);

  double volume() const;

  FArray min() const;
  FArray max() const;

  FArray random() const;

  //=== Other Functions =====================================================

  friend std::ostream& operator<<( std::ostream&, const FBoundingBox&);
  friend std::istream& operator>>( std::istream&, FBoundingBox&);

  //=== Protected Member Functions ============================================

protected:

private:

  bool fIsInside(const std::vector<double>& inR) const;

  //=== Variables ===========================================================
private:  
  double *range;
  positive dim;
  const double epsilon;
};

//===========================================================================
#ifndef OUTLINE
#include "FBoundingBox.icc"
#endif
//===========================================================================
#endif // __FBoundingBox_hh
 
