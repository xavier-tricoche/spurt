//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet3DRectilinear.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:10 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FPositionSet3DRectilinear_hh
#define __FPositionSet3DRectilinear_hh

#include "FPositionSet.hh"

class FIndex;

/** 
 *
 */

class FPositionSet3DRectilinear: public FPositionSet
{
  //==== Constructors ========================================================
public:

  /** 
   *\par Description:
   *Constructor
   *\pre
   *Xcoords and YCoords and ZCoords save the x/y/z -coordinate 
   *values in increasing order
   *\post 
   *Xcoords and YCoords and ZCoords are empty vectors 
   *(because of swap operation)
   *\exception
   *none
   *\param XCoords: vector with x coordinate values
   *\param YCoords: vector with y coordinate values
   *\param ZCoords: vector with z coordinate values
   */
  FPositionSet3DRectilinear( std::vector<double>& XCoords,
			     std::vector<double>& YCoords,
			     std::vector<double>& ZCoords );
  /**
  * Constructor
  * \pre dx, dy, dz are positive
  * \param nx number of points in x direction
  * \param ny number of points in y direction
  * \param nz number of points in z direction
  * \param dx x spacing
  * \param dy y spacing
  * \param dz z spacing
  */
  FPositionSet3DRectilinear( 
      unsigned int nx, unsigned int ny, unsigned int nz,
      double dx, double dy, double dz );
  
  //=== Member Functions ====================================================

  //=========================================================================
  //=== QUERY FUNCTIONS  ====================================================
  //=========================================================================

  /** 
   *\par Description:
   *Gets the class name.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  const FString& getClassName() const;

  /**
   *\par Description:
   *gets a reference (const) to the bounding box
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *none
   */
  const FBoundingBox& getBoundingBox() const;

  /** 
   *\par Description:
   *Returns the dimension of the positions belonging to this FPositionSet.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *none
   */
  
  void getPosition(FPosition& resultPos, const FIndex& pIndex) const;

  void getPosition( std::vector<double>& resultPos, const FIndex& pIndex) const;

  const std::vector<double>& getXCoords() const;
  const std::vector<double>& getYCoords() const;
  const std::vector<double>& getZCoords() const;

  positive getNbXCoords() const;
  positive getNbYCoords() const;
  positive getNbZCoords() const;

  // remove these three functions, they are no longer needed
  // the interface above is consistent with all other classes
  positive getDimensionX() const;
  positive getDimensionY() const;
  positive getDimensionZ() const;

  positive getNbPositions() const;

  positive memSize() const;

private:

  //position distribution for the different coordinates
  std::vector<double> pos[3];

  positive nbPositions;

  friend class FCellLocator3DRectilinear;
  friend class FGrid3DRectilinear;
};

#endif 
