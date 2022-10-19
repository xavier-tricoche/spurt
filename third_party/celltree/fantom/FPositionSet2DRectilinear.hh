//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet2DRectilinear.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:09 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FPositionSet2DRectilinear_hh
#define __FPositionSet2DRectilinear_hh

#include "FPositionSet.hh"
#include <vector>

/** 
    the class has the functionality its name says
 */
class FPositionSet2DRectilinear: 
  public FPositionSet
	      
{
  //  MAKE_SERIALIZABLE( FPositionSet2DRectilinear );

  //==== Constructors ========================================================
public:

  /** 
   *\par Description:
   *Constructor
   *\pre
   *Xcoords and YCoords save the x/y -coordinate 
   *values in increasing order
   *\post 
   *Xcoords and YCoords are empty vectors 
   *(because of swap operation)
   *\exception
   *none
   *\param
   *none
   *\param XCoords: vector with x coordinate values
   *\param YCoords: vector with y coordinate values
   */
  FPositionSet2DRectilinear( vector<double>& XCoords, vector<double>& YCoords );
 
  /**
   * Constructor
   * \pre nx and ny > 0
   * \pre dx and dy > 0
   * \param nx number of points in x direction
   * \param ny number of points in y direction
   * \param dx spacing in x direction
   * \param dy spacing in y direction
   */
  FPositionSet2DRectilinear( unsigned int nx, unsigned int ny, double dx, double dy );

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
   *Gets a position.
   *\pre
   *finalize() has been invoked
   *0 <= posId < nbPositions
   *\post
   *none
   *\exception
   *FIndexOutOfBoundsException
   *\param
   *resultPos FPosition to store the searched position (or vector to receive
   *the double coordinates).
   *\param
   *pIndex: index of Position to get
   */           
  void getPosition( FPosition& resultPos, const FIndex& pIndex ) const;

  void getPosition( std::vector<double>& resultPos, const FIndex& pIndex ) const;

  const std::vector<double>& getXCoords() const;

  const std::vector<double>& getYCoords() const;

  positive getNbPositions() const;

  positive memSize() const;

  //  void serialize( std::ostream& ) const;

  //  static FPositionSet2DRectilinear *rebuild( std::istream& );

  positive getNbXCoords() const;
  positive getNbYCoords() const;


protected:
  positive nbPositions;
  //position distribution for the different coordinates
  std::vector<double> pos[2];

private:

  friend class FCellLocator2DRectilinear;
  friend class FGrid2DRectilinear;
};

#endif // __FPositionSet_hh
