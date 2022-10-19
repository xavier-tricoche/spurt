//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet2DCurvilinear.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:09 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FPositionSet2DCurvilinear_hh
#define __FPositionSet2DCurvilinear_hh

#include "FPositionSet2DArbitrary.hh"

// forward-declaration
class FGrid2DCurvilinear;

/** 
 *  the class has the functionality its name says
 */
class FPositionSet2DCurvilinear: public FPositionSet2DArbitrary
{
  //==== Constructors ========================================================
public:

  /** 
   *\par Description:
   *Constructor
   *\pre
   *none
   *\post 
   *coords contains the x and y component for each position.
   *( pos0[x], pos0[y], pos1[x], pos1[y] ... )
   *\exception
   *none
   *\param
   *none
   *\param 
   */
  FPositionSet2DCurvilinear( vector<double>& coords,
			     positive newxCoords,
			     positive newyCoords );

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

  friend class FGrid2DCurvilinear;

  positive getNbXCoords() const;
  positive getNbYCoords() const;

private:

  positive xCoords, yCoords;

};

#endif // __FPositionSet2DCurvilinear_hh
