//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet3DCurvilinear.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:10 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FPositionSet3DCurvilinear_hh
#define __FPositionSet3DCurvilinear_hh

#include "FPositionSet3DArbitrary.hh"

// forward-declaration
class FGrid3DCurvilinear;

/** 
    the class has the functionality its name says
 */
class FPositionSet3DCurvilinear: public FPositionSet3DArbitrary
{
  //==== Constructors ========================================================
public:

  /** 
   *\par Description:
   *Constructor
   *\pre
   *coords contains the x,y and z component for each position.
   *( pos0[x], pos0[y], pos0[z], pos1[x] ... )
   *\post 
   *none
   *\exception
   *none
   *\param
   *none
   *\param 
   */
  FPositionSet3DCurvilinear( vector<double>& coords,
			     positive newxCoords,
			     positive newyCoords,
			     positive newzCoords );

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

  friend class FGrid3DCurvilinear;

  positive getNbXCoords() const;
  positive getNbYCoords() const;
  positive getNbZCoords() const;

private:
  
  positive xCoords, yCoords, zCoords;

};

#endif // __FPositionSet2DCurvilinear_hh
