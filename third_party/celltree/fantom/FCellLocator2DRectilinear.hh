//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellLocator2DRectilinear.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:03 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellLocator2DRectilinear_hh
#define __FCellLocator2DRectilinear_hh

#include "FCellLocator.hh"

/**
 */

class FCellLocator2DRectilinear : public FCellLocator
{
public:
  
  /**
   *\par Description:
   *Constructor: provides an empty FCellLocator pointing to the
   *given PositionSet.
   */
  FCellLocator2DRectilinear( const FPositionSet* positionSet,
			     bool triangulated=false);

  
  /**
   *\par Description:
   *Destructor overloaded from FCellLocator
   */
  ~FCellLocator2DRectilinear();

  //#####################################################################
  //  USER INTERFACE : REQUESTING CELL INDICES FOR POSITIONS
  //#####################################################################
  
  /**
   *\par Description:
   *Returns the index of the cell containing a given position.
   *Note that this might be impossible in some cases.
   *\pre
   *none
   *\post
   *aIndex contains the id of the found cell if any
   *\exception
   *none
   *\param
   *aPosition: the input position
   *\param
   *aIndex: the index to set to the found cell index, will be invalid
   *              if no cell containing the position was found.
   *\return
   *\c true if a cell was found
   *\c false if no cell was found containing the position
   */

  bool searchCell( FIndex& aIndex, const FPosition& aPosition ) const;


  positive memSize() const;
  
protected:

  bool isTriangulated;

};

#endif // __FCellLocator2DRectilinear_hh