//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellLocator.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:03 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellLocator_hh
#define __FCellLocator_hh

#include <vector>
#include <boost/shared_ptr.hpp>

#include "FPosition.hh"
#include "FObject.hh"

using namespace boost;

class FIndex;
class FPositionSet;

/**
 * FCellLocator is an abstract basis class for every CellLocators.
 * Those are designed to return one or more cells that (might) contain 
 * a given position, corresponding to a particular combination of
 * PositionSet and connectivity.
 */

class FCellLocator : public FObject
{
public:
  
  /**
   *\par Description:
   *Constructor: provides an empty FCellLocator pointing to the
   *given PositionSet.
   */
  FCellLocator( const FPositionSet* positionSet );

  /**
   *\par Description:
   *Destructor, has to be overloaded by subclasses
   */
  virtual ~FCellLocator();

  
  //#####################################################################
  //  USER INTERFACE : REQUESTING CELL INDICES FOR POSITIONS
  //#####################################################################
  
  /**
   *\par Description:
   *Returns the index of the (a) cell containing a given position.
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
  virtual 
  bool searchCell( FIndex& aIndex, const FPosition& aPosition ) const = 0;

  /** 
   *\par Description:
   *
   *\pre
   *
   *\post
   *
   *\exception
   *
   *\param
   *
   */
//  virtual void save (std::ostream& out) const =0;
  
  /**
   *\par Description:
   *
   *\pre
   *
   *\post
   *
   *\exception
   *
   *\param
   *
   */
//  virtual void load (std::istream& in)=0;

  
  /**
   *\par Description:
   *
   *\pre
   *
   *\post
   *
   *\exception
   *
   *\param
   *
   */
  virtual void info( std::ostream& tmp ) const;

  /*
   * \return
   * approximate size in bytes
   */
  virtual positive memSize() const = 0;
  
protected:
  
  positive nbCells;
  const FPositionSet* pSet;
};

#endif // __FCellLocator_hh
