//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:36:59 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGrid_hh
#define __FGrid_hh

class FCellLocator;
class FPositionSet;
class FCell;

#include <string>
#include <boost/shared_ptr.hpp>

#include "FArray.hh"
#include "FCellDefinitions.hh"

class FCellDefinitions;
class FPositionSet;
class FIndex;

using namespace boost;
/**
 * \ingroup DataSet
 *
 * FGrid is an abstract basis class for the management of a grid,
 * that is a geometric connectivity on one hand and the positions
 * of the vertices on the other hand. Therefore, it is in charge of
 *  - mapping positions to cells that include the position (therefore
 *    it owns a cellLocator that it builds according to its type)
 *  - mapping positionIndices to cells (getCellDefinitions()->getNeighborhoodInfo()-> etc.)
 *  - mapping an index to the cell with that index in the particular FGrid
 *  - providing the geometry information for cells
 *  - mapping positionIndices to positions (by getPositionSet()-> etc. )
 *
 * The management of the positions is entirely done by FPositionSet.
 * Specialized FGrid classes must be designed to optimize the memory
 * and access costs for typical or particular grid structures.
 *
 *
 */
class FGrid
{
public:

  /**
   *\par Description:
   *Constructor: provides an FGrid pointing on a particular FPositionSet.
   *FGrid does not own the FPositionSet but is able to consult its constant
   *reference as a friend class. A name can be given.
   *\pre
   *posSetPtr and refers to an existing FPositionSet object  and
   *cellDefPtr refers to an existing FCellDefinitions object
   *\post
   *none
   *\param
   *posSetPtr: pointer on a FPositionSet
   *\param
   *cellDefPtr: pointer on a set of FCellDefinitions ( connectivity ) 
   *\param newname:
   * name of the grid
   */
  FGrid( shared_ptr<FPositionSet> posSetPtr,
	 shared_ptr<FCellDefinitions> cellDefPtr,
	 const std::string& aName );
  
  /**
   *\par Description:
   *Destructor
   */
  virtual ~FGrid();


  //#####################################################################
  //  USER INTERFACE : REQUESTING INFORMATION ON THE GRID
  //#####################################################################




  
  const std::string& getName() const;  

  void setName(const std::string& newname);

  virtual const FString& getClassName() const;

  /**
   *\par Description:
   *Return a pointer on the grid cell containing the given position, if any.
   *\pre
   *The grid is not empty
   *\post
   *none
   *\exception
   *FEmptyObjectException
   *\param
   *aPosition: the request position
   *\return
   *\c searched cell or
   *\c NULL pointer if the position lies outside the grid
   */
  virtual shared_ptr<FCell> searchCell( const FPosition& aPosition ) const;

  /**
   *\par Description:
   *Returns the index of the grid cell containing the given position, if any.
   *\pre
   *The grid is not empty
   *\post
   *none
   *\exception
   *FEmptyObjectException
   *\param
   *aPosition: the request position
   *\return
   *\c searched cell or
   *\c NULL pointer if the position lies outside the grid
   */
  virtual FIndex searchCellIndex( const FPosition& aPosition ) const;

  /**
   *\par Description:
   *Return the FCell with given index in the grid, if any
   *\pre
   *The grid is not empty.
   *0 <= cellId < number of cells in the grid
   *\post
   *none
   *\exception
   *FEmptyObjectException
   *\exception
   *FIndexOutOfBoundsException
   *\param
   *cellId the index of the requested cell
   *\return
   *\c requested cell
   */
  virtual shared_ptr<FCell> getCell(const FIndex& cellId) const;

  /**
   *\par Description:
   *Return FPositionSet belonging to this grid.
   *\pre
   *\post
   *none
   *\exception
   *\param
   *\return
   *\c position set 
   */
  shared_ptr<FPositionSet> getPositionSet() const;

  /**
   *\par Description:
   *Return FCellDefinitions module belonging to this grid.
   *\pre
   *none
   *\post
   *none
   *\return
   *shared pointer on cell Definitions
   */
  shared_ptr<FCellDefinitions> getCellDefinitions() const;

  /**
   *\par Description:
   *Return FCellLocator module of this grid.
   *(for point location)
   *\pre
   *none
   *\post
   *none
   *\return
   *shared pointer on FCellLocator of this grid.
   */
  shared_ptr<const FCellLocator> getCellLocator() const;

  /**
   *\par Description:
   *gets the cell containing the position newPos
   * by following the line from oldpos to newpos
   * through the neighboring cells
   *
   *\pre
   *pos is Inside cell
   *
   *\post
   *pos is Inside cell
   *\c newPos is Inside cell or pos is place where the 
   *search ray hit the border
   *
   *\param cell
   * actual cell which is updated if possible 
   * to cell containing newPos
   *
   *\param cellId
   * id of cell in Grid
   *
   *\param indicesGone:
   *if not zero, all cell indices looked at are saved in here
   *
   *\param positionsGone:
   *if not zero, all positions gone are saved in here
   *
   *\return
   *\c true if newPos is Inside cell
   *\c false if cell containing newPos couldn't be found
   * and the border was hit.
   */
  bool getNextCell(shared_ptr<FCell> & cell , FIndex &cellId,
		   FPosition & pos, const FPosition& newPos,
		   vector<FIndex> * indicesGone=0,
		   vector<FPosition> * positionsGone=0) const;

  /*
   * 
   * \par Description:
   * Get the grid's size in memory
   * \return
   * Approximate size in bytes
   */
  virtual positive memSize() const;

  

  /*
   * \par Description
   * This function determines the correct type of grid that fits 
   * to the given FCellDefinitions and FPositionSet, creates an 
   * appropiated object and returns a pointer to the object.
   */
  static shared_ptr<FGrid> constructSuitableGrid( const std::string& name,
						  shared_ptr<FCellDefinitions> cellDef,
						  shared_ptr<FPositionSet> posSet );


protected:

  // position set
  shared_ptr<FPositionSet> posSet;

  // cell definitions module
  shared_ptr<FCellDefinitions> cellDef;

  // cell locator
  shared_ptr<FCellLocator> locator;

  // name of the geometry
  std::string name;
};

#endif // __FGrid_hh
