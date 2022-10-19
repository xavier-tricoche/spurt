//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGridStructureDescription.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:06 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGridStructureDescription_hh
#define __FGridStructureDescription_hh

#include "FObject.hh"
#include "FCell.hh"
#include <vector>

class FIndex;

/** 
 * This class is designed to hold information about the structure
 * of a grid, and compare the structuredness to propose a way
 * to efficiently use the structure
 */
class FGridStructureDescription
{

  /** 
   *\par Description:
   *default Constructor
   */
  FGridStructureDescription();
  
  /** 
   *\par Description:
   *Constructs a cell from the given cellType
   *\pre
   *myType must be a valid cell Type
   *\post
   *none
   *\exception
   *none
   *\param
   *myType: type of the cell we want to be created
   *\return
   *A pointer to a newly instantiated cell of type aType is 
   */
  static FCell* createCell(FCell::CellType myType);

  /** 
   *\par Description:
   *Returns the amount of FIndex values to fully describe a cell of
   *type myType.
   *\param
   *myType: a type from cellType
   *\return
   * number of needed FIndex.
   */
  static positive sizeOfCellType( FCell::CellType myType );

  /** 
   *\par Description:
   *Returns the dimension of the space that this cell lives in
   *\param
   *myType: a type from cellType
   *\return
   * dimension.
   */
  static positive dimensionOfCellType( FCell::CellType myType );

  /** 
   *\par Description:
   *Returns the type of the cell that a tensorproduct of
   * A and B produces.
   *\param
   *A: a type from cellType
   *\b B: a type from cellType
   *\return
   * tensorproduct cell of A and B.
   */
  static FCell::CellType cellTensorProduct( FCell::CellType A, FCell::CellType B );
  

};

#endif // __FGridStructureDescription_hh
