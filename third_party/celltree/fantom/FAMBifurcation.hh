//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMBifurcation.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:40 $
// Author:    $Author: garth $
// Version:   $Revision: 1.12 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMBifurcation_hh
#define __FAMBifurcation_hh

#include "FString.hh"
#include "FPosition.hh"
#include "FIndex.hh"

#include "FException.hh"
#include "FObject.hh"
#include "FVector.hh"
class FMatrix;
class FCell;


/** 
 *Structure to handle type and location of a bifurcation in a 2D unsteady
 *std::vector field.
 */
class FAMBifurcation : public FObject
{
public:

  //===========================================================================
  
  /**
   * Undocumented 
   */
  enum bifurcationType { CREATION, ANNIHILATION, HOPF, SWAP, 
			 WEDGE_SWAP, NONE_BIF };
  
//===========================================================================
  
  /** 
   *\par Description:
   *Constructor: 
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param pos1, pos2
   * edge vertices
   *\param type
   * edge type
   */

  /** 
   *\par Description:
   *Constructor: returns an empty FAMBifurcation.
   */
  FAMBifurcation();

  /** 
   *\par Description:
   *copy constructor
   *\param bif
   * Bifurcation to copy.
   */
  FAMBifurcation(const FAMBifurcation& bif);
  
  /** 
   *\par Description:
   *Constructor: returns a FAMBifurcation located at the given position.
   *\param pos
   * position where bifurcation takes place.
   */
  FAMBifurcation(const FPosition& _pos);

  /** 
   *\par Description:
   * Constructor: returns a FAMBifurcation located at the given position with 
   *the given type.
   *\param pos
   * Position where bifucation takes place.
   *\param type
   * Bifurcation type.
   *\param cellId
   * Containing cell index.
   */
  FAMBifurcation(const FPosition& _pos, const bifurcationType& _type,
		 const FIndex& _cellId);  

  /** 
   *\par Description:
   * Returns the position of the bifurcation.
   *\pre
   * A position has been set.
   *\post
   * none
   *\exception
   * FEmptyObjectException
   *\param result
   * Returned position of the bifurcation.
   */
  void getPosition(FPosition& result) const;

  /** 
   *\par Description:
   * Sets the position of the bifurcation.
   *\pre
   * none.
   *\post
   * The position of the bifurcation is set.
   *\exception
   * FEmptyObjectException
   *\param result
   * Returned position of the bifurcation.
   */
  void setPosition(const FPosition& result);

  /** 
   *\par Description:
   * Returns the type of the Bifurcation.
   *\pre
   *a type has been set.
   *\post
   *none
   *\exception
   *FException: no type has been set.
   *\param result
   * Returned type of the bifurcation.
   */
  void getType(bifurcationType& result) const;

  /** 
   *\par Description:
   * Sets the index of the cell containing bifurcation location
   *\pre
   * none
   *\post
   * The index of the including cell is set.
   *\exception
   * none
   *\param cellId
   * containing cell index
   */
  void setIncludingCellIndex(const FIndex& cellId);

  /** 
   *\par Description:
   * Returns the index of the cell containing bifurcation location
   *\pre
   * A cell index has been set.
   *\post
   * none
   *\exception
   * FException: no cell index has been set.
   *\param cellId
   * Returned cell index
   */
  void getIncludingCellIndex(FIndex& cellId) const;

  /** 
   *\par Description:
   * Sets the indices of the singularities involved in the bifurcation
   * (as singular paths indices).
   *\pre
   *none
   *\post
   * Pathindices set. Awful functionality.
   *\exception
   *none
   *\param singPathId1, singPathId2
   * given indices.
   */
  void setInvolvedSingPathsIndices(const FIndex& singPathId1,
				   const FIndex& singPathId2);

  /** 
   *\par Description:
   * Gets the indices of the singularities involved in the bifurcation
   * (as singular paths indices).
   *\pre
   * The indices have been set before.
   *\post
   * none
   *\param singPathId1, singPathId2
   * Returned indices.
   */
  void getInvolvedSingPathsIndices(FIndex& singPathId1,
				   FIndex& singPathId2);

  /** 
   *\par Description:
   *Destructor.
   */
  ~FAMBifurcation();

  /**
   *\par Description:
   *Returns the class name.
   */
  virtual const FString& getClassName() const;

  /**
   *\par Description:
   * The assignment operator. Function should be clear.
   */
  FAMBifurcation& operator=(const FAMBifurcation& bif);


  /** 
   *\par Description:
   * Prints the contents of the bifurcation, or merely numbers conveying
   * the content.
   *\param bif
   * bifurcation to print
   */
  friend std::ostream& operator<< (std::ostream &os,
				   const FAMBifurcation& bif);


private:
  
  FPosition pos;   //! bifurcation location

  bifurcationType type;   //! bifurcation type 

  FIndex singPathId1, singPathId2; //! involved singularities (as singular paths indices)

  FIndex cellId; //! index of the containing cell
};


/**
 *\par Description:
 * just transforms a bifurcationType into Text.
 */
std::ostream& operator<< (std::ostream &os,
		     const FAMBifurcation::bifurcationType& type);



//===========================================================================
#ifndef OUTLINE
#include "FAMBifurcation.icc"
#endif
//=========================================================================== 

#endif // __FAMBifurcation_hh
