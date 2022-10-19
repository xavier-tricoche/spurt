//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMGridBoundary.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:41 $
// Author:    $Author: garth $
// Version:   $Revision: 1.8 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMGridBoundary_hh
#define __FAMGridBoundary_hh
 
#include "FAMElement.hh"
#include "FIndex.hh"

#include <vector>
#include <map>  // for std::pair

/** 
 * Holds the index of a tensor field that represents the boundary.
 * You should set this property of a tensor field to reference another
 * tensor field which represents the (extracted) boundary.
 *
 * The original feature is the possibility to model the relation between
 * the indices of the initial grids 3D cells and the 2Din3D cells of the
 * boundary. This might come handy if you need to process the 3D field
 * depending on the information extracted from the boundary.
 */
class FAMGridBoundary : public FAMElement
{
public:
  /**
   * \par Description:
   *   Default Constructor.
   */
  FAMGridBoundary();
  /**
   * \par Description:
   *   Destructor.
   */
  ~FAMGridBoundary();
  
  /**
   * \par Description:
   *   Returns "FAMGridBoundary".
   * \return
   *   "FAMGridBoundary.
   */
  virtual const FString& getClassName() const;

  /**
   * \par Description:
   *   Save the objects state to one or more files using the given name base.
   */
  virtual void save( const FString& fileNameBase );

  /**
   * \par Description:
   *   Reads the objects state from one or more files using the given
   *   name base.
   */
  virtual void load( const FString& fileNameBase );

  /**
   * \par Description:
   *   Returns the index of the tensor field holding the boundary.
   * \return
   *   The index of the corresponding tensor field.
   */
  FIndex getGridBoundary( void ) const;

  /**
   * \par Description:
   *   Returns a structure containing the std::mapping of 3D to 2Din3D cells, in
   *   other words - which boundary cell belongs to which 3D cell?
   * \param relations
   *   A std::list of std::pairs of cellindices, the first index being the 3D cell,
   *   the second the index of the 2Din3D cell in the boundary field.
   */
  void getGridBoundaryCellRelations( std::vector< std::pair< FIndex, FIndex > > &relations ) const;


  /**
   * \par Description:
   *   Sets the index of the tensor field holding the boundary.
   * \param inBoundaryField
   *   The index.
   */
  void setGridBoundary( const FIndex& inBoundaryField );

  /**
   * \par Description:
   *   Sets the index of the tensor field holding the boundary.
   * \param relations
   *   A std::list of std::pairs of cellindices, the first index being the 3D cell,
   *   the second the index of the 2Din3D cell in the boundary field.
   */
  void setGridBoundaryCellRelations( const std::vector< std::pair< FIndex, FIndex > > &relations );

  /**
   * \par Description:
   *   Are there 3D to 2Din3D relations defined ?
   * \return
   *   True if, false if not.
   */
  bool hasRelations() const;
  
private:
  FIndex boundaryField; //! The Index of the tensor field holding the boundary
                        //! of the tensor field this analysis data belongs to.

  std::vector<std::pair <FIndex, FIndex> > relations; //! Mapping of indices of the
                                            //! 3D cells to 2Din3D cells.
};

//===========================================================================

#endif // __FAMGridBoundary_hh
