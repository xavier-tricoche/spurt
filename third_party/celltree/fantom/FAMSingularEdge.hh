//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularEdge.hh,v $
// Language:  C++
// Date:      $Date: 2003/01/09 17:30:37 $
// Author:    $Author: bobach $
// Version:   $Revision: 1.4 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMSingularEdge_hh
#define __FAMSingularEdge_hh

#include "FObject.hh"
#include "FIndex.hh"
#include "FPosition.hh"

/** 
 * Class to handle the singular edges that may be encountered on polyhedral
 * surfaces.
 */
class FAMSingularEdge : public FObject
{
public:

  //===========================================================================
  
  enum edgeType { SINK, SOURCE, NONE };  //! Main type enumeration

  //===========================================================================

  
  /** 
   *\par Description:
   *Constructor: returns an empty FAMSingularEdge.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  FAMSingularEdge();

  /** 
   *\par Description:
   *copy constructor
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *{\bf edge}: singularity to copy
   */
  FAMSingularEdge(const FAMSingularEdge& edge);

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
  FAMSingularEdge(const FPosition& pos1, const FPosition& pos2,
		  const FIndex& cellId1, const FIndex& cellId2,
		  const edgeType& type);  

  // undocumented
  void getEdgeVertices( FPosition& pos1, FPosition& pos2 );
  void getEdgeType( edgeType& type );
  void getEdgeCells( FIndex& cell1, FIndex& cell2 );

  virtual const FString& getClassName() const;
private:

  FPosition pos1, pos2;


  edgeType type;  //! singular point type


  FIndex cellId1, cellId2;  // indexes of the containing cells
};

#endif // __FAMSingularEdge_hh
