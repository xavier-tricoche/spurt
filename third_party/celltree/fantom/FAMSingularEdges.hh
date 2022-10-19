//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularEdges.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:42 $
// Author:    $Author: garth $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMSingularEdges_hh
#define __FAMSingularEdges_hh
 
#include "FAMElement.hh"
#include "FAMSingularEdge.hh"

#include <vector>

/**
 *  This class represents edges of 2Din3D cells that are considered
 *  singular due to the discontinuity of the std::vectorfield resulting from
 *  the projection of the three dimensional std::vector field on the polyhedral
 *  boundary.
 */
class FAMSingularEdges : public FAMElement
{
public:
  /**
   * \par Description:
   *  Standard constructor. Emtpy object.
   */
  FAMSingularEdges();

  /**
   * \par Description:
   *  Destructor. Cleans up.
   */
  ~FAMSingularEdges();
  
  /**
   * \par Description:
   *  Relicts from the inheritance of FObject.
   * \return
   *  "FAMSingularEdge"
   */
  virtual const FString& getClassName() const;

  
  /**
   * \par Description:
   *  Gets you the stored edges.
   */
  void getSingularEdges(std::vector<FAMSingularEdge>& outEdge) const;  

  /**
   * \par Description:
   *  Stores the edges you provide in inEdges.
   */
  void setSingularEdges(const std::vector<FAMSingularEdge>& inEdge);

  /**
   * \par Description:
   *  Inherited from FAMElement.
   */
  virtual void save( const FString& fileNameBase );

  /**
   * \par Description:
   *  Inherited from FAMElement.
   */
  virtual void load( const FString& fileNameBase );
  
private:
  std::vector<FAMSingularEdge> edges;
};

//===========================================================================

#endif // __FAMSingularEdges_hh
