//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularPoints.hh,v $
// Language:  C++
// Date:      $Date:  $
// Author:    $Author: jaenicke $
// Version:   $Revision:  $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMSingularPoints_hh
#define __FAMSingularPoints_hh

#include "FObject.hh"
#include "FAMElement.hh"
#include "FAMSingularPoint.hh"
#include <boost/shared_ptr.hpp>
  

class FAMSingularPoints : public FAMElement
{
  public:
  /** 
   *{\bf Description:}\\
   *Constructor: returns a FAMSingularPoints with given points.
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   */
  FAMSingularPoints (std::vector<FAMSingularPoint>& inSingularPoints);

  /** 
   *{\bf Description:}\\
   *Constructor: returns empty FAMSingularPoints object.
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   */
  FAMSingularPoints ();

  /** 
   *{\bf Destructor:}\\
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   */
  ~FAMSingularPoints ();
  

  /**
   * Undocumented 
   */
  const FString& getClassName () const;

  /** 
   *{\bf Description:}\\
   *returns the vector of the FAMSingularPoint.
   *\\{\bf Precondition:}\\
   *a vector of points has been set.
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *FEmptyObjectException
   *\param
   *{\bf result:} returned position of the singularity.
   */
  void getSingularPoints (std::vector<FAMSingularPoint>& outSingularPoints) const;

  /**
   * Undocumented 
   */
  void setSingularPoints (std::vector<FAMSingularPoint>& inSingularPoints);

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
  std::vector<FAMSingularPoint> points;

};

#endif // __FAMSingularPoints_hh
