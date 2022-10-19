//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSeparatrix.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:30 $
// Author:    $Author: garth $
// Version:   $Revision: 1.9 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMSeparatrix_hh
#define __FAMSeparatrix_hh

#include "FString.hh"
#include "FPosition.hh"
#include "FIndex.hh"

#include "FException.hh"

#include <list>
#include <ostream>

#include "FVector.hh"
#include "FObject.hh"

class FMatrix;
class FCell;

//===========================================================================
/** 
 * Class to handle the topological connections (separatrices) between 
 * singular points.
 */
class FAMSeparatrix : public FObject
{
public:

 
  /** 
   *{\bf Description:}\\
   *Constructor: returns an empty FAMSeparatrix.
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   */
  FAMSeparatrix();

  /** 
   *{\bf Description:}\\
   *copy constructor
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *\param
   *{\bf sing}: separatrix to copy
   */
  FAMSeparatrix(const FAMSeparatrix& sep);

  /** 
   *{\bf Description:}\\
   *Constructor: returns a separatrix connecting the given saddle point
   *with the given source or sink (reachedId) and sets the corresponding
   *curve length.
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *@param
   *{\bf saddleId}: index of a saddle point (separatrix' origin)
   *@param
   *{\bf reachedId}: index of a sink or source (separatrix' end)
   *@param
   *{\bf length}: separatrix' length
   */
  FAMSeparatrix(const FIndex& saddleId, const FIndex& reachedId, double length); 

  double getLength() const;
  void getSaddleIndex(FIndex& saddleId) const;
  void getSinkSourceIndex(FIndex& sinksourceId) const;

  /** 
   *{\bf Description:}\\
   *Destructor.
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   */
  virtual ~FAMSeparatrix();

  /**
   *{\bf Description:}\\
   *Returns the class name.
   */
  virtual const FString& getClassName() const;

  void setLastPoint(const FPosition& last);
  void getLastPoint(FPosition& last) const;

  /** 
   *{\bf Description:}\\
   *prints the contents of the separatrix
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *\param
   *{\bf singularity}: separatrix to print.
   */
  friend std::ostream& operator<< (std::ostream &os, const FAMSeparatrix& sep);

private:
  
  FIndex saddle, sinksource;
  double length;
  FPosition last;
};

#endif // __FAMSeparatrix_hh
