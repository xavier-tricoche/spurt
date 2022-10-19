//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMCycles.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:41 $
// Author:    $Author: garth $
// Version:   $Revision: 1.12 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMCycles_hh
#define __FAMCycles_hh
 
#include "FString.hh"
#include "FAMElement.hh"
#include "FAMCycle.hh"

#include <vector>

/** 
 * Container for FAMCycle objects. Nothing spectacular.
 */
class FAMCycles : public FAMElement
{
public:
  /// Undocumented.
  FAMCycles();
  /// Undocumented.
  ~FAMCycles();
  
  /// Undocumented.
  virtual const FString& getClassName() const;

  /// Undocumented.
  virtual void save( const FString& fileNameBase );
  /// Undocumented.
  virtual void load( const FString& fileNameBase );

  /// Undocumented.
  void getCycles(std::vector<FAMCycle>& outCyc) const;
  
  /// Undocumented.
  void setCycles(const std::vector<FAMCycle>& inCyc);

  /// Undocumented.
  void addCycle(const FAMCycle& inSing);
  /// Undocumented.
  void addCycles(const std::vector<FAMCycle>& inCyc);
  

private:
  std::vector<FAMCycle> points;
};

//===========================================================================

#endif // __FAMCycles_hh
