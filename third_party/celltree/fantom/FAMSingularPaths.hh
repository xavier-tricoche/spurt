//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularPaths.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:43 $
// Author:    $Author: garth $
// Version:   $Revision: 1.8 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMSingularPaths_hh
#define __FAMSingularPaths_hh
 
#include "FAMElement.hh"
#include "FAMSingularPath.hh"
#include <vector>

/** 
 * Undocumented ... bobach
 */
class FAMSingularPaths : public FAMElement
{
public:
  /// Undocumented.
  FAMSingularPaths();
  /// Undocumented.
  ~FAMSingularPaths();
  
  /// Undocumented.
  const FString& getClassName() const;

  /// Undocumented.
  virtual void save( const FString& fileNameBase );
  /// Undocumented.
  virtual void load( const FString& fileNameBase );

  /// Undocumented.
  void getSingularPaths(std::vector<FAMSingularPath>& outPath) const; 
  /// Undocumented. 
  void setSingularPaths(const std::vector<FAMSingularPath>& inPath);

  /// Undocumented.
  void addSingularPath(const FAMSingularPath& inPath);
  /// Undocumented.
  void addSingularPaths(const std::vector<FAMSingularPath>& inPath);
  
private:
  std::vector<FAMSingularPath> paths;
};

//===========================================================================

#endif // __FAMSingularPaths_hh
