//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMBifurcations.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:40 $
// Author:    $Author: garth $
// Version:   $Revision: 1.11 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMBifurcations_hh
#define __FAMBifurcations_hh
 
#include "FString.hh"
#include "FAMElement.hh"
#include "FAMBifurcation.hh"
#include <vector>

/** 
 * Container for FAMBifurcation objects. Nothing special.
 */
class FAMBifurcations : public FAMElement
{
public:
  /// Undocumented.
  FAMBifurcations();
  /// Undocumented.
  ~FAMBifurcations();
  
  /// Undocumented.
  const FString& getClassName() const;

  /// Undocumented.
  virtual void save( const FString& fileNameBase );
  /// Undocumented.
  virtual void load( const FString& fileNameBase );
  
  /// Undocumented.
  void getBifurcations(std::vector<FAMBifurcation>& outBif) const;  
  /// Undocumented.
  void setBifurcations(const std::vector<FAMBifurcation>& inBif);

  /// Undocumented.
  void addBifurcation(const FAMBifurcation& inBif);
  /// Undocumented.
  void addBifurcations(std::vector<FAMBifurcation>& inBif);
  
private:
  std::vector<FAMBifurcation> bifs;
};

//===========================================================================

#endif // __FAMBifurcations_hh
