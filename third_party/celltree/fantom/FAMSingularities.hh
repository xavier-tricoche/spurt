//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularities.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:42 $
// Author:    $Author: garth $
// Version:   $Revision: 1.14 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMSingularities_hh
#define __FAMSingularities_hh
 
#include "FAMElement.hh"
#include "FAMSingularPoint.hh"

#include <vector>
#include <map>

/** 
 * This class holds a set of singularities along with a std::map of cells to the
 * singularities they contain. The interesting information is stored in
 * the FAMSingularPoint class.
 * Singularities (represented by an FAMSingularPoint object) may be added as
 * a whole set or one by one.
 */
class FAMSingularities : public FAMElement
{
public:
  /// Undocumented.
  FAMSingularities();
  /// Undocumented.
  ~FAMSingularities();
  
  /// Undocumented.
  virtual const FString& getClassName() const;

  /// Undocumented.
  virtual void save( const FString& fileNameBase );
  /// Undocumented.
  virtual void load( const FString& fileNameBase );

  /// Undocumented.
  void getSingularities(std::vector<FAMSingularPoint>& outSing) const;  
  /// Undocumented.
  void setSingularities(const std::vector<FAMSingularPoint>& inSing);

  /// Undocumented.
  void getSingularityForCell(FIndex cellId, FAMSingularPoint& sing) const;
  
  void addSingularity(const FAMSingularPoint& inSing);
  /// Undocumented.
  void addSingularities(std::vector<FAMSingularPoint>& inSing);
  
private:
  std::vector<FAMSingularPoint> points;
  std::map<positive, positive> cellToSingularitiesIndex;
};

//===========================================================================

#endif // __FAMSingularities_hh
