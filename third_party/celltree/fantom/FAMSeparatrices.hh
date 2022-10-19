//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSeparatrices.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:41 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMSeparatrices_hh
#define __FAMSeparatrices_hh
 
#include "FAMElement.hh"
#include "FAMSeparatrix.hh"

#include <vector>

/** 
 * Undocumented ... tricoche
 */
class FAMSeparatrices : public FAMElement
{
public:
  FAMSeparatrices();
  ~FAMSeparatrices();
  
  virtual const FString& getClassName() const;

  virtual void save( const FString& fileNameBase );
  virtual void load( const FString& fileNameBase );

  void getSeparatrices(std::vector<FAMSeparatrix>& outSep) const;  
  void setSeparatrices(const std::vector<FAMSeparatrix>& inSep);
  
  void addSeparatrix(const FAMSeparatrix& inSep);
  void addSeparatrices(std::vector<FAMSeparatrix>& inSing);
  
private:
  std::vector<FAMSeparatrix> seps;
};

//===========================================================================

#endif // __FAMSeparatrices_hh
