//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid2DCurvilinear.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGrid2DCurvilinear_hh
#define __FGrid2DCurvilinear_hh

#include "FGrid.hh"

class FGrid2DCurvilinear : public FGrid
{
public:

  FGrid2DCurvilinear(  shared_ptr<FPositionSet> posSetPtr,
		       shared_ptr<FCellDefinitions> cellDefPtr,
		       const std::string& newname );

  virtual const FString& getClassName() const;
  virtual ~FGrid2DCurvilinear();
};

#endif
