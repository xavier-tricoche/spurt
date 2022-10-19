//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid3DCurvilinear.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGrid3DCurvilinear_hh
#define __FGrid3DCurvilinear_hh

#include "FGrid.hh"

class FGrid3DCurvilinear : public FGrid
{
public:
  FGrid3DCurvilinear();

  FGrid3DCurvilinear( shared_ptr<FPositionSet> posSetPtr,
		      shared_ptr<FCellDefinitions> cellDefPtr,
		      const std::string& newname );

  const FString& getClassName() const;

  ~FGrid3DCurvilinear();
};

#endif
