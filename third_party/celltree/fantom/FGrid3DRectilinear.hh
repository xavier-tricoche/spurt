//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid3DRectilinear.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGrid3DRectilinear_hh
#define __FGrid3DRectilinear_hh

#include "FGrid.hh"

class FGrid3DRectilinear: public FGrid
{
public:

  FGrid3DRectilinear( shared_ptr<FPositionSet> posSetPtr,
		      shared_ptr<FCellDefinitions> cellDefPtr,
		      const std::string& newname );

  ~FGrid3DRectilinear();

  const FString& getClassName() const;
};

#endif
