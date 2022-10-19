//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid2DRectilinear.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGrid2DRectilinear_hh
#define __FGrid2DRectilinear_hh

#include "FGrid.hh"
#include "FPositionSet.hh"
#include "FPositionSet2DRectilinear.hh"
#include "FCellLocator2DRectilinear.hh"
#include <vector>

class FCellDefinitions2DStructured;

class FGrid2DRectilinear: public FGrid
{
public:

  FGrid2DRectilinear( shared_ptr<FPositionSet> posSetPtr,
		      shared_ptr<FCellDefinitions> cellDefPtr,
		      const std::string& newname );

  const FString& getClassName() const;
  
  ~FGrid2DRectilinear();
};

#endif
