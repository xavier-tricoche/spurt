//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid2DArbitrary.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:36:59 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGrid2DArbitrary_hh
#define __FGrid2DArbitrary_hh

#include "FGrid.hh"

class FGrid2DArbitrary : public FGrid
{
public:
  FGrid2DArbitrary();

  FGrid2DArbitrary( shared_ptr<FPositionSet> posSetPtr,
		    shared_ptr<FCellDefinitions> cellDef,
		    const std::string& newname );

  const FString& getClassName() const;
  ~FGrid2DArbitrary();
};

#endif
