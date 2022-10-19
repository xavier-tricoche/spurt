//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid2Din3DArbitrary.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGrid2Din3DArbitrary_hh
#define __FGrid2Din3DArbitrary_hh

#include "FGrid.hh"

class FGrid2Din3DArbitrary : public FGrid
{
public:

  FGrid2Din3DArbitrary( shared_ptr<FPositionSet> posSetPtr,
			shared_ptr<FCellDefinitions> cellDef,
			const std::string& newname );
  
  ~FGrid2Din3DArbitrary();

  const FString& getClassName() const;

  /**
   *   not implemented, will always return false
   */
  FIndex searchCellIndex( const FPosition& aPosition ) const;
};

#endif

