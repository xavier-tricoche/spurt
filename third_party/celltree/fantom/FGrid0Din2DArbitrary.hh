//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid0Din2DArbitrary.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGrid0Din2DArbitrary_hh
#define __FGrid0Din2DArbitrary_hh

#include "FGrid.hh"

class FGrid0Din2DArbitrary : public FGrid
{
public:

  FGrid0Din2DArbitrary( shared_ptr<FPositionSet> posSetPtr,
			shared_ptr<FCellDefinitions> cellDef,
			const std::string& newname );
  
  ~FGrid0Din2DArbitrary();

  virtual const FString& getClassName() const;
  /**
   *   not implemented, will always return false
   */
  virtual FIndex searchCellIndex( const FPosition& aPosition ) const;
};

#endif

