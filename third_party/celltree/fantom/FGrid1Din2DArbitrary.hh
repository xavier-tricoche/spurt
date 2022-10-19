//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid1Din2DArbitrary.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGrid1Din2DArbitrary_hh
#define __FGrid1Din2DArbitrary_hh

#include "FGrid.hh"

class FGrid1Din2DArbitrary : public FGrid
{
public:

  FGrid1Din2DArbitrary( shared_ptr<FPositionSet> posSetPtr,
			shared_ptr<FCellDefinitions> cellDef,
			const std::string& newname );
  
  virtual ~FGrid1Din2DArbitrary();
  

  virtual const FString& getClassName() const;
  
  /**
   *   not implemented, will always return false
   */
  virtual FIndex searchCellIndex( const FPosition& aPosition ) const;
};

#endif

