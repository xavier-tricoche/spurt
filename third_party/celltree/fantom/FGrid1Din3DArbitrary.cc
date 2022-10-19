//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid1Din3DArbitrary.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FGrid1Din3DArbitrary.hh"
#include "FIndex.hh"

#include <iostream>

//--------------------------------------------------------------------------- 

FGrid1Din3DArbitrary::FGrid1Din3DArbitrary( shared_ptr<FPositionSet> posSetPtr,
					    shared_ptr<FCellDefinitions> cellDefPtr,
					    const std::string& newname ) 
  : FGrid( posSetPtr, cellDefPtr, newname )
{
}

//--------------------------------------------------------------------------- 
					       
FGrid1Din3DArbitrary::~FGrid1Din3DArbitrary()
{
    // locator freed implicitly via shared ptr
}

//--------------------------------------------------------------------------- 

FIndex FGrid1Din3DArbitrary::searchCellIndex( const FPosition& /*aPosition*/ ) const
{
    cerr << "searchCell() is not supported by FGrid1Din3D!" << endl;

    // functionality not implemented/possible for this grid type
    // -> return invalid index
    return FIndex();
}

const FString& FGrid1Din3DArbitrary::getClassName() const
{
  static FString name( "FGrid1Din3DArbitrary" );
  return name;
}
