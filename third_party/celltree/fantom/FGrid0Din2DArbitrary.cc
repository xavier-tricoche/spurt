//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid0Din2DArbitrary.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FGrid0Din2DArbitrary.hh"
#include "FIndex.hh"

#include <iostream>

//--------------------------------------------------------------------------- 

FGrid0Din2DArbitrary::FGrid0Din2DArbitrary( shared_ptr<FPositionSet> posSetPtr,
					    shared_ptr<FCellDefinitions> cellDefPtr,
					    const std::string& newname ) 
  : FGrid( posSetPtr, cellDefPtr, newname )
{
}

//--------------------------------------------------------------------------- 
					       
FGrid0Din2DArbitrary::~FGrid0Din2DArbitrary()
{
    // locator freed implicitly via shared ptr
}

//--------------------------------------------------------------------------- 

FIndex FGrid0Din2DArbitrary::searchCellIndex( const FPosition& /*aPosition*/ ) const
{
    cerr << "searchCell() is not supported by FGrid0Din2D!" << endl;

    // functionality not implemented/possible for this grid type
    // -> return invalid index
    return FIndex();
}

const FString& FGrid0Din2DArbitrary::getClassName() const
{
  static FString name( "FGrid0Din2DArbitrary" );
  return name;
}
