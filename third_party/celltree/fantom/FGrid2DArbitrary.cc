//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid2DArbitrary.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:36:59 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 


#include "FGrid2DArbitrary.hh"
#include "FCellLocatorNeighbors.hh"

using namespace std;

//--------------------------------------------------------------------------- 

FGrid2DArbitrary::FGrid2DArbitrary( shared_ptr<FPositionSet> posSetPtr,
				    shared_ptr<FCellDefinitions> cellDefPtr,
				    const std::string& newname ) 
    : FGrid( posSetPtr, cellDefPtr, newname )
{
    locator = shared_ptr<FCellLocator>( new FCellLocatorNeighbors( this ) ); 
}

//--------------------------------------------------------------------------- 

FGrid2DArbitrary::~FGrid2DArbitrary()
{
    // locator freed implicitly via shared ptr
}

//--------------------------------------------------------------------------- 

const FString& FGrid2DArbitrary::getClassName() const
{
  static FString name( "FGrid2DArbitrary" );
  return name;
}
