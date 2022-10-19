//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid3DArbitrary.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 


#include "FGrid3DArbitrary.hh"
#include "FCellLocatorNeighbors.hh"

using namespace std;

//--------------------------------------------------------------------------- 

FGrid3DArbitrary::FGrid3DArbitrary( shared_ptr<FPositionSet> posSetPtr,
				    shared_ptr<FCellDefinitions> cellDefPtr,
				    const std::string& newname ) 
    : FGrid( posSetPtr, cellDefPtr, newname )
{
}

//--------------------------------------------------------------------------- 

void FGrid3DArbitrary::buildKdTree()
{
    if( !locator )
        locator = shared_ptr<FCellLocator>( new FCellLocatorNeighbors( this ) );
}

//--------------------------------------------------------------------------- 

FGrid3DArbitrary::~FGrid3DArbitrary()
{
    // locator freed implicitly via shared ptr
}

//--------------------------------------------------------------------------- 
	
const FString& FGrid3DArbitrary::getClassName() const
{
  static FString name( "FGrid3DArbitrary" );
  return name;
}				       
