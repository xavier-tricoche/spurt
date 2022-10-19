//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid3DCurvilinear.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 


#include "FGrid3DCurvilinear.hh"
#include "FCellLocatorNeighbors.hh"

//--------------------------------------------------------------------------- 

FGrid3DCurvilinear::FGrid3DCurvilinear( shared_ptr<FPositionSet> posSetPtr,
					shared_ptr<FCellDefinitions> cellDefPtr,
					const std::string& newname ) 
  : FGrid( posSetPtr, cellDefPtr, newname )
{
    locator = shared_ptr<FCellLocator>( new FCellLocatorNeighbors(this) );
}

//--------------------------------------------------------------------------- 

FGrid3DCurvilinear::~FGrid3DCurvilinear()
{
    // locator freed implicitly via shared ptr
}

const FString& FGrid3DCurvilinear::getClassName() const
{
  static FString name( "FGrid3DCurvilinear" );
  return name;
}
