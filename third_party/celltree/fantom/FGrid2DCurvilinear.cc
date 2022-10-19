//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid2DCurvilinear.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 


#include "FGrid2DCurvilinear.hh"
#include "FCellLocatorNeighbors.hh"

using namespace std;

//--------------------------------------------------------------------------- 

FGrid2DCurvilinear::FGrid2DCurvilinear( shared_ptr<FPositionSet> posSetPtr,
					shared_ptr<FCellDefinitions> cellDefPtr,
					const std::string& newname ) 
    : FGrid( posSetPtr, cellDefPtr, newname )
{
    locator = shared_ptr<FCellLocator>( new FCellLocatorNeighbors(this) );
}

//--------------------------------------------------------------------------- 

FGrid2DCurvilinear::~FGrid2DCurvilinear()
{
    // locator freed implicitly via shared ptr
}

//--------------------------------------------------------------------------- 

const FString& FGrid2DCurvilinear::getClassName() const
{
  static FString name( "FGrid2DCurvilinear" );
  return name;
}
