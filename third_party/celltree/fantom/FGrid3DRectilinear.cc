//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid3DRectilinear.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FGrid3DRectilinear.hh"
#include "FCellLocator3DRectilinear.hh"
#include "FCellDefinitions3DStructured.hh"

#include <cassert>

using namespace std;

// -------------------------------------------------------------------------

FGrid3DRectilinear::FGrid3DRectilinear( shared_ptr<FPositionSet> posSetPtr,
					shared_ptr<FCellDefinitions> cellDefPtr,
					const std::string& newname ) 
    : FGrid( posSetPtr, cellDefPtr, newname )
{
    const FCellDefinitions3DStructured *cdS;
    
    cdS = dynamic_cast<const FCellDefinitions3DStructured*>( cellDef.get() );
    assert( cdS );

    locator = shared_ptr<FCellLocator>( new FCellLocator3DRectilinear( posSetPtr.get(), 
					cdS->getTriangulated() ) );
}

// --------------------------------------------------------------------------

FGrid3DRectilinear::~FGrid3DRectilinear()
{
    // locator freed implicitly via shared ptr
}

// --------------------------------------------------------------------------

const FString& FGrid3DRectilinear::getClassName() const
{
  static FString name( "FGrid3DRectilinear" );
  return name;
}
