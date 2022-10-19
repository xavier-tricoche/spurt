//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid2DRectilinear.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FGrid2DRectilinear.hh"
#include "FCellLocator2DRectilinear.hh"
#include "FCellDefinitions2DStructured.hh"

#include <cassert>

using namespace std;

// -------------------------------------------------------------------------

FGrid2DRectilinear::FGrid2DRectilinear( shared_ptr<FPositionSet> posSetPtr,
					shared_ptr<FCellDefinitions> cellDefPtr,
					const std::string& newname ) 
    : FGrid( posSetPtr, cellDefPtr, newname )
{
    const FCellDefinitions2DStructured *cdS;

    cdS = dynamic_cast<const FCellDefinitions2DStructured*>( cellDef.get() );
    assert( cdS );

    locator = shared_ptr<FCellLocator>( new FCellLocator2DRectilinear( posSetPtr.get(), 
					cdS->getTriangulated() ) );
}

// --------------------------------------------------------------------------

FGrid2DRectilinear::~FGrid2DRectilinear()
{
    // locator freed implicitly via shared ptr
}

// --------------------------------------------------------------------------

const FString& FGrid2DRectilinear::getClassName() const
{
  static FString name( "FGrid2DRectilinear" );
  return name;
}
