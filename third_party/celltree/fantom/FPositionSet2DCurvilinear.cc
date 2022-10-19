//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet2DCurvilinear.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:09 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 


#include "FPositionSet2DCurvilinear.hh"
#include "FIndex.hh"

#include "FException.hh"

FPositionSet2DCurvilinear::FPositionSet2DCurvilinear( vector<double>& coords,
						      positive newxCoords,
						      positive newyCoords )
    : FPositionSet2DArbitrary( coords )
{
    xCoords = newxCoords;
    yCoords = newyCoords;
}

//--------------------------------------------------------------------------- 

const FString& FPositionSet2DCurvilinear::getClassName() const
{
    static FString name("FPositionSet2DCurvilinear");
    return name;
}

//--------------------------------------------------------------------------- 

positive FPositionSet2DCurvilinear::getNbXCoords() const
{
    return xCoords;
}

//--------------------------------------------------------------------------- 

positive FPositionSet2DCurvilinear::getNbYCoords() const
{
    return yCoords;
}

