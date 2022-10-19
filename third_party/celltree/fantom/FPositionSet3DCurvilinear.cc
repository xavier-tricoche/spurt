//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet3DCurvilinear.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:10 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 


#include "FPositionSet3DCurvilinear.hh"
#include "FIndex.hh"

#include "FException.hh"

//--------------------------------------------------------------------------- 

FPositionSet3DCurvilinear::
FPositionSet3DCurvilinear( vector<double>& coords,
			   positive newxCoords,
			   positive newyCoords,
			   positive newzCoords )
    : FPositionSet3DArbitrary( coords )
{
    xCoords = newxCoords;
    yCoords = newyCoords;
    zCoords = newzCoords;
}

//--------------------------------------------------------------------------- 

const FString& FPositionSet3DCurvilinear::getClassName() const
{
    static FString name("FPositionSet3DCurvilinear");
    return name;
}

//--------------------------------------------------------------------------- 

positive FPositionSet3DCurvilinear::getNbXCoords() const
{
    return xCoords;
}

//--------------------------------------------------------------------------- 

positive FPositionSet3DCurvilinear::getNbYCoords() const
{
    return yCoords;
}

//--------------------------------------------------------------------------- 

positive FPositionSet3DCurvilinear::getNbZCoords() const
{
    return zCoords;
}
