//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile$
// Language:  C++
// Date:      $Date$
// Author:    $Author$
// Version:   $Revision$
//
//--------------------------------------------------------------------------- 

#include "FPositionSet2Din3DRectilinear.hh"
#include "FException.hh"
#include "FIndex.hh"

#include "binio.hh"

#include <cassert>

using namespace std;

//--------------------------------------------------------------------------- 

FPositionSet2Din3DRectilinear::
FPositionSet2Din3DRectilinear( vector<double>& xCoords,
			       vector<double>& yCoords,
			       const FMatrix& _trans )
    : FPositionSet2DRectilinear( xCoords, yCoords ), trans( _trans )
{
}

//--------------------------------------------------------------------------- 

const FString& FPositionSet2Din3DRectilinear::getClassName() const
{
    static FString ugly("FPositionSet2Din3DRectilinear");
    return ugly;
}

//--------------------------------------------------------------------------- 

void FPositionSet2Din3DRectilinear::get3DPosition( FArray& pos3D, 
						   const FIndex& id ) const
{
    positive ind = id.getIndex();
    assert( ind < nbPositions );

    pos3D.resize( 3 );

    for ( positive i=0 ; i<3 ; ++i )
    {
	pos3D[i] = 
	    trans(i,0)*pos[0][ ind%pos[0].size() ] +
	    trans(i,1)*pos[1][ ind/pos[0].size() ] +
	    trans(i,2);
    }
}

//--------------------------------------------------------------------------- 

void FPositionSet2Din3DRectilinear::get2DPosition( FArray& pos2D, 
						   const FArray& pos3D ) const
{
    pos2D.resize( 2 );
    
    for ( positive i=0 ; i<2 ; ++i )
	pos2D(i) = 
	    ( pos3D(0)-trans(0,2) )*trans(0,i) + 
	    ( pos3D(1)-trans(1,2) )*trans(1,i) +
	    ( pos3D(2)-trans(2,2) )*trans(2,i);
    
}

//--------------------------------------------------------------------------- 
