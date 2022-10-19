//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet2DArbitrary.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:09 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 


#include "FPositionSet2DArbitrary.hh"
#include "FIndex.hh"

#include "FException.hh"

#include <cassert>

//--------------------------------------------------------------------------- 

FPositionSet2DArbitrary::FPositionSet2DArbitrary( shared_ptr<FanyArray<double> > coords )
  : FPositionSet( 2 ) ,positions (coords)
{
    // some sanity checks on the submitted coordinates
    assert( positions->size() );
    assert( positions->size() % 2 == 0 );


    // determine bounding box
    double minX, minY, maxX, maxY;

    minX = maxX = (*positions)[0];
    minY = maxY = (*positions)[1];

    for( unsigned int i=1; i<positions->size()/2; i++ )
    {
	if( (*positions)[2*i] < minX )
	  minX = (*positions)[2*i];
	else if( (*positions)[2*i] > maxX )
	  maxX = (*positions)[2*i];

	if( (*positions)[2*i+1] < minY )
	    minY = (*positions)[2*i+1];
	else if( (*positions)[2*i+1] > maxY )
	    maxY = (*positions)[2*i+1];
    }

    bBox.setBoundingBox( (double)minX, (double)minY,
			 (double)maxX, (double)maxY );
}


FPositionSet2DArbitrary::FPositionSet2DArbitrary( vector<double>& coords )
  : FPositionSet( 2 ) ,positions (new FanyVector<double>(coords))
{
    // some sanity checks on the submitted coordinates
    assert( positions->size() );
    assert( positions->size() % 2 == 0 );


    // determine bounding box
    double minX, minY, maxX, maxY;

    minX = maxX = (*positions)[0];
    minY = maxY = (*positions)[1];

    for( unsigned int i=1; i<positions->size()/2; i++ )
    {
	if( (*positions)[2*i] < minX )
	  minX = (*positions)[2*i];
	else if( (*positions)[2*i] > maxX )
	  maxX = (*positions)[2*i];

	if( (*positions)[2*i+1] < minY )
	    minY = (*positions)[2*i+1];
	else if( (*positions)[2*i+1] > maxY )
	    maxY = (*positions)[2*i+1];
    }

    bBox = FBoundingBox( (double)minX, (double)minY,
			 (double)maxX, (double)maxY );
}

//--------------------------------------------------------------------------- 

const FString& FPositionSet2DArbitrary::getClassName() const
{
    static FString name("FPositionSet2DArbitrary");
    return name;
}

//--------------------------------------------------------------------------- 

void FPositionSet2DArbitrary::getPosition( FPosition& resultPos, 
					   const FIndex& pIndex ) const
{
    assert( (positive)pIndex < positions->size()/2 );

    resultPos.resize(2);
    resultPos[0] = (double)((*positions)[ 2*pIndex.getIndex() ] );
    resultPos[1] = (double)((*positions)[ 2*pIndex.getIndex()+1 ]);
}

//--------------------------------------------------------------------------- 

void FPositionSet2DArbitrary::getPosition( vector<double>& resultPos, 
					   const FIndex& pIndex ) const
{
    assert( (positive)pIndex < positions->size()/2 );

    resultPos.resize( 2 );
    resultPos[0] = (double)((*positions)[ 2*pIndex.getIndex() ]);
    resultPos[1] = (double)((*positions)[ 2*pIndex.getIndex()+1 ]);
}

//--------------------------------------------------------------------------- 

positive FPositionSet2DArbitrary::getNbPositions() const
{
    return positions->size()/2;
}

//--------------------------------------------------------------------------- 

positive FPositionSet2DArbitrary::memSize() const
{
    return positions->size() * sizeof( positions->front() );
}
