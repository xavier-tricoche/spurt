//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet3DArbitrary.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:09 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 


#include "FPositionSet3DArbitrary.hh"
#include "FanyMMappedArray.hh"

#include "FIndex.hh"

#include "FException.hh"

#include "eassert.hh"



FPositionSet3DArbitrary::FPositionSet3DArbitrary( shared_ptr< FanyArray< double > >  coords ,const FBoundingBox*bb )
  : FPositionSet( 3 ),sharedPositionPtr(coords),positions(*coords),
    bufvector(3)
{
  
    //set variables from superclass
    dimension   = 3;
  
    // some sanity checks on the submitted coordinates
    eassert( positions.size() );
    eassert( positions.size() % 3 == 0 );

    //if bounding box not given
    if(!bb){
      // determine bounding box
      double minX, minY, maxX, maxY, maxZ, minZ;
      minX = maxX = positions[0];
      minY = maxY = positions[1];
      minZ = maxZ = positions[2];

      for( positive i=1; i<positions.size()/3; i++ ) 
	{
	if( positions[3*i] < minX )
	  minX = positions[3*i];
	
	if( positions[3*i] > maxX )
	  maxX = positions[3*i];
	
	if( positions[3*i+1] < minY )
	  minY = positions[3*i+1];
	
	if( positions[3*i+1] > maxY )
	  maxY = positions[3*i+1];
	
	if( positions[3*i+2] < minZ )
	  minZ = positions[3*i+2];
	
	if( positions[3*i+2] > maxZ )
	  maxZ = positions[3*i+2];
	}
      
      bBox = FBoundingBox( minX, minY, minZ, maxX, maxY, maxZ );
    }
    else
      bBox=*bb;
}

void
FPositionSet3DArbitrary::setTree(shared_ptr<FkdTree> atree)
{
  tree = atree;
}


FPositionSet3DArbitrary::FPositionSet3DArbitrary( vector<double> & coords )
  : FPositionSet( 3 ),
    sharedPositionPtr(new FanyVector<double> (coords) ),
    positions(*sharedPositionPtr),
    bufvector(3)
{

    //set variables from superclass
    dimension   = 3;
  
    // some sanity checks on the submitted coordinates
    eassert( positions.size() );
    eassert( positions.size() % 3 == 0 );


    // determine bounding box
    double minX, minY, maxX, maxY, maxZ, minZ;
    minX = maxX = positions[0];
    minY = maxY = positions[1];
    minZ = maxZ = positions[2];

    for( positive i=1; i<positions.size()/3; i++ ) 
    {
	if( positions[3*i] < minX )
	    minX = positions[3*i];

	if( positions[3*i] > maxX )
	    maxX = positions[3*i];

	if( positions[3*i+1] < minY )
	    minY = positions[3*i+1];

	if( positions[3*i+1] > maxY )
	    maxY = positions[3*i+1];

	if( positions[3*i+2] < minZ )
	    minZ = positions[3*i+2];

	if( positions[3*i+2] > maxZ )
	    maxZ = positions[3*i+2];
    }

    bBox = FBoundingBox( minX, minY, minZ, maxX, maxY, maxZ );
}

//--------------------------------------------------------------------------- 

const FString& FPositionSet3DArbitrary::getClassName() const
{
    static FString name("FPositionSet3DArbitrary");
    return name;
}

//--------------------------------------------------------------------------- 

void FPositionSet3DArbitrary::getPosition( FPosition& resultPos, 
					   const FIndex& pIndex ) const
{
    eassert( (positive)pIndex < positions.size()/3 );

    positions.get_range(3*pIndex.getIndex(),3,&bufvector[0]);
    resultPos=bufvector;
}

//--------------------------------------------------------------------------- 

void FPositionSet3DArbitrary::getPosition( vector<double>& resultPos, 
					   const FIndex& pIndex ) const
{
#ifndef NODEBUG
  eassert( (positive)pIndex < positions.size()/3 );
#endif

  resultPos.resize(3);
  positions.get_range(3*pIndex.getIndex(),3,&resultPos[0]);
}

//--------------------------------------------------------------------------- 

positive FPositionSet3DArbitrary::memSize() const
{
    return positions.size() * sizeof( positions.front() ) ;
}

//--------------------------------------------------------------------------- 

positive FPositionSet3DArbitrary::getNbPositions() const
{
    return positions.size()/3;
}

//--------------------------------------------------------------------------- 

void FPositionSet3DArbitrary::
getDistribution(vector<positive>& sizes,
		vector<string> & names) const
{
  FanyMMappedArray<double>* x 
    =  dynamic_cast< FanyMMappedArray<double>*  >(sharedPositionPtr.get());

  if(!x)
    THROW_EXCEPTION(FNotImplementedException," this class is not distributed! ");

  x->arr->getFileNames(names);
  x->arr->getBlocksPerFile(sizes);

  for(unsigned i=0;i<sizes.size();i++)
    sizes[i]/=3;
}
