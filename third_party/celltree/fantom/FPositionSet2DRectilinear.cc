//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet2DRectilinear.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:09 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#include "FPositionSet2DRectilinear.hh"
#include "FException.hh"
#include "FIndex.hh"

#include "binio.hh"

#include <cassert>

using namespace std;

//--------------------------------------------------------------------------- 

FPositionSet2DRectilinear::FPositionSet2DRectilinear( vector<double>& xCoords,
						      vector<double>& yCoords )
    : FPositionSet( 2 )
{
    assert( xCoords.size() && yCoords.size() );

    pos[0].swap( xCoords );
    pos[1].swap( yCoords );

    //set variables from superclass
    bBox.setBoundingBox( pos[0].front(), pos[1].front(), 
			 pos[0].back(), pos[1].back() );

    nbPositions = pos[0].size() * pos[1].size();  
}

//---------------------------------------------------------------------------

namespace{
  template<typename T>
  struct increment_generator
  {
    increment_generator(const T& val) : val_ (val), i(0){}
    T operator()(){
      return val_*(i++);
    }
    const T val_;
    int i;
  };
}
/*
FPositionSet3DRectilinear::FPositionSet3DRectilinear( unsigned int nx, unsigned int ny, unsigned int nz, double dx, double dy, double dz )
  : FPositionSet( 3 )
{
  assert( dx > 0 && dy > 0 && dz > 0 );
  assert( nx > 0 && ny > 0 && nz > 0 );
  
    pos[0].resize( nx );
    std::generate( pos[0].begin(), pos[0].end(), increment_generator<double>( dx ));
    pos[1].resize( ny );
    std::generate( pos[1].begin(), pos[1].end(), increment_generator<double>( dy ));
    pos[2].resize( nz );
    std::generate( pos[2].begin(), pos[2].end(), increment_generator<double>( dz ));
  
    //set variables from superclass
    bBox.setBoundingBox( pos[0].front(), pos[1].front(), pos[2].front(),
			 pos[0].back(),  pos[1].back(),  pos[2].back() );

    nbPositions = pos[0].size() * pos[1].size() * pos[2].size();  
}
*/

FPositionSet2DRectilinear::FPositionSet2DRectilinear( unsigned int nx, unsigned int ny,
    double dx, double dy)
    : FPositionSet( 2 )
{
  assert( dx > 0 && dy > 0 );
  assert( nx > 0 && ny > 0 );
  
    pos[0].resize( nx );
    std::generate( pos[0].begin(), pos[0].end(), increment_generator<double>( dx ));
    pos[1].resize( ny );
    std::generate( pos[1].begin(), pos[1].end(), increment_generator<double>( dy ));

    //set variables from superclass
    bBox.setBoundingBox( pos[0].front(), pos[1].front(), 
			 pos[0].back(), pos[1].back() );

    nbPositions = pos[0].size() * pos[1].size();  
}

//--------------------------------------------------------------------------- 

const FString& FPositionSet2DRectilinear::getClassName() const
{
    static FString ugly("FPositionSet2DRectilinear");
    return ugly;
}

//--------------------------------------------------------------------------- 

void FPositionSet2DRectilinear::getPosition( FPosition& resultPos, 
					     const FIndex& pIndex ) const
{
    positive ind = pIndex.getIndex();
    assert( ind < nbPositions );

    resultPos.resize( 2 );

    resultPos[0] = pos[0][ ind%pos[0].size() ];
    resultPos[1] = pos[1][ ind/pos[0].size() ];
}

//--------------------------------------------------------------------------- 

void FPositionSet2DRectilinear::getPosition( std::vector<double>& resultPos, 
					     const FIndex& pIndex ) const
{
    positive ind = pIndex.getIndex();
    assert( ind < nbPositions );
    
    resultPos.resize(2);
    
    resultPos[0] = pos[0][ ind%pos[0].size() ];
    resultPos[1] = pos[1][ ind/pos[0].size() ];
}

//--------------------------------------------------------------------------- 

positive FPositionSet2DRectilinear::getNbPositions() const
{
    return nbPositions;
}

//--------------------------------------------------------------------------- 

positive FPositionSet2DRectilinear::memSize() const
{
    return (pos[0].size()+pos[1].size()) * sizeof( pos[0].front() );
}

//--------------------------------------------------------------------------- 

const vector<double>& FPositionSet2DRectilinear::getXCoords() const
{
    return pos[0];
}

//--------------------------------------------------------------------------- 

const vector<double>& FPositionSet2DRectilinear::getYCoords() const
{
    return pos[1];
}

positive FPositionSet2DRectilinear::getNbXCoords() const
{
  return pos[0].size();
}

positive FPositionSet2DRectilinear::getNbYCoords() const
{
  return pos[1].size();
}

// //--------------------------------------------------------------------------- 

// void FPositionSet2DRectilinear::serialize( ostream& out ) const
// {
//     binwrite( out, pos[0] );
//     binwrite( out, pos[1] );
// }

// //--------------------------------------------------------------------------- 

// FPositionSet2DRectilinear* FPositionSet2DRectilinear::rebuild( istream& in )
// {
//     vector<double> xc, yc;

//     binread( in, xc );
//     binread( in, yc );

//     return new FPositionSet2DRectilinear( xc, yc );
// }
// //--------------------------------------------------------------------------- 

