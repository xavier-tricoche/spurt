//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionSet3DRectilinear.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:10 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 


#include "FPositionSet3DRectilinear.hh"
#include "FException.hh"
#include "FIndex.hh"

using namespace std;

#include <cassert>

//--------------------------------------------------------------------------- 

FPositionSet3DRectilinear::FPositionSet3DRectilinear( std::vector<double>& XCoords,
						      std::vector<double>& YCoords,
						      std::vector<double>& ZCoords )
    : FPositionSet( 3 )
{
    assert( XCoords.size() && YCoords.size() && ZCoords.size() );

    pos[0].swap( XCoords );
    pos[1].swap( YCoords );
    pos[2].swap( ZCoords );
  
    //set variables from superclass
    bBox.setBoundingBox( pos[0].front(), pos[1].front(), pos[2].front(),
			 pos[0].back(),  pos[1].back(),  pos[2].back() );

    nbPositions = pos[0].size() * pos[1].size() * pos[2].size();  
}

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

//--------------------------------------------------------------------------- 

const FString& FPositionSet3DRectilinear::getClassName() const
{
  static FString name("FPositionSet3DRectilinear");
  return name;
}

//--------------------------------------------------------------------------- 

void FPositionSet3DRectilinear::getPosition( FPosition& resultPos, 
					     const FIndex& pIndex ) const
{
    positive ind = pIndex.getIndex();
    
    assert( ind < nbPositions );

    resultPos.resize(3);

    resultPos[0] = pos[0][ ind%pos[0].size() ];

    ind /= pos[0].size();

    resultPos[1] = pos[1][ ind%pos[1].size() ];
    resultPos[2] = pos[2][ ind/pos[1].size() ];
}

//--------------------------------------------------------------------------- 

void FPositionSet3DRectilinear::getPosition( std::vector<double>& resultPos, 
					     const FIndex& pIndex ) const
{
    positive ind = pIndex.getIndex();
    
    assert( ind < nbPositions );

    resultPos.resize(3);

    resultPos[0] = pos[0][ ind%pos[0].size() ];

    ind /= pos[0].size();

    resultPos[1] = pos[1][ ind%pos[1].size() ];
    resultPos[2] = pos[2][ ind/pos[1].size() ];
}

//--------------------------------------------------------------------------- 

positive FPositionSet3DRectilinear::getNbPositions() const
{
    return nbPositions;
}

//--------------------------------------------------------------------------- 

positive FPositionSet3DRectilinear::memSize() const
{
    return (pos[0].size()+pos[1].size()+pos[2].size()) * sizeof( pos[0].front() );
}

//--------------------------------------------------------------------------- 

const vector<double>& FPositionSet3DRectilinear::getXCoords() const
{
    return pos[0];
}

//--------------------------------------------------------------------------- 

const vector<double>& FPositionSet3DRectilinear::getYCoords() const
{
    return pos[1];
}

//--------------------------------------------------------------------------- 

const vector<double>& FPositionSet3DRectilinear::getZCoords() const
{
    return pos[2];
}

positive FPositionSet3DRectilinear::getNbXCoords() const
{
  return pos[0].size();
}

positive FPositionSet3DRectilinear::getNbYCoords() const
{
  return pos[1].size();
}

positive FPositionSet3DRectilinear::getNbZCoords() const
{
  return pos[2].size();
}

positive FPositionSet3DRectilinear::getDimensionX() const
{
  return pos[0].size();
}

positive FPositionSet3DRectilinear::getDimensionY() const
{
  return pos[1].size();
}

positive FPositionSet3DRectilinear::getDimensionZ() const
{
  return pos[2].size();
}
