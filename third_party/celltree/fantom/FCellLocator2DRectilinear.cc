//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellLocator2DRectilinear.cc,v $
// Language:  C++
// Date:      $Date: 2004/02/14 08:00:08 $
// Author:    $Author: hlawit $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FCellLocator2DRectilinear.hh"
#include "FPositionSet2DRectilinear.hh"
#include "FException.hh"
#include "FIndex.hh"

// --------------------------------------------------------------------------

FCellLocator2DRectilinear::
FCellLocator2DRectilinear( const FPositionSet* positionSet,
			   bool triangulated )
  : FCellLocator(positionSet)
{
    isTriangulated = triangulated;

    const FPositionSet2DRectilinear *pSet 
	= dynamic_cast<const FPositionSet2DRectilinear*>( FCellLocator::pSet );
   
    nbCells = (pSet->pos[0].size()-1) * (pSet->pos[1].size()-1);
    
    if (triangulated)
	nbCells*=2;
}

// --------------------------------------------------------------------------

FCellLocator2DRectilinear::~FCellLocator2DRectilinear()
{
}

// --------------------------------------------------------------------------


bool FCellLocator2DRectilinear::
searchCell( FIndex& aIndex, const FPosition& aPosition ) const
{
  try{

    const FPositionSet2DRectilinear *pSet 
      = dynamic_cast<const FPositionSet2DRectilinear*>( FCellLocator::pSet );

    if(!pSet)
      throw FException("wrong positionset type for this celllocator");

    positive indC[2]; //indexes in the two coordinates

    for(int i=0;i<2;i++){

      const vector<double> &p = pSet->pos[i];
      double ap = aPosition[i];

      if(ap > p[p.size()-1] != ap < p[0]) {
	aIndex.setToInvalid();
	return false;
      }

      //binary search direction of the actual coordinate axle
      //(search biggest p smaller than ap)

      positive h = p.size()-2, l = 0, m;

      while(h != l){

	m = (h+l+1) / 2; //aufgerundet
	if(ap < p[m])
	  h = m-1;
	else 
	  l = m;      
      }

      indC[i] = l;

    }

    positive ind = indC[1] * (pSet->pos[0].size()-1) + indC[0];

    if (isTriangulated) {

      double pmin0 = pSet->pos[0][indC[0]],   pmin1 = pSet->pos[1][indC[1]];
      double pmax0 = pSet->pos[0][indC[0]+1], pmax1 = pSet->pos[1][indC[1]+1];

      ind *= 2; 

      //if position is in upper left triangle, increase index by one
      ind += ( (aPosition[0]-pmin0) / (pmax0-pmin0) + 
	       (aPosition[1]-pmin1) / (pmax1-pmin1) ) > 1; // u+v > 1 ?

    }

    aIndex = ind;
    return true;

  }
  
  catch(FException e){
    e.addTraceMessage
      ("bool FCellLocator2DRectilinear::searchCell( FIndex& aIndex, const FPosition& aPosition ) const");
    throw e;
  }

}


positive FCellLocator2DRectilinear::memSize() const
{
  return sizeof(*this);
}
