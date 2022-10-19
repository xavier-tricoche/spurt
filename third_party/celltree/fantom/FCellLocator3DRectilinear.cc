#include "FCellLocator3DRectilinear.hh"
#include "FPositionSet3DRectilinear.hh"
#include "FException.hh"
#include "FIndex.hh"

#include <iostream>

FCellLocator3DRectilinear::~FCellLocator3DRectilinear()
{
}


FCellLocator3DRectilinear::
FCellLocator3DRectilinear( const FPositionSet* positionSet,
			   bool triangulated )
  : FCellLocator(positionSet)
{
   isTriangulated = triangulated;

   const FPositionSet3DRectilinear *pSet 
       = dynamic_cast<const FPositionSet3DRectilinear*>( FCellLocator::pSet );

   nbCells = (pSet->pos[0].size()-1) * (pSet->pos[1].size()-1)  * (pSet->pos[2].size()-1);

   if (triangulated)
     nbCells*=6;
}


bool FCellLocator3DRectilinear::
searchCell( FIndex& aIndex, const FPosition& aPosition ) const
{
  try{

    const FPositionSet3DRectilinear *pSet 
      = dynamic_cast<const FPositionSet3DRectilinear*>(FCellLocator::pSet);

    if(!pSet)
      throw FException("wrong positionset type for this celllocator");

    positive indC[3]; //indexes of cell in the 3 coordinates

    for(int i=0;i<3;i++){

      const vector<double>&p=pSet->pos[i];
      double ap=aPosition[i];

      
      if(ap>p[p.size()-1] != ap<p[0]){

	//1st exit point
	aIndex.setToInvalid();
	return false;
      }


      //binary search in direction of the actual coordinate axle

      positive h=p.size()-2,l=0,m;

      while(h!=l){

	m=(h+l+1)/2;
	if(ap<p[m])
	  h=m-1;
	else 
	  l=m;      
      }

      indC[i]=l;

    }


    positive ind = ( indC[2] *(pSet->pos[1].size()-1) + indC[1] ) * (pSet->pos[0].size()-1) + indC[0];

    //this case is modelled as in the obsolete class FCellSubLocatorStructured3D 
    if ( isTriangulated )
      {

	//koordinaten im jeweiligen voxel (im bereich 0..1)
	double voko[3];

	for(positive dir=0;dir<3;dir++)
	  {
	    const double* d = & (pSet->pos[dir][0]) + indC[dir];	 
	    voko[dir]= (aPosition[dir]-d[0]) / (d[1]-d[0]);
	  }

	//cout<<"voxelkoords:"<<voko[0]<<' '<<voko[1]<< ' '<<voko[2]<<endl;
	
	// "9999" steht fuer faelle, die nicht vorkommen
	static positive drk[8]={5, 3, 9999, 2, 4, 9999, 1, 0};

	//saving the results of the 3 following comparisons 
	//in three bits of idrk, the result of the 1st comparison as MSB
	positive idrk =    (voko[0]>voko[1]);
	idrk = (idrk<<1) | (voko[0]>voko[2]);
	idrk = (idrk<<1) | (voko[1]>voko[2]);

	//cout<<"idrk:"<<idrk<<endl;
	if(drk[idrk]==9999){cout<<"wrong idrk"<<endl;throw FException("error in implementation");}
      
	ind = 6*ind + drk[idrk];
      }


    //2nd exit point
    aIndex = ind;
    return true;

  }
  
  catch(FException e){
    e.addTraceMessage
      ("bool FCellLocator3DRectilinear::searchCell( FIndex& aIndex, const FPosition& aPosition ) const");
    throw e;
  }

}


positive FCellLocator3DRectilinear::memSize() const
{
  return sizeof(*this);
}
