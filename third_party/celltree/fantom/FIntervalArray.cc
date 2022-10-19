
#include "FIntervalArray.hh"
#ifdef OUTLINE
#define inline
#include "FIntervalArray.icc"
#endif

#include <numeric>
#include <iostream>
#include <cmath>

using namespace std;

/** 
 * constructor.
 * contructs a number of intervals
 * which are separated by the the entries of borders
 *
 *\param borders 
 *start indices of the arrays
 */

FIntervalArray::FIntervalArray(const std::vector<positive>& sizes)
{

  positive numIntervals= sizes.size();  

  if(!numIntervals)
    {sumSize=0;return;}

  //exclude terminating null entries
  while( ! sizes[numIntervals-1]  && numIntervals > 1 ) --numIntervals;
  
  if( numIntervals==1 ){
    table.resize(1);
    sumSize=sizes[0];

    //sumSize >> pMinSize is now 0:
    frexp( double(sumSize), &pMinSize);
    
    table[0].interval=0;
    table[0].border=sumSize;
    return;
  }
    
  //minsize is the smallest entry of the sizes vector (without the last entry)
  positive minSize = * std::min_element(sizes.begin(),sizes.begin() + numIntervals - 1 );

  frexp(double(minSize),&pMinSize);
  pMinSize--;

  minSize=1LU<<pMinSize;
  
  sumSize = std::accumulate(sizes.begin(),sizes.end(),positive(0));

  //arraygroesse ist aufgerundet sumSize/minSize
  positive tblSiz;

  if(pMinSize>=0)
    tblSiz = (sumSize+minSize-1)>>pMinSize;
  
  //if minSize was to small, and table size too big, use
  //other method (binary search instead of direct addressing)
  //will be indicated by pMinSize<0
  if( pMinSize<0 || tblSiz > sizes.size()*10000 ) {
    pMinSize = -1;
    cerr<<"FIntervalArray: distribution too bad, using binary search"<<endl;

    table.resize(numIntervals);

    table[0].interval = 0;
    table[0].border = sizes[0];

    for(unsigned i=1;i<numIntervals;i++)
      {
	table[i].interval = i;
	table[i].border   =  table[i-1].border + sizes[i];	
      }
      
    return;
  }



  table.resize(tblSiz);
  
  //  cerr<<"FintervalArray:size "<<table.size()*sizeof(table[0])<<" b"<<flush;
  
  positive border=sizes[0];
  positive interval=0;
  
  for(positive i=0 ; i<table.size() ; i++){
    
    if((i<<pMinSize)>=border){
      interval++;
      border+=sizes[interval];
    }
    
    table[i].interval = interval;
    table[i].border = border;
    
  }

}

