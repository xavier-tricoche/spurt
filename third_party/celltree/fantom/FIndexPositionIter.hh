
#ifndef __FIndexPositionIter_hh
#define __FIndexPositionIter_hh

#include <iterator>
#include "FanyArray.hh"

using namespace std;

struct FIndexPositionIter;
struct FIndexPositionRef{
  

  FAssignableAnyArray<double>  * positions;
  FAssignableAnyArray<positive> * indices;
  positive dimension,coordId;
  //coordinate in positions used for comparisons

  positive index;



  typedef struct value_type
  {
    double pos[3];
    positive ind;
  };

  typedef FIndexPositionIter pointer;
  typedef long difference_type;
  typedef FIndexPositionRef reference;


  
  
  reference & operator = (const reference & b)
  {
    positive indexA = dimension*index, indexB = dimension*b.index;
    
    for(positive i=0;i<dimension;i++){
      (*positions) [ indexA  + i ] = (*b.positions) [ indexB  + i ];
    }
    
    (*indices)[ index ] = (*b.indices) [ b.index ];

    return *this;
  }

  reference & operator = (const value_type & b)
  {
    positive indexA = dimension*index;
    
    for(positive i=0;i<dimension;i++){
      (*positions) [ indexA  + i ] = b.pos [ i ];
    }
    
    (*indices)[ index ] = b.ind;

    return *this;
  }

  operator value_type()
  {
    value_type x;
    positive indexB = dimension*index;
    
    for(positive i=0;i<dimension;i++){
      x.pos [ i ] = (*positions) [ indexB  + i ];
    }
    x.ind = (*indices)[index];
    
    return x;
  }

  operator double()
  {
    return (*positions)[ index*dimension + coordId ];
  }

};


struct FIndexPositionIter
{

  typedef random_access_iterator_tag iterator_category;

  typedef FIndexPositionRef::value_type value_type;
  typedef FIndexPositionIter pointer;
  typedef long difference_type;
  typedef FIndexPositionRef reference;

  reference ref;

  inline FIndexPositionIter(const FIndexPositionRef&x)
    :ref(x)
  {}

  inline FIndexPositionIter(){}


  inline reference & operator*(){return ref;}

  inline reference operator[](positive i)
  { 
    reference a=ref;
    a.index+=i;
    return a;
  }


  inline FIndexPositionIter operator + (positive i)
  {
    FIndexPositionIter a=*this;
    a.ref.index+=i;
    return a;
  }

  inline FIndexPositionIter operator - (positive i)
  {
    FIndexPositionIter a=*this;
    a.ref.index-=i;
    return a;
  }


  inline FIndexPositionIter operator ++ (int)
  {
    FIndexPositionIter a=*this;
    ++ref.index;
    return a;
  }

  inline FIndexPositionIter operator -- (int)
  {
    FIndexPositionIter a=*this;
    --ref.index;
    return a;
  }

  inline FIndexPositionIter& operator ++ ()
  {++ref.index;return *this;}

  inline FIndexPositionIter& operator -- ()
  {--ref.index;return *this;}

  inline FIndexPositionIter& operator += (positive i)
  {ref.index+=i;return *this;}

  inline FIndexPositionIter& operator -= (positive i)
  {ref.index-=i;return *this;}

  bool operator == (const FIndexPositionIter&b)
  {return ref.index == b.ref.index;}
  
  difference_type
  operator - (const FIndexPositionIter&b)
  {return ref.index-b.ref.index;}

};

#endif
