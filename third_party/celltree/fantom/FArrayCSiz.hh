#ifndef __FArrayCSiz_hh
#define __FArrayCSiz_hh

#include "FArray.hh"
#include"FInvalidDimensionException.hh"

template<int dim>
class FArrayCSiz:public FArray{

public:

  FArrayCSiz()
  {
    comp=dat;
    sizeOfArray=dim;
  }

  FArrayCSiz(const FArray&x)
  {
    #ifndef NODEBUG
    if(x.size()<dim)
      throw FInvalidDimensionException();
    #endif
    memcpy(dat,&x[0],sizeof(dat));
  }

  void FArrayCSiz::operator=(const FArray&x)
  {
    #ifndef NODEBUG
    if(x.size()<dim)
      throw FInvalidDimensionException();
    #endif
    memcpy(dat,&x[0],sizeof(dat));
  }


  
  ~FArrayCSiz()
  {
    //to prevent deletion of this array
    comp=0;
  }

private:
  void resize(unsigned int,bool){}

  double dat[dim];
};


#endif
