//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FRefTensor.cc,v $
// Language:  C++
// Date:      $Date: 2003/03/27 08:29:52 $
// Author:    $Author: garth $
// Version:   $Revision: 1.7 $
//
//--------------------------------------------------------------------------- 

#include "FRefTensor.hh"
#include "FException.hh"

FRefTensor::FRefTensor(FTensor&x)
{
  memcpy(this,&x,sizeof(x));
}

FRefTensor::FRefTensor()
{
}


FTensor& FRefTensor::operator = (const FTensor&x)
{
#ifndef NODEBUG
  if(sizeOfArray!=x.sizeOfArray) {
    FInvalidDimensionException e;
    e.addTraceMessage("FTensor& FRefTensor::operator = (const FTensor&x)");
    throw e;    
  }  
#endif
  memcpy(comp,x.comp,sizeOfArray*sizeof(*comp));
  return *this;
}


FRefTensor::FRefTensor(int lDimension,int lOrder,double*lComp=0)
{
  dimension=lDimension;
  order=lOrder;
  sizeOfArray=pow(dimension,order);
  comp=lComp;
}

FRefTensor::~FRefTensor()
{
  comp=0;//to force destructor of FTensor not to delete anyything
  dimension=0;
  sizeOfArray=0;order=0;
}

void FRefTensor::setPointerOnComponents(double *Comp)
{
  comp=Comp;
}

void FRefTensor::setDimensionOrderPointerAndSize(int d,int o,double *c,int s)
{
  dimension=d;
  order=o;
  if(!s)
    sizeOfArray=pow(d,o);
  else 
    sizeOfArray=s;
  comp=c;
}


#ifdef OUTLINE
#include "FRefTensor.icc"
#endif

//---------------------------------------------------------------------------
