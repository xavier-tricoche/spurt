//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTrilinearInterpolation.cc,v $
// Language:  C++
// Date:      $Date: 2003/03/27 08:29:53 $
// Author:    $Author: garth $
// Version:   $Revision: 1.6 $
//
//--------------------------------------------------------------------------- 

#include "FTrilinearInterpolation.hh"
#include "FException.hh"


void FTrilinearInterpolation::
buildParameters(double * p, const double * end)
{ 
  for(;p!=end;p+=8)
    buildParameters(p);
}



void FTrilinearInterpolation::
interpolateT( const double* p, const double * end,
	      const double*xx, double * ret)
{
     
  double x=xx[0],y=xx[1],z=xx[2];
  double xy=x*y,xz=x*z,yz=y*z,xyz=xy*z;
  for(;p!=end;p+=8,ret++)
    *ret =
      p[0] 
      + x * p[1] + y * p[3] + z * p[4]
      + xy * p[2]  +  xz * p[5]  +  yz * p[7]
      + xyz * p[6];
}

//-------------------------------------------------------------------

void FTrilinearInterpolation::
interpolatedTdx(const double * p, const double * end,
		const double*x, double * ret)
{
  double yz = x[1]*x[2];
  for(;p!=end;p+=8,ret++)
    *ret = p[1] + x[1] * p[2] + x[2] * p[5] + yz * p[6];
}   

//-------------------------------------------------------------------

void FTrilinearInterpolation::
interpolatedTdy(const double * p, const double * end,
		const double*x, double * ret)
{
  double xz = x[0]*x[2];
  for(;p!=end;p+=8,ret++)
    *ret=  p[3] + x[0] * p[2] + x[2] * p[7] + xz * p[6];      
}

void FTrilinearInterpolation::
interpolatedTdz(const double * p, const double * end,
		const double*x, double * ret)
{
  double xy = x[0]*x[1];
  for(;p!=end;p+=8,ret++)
    *ret= p[4] + x[0] * p[5] + x[1] * p[7] + xy * p[6];
}


//-------------------------------------------------------------------





void FTrilinearInterpolation::
buildParameters(vector<double>&p)
{ 
#ifndef NODEBUG
  if(p.size()%8!=0){
    FIndexOutOfBoundsException e("vector is not in the right size");
    e.addTraceMessage("void FTrilinearInterpolation::"
		      "buildParameters(vector<double>&p)");
    throw e;
  }
#endif

  params.swap(p);
  
  vector<double>::iterator pit;
  for(pit=params.begin();pit!=params.end();pit+=8)
    buildParameters(&(*pit));
}

//-------------------------------------------------------------------


void FTrilinearInterpolation::
interpolateT(const FVector&xx, FArray&retval) const
{

#ifndef NODEBUG
  if(params.size()/8!=retval.size()){
    FIndexOutOfBoundsException e("vector is not in the right size");
    e.addTraceMessage("void FTrilinearInterpolation::"
		      "interpolateT(const FVector&xx, FArray&ret) const");
    throw e;
  }
#endif
     
  double x=xx[0],y=xx[1],z=xx[2];
  double xy=x*y,xz=x*z,yz=y*z,xyz=xy*z;
  
  vector<double>::const_iterator p;
  double* ret = &retval[0];

  for(p=params.begin();p!=params.end();p+=8,ret++)
    *ret =
      p[0] 
      + x * p[1] + y * p[3] + z * p[4]
      + xy * p[2]  +  xz * p[5]  +  yz * p[7]
      + xyz * p[6];
}

//-------------------------------------------------------------------

void FTrilinearInterpolation::
interpolatedTdx(const FVector&xx, FArray&retval) const
{

#ifndef NODEBUG
  if(params.size()/8!=retval.size()){
    FIndexOutOfBoundsException e("vector is not in the right size");
    e.addTraceMessage("void FTrilinearInterpolation::"
		      "interpolateTdx(const FVector&xx, FArray&ret) const");
    throw e;
  }
#endif

  double y=xx[1],z=xx[2],yz = y*z;

  vector<double>::const_iterator p;
  double* ret = &retval[0];

  for(p=params.begin();p!=params.end();p+=8,ret++)
    *ret = p[1] + y * p[2] + z * p[5] + yz * p[6];
}   

//-------------------------------------------------------------------

void FTrilinearInterpolation::
interpolatedTdy(const FVector&xx, FArray&retval) const
{

#ifndef NODEBUG
  if(params.size()/8!=retval.size()){
    FIndexOutOfBoundsException e("vector is not in the right size");
    e.addTraceMessage("void FTrilinearInterpolation::"
		      "interpolateTdy(const FVector&xx, FArray&ret) const");
    throw e;
  }
#endif

  double x=xx[0],z=xx[2],xz = x*z;

  vector<double>::const_iterator p=params.begin();
  double* ret = &retval[0];

  for(p=params.begin();p!=params.end();p+=8,ret++)
    *ret=  p[3] + x * p[2] + z * p[7] + xz * p[6];      
}
//-------------------------------------------------------------------

void FTrilinearInterpolation::
interpolatedTdz(const FVector&xx, FArray&retval) const
{
#ifndef NODEBUG
  if(params.size()/8!=retval.size()){
    FIndexOutOfBoundsException e("vector is not in the right size");
    e.addTraceMessage("void FTrilinearInterpolation::"
		      "interpolateTdz(const FVector&xx, FArray&ret) const");
    throw e;
  }
#endif

  double x=xx[0],y=xx[1],xy = x*y;

  vector<double>::const_iterator p=params.begin();
  double* ret = &retval[0];

  for(p=params.begin();p!=params.end();p+=8,ret++)
    *ret= p[4] + x * p[5] + y * p[7] + xy * p[6];
}
