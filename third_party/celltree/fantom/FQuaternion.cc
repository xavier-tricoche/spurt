//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FQuaternion.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:54:48 $
// Author:    $Author: garth $
// Version:   $Revision: 1.10 $
//
//--------------------------------------------------------------------------- 

#include <fstream>
#include <iomanip>

#include "FQuaternion.hh"
#include "FException.hh"

using namespace std;

#ifdef OUTLINE
#include "FQuaternion.icc"
#endif

//--------------------------------------------------------------------------- 

FQuaternion::FQuaternion()
{
}

//---------------------------------------------------------------------------

FQuaternion::~FQuaternion()
{
}

//---------------------------------------------------------------------------

FQuaternion::FQuaternion(const double &value) : FObject()
{
  for(unsigned char i=0; i<4; i++)
    comp[i] = value;
}

//---------------------------------------------------------------------------

FQuaternion::FQuaternion(const FVector& v) : FObject()
{
#ifndef NODEBUG
  if (v.getDimension() != 4)
    throw FInvalidDimensionException("FQuaternion must use 4Dimensional vector for copyconstructor");
#endif
  for(unsigned int i=0; i<4; i++)
    comp[i] = v[i];
}

//---------------------------------------------------------------------------

FQuaternion::FQuaternion(const FQuaternion& v) : FObject()
{
  for(unsigned int i=0; i<4; i++)
    comp[i] = v[i];
}

//---------------------------------------------------------------------------

const FString& FQuaternion::getClassName() const
{
  static const FString classname("FQuaternion");

  return classname;
}

//---------------------------------------------------------------------------

void FQuaternion::toMatrix(FMatrix& m)
{ 
  double Nq = comp[0]*comp[0] + comp[1]*comp[1] + 
              comp[2]*comp[2] + comp[3]*comp[3];
  double s = (Nq > 0.0) ? (2.0 / Nq) : 0.0;
  double xs = comp[0]*s,      ys = comp[1]*s,	   zs = comp[2]*s;
  double wx = comp[3]*xs,     wy = comp[3]*ys,     wz = comp[3]*zs;
  double xx = comp[0]*xs,     xy = comp[0]*ys,     xz = comp[0]*zs;
  double yy = comp[1]*ys,     yz = comp[1]*zs,     zz = comp[2]*zs;

  m(0,0) = 1.0 - (yy + zz); 
  m(1,0) = xy + wz; 
  m(2,0) = xz - wy;
  m(0,1) = xy - wz; 
  m(1,1) = 1.0 - (xx + zz); 
  m(2,1) = yz + wx;
  m(0,2) = xz + wy; 
  m(1,2) = yz - wx; 
  m(2,2) = 1.0 - (xx + yy);
  m(0,3) = m(1,3) = m(2,3) = m(3,0) = m(3,1) = m(3,2) = 0.0;
  m(3,3) = 1.0;
}

//---------------------------------------------------------------------------

const FQuaternion FQuaternion::qOne()
{
  FQuaternion qone;
  qone[0] = qone[1] = qone[2] = 0.0 ;
  qone[3] = 1.0;
  return qone;
}

//---------------------------------------------------------------------------

void FQuaternion::mult(const FQuaternion &qR, FQuaternion &dest)
{
  dest[3] = comp[3]*qR[3] - comp[0]*qR[0] - comp[1]*qR[1] - comp[2]*qR[2];
  dest[0] = comp[3]*qR[0] + comp[0]*qR[3] + comp[1]*qR[2] - comp[2]*qR[1];
  dest[1] = comp[3]*qR[1] + comp[1]*qR[3] + comp[2]*qR[0] - comp[0]*qR[2];
  dest[2] = comp[3]*qR[2] + comp[2]*qR[3] + comp[0]*qR[1] - comp[1]*qR[0];
}

//---------------------------------------------------------------------------

ostream& operator<<(ostream& os, const FQuaternion& a)
{
  os << "[ ";
  for(unsigned int i=0; i<4; i++)
    {
        os << setw(8) << setprecision(4) << a.comp[i];
        if (i<3)
        os << ",";
    }
  os << " ]";
  return os;
}     

//===========================================================================
