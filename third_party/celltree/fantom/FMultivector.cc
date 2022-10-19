//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FMultivector.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:27 $
// Author:    $Author: garth $
// Version:   $Revision: 1.7 $
//
//--------------------------------------------------------------------------- 

#include "FMultivector.hh"
#include <cmath>
#include <iostream>

using namespace std;

//---------------------------------------------------------------------------

FMultivector::~FMultivector()
{
}

//---------------------------------------------------------------------------

ostream& operator<< (ostream& os, const FMultivector& t)
{
  os <<"[ ";
  for (unsigned char i=0; i<(unsigned char)FMultivector::pow (t.dimension,
							 t.order); i++)
    os << t.comp[i] << " ";
  os <<" ]";
  return os;
}

//---------------------------------------------------------------------------


//===========================================================================

#ifdef OUTLINE
#include "FMultivector.icc"
#endif

//---------------------------------------------------------------------------
