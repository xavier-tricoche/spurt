//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FBasicODE.cc,v $
// Language:  C++
// Date:      $Date: 2003/02/06 10:22:32 $
// Author:    $Author: garth $
// Version:   $Revision: 1.4 $
//
//--------------------------------------------------------------------------- 

#include "FBasicODE.hh"

#ifdef OUTLINE
#include "FBasicODE.icc"
#endif

//===========================================================================

FBasicODE::~FBasicODE()
{
}

//---------------------------------------------------------------------------

FBasicODE::FBasicODE(int d)
{
  dim = d;

  eps  = 1.0E-7;
  h1   = 0.01;
  hmin = 1.0E-5;

  nok   = -1;
  nbad  = -1;

  kmax  = 2000;
  kount = -1;

  dxsav = 0.1;

  nstpmax = 5000;
}

//===========================================================================
