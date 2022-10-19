//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSimpleBasicODE.cc,v $
// Language:  C++
// Date:      $Date: 2003/03/27 08:29:52 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FSimpleBasicODE.hh"
#include "FSimpleODEController.hh"
#include "FVector.hh"
#include "math.h"
#include "FException.hh"

//===========================================================================

FSimpleBasicODE::~FSimpleBasicODE()
{
}

//---------------------------------------------------------------------------

FSimpleBasicODE::FSimpleBasicODE(int d, FSimpleODEController* controller) 
    : FObject (), dim(d)
      
{
  this->controller=controller;
}

//---------------------------------------------------------------------------

FSimpleBasicODE::FSimpleBasicODE()
{
}
			      
//---------------------------------------------------------------------------

int FSimpleBasicODE::stepsTaken() const
{ return currentStepNumber; }

int FSimpleBasicODE::goodSteps() const
{ return currentGoodSteps; }

int FSimpleBasicODE::badSteps() const
{ return currentStepNumber - currentGoodSteps; }

//===========================================================================
