//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSimpleODEController.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:54:48 $
// Author:    $Author: garth $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#include "FSimpleODEController.hh"
#include "FArray.hh"



FSimpleODEController::FSimpleODEController ()
{
}

  
FSimpleODEController::~FSimpleODEController ()
{
}

const FString& FSimpleODEController::getClassName() const
{
  static const FString className("FSimpleODE");
  return className;
}

int FSimpleODEController::ODECheck ( int /*stepNumber*/, 
                                     const double /*currentParam*/,
                                     const FArray& /*currentVal*/,
                                     const FArray& dydx)
{
  // if we have extraordinary small values we assume a deadlock and
  // exit...
  if (dydx.norm() < 1.0e-9)
    return( -1 );
  return (0);
}

double FSimpleODEController::ODENextStepsize ( const double& /*currentStepsize*/,
                                               const double& predictedStepsize,
                                               int /*stepNumber*/,
                                               const double &/*currentParam*/,
                                               const FArray& /*currentVal*/,
                                               const FArray& /*dydx*/ )
{
  // the simplest form, just thrust the ode integrator...
  return (predictedStepsize);
}

bool FSimpleODEController::ODEStepOffered ( int /*stepNumber*/, 
                                            const double& /*currentParam*/,
                                            const FArray& /*currentVal*/,
                                            const FArray& /*dydx*/ )
{
  // default is to judge the step as good.
  return true;
}
