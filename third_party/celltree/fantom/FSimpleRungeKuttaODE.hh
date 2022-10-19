//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSimpleRungeKuttaODE.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:54:49 $
// Author:    $Author: garth $
// Version:   $Revision: 1.7 $
//
//--------------------------------------------------------------------------- 

#ifndef __FSimpleRungeKuttaODE_hh
#define __FSimpleRungeKuttaODE_hh

//===========================================================================

#include "FSimpleBasicODE.hh"
#include "FArray.hh"

#include <cmath>
#include <vector>

class FSimpleODEController;


/** 
 * This is a base class for all ordinary differential equations solvers.
 */
class FSimpleRungeKuttaODE : public FSimpleBasicODE
{
public:

  /**
   * \brief
   * Constructor.
   * \pre
   * none
   * \post
   * ODE solver created.
   * \exception
   * none
   * \param d Dimension of the desired function.
   */
  FSimpleRungeKuttaODE (int d, FSimpleODEController* controller);

  /**
   * \brief
   * Destructor.
   * \pre
   * none
   * \post
   * ODE solver data deleted.
   * \exception
   * none
   */
  virtual ~FSimpleRungeKuttaODE ();

  /** Returns the class name as string.
   *  \return class name   
   */
  virtual const FString& getClassName() const;

  /**
   * \brief
   * integrates the function until controller->check() tells us to stop.
   * the function evaluation is done by controller->derivatives()
   * the stepsize is piped through controller->nextStepSize() to give the
   * user the possibility to modify it
   * each step taken is reported to controller->ODEStepOffered()
   * \pre
   * none
   * \post
   * ystart contains last calculated Position.
   * \param ystart starting values of dependent variables
   * \param steps the number of steps that should be taken,
   * \param t the t parameter until which integration should go,
   *     zero means no limit
   * \exception
   * none
   * \return the value (FNewStreamlineBase::INTEGRATION_RESULT) that  FSimpleODEController::ODECheck() returned
   *         or -1 if the integrator failed to take a step
   */
  virtual int integrate( FVector& ystart,
			 double initialStepSize = 1.0e-5,
			 int steps = 0, 
			 double t=0.0,
			 double startT=0.);

  /**
   * \brief
   * integrates the function a number of steps or .
   * \pre
   * none
   * \post
   * resultValues contains the values calculated in the steps.
   * resultParameters contains the t parameters of the steps.
   * start contains last calculated Position.
   * \param
   * t - the t parameter distance from the current t to integrate
   *     zero means no limit
   * \exception
   * none
   * \return the value  (FNewStreamlineBase::INTEGRATION_RESULT) that FSimpleODEController::ODECheck() returned
   *         or -1, if the integrator failed to take the step
   */
  virtual int takeSteps( FVector& result, 
			 int steps = 1, 
			 double t=0.0);

  virtual double getCurrentStepSize( void );
  
  // numerical tolerance for integration
  static double epsilon;

  // numerical tolerance to monitor integration progression in physical space
  static double minStepSize;
  
protected:

  /**
   * Default constructor.
   */
  FSimpleRungeKuttaODE ();


  /**
   * \brief
   * Fifth-order <b>R</b>unge-<b>K</b>utta <b>Q</b>uality-controlled <b>S</b>tep (<b>rkqs</b>) 
   * with monitoring of local truncation
   * error to ensure accuracy and adjust stepsize.
   * This routine was taken from Press/Teukolsky/Vetterling/Flannery:
   * 'Numerical Recipes in C', Cambridge Press, p.719.
   *
   * \pre
   * none
   * \post
   * ode solver created
   * \exception
   * none
   * \param y     dependent variable vector y[0..N-1], will be replaced
   *              by their new values
   * \param dydx  derivative dydx[0..N-1] at the starting value of the
   * \param x     independent variable x, will be replaced by its new value
   * \param htry  stepsize to be attempted
   * \param yscal the vector yscal[0..N-1] against which the error is
   *              scaled
   * \param hdid  stepsize that was actually accomplished
   * \param hnext extimated next stepsize
   */
  bool rkqs (FVector& y, double& x, const double& htry,
                                   double& hdid, double& hnext );
  
  /**
   * \brief
   * <b>R</b>unge-<b>K</b>utta <b>C</b>ash-<b>K</b>arp (<b>rkck</b>) step
   * Given values for n variables y[0..N-1] and their derivatives
   * dydx[0..N-1] known at x, use the fifth-order Cash-Karp
   * Runge-Kutta method to advance the solution over an interval h and
   * return the incremented variables as yout[0..N-1]. Also return an
   * estimate of the local truncation error in yerr[0..N-1] using the
   * embedded fourth-order method.
   * This routine was taken from Press/Teukolsky/Vetterling/Flannery:
   * 'Numerical Recipes in C', Cambridge Press, p.719f.
   *
   * \pre
   * none
   * \post
   * ode solver created
   * \exception
   * none
   * \param y    dependent variable vector y[0..N-1]
   * \param dydx derivative dydx[0..N-1] at the starting value of the
   * \param x    independent variable x.
   * \param h    stepsize
   * \param yout incremented variables (output)
   * \param yerr estimate of the local truncation error (output)
   */
  bool rkck (FVector& y, FVector& dydx, const double& x,
             const double& h, FVector& yout, FVector& yerr );
  
  
//===========================================================================
//===========================================================================
//================ V A R I A B L E S   ! ! ==================================
//===========================================================================
//===========================================================================
  
  /// Current value vector.
  FArray currentVal;
  // used to compute the step-equation
  FArray tmpVal;
  double currentParam;
  double currentStepsize;
  FArray dydx;


  /// Current value vector.
  FArray currentValSav;
  // used to compute the step-equation
  FArray tmpValSav;
  double currentParamSav;
  double currentStepsizeSav;
  int currentStepNumberSav;
  int currentGoodStepsSav;
  FArray dydxSav;
  
  // could (should?) be taken into account too
  // variables to use if rkck is called with the same values, so
  // the many multiplications are not made again...
  double rkck_x;
  double rkck_h;
  FVector rkck_y;
  FVector rkck_out;
  FVector rkck_err;
  // a flag indicating that the variables have been set in this run..
  bool reuseValuesValid;

  FVector rkck_ak2, rkck_ak3, rkck_ak4,
    rkck_ak5,rkck_ak6,
    ytempRKCK, ytempRKQS;

  // values that are made global in order to support the "suspend" mode..
  FVector yScale;
  // the old y value
  FVector ySave;
  FVector yError;


  
};

//===========================================================================

#endif // __FSimpleRungeKuttaODE_hh
