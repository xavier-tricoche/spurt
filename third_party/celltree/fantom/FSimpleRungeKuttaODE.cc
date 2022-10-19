//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSimpleRungeKuttaODE.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/26 13:36:01 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.14 $
//
//--------------------------------------------------------------------------- 

#include "FSimpleRungeKuttaODE.hh"
#include "FSimpleODEController.hh"
#include "FVector.hh"
#include "math.h"
#include "FException.hh"

#include <iostream>

static double SAFETY=0.9;
static double PGROW= -0.2;
static double PSHRINK=-0.25;
static double ERRCON =  1.89e-4;	// (5/SAFETY)^(1/PGROW)
static double TINY=1.0e-30;

double FSimpleRungeKuttaODE::minStepSize = 1.0E-10;
double FSimpleRungeKuttaODE::epsilon = 1.0E-9;

inline double FMAX(double a, double b)
{
  return (a>b ? a : b);
}

//---------------------------------------------------------------------------

inline double FMIN(double a, double b)
{
 return (a<b ? a : b);
}

//===========================================================================

FSimpleRungeKuttaODE::~FSimpleRungeKuttaODE()
{
}

//---------------------------------------------------------------------------

FSimpleRungeKuttaODE::FSimpleRungeKuttaODE(int d,
                                           FSimpleODEController* controller)
  : FSimpleBasicODE(d, controller), currentVal(d), tmpVal(d), dydx(d),
    currentValSav(d), tmpValSav(d), dydxSav(d), rkck_y(d),
    rkck_out(d), rkck_err(d), rkck_ak2(d), rkck_ak3(d), rkck_ak4(d),
    rkck_ak5(d), rkck_ak6(d), ytempRKCK(d), ytempRKQS(d),
    yScale(d), ySave(d), yError(d)
{
  currentStepsize = 1.0e-5;
}

//---------------------------------------------------------------------------

FSimpleRungeKuttaODE::FSimpleRungeKuttaODE()
{
}

//---------------------------------------------------------------------------
// very simple polygon-shoot method ( but may be fast indeed :) )
//---------------------------------------------------------------------------
int FSimpleRungeKuttaODE::integrate ( FVector& ystart,
                                      double initialStepSize,
                                      int steps, 
				      double t, 
				      double startT)
{
  
  try {

    // we initialize the integrator.
    currentVal=ystart;

    // set the starting parameter to zero
    currentParam = startT;
 
    // currentStepsize
    currentStepsize = initialStepSize;

    controller->ODEDerivatives( dydx,
				currentStepNumber,
                                currentParam,
                                ystart );

    // check the first value, maybe it is out of grid
    int result = controller->ODECheck( 0, currentParam, ystart, dydx);
    if(result > 0)
      return result;

    // offer the current state to the controller
    controller->ODEStepOffered( currentStepNumber, currentParam,
                                ystart, dydx);
    
    // we take as many steps as we are allowed...
    currentStepNumber = 0;
    currentGoodSteps = 0;

    // we use the step-function that does the actual work
    int tmpInt = takeSteps( ystart, steps, t);
    return tmpInt;
  }
  CATCH_N_RETHROW(FException);
}

//---------------------------------------------------------------------------

int FSimpleRungeKuttaODE::takeSteps ( FVector& ystart, int steps, double/* t */ )
{
  double hnext, hdid;
  bool rkqsResult;
  
  try {
    // tell the integrator how far it should integrate by steps
    int stopStepCount = currentStepNumber+steps;

    // tell the integrator how far it should integrate by parameter space
//    double stopParam=currentParam+t;
//    bool finishForReachedTParameter=false;
    while ( (steps == 0) ||
            (currentStepNumber < stopStepCount) )  {

      for( int i=0; i<dim; i++ ) {
        yScale[i] = fabs(currentVal[i])+fabs(dydx[i]*currentStepsize)+TINY;
        // scaling used to monitor accuracy. This general-purpose
        // choice can be modified, if need be.
      }
      
      // save the current config if we encounter a bad step and
      // have to repeat the last...
      currentValSav=currentVal;
      currentParamSav=currentParam;
      currentStepsizeSav=currentStepsize;
      currentStepNumberSav=currentStepNumber;
      currentGoodStepsSav=currentGoodSteps;
      dydxSav=dydx;

      //----------------------------------------------------------
      //---------- TAKE THE RUNGE KUTTA STEP !!!  --------------
      //----------------------------------------------------------

      do {
	rkqsResult = rkqs( currentVal, currentParam,
			   currentStepsize, hdid, hnext );

	// if taking the step failed, me must react.
	// that is we return -1 to signalize "hey, internal integrator failure."
	if ( !rkqsResult )
	  {
	    currentStepsize /= 2;

	    if ( fabs( currentStepsize ) <= FSimpleRungeKuttaODE::minStepSize )
        {
	  // INTEGRATOR_COLLAPSED
	  return -1; // was 4 but -1 is integrator collapsed and that
                     // is what we want here.
        }
	  }
      } while ( !rkqsResult );
      
	

      // update the dydx value, since it is needed (in the next step, too)
      /*bool odeDerivativesResult =*/
      controller->ODEDerivatives( dydx,
                				  currentStepNumber,
                                  currentParam,
                                  currentVal );
      
      // If the step was good we even increment goodsteps :)
      if( hdid == currentStepsize )
        currentGoodSteps++;
      
      // we ask the controller if we should stop the integration
      // if -1 is returned we repeat the next step with modified
      // stepsize
      int odeCheckResult = controller->ODECheck( currentStepNumber, currentParam,
                                         currentVal, dydx );

      
      //----------------------------------------------------------
      //---------- OFFER A GOOD STEP TO THE CLIENT  --------------
      //---------- OR TAKE LAST STEP AGAIN WITH NEW STEPSIZE -----
      //----------------------------------------------------------
      
      // offer the current state to the controller
      if ( odeCheckResult != -1 )
        {
          // if this step was good we offer it to the client
          controller->ODEStepOffered( currentStepNumber,
                                      currentParam,
                                      currentVal,
                                      dydx );
        }
      else {
        // restore the saved status, have to repeat the last...
        currentVal=currentValSav;
        currentParam=currentParamSav;
        currentStepNumber=currentStepNumberSav;
        currentGoodSteps=currentGoodStepsSav;
        dydx=dydxSav;
      }
      
 
      currentStepsize = controller->ODENextStepsize( hdid,
                                                     hnext,
                                                     currentStepNumber,
                                                     currentParam,
                                                     currentVal,
                                                     dydx );
      
      if( fabs(currentStepsize) < FSimpleRungeKuttaODE::minStepSize ) {
	// INTEGRATOR_COLLAPSED
        return -1;
      }
      
      // values < 0 are reserved error values of the ode conrollers ...
      if ( odeCheckResult > 0 ) 
      {
        ystart = currentVal;
        return odeCheckResult;
      }
     
      currentStepNumber ++; 
    } // while

    ystart = currentVal;

    // if steps == 0, we must have left the loop after a step size underflow
    // that has not been caused by a vanishing vector norm (cf. ODECheck):
    // We have left the grid!
//    if (!steps)
//      return 4;

    return 0;
  }
  CATCH_N_RETHROW(FException);
}
			      
//---------------------------------------------------------------------------

const FString& FSimpleRungeKuttaODE::getClassName() const
{
  static const FString className("FSimpleRungeKuttaODE");

  return className;
} 
     
//---------------------------------------------------------------------------

bool FSimpleRungeKuttaODE::rkqs (FVector& y,
                                 double& x,
                                 const double& htry,
                                 double& hdid,
                                 double& hnext )
{

  double errmax, h, htemp, xnew;
  
  h = htry;


  //   cout << "entering rkqs" << endl;
  try {
    bool rkckResult;
    for (;;) {

      rkckResult = rkck( y, dydx, x, h, ytempRKQS, yError );	// take a step


      if (!rkckResult) {
//	cout << "." << flush;
	return false;
      }

      
      //------------------------------------------------------------------
      // evaluate accuracy
      //------------------------------------------------------------------
      errmax = 0.0;
      for ( int i=0; i<dim; i++ )
        errmax = FMAX( errmax, fabs(yError[i]/yScale[i]) );
      
      //------------------------------------------------------------------
      // scale relative to required tolerance
      //------------------------------------------------------------------
      errmax /= FSimpleRungeKuttaODE::epsilon;
      
      if (errmax > 1.0) {

        // truncation error too large, reduce stepsize
        
        htemp = h * SAFETY  * pow(errmax, PSHRINK);
        
        h =
          (h>=0.0)
          ?
          FMAX(htemp, (double)0.1*h)
          :
          FMIN(htemp, (double)0.1*h) ;
          
        // no reduction bigger than factor 10

        xnew = x + h;

        if( xnew == x )
          std::cerr << "rkqs : stepsize underflow in rkqs" << std::endl;
        
        continue; // another try
      }
      else { // errmax <= 0
        
        // step succeeded. Compute size of next step
        // reset error values
        if( errmax > ERRCON ) 
          hnext = SAFETY*h*pow( errmax, PGROW );
        else    
          hnext = 5.*h;
        
        // save the stepsize that has been taken
        hdid = h;
        // and accumulate x
        x += h;
        for( int i=0; i<dim; i++ )
          y[i] = ytempRKQS[i];
        
        // we'll check now, wether that lies outside the boundary
        FVector tmp1(dim), tmp2(dim);
        FVector dydxtmp(dydx);
        controller->ODEDerivatives( dydxtmp,
				    currentStepNumber,
                                    x,
                                    ytempRKQS );
        while (1) {
          // reset error values;
	  rkckResult = rkck(ytempRKQS, dydxtmp, x, hnext, tmp1, tmp2); 
		  
	  // If the step could not be computed ( out of definition range )
	  // we lower the stepsize, otherwise we can go on (break)
	  if (rkckResult)
	    break;
		  
	  hnext = 0.5 * hnext;
	  if (hnext < FSimpleRungeKuttaODE::minStepSize)
	    {
//	      cout << "step size underflow in rkqs" << endl;
	      return false;
	      // no way to set step size without leaving the grid
	    }
        }

        if (hnext < FSimpleRungeKuttaODE::minStepSize) {
#ifndef NODEBUG
          std::cerr << "step size too small" << std::endl;
#endif
          hnext = FSimpleRungeKuttaODE::minStepSize;
        }

	break;
      }
    }
  }
  CATCH_N_RETHROW(FException);

  return true;
}

//---------------------------------------------------------------------------

bool FSimpleRungeKuttaODE::rkck ( FVector& y,
				  FVector& dydx,
				  const double& x,
				  const double& h,
				  FVector& yout,
				  FVector& yerr )
{

  bool success;
  
  try {
    
    if (rkck_h == h && rkck_y == y) {
      yerr = rkck_err;
      yout = rkck_out;
      
      //     cout << "reuse" << endl;
      
      return true;
    }
    
    rkck_h = h;
    rkck_x = x;
	
    //   cout << "rkck_y is set to y: ";
    rkck_y = y;
    //   cout << "new value is " << rkck_y << endl;
    
    //   cout << endl << "in rkck: h = " << h << ", y = " << y << ", x = " << x << ", dydx = " << dydx << endl;
    
    int i;
    
    static double a2=0.2, a3=0.3, a4=0.6, a5=1.0, a6=0.875, 
      b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42=-0.9, 
      b43=1.2, b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0, 
      b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0, 
      b64=44275.0/110592.0, b65=253.0/4096.0, 
      c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0, 
      dc5=-277.0/14336.0;
    
    double dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, 
      dc4=c4-13525.0/55296.0, dc6=c6-0.25;
    
    
    
    // first step
    for( i=0; i<dim; i++ )
      ytempRKCK[i] = y[i]+b21*h*dydx[i];
    
    // second step
    //    computeDerivatives (rkck_ak2, x+a2*h, ytempRKCK);
	
    success = controller->ODEDerivatives( rkck_ak2, currentStepNumber,
					  x+a2*h, ytempRKCK );
    if ( !success )
      {
	rkck_h = -123456789.987654321;

	return false;
      }
	
	
    
    for( i=0; i<dim; i++ )
      ytempRKCK[i] = y[i]+h*(b31*dydx[i]+b32*rkck_ak2[i]);
    
    // third step
    //    computeDerivatives (rkck_ak3, x+a3*h, ytempRKCK);
    success = controller->ODEDerivatives( rkck_ak3, currentStepNumber,
					  x+a3*h, ytempRKCK );
    if ( !success )
      {
	rkck_h = -123456789.987654321;

	return false;
      }
    
    for( i=0; i<dim; i++ )
      ytempRKCK[i] = y[i]+h*(b41*dydx[i]+b42*rkck_ak2[i]+b43*rkck_ak3[i]);
    
    // fourth step
    //    computeDerivatives (rkck_ak4, x+a4*h, ytempRKCK);
    success = controller->ODEDerivatives( rkck_ak4, currentStepNumber,
					  x+a4*h, ytempRKCK );
    if ( !success )
      {
	rkck_h = -123456789.987654321;
	return false;
      }
    
    for( i=0; i<dim; i++ )
      ytempRKCK[i] = y[i]+h*(b51*dydx[i]+b52*rkck_ak2[i]+b53*rkck_ak3[i]
                             +b54*rkck_ak4[i]);
    
    // fifth step
    //    computeDerivatives (rkck_ak5, x+a5*h, ytempRKCK);
    success = controller->ODEDerivatives( rkck_ak5, currentStepNumber,
					  x+a5*h, ytempRKCK );
    if ( !success )
      {
	rkck_h = -123456789.987654321;

	return false;
      }
    
    for( i=0; i<dim; i++ )
      ytempRKCK[i] = y[i]+h*(b61*dydx[i]+b62*rkck_ak2[i]+b63*rkck_ak3[i]
                             +b64*rkck_ak4[i]+b65*rkck_ak5[i]);
    
    // sixth step
    //    computeDerivatives (rkck_ak6, x+a6*h, ytempRKCK);
    success = controller->ODEDerivatives( rkck_ak6, currentStepNumber,
					  x+a6*h, ytempRKCK );
    if ( !success )
      {
	rkck_h = -123456789.987654321;

	return false;
      }
    
    // get yout via 5th order
    for( i=0; i<dim; i++ )
      yout[i]=y[i]+h*(c1*dydx[i]+c3*rkck_ak3[i]+c4*rkck_ak4[i]
                      +c6*rkck_ak6[i]);
    
    // estimate error by difference between 4th and 5th order
    for( i=0; i<dim; i++ )
      yerr[i]=h*(dc1*dydx[i]+dc3*rkck_ak3[i]+dc4*rkck_ak4[i]
                 +dc5*rkck_ak5[i]+dc6*rkck_ak6[i]);
    
    // save values for reuse
    rkck_err = yerr;
    rkck_out = yout;


    return true;
  }
  catch (FException& e) {
	
    // no valid error value could be computed for reuse
    // set rkck_h to dummy value
    rkck_h = -123456789.987654321;
	
    e.addTraceMessage("void FSimpleRungeKuttaODE::rkck (FVector& y, FVector& dydx, const double& x, const double& h, FVector& yout, FVector& yerr )");

    std::cout << e << std::endl;

    throw;
  }
} 

//===========================================================================

double FSimpleRungeKuttaODE::getCurrentStepSize( void )
{
  return currentStepsize;
}

