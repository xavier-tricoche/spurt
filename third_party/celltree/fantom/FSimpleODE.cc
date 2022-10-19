//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSimpleODE.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:54:48 $
// Author:    $Author: garth $
// Version:   $Revision: 1.7 $
//
//--------------------------------------------------------------------------- 

#include "FSimpleODE.hh"
#include "FSimpleODEController.hh"
#include "FVector.hh"
#include "math.h"
#include "FException.hh"

//===========================================================================

FSimpleODE::~FSimpleODE()
{
}

//---------------------------------------------------------------------------

FSimpleODE::FSimpleODE(int d, FSimpleODEController* controller) 
    : FSimpleBasicODE(d, controller), currentVal(d), dydx(d),
      tmpVal(d), currentValSav(d), tmpValSav(d), dydxSav(d), stepSizeControl(false)
{
}

//---------------------------------------------------------------------------

FSimpleODE::FSimpleODE()
{
}

//---------------------------------------------------------------------------
// very simple polygon-shoot method ( but may be fast indeed :) )
//---------------------------------------------------------------------------
int FSimpleODE::integrate ( FVector& ystart, double initialStepSize, int steps, double t, double startT)
  
  // integrate starting values ystart[0..N-1] 
  //
  // INPUT
  // =====
  //   ystart[0..N-1]: starting values of dependent variables
  //   minDist: minimal distance between two consecutive points along
  //            the integral curve.
  //
  // OUTPUT
  // ======
  //   ystart              : values at the end of integration interval
  //   yp[0..N-1][0..kmax] : intermediate values
  //
{
  try {
    // we initialize the integrator.
    currentVal=ystart;

    // set the starting parameter to zero
    currentParam = startT;
    
    // currentStepsize
    currentStepsize = initialStepSize;
    
    bool derivativeResult;
    derivativeResult  = controller->ODEDerivatives( dydx,
							currentStepNumber,
							currentParam,
							currentVal );
    if(!derivativeResult) return 5;

    // offer the current state to the controller
    controller->ODEStepOffered( currentStepNumber, currentParam,
                                ystart, dydx);
    // we take as many steps as we are allowed...
    currentStepNumber = 0;
    currentGoodSteps = 0;

    // we use the step-function that does the actual work
    int tmpInt = takeSteps( ystart, steps, t );
    return tmpInt;
  }
  catch (FException& e) {
    e.addTraceMessage("int FSimpleODE::integrate ( FVector& ystart, double initialStepSize, int steps)");
    throw;
  }
}
			      
//---------------------------------------------------------------------------

int FSimpleODE::takeSteps ( FVector& ystart, int steps, double /* t */)
  
  // integrate starting values ystart[0..N-1] 
  //
  // INPUT
  // =====
  //   ystart[0..N-1]: starting values of dependent variables
  //   minDist: minimal distance between two consecutive points along
  //            the integral curve.
  //
  // OUTPUT
  // ======
  //   ystart              : values at the end of integration interval
  //   yp[0..N-1][0..kmax] : intermediate values
  //
{

  try {
    int stopStepCount = currentStepNumber+steps;
    bool derivativeResult;

    while ( (steps == 0) ||
            (currentStepNumber < stopStepCount) )
      {

        // save the current config if we encounter a bad step and
        // have to repeat the last...
        currentValSav=currentVal;
        currentParamSav=currentParam;
        currentStepsizeSav=currentStepsize;
        currentStepNumberSav=currentStepNumber;
        currentGoodStepsSav=currentGoodSteps;
        dydxSav=dydx;

        //----------------------------------------------------------
        //------------ THE STEP ITSELF -----------------------------
        //----------------------------------------------------------
        currentVal.plus( dydx * (0.5*currentStepsize), tmpVal);
        derivativeResult=controller->ODEDerivatives( dydx,
				    currentStepNumber,
                                    currentParam + 0.5*currentStepsize,
                                    tmpVal );
	if(!derivativeResult&&!stepSizeControl) return 5;

        currentVal+=( dydx * currentStepsize );
        currentParam += currentStepsize;
        currentStepNumber ++;

        //----------------------------------------------------------


        //----------------------------------------------------------
        //------------ INITIALIZE THE NEXT STEP --------------------
        //----------------------------------------------------------

        // we have finished our step, now get the derivatives information
        // for the current position

        derivativeResult=controller->ODEDerivatives( dydx,
				    currentStepNumber,
				    currentParam,
				    currentVal );
	if(!derivativeResult&&!stepSizeControl) return 5;
	
        //----------------------------------------------------------


        //----------------------------------------------------------
        //----- CHECK IF THE TERMINATION CRITERION HOLDS  ----------
        //----------------------------------------------------------
        
        // we ask the controller if we should stop the integration
        // if -1 is returned we repeat the next step with modified
        // stepsize
        int odeCheckResult = controller->ODECheck( currentStepNumber, currentParam,
                                           currentVal, dydx );

        //----------------------------------------------------------


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
                                        dydx);
            currentGoodSteps++;
          }
        else {
          // restore the saved status, have to repeat the last...
          currentVal=currentValSav;
          currentParam=currentParamSav;
          currentStepNumber=currentStepNumberSav;
          currentGoodSteps=currentGoodStepsSav;
          dydx=dydxSav;
        }

	// give the user the opportunity to set a new stepsize himself,
	// regardless of the need for it.
	 double currentStepsizeTMP = controller->ODENextStepsize( currentStepsize,
							 currentStepsize,
							 currentStepNumber,
							 currentParam,
							 currentVal,
							 dydx );
	 if(stepSizeControl) currentStepsize=currentStepsizeTMP;
        // if we had a bad step: check if ODENextStepSize has already
        // taken care of the stepsize adjustment, if not, do it ourselves
        if ( currentStepNumber == currentStepNumberSav &&
             currentStepsize   >= currentStepsizeSav )
	  {
	    if(stepSizeControl)currentStepsize= 0.5 * currentStepsizeSav;

	    // here we take care of a pathetic case where we reach
	    // the smallest stepsize that a double can offer.
	    // it is likely that someting BAD happened. stop the process
	    // with -1.
	    if (currentStepsize == currentStepsizeSav)
	      return -1;
	  }
	
      
//////////
//  FIXME mario und alex: t parameter abfrage auskommentiert, soll demnächst nach ODECheck ausgelagert werden
/////////	    
//         if ( finishForReachedTParameter )
//           break;
        
//         // if we wish a certain t parameter to be reached-> make it so !
//         // we make sure we do not miss the point where we should stop
//         // take another step and leave since then we are exactly where
//         // we want to be
//         if (t!=0.0)
//           if (stopParam-currentParam < currentStepsize) {
//             currentStepsize = stopParam-currentParam;
//             finishForReachedTParameter=true;
//           }
        
        //----------------------------------------------------------

        //----------------------------------------------------------
        //---------- LEAVE IF WE GOT A TERMINATION -----------------
        //---------- INDICATION FROM CHECK ()      -----------------
        //----------------------------------------------------------
        
        // values < 0 are reserved...
        if ( odeCheckResult > 0 ) {
          ystart = currentVal;
          return odeCheckResult;
        }

        //----------------------------------------------------------
      }
    ystart = currentVal;
    return 0;
  }
  catch (FException& e) {
    e.addTraceMessage("int FSimpleODE::takeSteps ( FVector& ystart, int steps)");
    throw;
  }
}
			      
//---------------------------------------------------------------------------

const FString& FSimpleODE::getClassName() const
{
  static const FString className("FSimpleODE");

  return className;
} 

//===========================================================================

double FSimpleODE::getCurrentStepSize( void )
{
  return currentStepsize;
}


//===========================================================================
void FSimpleODE::setUseStepSizeControl(bool stepSizeControl )
{
  this->stepSizeControl=stepSizeControl;
}
