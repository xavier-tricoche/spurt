//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSimpleODE.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:54:49 $
// Author:    $Author: garth $
// Version:   $Revision: 1.5 $
//
//--------------------------------------------------------------------------- 

#ifndef __FSimpleODE_hh
#define __FSimpleODE_hh

//===========================================================================

#include "FSimpleBasicODE.hh"
#include "FArray.hh"

#include <cmath>
#include <vector>

class FSimpleODEController;


/** 
 * This is a base class for all ordinary differential equations solvers.
 */
class FSimpleODE : public FSimpleBasicODE
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
  FSimpleODE (int d, FSimpleODEController* controller);

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
  virtual ~FSimpleODE ();

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
   * start contains last calculated Position.
   * \param
   * steps - the number of steps that should be taken,
   * t - the t parameter until which integration should go,
   *     zero means no limit
   * \exception
   * none
   * \return the value that controller->check() returned
   *         or -1 if the step could not be computed.
   */
  virtual int integrate( FVector& start,
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
   * \return the value that controller->check() returned
   *         or -1 if the step could not be computed.
   */
  virtual int takeSteps( FVector& result,
			 int steps = 1,
			 double t=0.0 );

  virtual double getCurrentStepSize( void );
  virtual void setUseStepSizeControl( bool stepSizeControl );

protected:

  /**
   * Default constructor.
   */
  FSimpleODE ();

  /// Current value vector.
  FArray currentVal;
  FArray dydx;
  // used to compute the step-equation
  FArray tmpVal;
  double currentParam;
  double currentStepsize;

  /// Current value vector.
  FArray currentValSav;
  // used to compute the step-equation
  FArray tmpValSav;
  double currentParamSav;
  double currentStepsizeSav;
  int currentStepNumberSav;
  int currentGoodStepsSav;
  FArray dydxSav;

  bool stepSizeControl;

};

//===========================================================================

#endif // __FSimpleODE_hh
