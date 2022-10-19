//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSimpleBasicODE.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:54:48 $
// Author:    $Author: garth $
// Version:   $Revision: 1.5 $
//
//--------------------------------------------------------------------------- 

#ifndef __FSimpleBasicODE_hh
#define __FSimpleBasicODE_hh

//===========================================================================

#include "FObject.hh"
#include "FArray.hh"

#include <cmath>
#include <vector>

class FSimpleODEController;


/** 
 * This is a base class for all ordinary differential equations solvers.
 */
class FSimpleBasicODE : public FObject
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
  FSimpleBasicODE (int d, FSimpleODEController* controller);

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
  virtual ~FSimpleBasicODE ();

  /** Returns the class name as string.
   *  \return class name   
   */
  virtual const FString& getClassName() const = 0;

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
   * \param start
   *    The function value to begin the integration at.
   * \param initialStepSize
   *    A first guess at how big the steps should be.
   *    In case of a second order integrator, this value may be
   *    the one used througout the integration.
   * \param steps
   *    The number of steps that may be taken at most. The default
   *    value of 0 indicates that there is no limit.
   * \param t
   *    The t parameter until which integration should go,
   *    zero means no limit. (default)
   * \exception
   *    none
   * \return
   *    the value that controller->check() returned
   *    or -1 if the step could not be computed.
   */
  virtual int integrate( FArray& start,
                         double initialStepSize = 1.0e-5,
                         int steps = 0, 
			 double t=0.0,
			 double startT = 0.) = 0;

  /**
   * \brief
   * Continues the work that integrate() began.
   * \pre
   * none
   * \post
   * resultValues contains the values calculated in the steps.
   * resultParameters contains the t parameters of the steps.
   * result contains last calculated Position.
   * \param steps
   *  The number of steps that shall be taken. Default is 1, which means
   *  the function returns after each step.
   * \param t
   *  the t parameter distance from the current t to integrate
   *  zero means no limit
   * \exception
   * none
   * \return the value that controller->check() returned
   *         or -1 if the step could not be computed.
   */
  virtual int takeSteps( FArray& result,
			 int steps = 1, 
			 double t=0.0) = 0;

  /**
   * \brief
   * Returns the number of calculated intermediate results.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \return Number of steps of last ODE solving or -1 if no
   * solving took place so far
   */
  int stepsTaken() const;

  /**
   * \brief
   * Get the number of good steps.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \return Number of 'good' steps of last ODE solving or -1 if no
   * solving took place so far
   */
  int goodSteps() const;

  /**
   * \brief
   * Get the number of bad steps.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \return Number of 'bad' (but retried and fixed) steps of last ODE
   * solving or -1 if no solving took place so far 
   */
  int badSteps() const;

  virtual double getCurrentStepSize( void ) = 0;
  
protected:

  /**
   * Default constructor.
   */
  FSimpleBasicODE ();

  ///  dimension of the function do be computed
  int dim;

  ///  pointer to the class that controls the ode integration
  FSimpleODEController* controller;

  int currentStepNumber;
  int currentGoodSteps;
  
};

//===========================================================================

#endif // __FSimpleBasicODE_hh
