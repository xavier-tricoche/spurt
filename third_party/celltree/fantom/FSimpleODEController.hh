//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSimpleODEController.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:54:48 $
// Author:    $Author: garth $
// Version:   $Revision: 1.5 $
//
//--------------------------------------------------------------------------- 

#ifndef __FSimpleODEController_hh
#define __FSimpleODEController_hh

//===========================================================================

#include "FObject.hh"
#include "FException.hh"

#include <cmath>
#include <vector>

#include "FArray.hh"

/** 
 * This is a base class for all ordinary differential equations solvers.
 */
class FSimpleODEController : public FObject
{
public:
  
  
  /**
   * \brief
   * Default constructor.
   */
  FSimpleODEController ();
  
  /**
   * \brief
   * Destructor.
   */
  virtual ~FSimpleODEController ();
  
  /** Returns the class name as string.
   *  \return class name   
   */
  virtual const FString& getClassName() const;
  
  /**
   * \brief
   * Check if we shall go any further.
   * \pre
   * ODE solver created.
   * \post
   * None.
   * \exception
   * None.
   * \param stepNumber Current step number
   * \param currentParam current integration parameter
   * \param currentVal Current functional values of ODE.
   * \param dydx Current derivatives.
   * \return any value greater than 0 indicates that the integration
   * should stop.
   *  -1 indicates that we should step back.
   */
  virtual int ODECheck ( int stepNumber, 
                         const double currentParam,
                         const FArray& currentVal,
                         const FArray& dydx);
  
  /**
   * \brief
   * compute derivatives
   * \pre
   * none
   * \post
   * None.
   * \exception
   * None.
   * \param dydx Current derivatives.
   * \param stepNumber current step number
   * \param currentParam
   * \param currentVal Current functional values of ODE.
   * \return Flag wether calculation should go on or not.
   */
  virtual bool ODEDerivatives ( FArray& dydx,
				int stepNumber,
                                const double& currentParam,
                                const FArray& currentVal ) = 0;

  /**
   * \brief
   * controls the stepsize (if no control is wanted simply return pred.)
   * \param currentStepSize
   * \param predictedStepsize Step size predicted for next step
   * \param stepNumber Current step number
   * \param currentParam Current integration parameter
   * \param currentVal Current functional values of ODE.
   * \param dydx Current derivatives.
   * \return corrected step size
   */
  virtual double ODENextStepsize ( const double& currentStepsize,
                                   const double& predictedStepsize,
                                   int stepNumber,
                                   const double &currentParam, 
                                   const FArray& currentVal,
                                   const FArray& dydx );

  /**
   * \brief
   * this function is called everytime that the ode solver has finished
   * a step. so we can decide whether we take it or not.
   * \param stepNumber Current step number
   * \param currentParam current integration parameter
   * \param currentVal Current functional values of ODE.
   * \param dydx Current derivatives.
   * \return
   * none.
   */
  virtual bool ODEStepOffered ( int stepNumber, 
                                const double& currentParam,
                                const FArray& currentVal,
                                const FArray& dydx );

protected:
  

};

//===========================================================================

#endif // __FSimpleODEController_hh
 
