//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FBasicODE.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:54:46 $
// Author:    $Author: garth $
// Version:   $Revision: 1.13 $
//
//--------------------------------------------------------------------------- 

#ifndef __FBasicODE_hh
#define __FBasicODE_hh

//===========================================================================

#include "FObject.hh"
#include "FArray.hh"
#include "FException.hh"

// set some constant for algorithms

/// FBasicODE parameter.
#define SAFETY   0.9
/// FBasicODE parameter.
#define PGROW   -0.2
/// FBasicODE parameter.
#define PSHRINK -0.25
/// FBasicODE parameter.
#define ERRCON   1.89e-4	// (5/SAFETY)^(1/PGROW)
/// FBasicODE parameter.
#define TINY     1.0e-30


/** 
 * This is a base class for all ordinary differential equations solvers.
 */
class FBasicODE : public FObject
{
public:
  
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
  virtual ~FBasicODE ();

  /** Returns the class name as string.
   *  \return class name   
   */
  virtual const FString& getClassName() const;

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
  int steps();

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
  int goodSteps();

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
  int badSteps();
  
  /**
   * \brief
   * Set minimal step size (can be zero).
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \param val Minimal step size
   */
  void setMinimalStepSize( const double val );

  /**
   * \brief
   * Get minimal step size.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \return Minimal stepsize
   */
  double getMinimalStepSize() const;

  /**
   * \brief
   * Set the desired accuracy.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \param value Accuracy value.
   */
  void setAccuracy( const double value );

  /**
   * \brief
   * Get the desired accuracy.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \return The desired accuracy
   */
  double getAccuracy();

  /**
   * \brief
   * Sets max number of intermediate steps; saving of intermediate
   * results is switched off when num=0.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \param num Maximal number of intermediate steps.
   */
  void setMaxIntermediateSteps( int num );

  /**
   * \brief
   * Gets max number of intermediate steps.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \return Maximum number of intermediate steps.
   */
  int maxIntermediateSteps();

  /**
   * \brief
   * Sets minimum distance between intermediate steps. Only relevant
   * if maxIntermediateSteps > 0.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \param value Minimal distance.
   */
  void setMinDistanceBetweenIntermediateSteps( double value );

protected:

  /**
   * \brief
   * Check if we shall go any further.
   * \pre
   * ODE solver created.
   * \post
   * None.
   * \exception
   * None.
   * \param x Current location.
   * \param y Current functional values of ODE.
   * \param dydx Current derivatives.
   * \return Flag wether calculation should go on or not.
   */
  virtual int check (const double x,
		     const FVector& y,
		     const FVector& dydx) = 0;

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
  FBasicODE (int d);

  /**
   * Default constructor.
   */
  FBasicODE ();

  ///  dimension of the function do be computed
  int dim; 
 
  /// the desired accurancy
  double eps;

  /// the guessed first stepsize
  double h1; 

  /// minimum allowed stepsize
  double hmin; 

  /// number of 'good' steps
  int nok;
  
  /// number of 'bad' (but retried and fixed) steps
  int nbad;

  /** maximal number of intermediate results; 
   * results only get stored if kmax != 0
   */
  int kmax;

  /// actual number of intermediate steps taken
  int kount; 

  /** because steps occur at unequal intervals results are 
   * only stored at intervals greater than dxsav. */
  double dxsav; 

  /// maximal number of integration steps
  int nstpmax; 

  /// Current error vector.
  FVector currentErr;

  /// Current value vector.
  FVector currentVal;
};

/// Double valued absolute function.
double ABS(double a); 

/// Double valued sign function.
double SIGN( double a, double b );

/// Double valued max function.
double FMAX(double a, double b);

/// Double valued min function.
double FMIN(double a, double b);

//===========================================================================

#ifndef OUTLINE
#include "FBasicODE.icc"
#endif

//===========================================================================

#endif // __FBasicODE_hh
 
