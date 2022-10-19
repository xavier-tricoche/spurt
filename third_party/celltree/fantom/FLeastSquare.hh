//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FLeastSquare.hh,v $
// Language:  C++
// Date:      $Date: 2001/08/02 14:54:34 $
// Author:    $Author: bessen $
// Version:   $Revision: 1.5 $
//
//--------------------------------------------------------------------------- 

#ifndef __FLeastSquare_hh
#define __FLeastSquare_hh

#include <cmath>

#include "FVector.hh"
#include "FMatrix.hh"
#include "FMath.hh"

//===========================================================================

/** 
 * The FLeastSquare class provides implementations of methods solving
 * least-square problems.
 */

class FLeastSquare
{
public:

  /**
   * \brief
   * Constructor, initialize object with reasonable values
   * \pre
   * None.
   * \post
   * Problem set up for solving.
   * \param
   * x Set of data points.
   * \param
   * y Set of data points.
   * \param
   * sig Individual standard deviations for data points.
   * \param
   * funcs Routine funcs(x,a,yfit,dyda) which evaluates the fitting function
   * yfit and its derivatives dyda with respect to the fitting parameters
   * a at x.
   */
  FLeastSquare(const FVector& x, const FVector& y, 
	       const FVector& sig, 
 	void (*funcs)(double, const FVector&, double&, FVector&));

  /** Destructor
   */
  ~FLeastSquare();

  /** Solver, uses Levenberg-Marquardt method.
   * \\\par Description:
   * Calculate best-fit values for the paramters with Levenberg-Marquardt. 
   * On the first call provide an initial guess for the parameters a and
   * Chi square. Set alamda<0 for initialization (sets alamda=.001).
   * If a step succeeds chisq becomes smaller and alamda decreases by 
   * a factor of 10. If a step fails alamda grows by a factor of 10.
   * You must call this routine repeatedly until convergence is achieved.
   * [ Then make one final call with alamda=0 so that covar returns the
   * covariance matrix and alpha the curvature matrix. ** NOT IMPLEMENTED **]
   * \pre
   * FLeastSquare initialized.
   * \post
   * \param
   * ia Input array,non-zero values indicate components that schould be fitted.
   * \param
   * a Parameter array for which best-fit values are to be found.
   * \param
   * chisq Chi square.
   * \param
   * alamda Value of alamda. 
   */
  void levenbergMarquardt(FVector& a, double& chisq, double& alamda);

private:

  /** Default constructor. (shouldn´t be used...)
   */

  FLeastSquare();
  
  void mrqcof(const FVector& x, const FVector& y, const FVector& sig,
	      const FVector& a,  
	      FMatrix& alpha, FVector& beta, double& chisq,
	      void (*funcs)(double, const FVector&, double&, FVector&));

  void covsrt(FMatrix& covar);

  FVector x,y,sig;
  void (*funcs)(double, const FVector&, double&, FVector&);
  FVector atry, beta, da;
  FMatrix oneda, covar, alpha;
  double ochisq;
};

//===========================================================================

#ifndef OUTLINE
#include "FLeastSquare.icc"
#endif

#endif 
