//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FLeastSquare.cc,v $
// Language:  C++
// Date:      $Date: 2000/05/31 12:51:38 $
// Author:    $Author: jfrey $
// Version:   $Revision: 1.1 $
//
//---------------------------------------------------------------------------

#include "FLeastSquare.hh"

#ifdef OUTLINE
#include "FLeastSquare.icc"
#endif

//---------------------------------------------------------------------------

void FLeastSquare::levenbergMarquardt(FVector& a,
				 double& chisq, double& alamda)
{
  int ma = a.getDimension();

  if(alamda < 0.0)     // Init
    { 
      atry = FVector(ma);
      beta = FVector(ma);
      da   = FVector(ma);

      oneda = FMatrix(ma, 1);

      alpha = FMatrix(ma, ma);
      covar = FMatrix(ma, ma);

      alamda=0.001;

      mrqcof(x, y, sig, a, alpha, beta, chisq, funcs);

      ochisq=chisq;

      for (int j=0; j<ma; j++) atry[j]=a[j]; 

    }

  for(int j=0; j<ma; j++)    // Alter linearized fitting matrix
                             // by augmenting diagonal elements. 
    {
      for(int k=0; k<ma; k++) covar(j,k)=alpha(j,k);

      covar(j,j)=alpha(j,j) * (1.0 + alamda);
      oneda(j,0)=beta[j];
    }

  FMath::FGaussJ(covar, oneda); // Matrix solution

  for(int j=0; j<ma; j++) 
    da[j]=oneda(j,0);

  if(alamda == 0.0)         // Once converged, evaluate covariance 
    {
      covsrt(covar);
      covsrt(alpha); // Spread out alpha to its full size too.
    }
  else
    {
      for(int j=0, l=0; l<ma; l++)                //Did the trial succeed?
	{
	  atry[l]=a[l]+da[j];
	  j++;
	}

      mrqcof(x, y, sig, atry, covar, da, chisq, funcs);
      
      if(chisq < ochisq)    // Success, accept the new solution.
	{
	  alamda *= 0.1;
	  ochisq = chisq;
	  for(int j=0; j<ma; j++) 
	    {
	      for(int k=0; k<ma; k++) alpha(j,k)=covar(j,k);
	      beta[j]=da[j];
	    }
	  
	  for(int l=0; l<ma; l++) a[l]=atry[l];
	} 
      else // Failure, increase alamda and return.
	{
	  alamda *= 10.0;
	  chisq=ochisq;
	}
    }
}

//---------------------------------------------------------------------------

void FLeastSquare::mrqcof(const FVector& x, const FVector& y, 
			   const FVector& sig, const FVector& a, 
			   FMatrix& alpha, 
			   FVector& beta, double& chisq,
	    void (*funcs)(double, const FVector&, double&, FVector&))
{
  int ma = a.getDimension();
  double ymod,wt,sig2i,dy;

  FVector dyda(ma);
 
  for(int j=0; j<ma; j++)   // Initialize (symmetric) alpha, beta.
    {
      for (int k=0; k<j; k++) alpha(j,k)=0.0;
      beta[j]=0.0;
    }

  chisq=0.0;

  for(unsigned int i=0; i<x.getDimension(); i++) // Sum loop over all data.
    {
      (*funcs)(x[i], a, ymod, dyda);

      sig2i=1.0/(sig[i]*sig[i]);
      dy=y[i]-ymod;

      int j,k,m,l;

      for(j=0,l=1; l<=ma; l++) 
	{
	  wt=dyda[l-1]*sig2i;

	  for(j++, m=1, k=0; m<=l; m++)
	    alpha(j-1,(++k)-1) += wt*dyda[m-1];
	  beta[j-1] += dy*wt;
	}
      chisq += dy*dy*sig2i;  // And find chi².
    }
  
  for (int j=1; j<ma; j++)      // Fill in the symmetric side.
    for (int k=0; k<j; k++) alpha(k,j)=alpha(j,k);
}

//---------------------------------------------------------------------------

void FLeastSquare::covsrt(FMatrix& covar)
{
  int ma = covar.getDimensionX();
  int k=ma;
  double c;

  for(int j=ma-1; j>=0; j--)
    {
      for(int i=0; i<ma; i++) MM_SWAP(covar(i,k), covar(i,j), c);
      for(int i=0; i<ma; i++) MM_SWAP(covar(k,i), covar(j,i), c);
      k--;
    }
}
//---------------------------------------------------------------------------
