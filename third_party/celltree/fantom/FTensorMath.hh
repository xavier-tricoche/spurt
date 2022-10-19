#ifndef FTensorMath_hh
#define FTensorMath_hh

#include "FVector.hh"
#include "FTensor.hh"

#include <cmath>

#ifndef NODEBUG
#include "eassert.hh"
#endif

#include "FMathValues.hh"

namespace FMath
{
  /**
   * Calculate the Fractional Anisotropy used in medical visualization
   * defined by
   * \f[ \hat{\lambda} = \frac{1}{3}\sum_i^N {\lambda_i} \f]
   * \f[ FA = \sqrt{
   *             \frac{3}{2}
   *             \frac{\sum_i (\lambda_i-\hat{\lambda})^2}
   *                  {\sum_i (\lambda_i)^2}
   *          }
   *  \f]
   * \param eigenvalues 
   * Eigenvalues of tensor
   */
  inline double fractionalAnisotropy( const FVector &eigenvalues )
  {
    double normSquare  = eigenvalues.normSquare();
    if(normSquare < F_SQRT_DBL_EPSILON) return 0;
    const double trace = (eigenvalues[0] + eigenvalues[1] + eigenvalues[2])/3.;
    double fa = sqrt(3./2.*(eigenvalues - FArray(trace, trace, trace)).normSquare() / normSquare);

    if(fa > 1.) { fa =1.;}
    if(fa < 0.) { fa = 0.;}

    return fa;
  }

  inline double anisotropyHelper( const FTensor &t )
  {
    // a1-a2+3.*a3
    // = Dxx^2+Dyy^2+Dzz^2
    //   -(DxxDyy + DxxDzz + DyyDzz)
    //   +3(Dxy^2+Dxz^2+Dyz^2)
    return 
       t(0,0)*t(0,0)+t(1,1)*t(1,1)+t(2,2)*t(2,2)
      -t(0,0)*t(1,1)-t(1,1)*t(2,2)-t(2,2)*t(0,0)
      +(t(0,1)*t(0,1)+t(1,2)*t(1,2)+t(2,0)*t(2,0))*3.;
   }

  inline double fractionalAnisotropy( const FTensor &t )
  {
#if 0
    static double normalization = sqrt(3./2.);
    const double upt = t(0,1) + t(1,2) + t(2,0);
    if(upt < 0.) THROW_EXCEPTION( FException, "Tensor not positive definite" );
    const double up = upt*upt; // square

    double denom = t(0,0)*t(0,0) + t(1,1)*t(1,1)+t(2,2)*t(2,2)+2.*(t(0,1)*t(0,1)+t(0,2)*t(0,2)+t(1,2)*t(1,2));
    if(denom < 1e-9) THROW_EXCEPTION( FException, "Tensor not really positive definite");
    return normalization * sqrt( 1. - ( up/denom));
//#elif 0
    static double normalization = sqrt(3./2.);
    const double upt = anisotropyHelper(t);
    const double denom = t.normSquare();
    if(denom < 1e-5) THROW_EXCEPTION( FException, "Tensor not really positive definite");
    return normalization*sqrt(upt/denom);
#else
    // as used by Xavier Tricoche and Gordon Kindlmann
//#warning check if scaling is right
//  static double normalization = sqrt(2./3.);
  const double cross = t(0,1)*t(0,1)+t(0,2)*t(0,2)+t(1,2)*t(1,2);
  const double j2=t(0,0)*t(1,1)+t(0,0)*t(2,2)+t(1,1)*t(2,2) -cross;
  const double j4=t(0,0)*t(0,0)+t(1,1)*t(1,1)+t(2,2)*t(2,2)+2.*cross;
    if(j4< 1e-5) THROW_EXCEPTION( FException, "Tensor not really positive definite");
    if(j4-j2< 1e-9) THROW_EXCEPTION( FException, "Tensor not really positive definite");
  return sqrt((j4-j2)/j4);
}
#endif
    inline double relativeAnisotropy( const FTensor &t )
    {
    static double normalization = sqrt(3./2.);
    const double upt = anisotropyHelper(t);
    const double denom = t(0,0)+t(1,1)+t(2,2); // t.trace()
    if(denom < 1e-5) THROW_EXCEPTION( FException, "Tensor not really positive definite");
    return normalization*sqrt(upt)/denom;
  }


  /**
   * Function is used by Higher Order Tensor methods.
   * \pre t is a symmetric tensor
   * get the scalar value of t*dir*...*dir where dir is multiplied t.getOrder() times.
   */
  inline double evaluateTensor( const FTensor &t, const FVector& dir )
  {
#ifndef NODEBUG
    eassert( t.getDimension() == dir.getDimension() );
    try{
#endif
      FTensor value(t);
      while( value.getOrder() > 0 )
      {
        value.multAndContract(dir, 0);
      }
      return value();
#ifndef NODEBUG
    }CATCH_N_RETHROW(FException);
    return 0;
#endif
  }

  /**
   * Function is used by Higher Order Tensor methods.
   * \pre t is a symmetric tensor
   * \param value same value as evaluateTensor returns
   * \param gradient if dir is normalized the gradient vector on the unit sphere
   */
  inline void evaluateTensorAndGradient( double &value, FVector &gradient, const FTensor &t, const FVector &dir)
  {
#ifndef NODEBUG
    eassert( t.getDimension() == dir.getDimension() );
    try{
#endif
      FTensor t1(t);
      while(t1.getOrder() > 1)
        t1.multAndContract(dir, 0);
      FTensor t0(t1);
      t0.multAndContract(dir, 0);
      value = t0();
      gradient = ( FVector(t1)-(dir*value)) *(double)t.getOrder();
#ifndef NODEBUG
    }CATCH_N_RETHROW( FException );
#endif
  }



  /**
   * Westin's linear anisotropy \$[c_l\$]
   * vals holds eigenvalues which have to be sorted already
   * vals[0] >= vals[1] >= vals[2] */
  inline double linearAnisotropy( const FArray &vals )
  {
#ifndef NODEBUG
    eassert(vals[0] >= vals[1] && vals[1] >= vals[2]);
#endif
    // basser wrote vals[0]-vals[2] which is wrong!
    // these are westins values that c_l + c_p + c_s = 1.0
    return (vals[0]-vals[1])/(vals[0]+vals[1]+vals[2]);
  }

  /** 
   * Westin's planar anisotropy \$[c_p\$]
   * vals holds eigenvalues which have to be sorted already
   * vals[0] >= vals[1] >= vals[2] */
  inline double planarAnisotropy( const FArray& vals )
  {
#ifndef NODEBUG
    eassert(vals[0] >= vals[1] && vals[1] >= vals[2]);
#endif
    return 2.0*(vals[1]-vals[2])/(vals[0]+vals[1]+vals[2]);
  }

  /** 
   * Westin's spherical anisotropy \$[c_s\$]
   * vals holds eigenvalues which have to be sorted already
   * vals[0] >= vals[1] >= vals[2] */
  inline double sphericalAnisotropy( const FArray& vals )
  {
#ifndef NODEBUG
    eassert(vals[0] >= vals[1] && vals[1] >= vals[2]);
#endif
    return 3.0*vals[2]/(vals[0]+vals[1]+vals[2]);
  }
    
};

#endif
