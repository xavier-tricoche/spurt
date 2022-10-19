#ifndef FSpecialFunctions_hh
#define FSpecialFunctions_hh

#include "FException.hh"

#include <cmath>
#include <vector>
#include "FMathValues.hh"
#include <complex>
namespace FMath
{
  DEFINE_DEFAULTMSG_EXCEPTION_TYPE(FDomainError,"Argument outside function domain.");
  DEFINE_DEFAULTMSG_EXCEPTION_TYPE(FOverflowError,"Argument outside function domain.");


  
#define IS_ODD(a) (a%2==0 ? false : true )

  /** the sinc function sinc(x) = 1.0 iff x==0 , sin(M_PI*x)/M_PI/x else
   */
  inline double sinc(const double x)
  {
    return (x==0) ? 1. : sin(M_PI*x)/M_PI/x;
  }

  namespace laguerre
  {
    //double laguerre1(const double a, const double x);
    //double laguerre2(const double a, const double x);
    //double laguerre3(const double a, const double x);
    //double laguerreN(const int n, const double a, const double x);

    /** Laguerre function of nu=1/2
     *
     * source http://en.wikipedia.org/wiki/Rice_distribution
     */
    double laguerre1_2(const double x );
  }

  
  //! error function namespace
  namespace erf
  {
    std::complex<double> cerfc_continued_fraction(const std::complex<double>& z);
    std::complex<double> cerf_series(const std::complex<double>& z);
    std::complex<double> cerf_rybicki(const std::complex<double>& z);
    std::complex<double> cerf_continued_fraction(const std::complex<double>& z);

    /** compute the complex error function by calling one of the functions
        above
        **/
    std::complex<double> erf(const std::complex<double>& z);

    // the real valued error function erf(double x)
    // should be part of C99 standard
  }
  
  namespace zeta
  {
    double hzeta(const double s, const double q);
  }

  namespace psi
  {
    double psi_int(const int n);
    double psi_x(const double x);
    double psi_n(const int n, const double x);
  }

  namespace gamma
  {
    double lngamma_sgn(double& sgn, const double x);
    //static double gamma_xgthalf(const double x);
    double gammastar_ser( const double x);
    double gamma(const double x); // simple implementation may be unprecise
  }

  namespace bessel
  {
    /** zeroth order bessel function
     *  only implemented up to fabs(x) < 4.0
     *
     *  untested
     */
    double bessel_J0(const double x);
    double bessel_I0(const double x);
    double bessel_I1(const double x);
  }

  /**
   * Compute P_m^m(x) from analytic result:
   * P_m^m(x) = (-1)^m (2m-1)!! (1-x^2)^(m/2), m>0
   * P_0^0(x) = 1
   */
  double legendre_Pmm(int m, double x);

  inline double legendre_P0(double /*x*/)
  {
    return 1.0;
  }

  inline double legendre_P1(double x )
  {
    return x;
  }

  inline double legendre_P2(double x )
  {
    return 0.5*(3.0*x*x-1.0);
  }

  inline double legendre_P3( double x )
  {
    return 0.5*x*(5.0*x*x-3.0);
  }

  inline double legendre_Pl(const int l, const double x )
  {
    if(l<0||x<-1.0|| x>1.0)
    {
      THROW_DEFAULT_EXCEPTION( FDomainError );
    }
    else if( l == 0)
    {
      return legendre_P0(x);
    }
    else if( l == 1 )
    {
      return legendre_P1(x);
    }
    else if( l == 2 )
    {
      return legendre_P2(x);
    }
    else if( x == 1.0 )
    {
      return 1.0;
    }
    else if( x ==-1.0 )
    {
      return (IS_ODD(l) ? -1.0: 1.0);
    }
    else if( l < 100000) // magic highest number to compute
    {
      // recurrence l P_l = (2l-1) z P_{l-1}-(l-1) P_{l-2}
      double pellm2 = legendre_P0(x);
      double pellm1 = legendre_P1(x);
      double pell = pellm1;
      int ell;

      for(ell=2; ell <= l; ++ell)
      {
        pell = (x*(2*ell-1)*pellm1-(ell-1)*pellm2)/ell;
        pellm2=pellm1;
        pellm1=pell;
      }
      return pell;
    }
    else
    {
      THROW_DEFAULT_EXCEPTION(FNotImplementedException);
#if 0
      // Asymptotic expansion.
      // by Olver, p. 473
      double u = l+0.5;
      double th = acos(x);
      double J0;
      double Jm1;
      J0 = bessel_J0(u*th);
      Jm1= bessel_Jn(-1, u*th );

      double pre;
      double B00;
      double c1;

      /* B00 = 1/8 (1-th cot(th) / th^2
       * pre = sqrt( th/sin(th) )
       */
      if(th < ROOT4_DBL_EPSILON )
      {
        B00 = (1.0+th*th/15.0)/24.0;
        pre = 1.0+th*th/12.0;
      }
      else
      {
        double sin_th = sqrt(1.0 - x*x);
        double cot_th = x/sin_th;
        B00 = 1.0/8.0* (1.0-th*cot_th)/(th*th);
        pre = sqrt(th/sin_th);
      }

      c1 = th/u*B00;
      return pre*(J0+ c1*Jm1);
#endif
    }
  }

  namespace poch
  {
    double lnpoch_pos( const double a, const double x );
    double lnpoch(const double a, const double x);
  }

  double legendre_Plm(const int l, const int m, const double x );



  /** size of array returned by legendre array functions */
  inline int legendreArraySize(const int lmax, const int m)
  {
    return lmax-m+1;
  }

  /**
   * calculate P_m^m(x) from analytic result:
   * P_m^m(x) = (-1)^m(2m-1)!! (1-x^2)^(m/2), m > 0
   *          = 1, m=0
   */
  //static double legendrePmm(const int m, const double x);

  /**
   * legendre polynomial as an array of all l up to lmax
   * \param result has already to be scaled to required size.
   * array must have at least size resturned by legendreArraySize(...) as it is not resized
   */
  void legendrePlArray(std::vector<double> &result, const int lmax, const double x);

  /**
   * derivatives of legendre polynomial.
   * \param result legendre Polynomial values as defined by legendrePlArray
   * array must have at least size resturned by legendreArraySize(...) as it is not resized
   * \param result_deriv derivatives of legendre polynomials
   */
  void legendrePlDerivArray(std::vector<double>& result, std::vector<double>& result_deriv, const int lmax, const double x);

  /**
   * Associated legendre polynomials  \f[P_l^m\f]
   */
  void legendrePlmArray(std::vector<double> &results, const int lmax, const int m, const double x);

  /**
   * Derivatives of associated legendre polynomials \f[P_l^m\f]
   */
  void legendrePlmDerivArray(std::vector<double> &results, std::vector<double> &deriv, const int lmax, const int m, const double x);


  /**
   * Compute a simple associated legendre polynomial already normalized for spherical harmonics computation 
   */
  double legendreSHPlm(const int l, int m, const double x);
    
  /**
   * Special kind of associated legendre polynomials already normalized for spherical harmonics
   * \f[\frac{\sqrt{2l+1}}{4\pi}\sqrt{\frac{(l-m)!}{(l+m)!}}P_l^m(x)\f]
   * for m>=0 and stores it in an array containing P_m^m ... P_lmax^m
   */
  void legendreSHPlmArray(std::vector<double> &results, const int lmax, const int m, const double x);

  /**
   * Derivatives of normalized legendre polynomials
   * for m>=0 and stores it in an array containing P_m^m ... P_lmax^m
   */
  void legendreSHPlmDerivArray(std::vector<double> &results, std::vector<double> &deriv, const int lmax, const int m, const double x);

  // my helper, should be removed one day
  double legendre(const unsigned int l, const unsigned int m, const double x);
  /*
   * old code, should be removed
   * Implementation taken from "Spherical Harmonics" by David Eberly
   * Magic Software, In.
   * http://www.magic-software.com
   *
   * This implementation is highly unstable near x=+/-1.0 for m=+/-1
   */
  // static
  inline double legendreDerivate( unsigned int l, unsigned int m, double x)
  {
    //eassert(l>=0);
    //eassert(m>=0);

    if(l==0) // that's easy
      return 0;

    double denom = 1 - x*x;
    double denomsqrt = sqrt(denom);

    if(m==0)
    {
      return FMath::legendre(l,1,x)/denomsqrt;
    }

    // hard core implementation, not very fast, but should work,
    // we could buffer a lot of information when evaluation complete SH here.
    double f=1;
    //  if( m %2 == 0) f = 1;
    double s1 = -f*m *x*FMath::legendre(l,m,x);
    double s2 = f*(l+m)*(l-m+1)*FMath::legendre(l,m-1,x);
    return s1/denom -s2/denomsqrt;
  }

  inline double func_log_1plusx(double x)
  {
    // FIXME: there are better implementations
    return log(1.0+x);
  }
}
#endif

