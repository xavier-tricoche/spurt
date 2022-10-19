#ifndef FRandom_hh
#define FRandom_hh

// Helpers for programs using random variables
//
// All implementations are not thread save!

#include "FFixArray.hh"
#include "FSpecialFunctions.hh"
#include "eassert.hh"
#include <iostream>

/**
 * Generate a single value
 */
class FDiracDistribution
{
  public:
    FDiracDistribution( double mean = 0. ) : mean(mean){}

    double random()
    {
      return mean;
    }

    double variance()
    {
      return 0.;
    }

    double pdf(const double x) const
    {
      return x==mean ? 1.: 0.;
    }

    double cdf(const double x) const
    {
      return x<mean? 0. : 1.;
    }

  private:
    double mean;
};

class FNormalRandomTraits
{
  public:
    FNormalRandomTraits(){}

    double pdf(double x) const
    {
      static const double FACT = 1./sqrt(2.*M_PI);
      return FACT*exp(-x*x/2);
    }

    double cdf(double x) const
    {
      return 0.5+0.5*erf(x*M_SQRT1_2 );
    }

    double mean() const
    {
      return 0;
    }

    double variance() const
    {
      return 1.;
    }

    double skewness() const
    {
      return 0.;
    }

    double kurtosis() const
    {
      return 0.;
    }

    double moments(int i) const
    {
      if(i==0) return mean();
      if(i==1) return variance();
      if( i%2 ==0)
        //return fact(i)/(pow(2, i/2) *fact(i/2))
        eassert(false); // TODO: implement fact
      else
        return 0.;
    }

    double entropy() const
    {
      return log(sqrt(2*M_PI*M_E));
    }

    // get distribution parameters
    double getSigma() const { return 1.0; }

};

/**
 * Generate gaussian random variables.
 * Random has to return values in [0,1) (0 is valid, 1 is not allowed!)
 * Valus have to be calculated by Random::random()
 *
 * it evaluates:
 * \f[ \eta_1=\sqrt{-2 \ln(1-\chi_2)}cos(2 \pi \chi_1) \f]
 * \f[ \eta_2=\sqrt{-2 \ln(1-\chi_2)}sin(2 \pi \chi_1) \f]
 */
template < class UniformRandom >
class FNormalRandom : public FNormalRandomTraits
{
    public:
    FNormalRandom()
    {
        pos = 1; //< next time, calculate a new value and return eta[0]
    }

    double random()
    {
        (++pos) %= 2;
        if(pos == 0) init_rand();
        return eta[pos];
    }
   protected:
    // calculate two independant gaussian random variables
    inline void init_rand(void)
    {
        double tau = rand.random()*2.*M_PI;
        double rho = sqrt(-2.*log(1.-rand.random()));

        eta[0] = rho *cos(tau);
        eta[1] = rho *sin(tau);
    }
        // current value to return
    int pos;
    double eta[2];

    private:
    UniformRandom rand;
};

class FNormalRandom2Traits
{
  public:
    typedef FNormalRandom2Traits self_type;
    FNormalRandom2Traits( double mu, double sigma )
    : mu(mu), sigma(sigma) {}

    // shift the distribution by shift
    self_type& shift( double shift )
    {
      mu += shift;
      return *this;
    }

    self_type& operator*=( double scale )
    {
      mu *= scale;
      sigma *= fabs( scale );
      return *this;
    }

    self_type& operator+=( const self_type& rhs )
    {
      mu += rhs.mu;
      sigma = sqrt( sigma*sigma + rhs.sigma*rhs.sigma );
      return *this;
    }

    void setParameters( const double mu, const double sigma)
    {
      this->mu = mu;
      this->sigma = sigma;
    }

    double pdf(double x)
    {
      static const double FACT = 0.25*M_2_SQRTPI;
      return FACT/sigma*exp(-(x-mu)*(x-mu)/(2.*sigma*sigma));
    }

    double cdf(double x) const
    {
      return 0.5+ 0.5*erf(( x-mu )/sigma*M_SQRT1_2);
    }

    double mean() const
    {
      return mu;
    }

    double variance() const
    {
      return sigma*sigma;
    }

    double entropy() const
    {
      return log(sigma*sqrt(2*M_PI*M_E));
    }


    double getSigma() const { return sigma; }
    // parameters
    double mu;
    double sigma;
};

template < class Random >
class FNormalRandom2 : public FNormalRandom2Traits
{
  public:
    typedef FNormalRandom2<Random> self_type;
    typedef FNormalRandom2Traits traits_type;

    FNormalRandom2( double mu, double sigma )
    : FNormalRandom2Traits( mu, sigma ), rand(){}

    FNormalRandom2( const traits_type &bt )
    : traits_type( bt ), rand() {}

    double random()
    {
      return rand.random()*sigma + mu;
    }

    FNormalRandom<Random> rand;
};


class FExponentialDistributionTraits
{
  public:
        FExponentialDistributionTraits( double lambda = 1. )
        : lambda(lambda)
        {}

        void setParameters( const double lambda )
        {
          this->lambda = lambda;
        }

        double mean() const
        {
          return 1./lambda;
        }

        double median() const
        {
          return log(2.)/lambda;
        }

        double variance() const
        {
          return sqrt(1./lambda);
        }

        double skewness() const
        {
          return 2.;
        }

        double kurtosis() const
        {
          return 6.;
        }

        double entropy() const
        {
          return 1.-log(lambda);
        }
        
    protected:
        double lambda;
 
};

/**
 * create a exponential distribution
 * \f[ f(x)= \lambda \exp{\left(-\lambda x\right) } \f]
 * by evaluating:
 * \f[ \chi = -1/\lambda \ln (1-\lambda) \f]
 *
 * Random has to return uniformly distributed values in [0,1) (0 is valid, 1 is not allowed!)
 */
template < class UniformRandom >
class FExponentialDistribution : public FExponentialDistributionTraits
{
    public:
        FExponentialDistribution( double lambda = 1. )
        : FExponentialDistributionTraits( lambda )
        {}

        FExponentialDistribution( const FExponentialDistributionTraits& t )
        : FExponentialDistributionTraits( t ){}

        double random()
        {
            return -1./lambda * log(1.-rand.random());
        }
   private:
        UniformRandom rand;
};

#if 0
/**
 * Computes random variables on the unit sphere i.e.
 * \f[ \eta =(x1,...,xn) \in S_{n-1} \subset R^n\f]
 * with
 * \f[ x_1^2+\ldots+x_n^2=1\f]
 * by evaluating
 * \f[ \eta =\frac{(\gamma_1,\ldots,\gamma_n)}{\sqrt{\gamma_1^2+\ldots+\gamma_n^2}}\f]
 */
template < class NormalRandom , int N = 3 >
class FRandomUnitSphere : private NormalRandom
{
    public:
    FRandomUnitSphere()
    {
    }

    FFixArray<double, N>& random()
    {
        for(int i=0; i<values; ++i)
        {
            values[i] = NormalRandom::random();
        }
        values.normalize();
        return values;
    }
    
    private:
    FFixArray<double, N> values;
};
#endif

template <class UniformRandom>
class FLaplaceRandom
{
  public:
    FLaplaceRandom( double mu=0.0, double b=1.0 )
      : mu(mu), b(b)
      {}

    void setParameters( const double mu, const double b )
    {
      this->mu = mu;
      this->b  = b;
    }
    
    double random()
    {
      double u=randomg.random();
      double l;
      if(u<0.5)
        l=std::log(2.0*u);
      else
        l=-std::log(2.0*(1.-u));
      l*=sqrt(b/2.0);
      l+=mu;
      return l;
    }

    double pdf(double x) const
    {
      return 1./(2.*b)*exp(-fabs(x-mu)/b);
    }

    double cdf(double x) const
    {
      if(x<mu)
        return 0.5*exp(-(mu-x)/b);
      else
        return 1.-0.5*exp(-(x-mu)/b);
    }

    double mean() const
    {
      return mu;
    }

    double median() const
    {
      return mu;
    }

    double variance() const
    {
      return 2*b*b;
    }

    double entropy() const
    {
      return 1.+log(2*b);
    }

  private:
    double mu, b;
    UniformRandom randomg;
};

template <class Random>
class FSumOfNRandoms
{
  public:
    FSumOfNRandoms( unsigned int additions, const Random& random )
      : nb(additions), randomg(random)
    {}

    double random()
    {
      double val=0.0;
      for(unsigned int i=0; i< nb; ++i)
        val += randomg.random();
      return val/nb;
    }
  protected:
    Random randomg;
    unsigned int nb;
};

template <class NormalRandom>
class FRayleighRandom
{
  /**
   * \todo FIXME: untested
   */
  public:
    FRayleighRandom( double sigma )
    : rand( 0.0, sigma )
    {}

    void setParameters( const double sigma )
    {
      rand.setParameters( 0.0, sigma );
    }
    
    double random()
    {
      double r1= rand.random();
      double r2= rand.random();

      return sqrt(r1*r1+r2*r2);
    }

    double pdf(const double x) const
    {
      return x*exp(-x*x/2./rand.variance())/rand.variance();
    }

    double cdf(const double x) const
    {
      return 1.-exp(-x*x/2./rand.variance());
    }

    double mean() const
    {
      return rand.getSigma()*sqrt(M_PI_2);
    }

    double median() const
    {
      return rand.getSigma()*sqrt(log(4.));
    }

    double moment2() const
    {
      return 2.*rand.getSigma()*rand.getSigma();
    }

    double moment3() const
    {
      double s = rand.getSigma();
      s=s*s*s;
      return 3.*s*sqrt(M_PI_2);
    }

    double moment4() const
    {
      double s = rand.getSigma();
      s*=s;
      s*=s;
      return 8.*s;
    }
    
    // mathworld:
    // http://mathworld.wolfram.com/RayleighDistribution.html
    double moment(int k) const
    {
      if(k==1) return mean();
      if(k==2) return moment2();
      if(k==3) return moment3();
      if(k==4) return moment4();
      
      const double sigma = rand.getSigma();
      return pow(2., k/2)* pow(sigma, k)* FMath::gamma::gamma(1+0.5*k); 
    }

    double variance() const
    {
      return (2.-M_PI_2)*rand.variance();
    }

  /**
   * \todo FIXME: on Wikipedia this has a small gamma in it, don't know
   * where it comes from.
   */
    //double entropy() const
    //{
    //  return 1.+log( 1./M_SQRT2/(rand.variance()*rand.getSigma()))+gamma/2.;
    //}

  private:
      NormalRandom rand;
};


class FRiceRandomTraits
{
  public:
  public:
    FRiceRandomTraits( double A, double sigma )
      : A( A ), sigma( sigma )
    {
      eassert(sigma >= 0);
      eassert(A     >= 0);
      //eassert(rand.mean() == A/M_SQRT2 );
#ifndef NODEBUG
//      std::cout << "RiceRandom: " << rand.mean() << std::endl;
#endif
    }

    void setParameters( const double A, const double sigma )
    {
      this->A = A;
      this->sigma = sigma;
    }

    /**
     * \f[ 
     *    P(x | A, \sigma) = \frac{A}{\sigma^2} \exp{\left( -\frac{A^2 + x^2}{2\sigma^2}\right )}I_0\left(Ax/\sigma^2\right) 
     * \f]
     */
    double pdf(const double x) const
    {
      //eassert(false);
      // this is currently not working
      using FMath::bessel::bessel_I0;
      if(x<=0) return 0.0;

      double V=A;
      double Z=x;
      //std::cout << Z << " " << V << std::endl; 

      double sigsqr = 1./sigma/sigma;
      double f1 = Z*sigsqr;
      double exp_param =-0.5*(Z*Z+V*V)*sigsqr;
      double bessel_param = Z*V*sigsqr;
      double f3;
      if( bessel_param > 3.75)
       //this may have problems for small sigma or large Z*V
        f3 = exp(exp_param+bessel_param)*large_bessel(bessel_param);
      else
        f3 = exp(exp_param)*bessel_I0(bessel_param);

        //std::cout << sigsqr << " " << f1 << " " << exp_param << " " << bessel_param << "     " << f3 << std::endl;
      return f1*f3;
    }

    double pdftest(const double x) const
    {
      //eassert(false);
      // this is currently not working
      // should be the same as above as pdf(x) as soon as
      // I figure out how to set the parameters e and sigma
      if(x<=0.0) return 0.0;
      const double e= 0.5;
      double s1 = e/sigma/sqrt(2.*M_PI);
      double s2 = sqrt(1-e*e)* erf( sqrt(1-e*e)/e*x/sigma)* x/sigma/sigma;
      double ex1 = exp(-0.5*( x/e/sigma)*(x/e/sigma));
      double ex2 = exp(-0.5*( x/  sigma)*(x/  sigma));
      return (s1*ex1+s2*ex2);
    }

    double cdftest(const double x ) const
    {
      eassert(false);
      // currently not working.
      // gives similar results like pdftest thus seems to be ok,
      // but I do not know how to set the parameters e and sigma
      double e=0.5;
      const double sq = sqrt(1.-e*e);
      double s1 = erf(x/e/sigma);
      double s2 = sq*erf(sq/e*x/sigma)*exp(-0.5*x*x/sigma/sigma);
      return s1-s2;
    }

    double variance() const
    {
      const double mean = this->mean();
      return this->moment2() - mean*mean;;
    }


    /**
     * \f[
     *    \mu_1 = \sigma\sqrt{\frac{\pi}{2}}L_{1/2}(\frac{A^2}{2\sigma^2})
     * \f]
     */
    double mean() const
    {
      // ok, working
      using FMath::laguerre::laguerre1_2;
      const double Asqr   = A*A;
      const double l= laguerre1_2(-0.5*Asqr/sigma/sigma);
      return sigma*sqrt(M_PI_2)*l;
    }

    // 1/N sum x*x when x is the random variable
    // (no mean subtracted here!)
    double moment2() const
    {
      // tested, working
      // source: wikipedia Rice_distribution
      // same as in moment(2) with: gamma(2) = 1! = 1
      const double sigsqr= sigma*sigma;
      const double Asqr  = A*A;
      /*
      using FMath::laguerre::laguerre1_2;
      const double l = laguerre1_2(-0.5*Asqr/sigsqr);
      return 2.*sigsqr+Asqr-M_PI_2*sigsqr*l*l;
      */
      return 2.*sigsqr+Asqr;
    }

    double moment4() const
    {
      const double sig2 = sigma*sigma;
      const double sig4 = sig2*sig2;
      const double Asqr  = A*A;
      return 8.*(sig4+sig2*Asqr)+Asqr*Asqr;
    }

    double moment(int k) const
    {
      if(k==1) return mean();
      if(k==2) return moment2();
      if(k==4) return moment4();
      
      const double Asqr=A*A;
      const double s= sqrt(sigma);
      return pow(s,k)*pow(2., k/2.)*FMath::gamma::gamma(1+k/2.)*
        FMath::laguerre::laguerre1_2(-0.5*Asqr/sigma/sigma);
      
    }

    // rician distribution parameters
    
    /**
     * Rice factor K of the faded envelope is a measure of the power 
     * in the faded envelope that has been producedd by the means of
     * of the two random variables.
     *
     * This can be interpreted as a Signal to Noise Ratio (SNR), too.
     * (e.g. as done in some papers on Magnetic Resonance Imaging (MRT)
     */
    double getK() const
    {
      const double Asqr   = A*A; // FIXME: do we need another factor? check!
      const double sigsqr = sigma*sigma; 
      return Asqr/sigsqr;
    }
  

  protected:
    // because of numerical instabilities, we need this here
    double large_bessel( const double x) const
    {
      double ax, ans;
      double y;

      ax=fabs(x);
      y=3.75/ax;
      ans=(1./sqrt(ax))*(0.3984228+y*(0.1328592e-1
            +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
                  +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
                        +y*0.392377e-2))))))));
      return ans;
    }

  protected:
    double A, sigma;
};

// brude force computation of a rice distribution
// random has to be gaussian distributed

/* Rician distribution arises in communications theory as the magnitude
 * of the fourier transform of a consine wave of amplitude A corrupted
 * by additive white noise.
 *
 * \f[\xi = x_r*x_r + x_c*x_c \f]
 * where x is a complex signal and its real and complex part is a random
 * process with mean A/M_SQRT2 and variance sigma.
 */
template <class NormalRandom>
class FRiceRandom : public FRiceRandomTraits
{
  public:
  /**
   * \todo FIXME: mu and sigma are not the right parameters to be 
   * used here, fix this to mach reasonable parameters of
   * the rician distribution
   *
   */
  FRiceRandom( double A, double sigma )
      : FRiceRandomTraits( A, sigma ), rand(A/M_SQRT2, sigma)
  {}

  FRiceRandom( const FRiceRandomTraits& traits )
    :FRiceRandomTraits( traits ), rand( traits.A/M_SQRT2, traits.sigma )
  {}

  double random()
  {
    double r1 = rand.random();
    double r2 = rand.random();

    return sqrt(r1*r1+r2*r2);
  }

  private:
    NormalRandom rand;
};

#if 0
template<class NormalDistribution>
class ChiSquareRandom
{
  // untested
  public:
  ChiSquareRandom(int k) : k(k){}
  ChiSquareRandom(int k, const NormalDistribution& n) : k(k), normal(n){}
  
  double random()
  {
    double sum=0;
    for(int i=0; i<k; ++i)
    {
      double r= normal.random();
      sum += r*r;
    }
    return sum;
  }

  double pdf(const double x) const
  {
    if(x<=0) return 0.0;
    return pow(k,n/2.-1.)*exp(-x/2.) /  (pow(2.,k/2.)*gamma(k/2));
  }

  double mean() const
  {
    return k;
  }

  double variance() const
  {
    return 2*k;
  }

  int k;
  NormalDistribution normal;
};
#endif


/** Bronstein Semedjajew Musil Muehlig, Taschenbuch der Mathematik
 * P. 780, Ch 16.2.4.3 Logarithmische Normalverteilung
 *
 * Produces a distribution for X where Y=log X and Y is a normal distribution
 */
template< class NormalDistribution >
class LogNormalRandom
{
  // untested
  public:
    LogNormalRandom( double muL, double sigmaL ) : normal(muL, sigmaL ), muL(muL), sigmaLsqr(sigmaL*sigmaL){}

    double random()
    {
      return exp( normal.random());
    }


    //! http://en.wikipedia.org/wiki/Log-normal_distribution
    double pdf(double t) const
    {
      double p= (log(t)-muL);
      p = p*p;
      p/= sigmaLsqr;
      p*= 2.;
      return exp(-p)/(t*sqrt(sigmaLsqr*M_2_PI));
    }

    //! http://en.wikipedia.org/wiki/Log-normal_distribution
    double cdf(double t) const
    {
      return 0.5+ 0.5*erf( (log(t)-muL)/sqrt( 2.*sigmaLsqr));
    }

    double median() const
    {
      return exp(muL);
    }

    double mode() const
    {
      return exp(muL-sigmaLsqr);
    }

    double skewness() const
    {
      return exp(-muL-sigmaLsqr/2)*(exp(sigmaLsqr + 2.))*sqrt( exp(sigmaLsqr - 1.));
    }

    double kurtosis() const
    {
      return exp( 4.*sigmaLsqr ) + 2.*exp( 3*sigmaLsqr ) + 3.*exp( 2.*sigmaLsqr) - 6.;
    }

    double entropy() const
    {
      return 0.5 + 0.5*log(M_2_PI*sigmaLsqr);
    }
    
    double mean()
    {
      return exp( muL+0.5*sigmaLsqr);
    }

    double variance()
    {
      return (exp(sigmaLsqr)-1.)*exp(2.*muL+sigmaLsqr);
    }

    void setParameters( const double muL, const double sigmaL )
    {
      normal.setParameters( muL, sigmaL );
      this->muL = muL;
      this->sigmaL = sigmaL;
    }

    // http://de.wikipedia.org/wiki/Lognormalverteilung
    void setParametersFromMeanAndVariance( const double mean, const double var )
    {
      double sigmaL = sqrt( log( 1.+ var/mean/mean));
      double muL    = log(mean)-var/2.;
      setParameters( muL, sigmaL );
    }

    NormalDistribution normal;
    double muL;
    double sigmaLsqr;
};

#endif
