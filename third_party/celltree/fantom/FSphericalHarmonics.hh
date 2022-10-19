#ifndef FSphericalHarmonics_HH
#define FSphericalHarmonics_HH

#include "stdAliases.hh"
#include "FArray.hh"
#include "FTensor.hh"

#include <complex>
#include <iostream>
#include "eassert.hh"
#include "FMatrix.hh"

#define NEW_SH

/**
 * \ingroup math
 *
 * Spherical harmonics
 * encapsulates both functions to calculate (static things below)
 * and objects to manage parameters calculated.
 *
 * \theta denotes the colatitudal coordinate \theta \in [0, \pi]
 * \phi denotes the longitudal coordinate \phi \in [0,2\pi)
 *
 * ATTENTION: This is heavy work in progress. Legendre functions
 * and evaluations should work, least-squares-matching gives
 * some funny results I have to evaluate first.
 *
 * These spherical harmonics are defined as
 * \f[ Y_l^m(\theta,\Phi)=X_l^m(\theta)e^{\Im m\Phi} \f]
 * where
 * \f[ X_l^m(\theta)=
 *      \left[ \frac{2l+1}{4\pi}\frac{(l-m)!}{l+m)!}\right]^{1/2}
 *      P_l^m(cos(\theta)) 
 * \f]
 * and
 * \f[
 *    Y_l^{-m}(\theta, \phi) = (-1)^m {Y_l^{m}}^*(\theta, \phi)
 * \f]
 * The Condon-Shortley Phase \f( (-1)^m\f) is already included in the Legendre polynomials.
 *
 * This definition gives \f$Y_l^m\f$s which are fully normalized in the sense
 * \f[
 *  \int_S Y_{l'}^{m'}* Y_l^m dS = \delta_{mm'}\delta_{ll'}
 * \f]
 * where \f$ dS = sin(\theta d\theta d\Phi)\f$.
 *
 * As for pure real spherical harmonics, complex parts are not needed,
 * since \f$X_l^{-m}=(-1)^mX_l^m \f$
 * coefficients of negative m relate to the complex part when evaluating
 * the function therealX(...)
 */
class FSphericalHarmonics 
{
  public: // functions to manage Spherical harmonics object

    /**
     * constructs an empty spherical harmonic object
     */
    FSphericalHarmonics();

    /**
     * Constructor that uses array data as coefficients of Spherical harmonic.
     * Storage is allocated, that the SH has at least data.size() coefficients.
     */
    FSphericalHarmonics(const FArray& data) : comp(data)
    {
      int l; int m;
      revindex(l,m,data.size());
      //assert(m == 0);
      nbBands=l+1;
      //            std::cout <<" constructing from array: " <<  comp << " bands = " << l << " " << comp.size() << data.size()<< std::endl;
      shrot = 0;
      std::cout << "SH: " << *this << std::endl;
    }

    /**
     * Constructor that creates SH representing the tensor specified.
     *
     * \warning Only works for tensors of order 0 and 2 atm.
     */
    FSphericalHarmonics(const FTensor& t);

    /**
     * Resize the SH to contain nbBands bands
     */
    void resizeBands(unsigned int nbBands, bool clear=true);

    /**
     * Sets (*this) to represent the tensor specified.
     * 
     * \warning Only works for tensors of order 0 and 2 atm.
     * \post
     * Spherical harmonics components of l=0 and l=2 represent
     * the second order tensor t
     * \param
     * t Tensor to match
     * \param resize
     * if true, the SH will be resized, otherwise, current size is kept
     * \param clear
     * values that are not set are kept
     */
    void setFromTensor(const FTensor& t, bool resize = true, bool clear=true);

    /**
     * \pre
     * SH has rank 0 or 2 (== 1 or 3 bands)
     * and is symmetric.
     * \returns
     * FTensor representing the symmetric part of the spherical harmonic
     */
    FTensor toTensor() const;

    /**
     * \returns
     * FArray holding the components of this spherical harmonic
     */
    const FArray& toArray() const;

    /**
     * Copy constructor
     */
    FSphericalHarmonics(const FSphericalHarmonics& rhs)
      : nbBands(rhs.nbBands), comp( rhs.comp ), shrot(new FSphericalHarmonics(*shrot)) 
      {}

    /**
     * Conversion of coefficients to FArray
     */
    operator FArray& ()
    {
      return comp;
    }

    operator const FArray& () const
    {
      return comp;
    }

    /**
     * Destructor
     */
    virtual ~FSphericalHarmonics();

    /**
     * constructs an spherical harmonic object with 
     * bands bands/degrees/rank/mode (however it is called) 
     */
    FSphericalHarmonics(unsigned int bands);

    /**
     * calculate value of real spherical harmonic at
     * position defined by theta and phi
     */
    virtual double evaluate(double theta, double phi) const;

    /**
     * Evaluate a complex spherical harmonic with coefficients cdata at position (theta, phi)
     * \param cdata
     * complex valued coefficients
     * \param theta
     * \param phi
     */
    static std::complex< double > evaluateComplex( const std::vector< std::complex< double > > &cdata, double theta, double phi);

    /**
     * access constant (l,m) of this spherical harmonic
     */
    const double & operator()(unsigned int l, int m) const;

    /**
     * access (l,m) of this spherical harmonic
     */
    double& operator()(unsigned int l, int m);

    /**
     * Access coefficients in linear array.
     */
    const double& operator()(unsigned int i) const
    {
      return comp[i];
    }

    /**
     * Access coefficients in linear array.
     */
    double operator[](unsigned int i) const
    {
      return comp[i];
    }

    /**
     * Access coefficients in linear form (non-const)
     */
    double& operator[](unsigned int i)
    {
      return comp[i];
    }

    /**
     * assign a spherical harmonic that matches best the 
     * values defined at the directions vecs and radius vals
     * matching is done by least square fitting of the given
     * points
     */
    void precomputeMatch(FMatrix &BI, const std::vector<FVector> &vecs) const;
    static void precomputeMatch(FMatrix &BI, const std::vector<FVector> &vecs, int order);
    virtual void match(const std::vector<FVector> &vecs, const FArray& vals);
    void matchSymmetric( const std::vector<FVector> &vecs, const FArray &b);

    /** 
     * BI is the same as in @see precomputeMatch or @see precomputeSymmtericMatch 
     */
    virtual void match(const FMatrix &BI, const FArray& vals);

    /**
     * Same as @see precomputeMatch but forces BI to compute a symmetric spherical harmonic
     * i.e. all coefficients of odd orders are zero.
     */
    void precomputeSymmetricMatch(FMatrix &BI, const std::vector<FVector> &vecs) const;

    /**
     * @see precomputeMatch
     * @see precomputeSymmetricMatch
     */
    static void precomputeSymmetricMatch(FMatrix &BI, const std::vector<FVector> &vecs, int order);

  public: 

    // static functions for calculation with Spherical harmonics
    // use these as a basic library.

    /*!
     *! Associated Legendre symbol
     * \f$P_l^m(x)\f$
     * \param x in [-1,1]
     * \param l
     * \param m mode where 0 <=m <= l
     *
     * \f[
     *   P_l^m(x) = (-1)^m(1-x^2)^{m/2}\frac{d^m}{dx^m}P_l(x)
     * \f]
     *
     * Reulting in the terms:
     * \f{eqnarray}
     *   P_0^0(x) = 1\\
     *   P_1^1(x) = -(1-x^2)^\frac{1}{2} \\
     *   P_1^0(x) = x\\
     *   P_2^2(x) = 3(1-x^2)\\
     *   P_2^1(x) = -3(1-x^2)^\frac{1}{2}x\\
     *   P_2^0(x) = \frac{1}{2}(3x^2-1)
     * \f}
     */
  //static double legendre( unsigned int l, unsigned int m, double x);

  /*
   * Weights for spherical harmonics to make it an orthonormal basis:
   * \f[\sqrt{ \frac{ 2l+1 }{ 2\pi } \frac{ (l-m)! }{ (l+m)! } }\f]
   */
  static double funcK(unsigned int l, int m);


  /**
   * Calculates spherical harmonics of mode (l,m).
   * \f$\theta \in [0..\pi]\f$
   * \f$\phi \in [0..2\pi]\f$
   */
  static std::complex<double> complexY(unsigned int l, int m, double theta, double phi);

  //static double sphericalHarmonics(unsigned int l, int m, double theta, double phi);

  //! compute \f$ \langle Y(\theta,\phi), Y(\theta,\phi)\rangle \f$
  static double squaredY(unsigned int l, int m, double theta, double phi);
#if 0 
  //! \f$ \Re\left[ Y_l^m(\theta,\phi)\right] \f$
  static double realY(unsigned int l, int m, double theta, double phi);

  //! \f$ \Im\left[ Y_l^m(\theta, \phi)\right] \f$
  static double imagY(unsigned int l, int m, double theta, double phi);
#endif
  /** SH as used i.e. in SH-lighting
   * \f$ Y_l^m(\theta,\phi) + {Y_l^m(\theta,\phi)}^{*} \f$
   * thus interpreting the coefficients the following way:
   * Let \f$[ a^{r}_lm\f$ be our coefficients of a real Spherical Harmonic and \f$a^{c}_lm\f$
   * be the complex coefficients of a complex spherical harmonic function.
   * \f{eqnarray}
   *    a^{r}_{l,m} = \Re \left[ a^{c}_{l,m} \right] + \Re \left[ a^{c}_{l,-m} \right] if m > 0 \\
   *    a^{r}_{l,m} = \Im \left[ a^{c}_{l,m} \right] - \Im \left[ a^{c}_{l,-m} \right] if m < 0 \\
   *    a^{r}_{l,m} = \Re \left[ a^{c}_{l,m} \right] if m=0
   * \f}
   */
  static double therealY(unsigned int l, int m, double theta, double phi);

  /*
   * Same as above, FVector has to be normalized
   * Using cartesian coordinates v
   * \f[
   * e^{i\phi} = (x+iy)/\sqrt{x^2+y^2}
   * \theta = sin^{-1}(\sqrt{\frac{ x^2+y^2 }{ x^2+y^2+z^2 }})
   *        = cos^{-1}(\frac{ z }{ \sqrt{x^2+y^2+z^2} }
   * \f]
   */
  //        static double sphericalHarmonics(unsigned int l, int m, const FVector &v);

  /*
   * Calculates \f$| Y_l^m(\theta, \phi)|^2\f$
   */
  //        static double sphericalHarmonicsEnergy(unsigned int l, int m, double theta, double phi);

  /**
   * \returns Number of components in Spherical Harmonic
   */
  inline size_t size()const{return comp.size();}
  inline size_t nbComponents() const { return size(); }

  /**
   * Sets all components x with fabs(x) < th to zero
   * \param th threshold
   * \post
   * Spherical harmonic with coefficients smaller than th set to zero
   */
void threshold(double th)
{
  //             std::cerr << "WARNING: " << " thresholding Spherical harmonic to " << th << "." << std::endl;
  for(unsigned int i=0; i< comp.size(); ++i)
    if(fabs(comp[i]) < th) comp[i]=0.0;
}

void clear()
{
  for(unsigned int i=0; i< comp.size(); ++i)
  {
    comp[i] = 0.;
  }
}

friend std::ostream& operator <<(std::ostream& o,const FSphericalHarmonics& sh);
public:
/**
 * helper function for calculation of normalization factors funcK
 * that computes the parts of the factorial fraction that is needed
 * \f[(l-m)!/(l+m)! = \frac{1}{(l-m+1)\cdot \ldots \cdot (l+m)}\f]
 */
static double calcdenomfactorial( unsigned int from, unsigned int to);

/**
 * convert l,m to the index in FArray
 * You want to use this if you extend functionality!
 * You won't calculate your own indices, that will go awfully wrong!
 */
static unsigned int index(unsigned int l, int m);
static void revindex(int &l, int &m, int size);
#ifndef NEW_SH
#if 0
static void derivativeY(FVector& deriv, int l, int m, double phi, double theta);
static double derivativeY(int l, int m, double theta, double phi, bool dir_is_theta);
#endif
#endif
//! Derivatives of the associated legendre functions \ref legendre.
/**
 * \warning Derivatives currently have low precision near the z-Axis, i.e. theta=0 and theta=2 M_PI
 *
 * \f[
 *    \frac{dP_n^0}{dx} = \frac{dP_n}{dx} = \frac{-n(n+1) P_n^{-1}}{\sqrt{1-x^2}} = \frac{ P_n^1}{\sqrt{1-x^2}}
 * \f]
 * \f[
 *    \frac{P_n^m}{dx} = \frac{mxP_n^m -(n+m)(n-m+1)\sqrt{1-x^2}P_n^{m-1}}{1-x^2}
 * \f]
 *
 * And hould result in
 * \f{eqnarray}
 *   DP_0^0 = 0 \\
 *   DP_1^0 = 1 \\
 *   DP_1^1 = \frac{-x}{\sqrt{1-x^2}} \\
 *   DP_2^0 = 3x \\
 *   DP_2^1 = 3\frac{1-2x^2}{\sqrt{1-x^2}} \\
 *   DP_2^2 = -6x \\
 *   DP_3^0 = \frac{3}{2}(5x^2-1) \\
 *   DP_3^1 = \frac{3}{2}\frac{x(11-15x^2)}{\sqrt{1-x^2}} \\
 *   DP_3^2 = 15(1-3x^2) \\
 *   DP_3^3 = -45\sqrt{1-x^2}
 * \f}
 */ 
//static double legendreDerivate( unsigned int l, unsigned int m, double x);

void derivative(FVector& deriv, double theta, double phi) const;

#ifndef NEW_SH
double derivate( double theta, double phi, bool use_theta ) const;
#endif

//! Get the coefficients of a complex SH matching this SH
void getComplexCoefficients( std::vector< std::complex< double > >& data ) const;

//! set this SH to the real part of the specified complex data
//! this produces a SH that matches the real part exactly but ignores
//! the complex part!!!
void setFromComplexCoefficients( std::vector< std::complex< double > > cdata );

//! set this SH to the complex part of the specified complex data
//! \warning: untested
void setComplexFromComplexCoefficients( std::vector< std::complex< double > > cdata);

public:
/**
 * Number of bands in SH.
 *
 * 1 = scalar only, completely invariant
 * 2 = asymmetric part, 3 components (+ scalar one above)
 * 3 = symmetric part, 5 components (+ both above) == 2nd Order Tensor
 */
inline unsigned int bands() const
{
  return nbBands;
}

/**
 * Returns the rank of Spherical harmonics i.e. bands()-1
 */
inline unsigned int rank() const
{
  return nbBands -1;
}

inline unsigned int getOrder() const { return rank(); }

/**
 * compute the eigensystem of this spherical harmonic interpreted
 * as symmetric tensor of rank 2.
 * Higher orders and antisymmetric part are ignored!
 * \pre 
 * Spherical harmonics is valid and of rank >= 2
 * \post
 * vals store three real eigenvalues
 * v stores the eigenvectors
 */
void getEigenSystem(FVector& vals, FVector v[3]) const;


static void test();
#if 0
void getSinCoefficients(std::vector<double>& s)
{
  for(int l=0; l< bands(); ++l)
    for(int m=-1; m>= -l; --m)
    {
      s.push_back( (*this)(l,m) );
    }
}

void getCosCoefficients( std::vector<double>& c)
{
  for(int l=0; l< bands(); ++l)
    for(int m=0; m<=l; ++m)
    {
      s.push_back( (*this)(l,m) );
    }
} 
void rotateZ()
{
  std::vector<double> cp; getCosCoefficients(cp);
  std::vector<double> sp; getSinCoefficients(sp);

  int L = cp.size();

  std::vector<double> cangl(L);
  std::vector<double> sangl(L);

  for(int i=0; i<L; ++i)
    cangl[i] = cos(i*angl);

  for(int i=0; i<L; ++i)
    sangl[i] = sin(i*angl);

  for(int l=0; l<=L; ++l)
  {
    Sl=shsin(S,l);
    ...
  }
#endif

  /**
   * returns a matrix r that rotates the Spherical harmonic around the z axis by the angle theta
   */
  static void getZRotationMatrix( FMatrix& r, double theta, int bands );


  protected: // member variables
  unsigned int nbBands; //< number of bands in spherical harmonic (i.e. L+1, bands, rank+1...)
  //        std::vector< std::complex< double > > comp;
  FArray comp; //< components of real valued SH 

  /**
   * The implementation of derivatives has trouble when being close to the
   * z-axis, thus we have a copy of this spherical harmonics dynamically 
   * created when needed. This points to zero if nothing is calculated,
   * otherwise to a rotated version of *this.
   *
   * WARNING: Currently not all functions invalidate this if the spherical 
   * harmonic changes.
   */
  mutable FSphericalHarmonics* shrot;
};
/*
   class FRealSphericalHarmonics : public FSphericalHarmonics
   {
   public:
   FRealSphericalHarmonics(unsigned int bands);
   virtual unsigned int index(unsigned int l, int m) const;
   virtual double evaluate(double theta, double phi) const;
   virtual bool match(const std::vector<FVector> &vecs, const std::vector<double> vals);
   };

   class FSymmetricSphericalHarmonics : public FSphericalHarmonics
   {
   public:
   FSymmetricSphericalHarmonics(unsigned int bands);
   virtual unsigned int index(unsigned int l, int m) const;
   virtual double evaluate(double theta, double phi) const;
   virtual bool match(const std::vector<FVector> &vecs, const std::vector<double> vals);
   };
   */



namespace FMath
{

  /**
   * Conversion from spherical coordinates to 3D Vector
   *
   * computes \f[ \left(\begin{array}{c}\sin(\theta)\cos(\phi)\\ 
   * \sin(\theta)\sin(\phi)\\
   * \cos(\theta)\end{array}\right)r
   * \f]
   */
  extern void toVector(FVector &v, double theta, double phi, double r);

  /**
   * Inverse computation to \ref toVector
   */
  extern void toThetaPhi(double &theta, double &phi, double &r, const FVector& v);

  /**
   * Tangent along constant \f[\phi\f] value on the sphere
   */
  extern void toTangentPhi(FVector &v, double theta, double phi);

  /**
   * Tangent vector along constant \f[\theta\f] direction on the sphere
   */
  extern void toTangentTheta(FVector &v, double theta, double phi);

  /** 
   * Create rotation matrices for real valued spherical harmonics for band l
   */
  void getMatrices(std::vector<FMatrix> &mats, const FMatrix& ROT, const int l);

  /**
   * Get a rotation matrix to rotate real valued spherical harmonics
   * up to order l ( including l )
   */
  void getRotationMatrix( FMatrix& mat, const FMatrix& rotation, const int l);


  /**
   * Computes a rotation matrix described by 
   * \f[ R_sh(\alpha, \beta, \gamma) = Z_{\gamma} X_{-90} Z_{\beta} X_{+90} Z_{\alpha} \f]
   * 
   * Warning, this is SH coordinate system, which has Z up and vectors like
   * (x, z, y)
   * alpha rotates around z ( 1,0,0 ) by alpha=M_PI_2 via (0,0, -1 ) etc...
   * beta around x, thus ( 0,1,0 ) by beta = M_PI_2 becomes (0,0, 1) etc...
   * gamma again around z
   */
  void getROTMatrix( FMatrix &rot, double alpha, double beta, double gamma );

} // namespace FMath

#endif
