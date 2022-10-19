#include "FSpecialFunctions.hh"

#include <math.h>

#include <iostream>
#include "FException.hh"

#include "FMatrix.hh"
#include "eassert.hh"

#include "FTensor.hh"

#include <limits>

#define VERBOSE
//using namespace std;

#define EPSILON 1e-9
#define SQRT_EPSILON 1e-8

#include "FMathValues.hh"

using namespace FMath;

namespace FMath
{

  namespace laguerre{
    
    double laguerre1_2( const double x )
    {
      using namespace bessel;
      double x_2 = 0.5*x;
      //WARNING: this may be split up as the other stuff where
      //an exp gets multiplied by bessel function to extract the
      //exponential in the approximation of the bessel function
      //to get better results, if needed
      return exp(x_2)*((1.-x)*bessel_I0(-x_2)-x*bessel_I1(-x_2));
    }
    
  }
  
  namespace erf{
    /*
     * Abramowitz and Stegun: Eq. (7.1.14) gives this continued fraction
     * for erfc(z)
     * erfc(z) = sqrt(pi).exp(-z^2).  1   1/2   1   3/2   2   5/2  
     *                               ---  ---  ---  ---  ---  --- ...
     *                               z +  z +  z +  z +  z +  z +
     *
     * This is evaluated using Lentz's method, as described in the
     * narative of Numerical Recipes in C.
     *
     * The continued fraction is true providing real(z) > 0. In practice
     * we like real(z) to be significantly greater than 0, say greater
     * than 0.5.
     */
    std::complex<double> cerfc_continued_fraction(const std::complex<double>& z)
    {
      const double tiny = std::numeric_limits<double>::min();
      const double eps = std::numeric_limits<double>::epsilon();
      // first calculate z+ 1/2   1 
      //                    ---  --- ...
      //                    z +  z + 
      std::complex<double> f(z);
      std::complex<double> C(f);
      std::complex<double> D(0.0);
      std::complex<double> delta;
      double a;

      a = 0.0;
      do {
        a += 0.5;
        D = z + a * D;
        C = z + a / C;
        if ((D.real() == 0.0) && (D.imag() == 0.0))
          D = tiny;
        D = 1.0 / D;
        delta = C * D;
        f = f * delta;
      } while (abs(1.0 - delta) > eps);

      // Do the first term of the continued fraction
      f = 1.0 / f;

      // and do the final scaling
      f = f * exp(-z * z) / sqrt(M_PI);

      return f;
    }

    std::complex<double> cerf_continued_fraction(const std::complex<double>& z)
    {
      if (z.real() > 0)
        return 1.0 - cerfc_continued_fraction(z);
      else
        return -1.0 + cerfc_continued_fraction(-z);
    }



    /*
     * Abramawitz and Stegun: Eq. (7.1.5) gives a series for erf(z) good
     * for all z, but converges faster for smallish abs(z), say abs(z) < 2.
     */
    std::complex<double> cerf_series(const std::complex<double>& z)
    {
      const double tiny = std::numeric_limits<double>::min();
      std::complex<double> sum(0.0);
      std::complex<double> term(z);
      std::complex<double> z2(z*z);

      for (int n = 0; (n < 3) || (abs(term) > abs(sum) * tiny); n++) {
        sum += term / static_cast<double>(2 * n + 1);
        term *= -z2 / static_cast<double>(n + 1);
      }

      return sum * M_2_SQRTPI; //2.0 / sqrt(M_PI);
    }

    /*
     * Numerical Recipes quotes a formula due to Rybicki for evaluating
     * Dawson's Integral:
     *
     * exp(-x^2) integral exp(t^2).dt = 1/sqrt(pi) lim  sum  exp(-(z-n.h)^2) / n
     *            0 to x                           h->0 n odd
     *
     * This can be adapted to erf(z).
     */
    std::complex<double> cerf_rybicki(const std::complex<double>& z)
    {
      double h = 0.2; // numerical experiment suggests this is small enough

      // choose an even n0, and then shift z->z-n0.h and n->n-h. 
      // n0 is chosen so that real((z-n0.h)^2) is as small as possible. 
      int n0 = 2 * static_cast<int>(z.imag() / (2 * h) + 0.5);

      std::complex<double> z0(0.0, n0 * h);
      std::complex<double> zp(z - z0);
      std::complex<double> sum(0.0, 0.0);

      // limits of sum chosen so that the end sums of the sum are
      // fairly small. In this case exp(-(35.h)^2)=5e-22 
      for (int np = -35; np <= 35; np += 2) {
        std::complex<double> t(zp.real(), zp.imag() - np * h);
        std::complex<double> b(exp(t * t) / static_cast<double>(np + n0));
        sum += b; 
      }

      sum *= M_2_PI * exp(-z*z); // 2.0 * exp(-z * z) / pi;

      return std::complex<double>(-sum.imag(), sum.real());
    }



    
    /*
     * This function calculates a well known error function erf(z) for
     * complex z. Three methods are implemented. Which one is used
     * depends on z. 
     */
    std::complex<double> erf(const std::complex<double>& z)
    {
      // Use the method appropriate to size of z - 
      // there probably ought to be an extra option for NaN z, or infinite z
      if (abs(z) < 2.0)
        return cerf_series(z);
      else {
        if (std::abs(z.real()) < 0.5)
          return cerf_rybicki(z);
        else
          return cerf_continued_fraction(z);
      }
    }

  }


  
  namespace gamma{
    double lngamma(const double x);
  }

  double exp_mult( const double x, const double y)
  {
    const double ay = fabs(y);

    if(y == 0.0)
      return 0.0;
    else if((x<0.5*F_LOG_DBL_MAX && x > 0.5*F_LOG_DBL_MIN)
        && (ay < 0.8*F_SQRT_DBL_MAX && ay > 1.2*F_SQRT_DBL_MIN)
        )
    {
      const double ex = exp(x);
      return y*ex;
    }
    else
    {
      const double ly = log(ay);
      const double lnr = x+ly;

      if(lnr > F_LOG_DBL_MAX -0.01)
      {
        THROW_EXCEPTION( FException, "Overflow" );
      }
      else if (lnr < F_LOG_DBL_MIN +0.01)
      {
        THROW_EXCEPTION( FException, "Underflow" );
      }
      else
      {
        const double sy = F_SIGN(y);
        const double M = floor(x);
        const double N = floor(ly);
        const double a = x -M;
        const double b = ly -N;
        return sy * exp(M+N)*exp(a+b);
      }
    }
  }

  //namespace detail{


  /*-*-*-*-*-*-*-*-*-*-*-* Private Section *-*-*-*-*-*-*-*-*-*-*-*/

#define FACT_TABLE_MAX  170
#define FACT_TABLE_SIZE (FACT_TABLE_MAX+1)
  static struct {int n; double f; long i; } fact_table[FACT_TABLE_SIZE] = {
    { 0,  1.0,     1L     },
    { 1,  1.0,     1L     },
    { 2,  2.0,     2L     },
    { 3,  6.0,     6L     },
    { 4,  24.0,    24L    },
    { 5,  120.0,   120L   },
    { 6,  720.0,   720L   },
    { 7,  5040.0,  5040L  },
    { 8,  40320.0, 40320L },

    { 9,  362880.0,     362880L    },
    { 10, 3628800.0,    3628800L   },
    { 11, 39916800.0,   39916800L  },
    { 12, 479001600.0,  479001600L },

    { 13, 6227020800.0,                               0 },
    { 14, 87178291200.0,                              0 },
    { 15, 1307674368000.0,                            0 },
    { 16, 20922789888000.0,                           0 },
    { 17, 355687428096000.0,                          0 },
    { 18, 6402373705728000.0,                         0 },
    { 19, 121645100408832000.0,                       0 },
    { 20, 2432902008176640000.0,                      0 },
    { 21, 51090942171709440000.0,                     0 },
    { 22, 1124000727777607680000.0,                   0 },
    { 23, 25852016738884976640000.0,                  0 },
    { 24, 620448401733239439360000.0,                 0 },
    { 25, 15511210043330985984000000.0,               0 },
    { 26, 403291461126605635584000000.0,              0 },
    { 27, 10888869450418352160768000000.0,            0 },
    { 28, 304888344611713860501504000000.0,           0 },
    { 29, 8841761993739701954543616000000.0,          0 },
    { 30, 265252859812191058636308480000000.0,        0 },
    { 31, 8222838654177922817725562880000000.0,       0 },
    { 32, 263130836933693530167218012160000000.0,     0 },
    { 33, 8683317618811886495518194401280000000.0,    0 },
    { 34, 2.95232799039604140847618609644e38,  0 },
    { 35, 1.03331479663861449296666513375e40,  0 },
    { 36, 3.71993326789901217467999448151e41,  0 },
    { 37, 1.37637530912263450463159795816e43,  0 },
    { 38, 5.23022617466601111760007224100e44,  0 },
    { 39, 2.03978820811974433586402817399e46,  0 },
    { 40, 8.15915283247897734345611269600e47,  0 },
    { 41, 3.34525266131638071081700620534e49,  0 },
    { 42, 1.40500611775287989854314260624e51,  0 },
    { 43, 6.04152630633738356373551320685e52,  0 },
    { 44, 2.65827157478844876804362581101e54,  0 },
    { 45, 1.19622220865480194561963161496e56,  0 },
    { 46, 5.50262215981208894985030542880e57,  0 },
    { 47, 2.58623241511168180642964355154e59,  0 },
    { 48, 1.24139155925360726708622890474e61,  0 },
    { 49, 6.08281864034267560872252163321e62,  0 },
    { 50, 3.04140932017133780436126081661e64,  0 },
    { 51, 1.55111875328738228022424301647e66,  0 },
    { 52, 8.06581751709438785716606368564e67,  0 },
    { 53, 4.27488328406002556429801375339e69,  0 },
    { 54, 2.30843697339241380472092742683e71,  0 },
    { 55, 1.26964033536582759259651008476e73,  0 },
    { 56, 7.10998587804863451854045647464e74,  0 },
    { 57, 4.05269195048772167556806019054e76,  0 },
    { 58, 2.35056133128287857182947491052e78,  0 },
    { 59, 1.38683118545689835737939019720e80,  0 },
    { 60, 8.32098711274139014427634118320e81,  0 },
    { 61, 5.07580213877224798800856812177e83,  0 },
    { 62, 3.14699732603879375256531223550e85,  0 },
    { 63, 1.982608315404440064116146708360e87,  0 },
    { 64, 1.268869321858841641034333893350e89,  0 },
    { 65, 8.247650592082470666723170306800e90,  0 },
    { 66, 5.443449390774430640037292402480e92,  0 },
    { 67, 3.647111091818868528824985909660e94,  0 },
    { 68, 2.480035542436830599600990418570e96,  0 },
    { 69, 1.711224524281413113724683388810e98,  0 },
    { 70, 1.197857166996989179607278372170e100,  0 },
    { 71, 8.504785885678623175211676442400e101,  0 },
    { 72, 6.123445837688608686152407038530e103,  0 },
    { 73, 4.470115461512684340891257138130e105,  0 },
    { 74, 3.307885441519386412259530282210e107,  0 },
    { 75, 2.480914081139539809194647711660e109,  0 },
    { 76, 1.885494701666050254987932260860e111,  0 },
    { 77, 1.451830920282858696340707840860e113,  0 },
    { 78, 1.132428117820629783145752115870e115,  0 },
    { 79, 8.946182130782975286851441715400e116,  0 },
    { 80, 7.156945704626380229481153372320e118,  0 },
    { 81, 5.797126020747367985879734231580e120,  0 },
    { 82, 4.753643337012841748421382069890e122,  0 },
    { 83, 3.945523969720658651189747118010e124,  0 },
    { 84, 3.314240134565353266999387579130e126,  0 },
    { 85, 2.817104114380550276949479442260e128,  0 },
    { 86, 2.422709538367273238176552320340e130,  0 },
    { 87, 2.107757298379527717213600518700e132,  0 },
    { 88, 1.854826422573984391147968456460e134,  0 },
    { 89, 1.650795516090846108121691926250e136,  0 },
    { 90, 1.485715964481761497309522733620e138,  0 },
    { 91, 1.352001527678402962551665687590e140,  0 },
    { 92, 1.243841405464130725547532432590e142,  0 },
    { 93, 1.156772507081641574759205162310e144,  0 },
    { 94, 1.087366156656743080273652852570e146,  0 },
    { 95, 1.032997848823905926259970209940e148,  0 },
    { 96, 9.916779348709496892095714015400e149,  0 },
    { 97, 9.619275968248211985332842594960e151,  0 },
    { 98, 9.426890448883247745626185743100e153,  0 },
    { 99, 9.332621544394415268169923885600e155,  0 },
    { 100, 9.33262154439441526816992388563e157,  0 },
    { 101, 9.42594775983835942085162312450e159,  0 },
    { 102, 9.61446671503512660926865558700e161,  0 },
    { 103, 9.90290071648618040754671525458e163,  0 },
    { 104, 1.02990167451456276238485838648e166,  0 },
    { 105, 1.08139675824029090050410130580e168,  0 },
    { 106, 1.146280563734708354534347384148e170,  0 },
    { 107, 1.226520203196137939351751701040e172,  0 },
    { 108, 1.324641819451828974499891837120e174,  0 },
    { 109, 1.443859583202493582204882102460e176,  0 },
    { 110, 1.588245541522742940425370312710e178,  0 },
    { 111, 1.762952551090244663872161047110e180,  0 },
    { 112, 1.974506857221074023536820372760e182,  0 },
    { 113, 2.231192748659813646596607021220e184,  0 },
    { 114, 2.543559733472187557120132004190e186,  0 },
    { 115, 2.925093693493015690688151804820e188,  0 },
    { 116, 3.393108684451898201198256093590e190,  0 },
    { 117, 3.96993716080872089540195962950e192,  0 },
    { 118, 4.68452584975429065657431236281e194,  0 },
    { 119, 5.57458576120760588132343171174e196,  0 },
    { 120, 6.68950291344912705758811805409e198,  0 },
    { 121, 8.09429852527344373968162284545e200,  0 },
    { 122, 9.87504420083360136241157987140e202,  0 },
    { 123, 1.21463043670253296757662432419e205,  0 },
    { 124, 1.50614174151114087979501416199e207,  0 },
    { 125, 1.88267717688892609974376770249e209,  0 },
    { 126, 2.37217324288004688567714730514e211,  0 },
    { 127, 3.01266001845765954480997707753e213,  0 },
    { 128, 3.85620482362580421735677065923e215,  0 },
    { 129, 4.97450422247728744039023415041e217,  0 },
    { 130, 6.46685548922047367250730439554e219,  0 },
    { 131, 8.47158069087882051098456875820e221,  0 },
    { 132, 1.11824865119600430744996307608e224,  0 },
    { 133, 1.48727070609068572890845089118e226,  0 },
    { 134, 1.99294274616151887673732419418e228,  0 },
    { 135, 2.69047270731805048359538766215e230,  0 },
    { 136, 3.65904288195254865768972722052e232,  0 },
    { 137, 5.01288874827499166103492629211e234,  0 },
    { 138, 6.91778647261948849222819828311e236,  0 },
    { 139, 9.61572319694108900419719561353e238,  0 },
    { 140, 1.34620124757175246058760738589e241,  0 },
    { 141, 1.89814375907617096942852641411e243,  0 },
    { 142, 2.69536413788816277658850750804e245,  0 },
    { 143, 3.85437071718007277052156573649e247,  0 },
    { 144, 5.55029383273930478955105466055e249,  0 },
    { 145, 8.04792605747199194484902925780e251,  0 },
    { 146, 1.17499720439091082394795827164e254,  0 },
    { 147, 1.72724589045463891120349865931e256,  0 },
    { 148, 2.55632391787286558858117801578e258,  0 },
    { 149, 3.80892263763056972698595524351e260,  0 },
    { 150, 5.71338395644585459047893286526e262,  0 },
    { 151, 8.62720977423324043162318862650e264,  0 },
    { 152, 1.31133588568345254560672467123e267,  0 },
    { 153, 2.00634390509568239477828874699e269,  0 },
    { 154, 3.08976961384735088795856467036e271,  0 },
    { 155, 4.78914290146339387633577523906e273,  0 },
    { 156, 7.47106292628289444708380937294e275,  0 },
    { 157, 1.17295687942641442819215807155e278,  0 },
    { 158, 1.85327186949373479654360975305e280,  0 },
    { 159, 2.94670227249503832650433950735e282,  0 },
    { 160, 4.71472363599206132240694321176e284,  0 },
    { 161, 7.59070505394721872907517857094e286,  0 },
    { 162, 1.22969421873944943411017892849e289,  0 },
    { 163, 2.00440157654530257759959165344e291,  0 },
    { 164, 3.28721858553429622726333031164e293,  0 },
    { 165, 5.42391066613158877498449501421e295,  0 },
    { 166, 9.00369170577843736647426172359e297,  0 },
    { 167, 1.50361651486499904020120170784e300,  0 },
    { 168, 2.52607574497319838753801886917e302,  0 },
    { 169, 4.26906800900470527493925188890e304,  0 },
    { 170, 7.25741561530799896739672821113e306,  0 }
  };

  double lnfact(const unsigned int n)
  {
    if(n <= FACT_TABLE_MAX)
      return log(fact_table[n].f);
    else
      return FMath::gamma::lngamma( n+1.0);
  }



  namespace zeta
  {


    /* coefficients for Maclaurin summation in hzeta()
     * B_{2j}/(2j)!
     */
    static double hzeta_c[15] = {
      1.00000000000000000000000000000,
      0.083333333333333333333333333333,
      -0.00138888888888888888888888888889,
      0.000033068783068783068783068783069,
      -8.2671957671957671957671957672e-07,
      2.0876756987868098979210090321e-08,
      -5.2841901386874931848476822022e-10,
      1.3382536530684678832826980975e-11,
      -3.3896802963225828668301953912e-13,
      8.5860620562778445641359054504e-15,
      -2.1748686985580618730415164239e-16,
      5.5090028283602295152026526089e-18,
      -1.3954464685812523340707686264e-19,
      3.5347070396294674716932299778e-21,
      -8.9535174270375468504026113181e-23
    };



    double hzeta(const double s, const double q)
    {
      if(s <= 1.0 || q <= 0.0)
        THROW_DEFAULT_EXCEPTION( FDomainError );

      const double max_bits = 54.0;
      const double ln_term0 = -s *log(q);

      if(ln_term0 < F_LOG_DBL_MIN +1.0)
      {
        THROW_EXCEPTION( FException, "Underflow" );
      }
      else if(ln_term0 > F_LOG_DBL_MAX -1.0)
      {
        THROW_EXCEPTION( FException, "Overflow" );
      }
      else if((s > max_bits && q < 1.0) || ( s > 0.5*max_bits && q < 0.25))
      {
        return pow(q,-s);
      }
      else if( s > 0.5*max_bits && q < 1.0)
      {
        const double p1 = pow(q,-s);
        const double p2 = pow(q/(1.0+q), s);
        const double p3 = pow(q/(2.0+q), s);
        return p1*(1.0+p2+p3);
      }
      else
      {
        // Euler-Maclaurin summation formula
        // Moshier, p 400 with several typo corrections
        const int jmax = 12;
        const int kmax = 10;

        int j, k;
        const double pmax = pow(kmax+q, -s);
        double scp = s;
        double pcp = pmax/(kmax+q);
        double ans = pmax*((kmax+q)/(s-1.0)+0.5);

        for(k=0; k<kmax; ++k)
          ans += pow(k+q, -s);

        for(j=0; j<=jmax; ++j)
        {
          double delta = hzeta_c[j+1] *scp*pcp;
          ans += delta;
          if(fabs(delta/ans) < 0.5*F_DBL_EPSILON) break;
          scp *=(s+2*j+1)*(s+2*j+2);
          pcp /= (kmax+q)*(kmax+q);
        }
        return ans;
      }
    }

  };

  //} // ns zeta

  namespace{

    struct cheb_series_struct {
      double *c;        // coefficients
      int order;        // order of expansion
      double a;         // lower interval point
      double b;         // upper interval point
      int order_sp;     // effective single precision order
    };
    typedef struct cheb_series_struct cheb_series;

    static inline double cheb_eval(const cheb_series * cs, const double x)
    {
      int j;
      double d  = 0.0;
      double dd = 0.0;

      double y  = (2.0*x - cs->a - cs->b) / (cs->b - cs->a);
      double y2 = 2.0 * y;

      double e = 0.0;

      for(j = cs->order; j>=1; j--) {
        double temp = d;
        d = y2*d - dd + cs->c[j];
        e += fabs(y2*temp) + fabs(dd) + fabs(cs->c[j]);
        dd = temp;
      }

      {
        double temp = d;
        d = y*d - dd + 0.5 * cs->c[0];
        e += fabs(y*temp) + fabs(dd) + 0.5 * fabs(cs->c[0]);
      }

      return d;
      //result->err = GSL_DBL_EPSILON * e + fabs(cs->c[cs->order]);

      //return GSL_SUCCESS;
    }
  } // anonymous ns


  namespace psi
  {

#define PSI_TABLE_NMAX 100
    static double psi_table[PSI_TABLE_NMAX+1] = {
      0.0,  /* Infinity */              /* psi(0) */
      -M_EULER,                          /* psi(1) */
      0.42278433509846713939348790992,  /* ...    */
      0.92278433509846713939348790992,
      1.25611766843180047272682124325,
      1.50611766843180047272682124325,
      1.70611766843180047272682124325,
      1.87278433509846713939348790992,
      2.01564147795560999653634505277,
      2.14064147795560999653634505277,
      2.25175258906672110764745616389,
      2.35175258906672110764745616389,
      2.44266167997581201673836525479,
      2.52599501330914535007169858813,
      2.60291809023222227314862166505,
      2.67434666166079370172005023648,
      2.74101332832746036838671690315,
      2.80351332832746036838671690315,
      2.86233685773922507426906984432,
      2.91789241329478062982462539988,
      2.97052399224214905087725697883,
      3.02052399224214905087725697883,
      3.06814303986119666992487602645,
      3.11359758531574212447033057190,
      3.15707584618530734186163491973,
      3.1987425128519740085283015864,
      3.2387425128519740085283015864,
      3.2772040513135124700667631249,
      3.3142410883505495071038001619,
      3.3499553740648352213895144476,
      3.3844381326855248765619282407,
      3.4177714660188582098952615740,
      3.4500295305349872421533260902,
      3.4812795305349872421533260902,
      3.5115825608380175451836291205,
      3.5409943255438998981248055911,
      3.5695657541153284695533770196,
      3.5973435318931062473311547974,
      3.6243705589201332743581818244,
      3.6506863483938174848844976139,
      3.6763273740348431259101386396,
      3.7013273740348431259101386396,
      3.7257176179372821503003825420,
      3.7495271417468059598241920658,
      3.7727829557002943319172153216,
      3.7955102284275670591899425943,
      3.8177324506497892814121648166,
      3.8394715810845718901078169905,
      3.8607481768292527411716467777,
      3.8815815101625860745049801110,
      3.9019896734278921969539597029,
      3.9219896734278921969539597029,
      3.9415975165651470989147440166,
      3.9608282857959163296839747858,
      3.9796962103242182164764276160,
      3.9982147288427367349949461345,
      4.0163965470245549168131279527,
      4.0342536898816977739559850956,
      4.0517975495308205809735289552,
      4.0690389288411654085597358518,
      4.0859880813835382899156680552,
      4.1026547480502049565823347218,
      4.1190481906731557762544658694,
      4.1351772229312202923834981274,
      4.1510502388042361653993711433,
      4.1666752388042361653993711433,
      4.1820598541888515500147557587,
      4.1972113693403667015299072739,
      4.2121367424746950597388624977,
      4.2268426248276362362094507330,
      4.2413353784508246420065521823,
      4.2556210927365389277208378966,
      4.2697055997787924488475984600,
      4.2835944886676813377364873489,
      4.2972931188046676391063503626,
      4.3108066323181811526198638761,
      4.3241399656515144859531972094,
      4.3372978603883565912163551041,
      4.3502848733753695782293421171,
      4.3631053861958823987421626300,
      4.3757636140439836645649474401,
      4.3882636140439836645649474401,
      4.4006092930563293435772931191,
      4.4128044150075488557724150703,
      4.4248526077786331931218126607,
      4.4367573696833950978837174226,
      4.4485220755657480390601880108,
      4.4601499825424922251066996387,
      4.4716442354160554434975042364,
      4.4830078717796918071338678728,
      4.4942438268358715824147667492,
      4.5053549379469826935258778603,
      4.5163439489359936825368668713,
      4.5272135141533849868846929582,
      4.5379662023254279976373811303,
      4.5486045001977684231692960239,
      4.5591308159872421073798223397,
      4.5695474826539087740464890064,
      4.5798567610044242379640147796,
      4.5900608426370772991885045755,
      4.6001618527380874001986055856
    };

#define PSI_1_TABLE_NMAX 100
    static double psi_1_table[PSI_1_TABLE_NMAX+1] = {
      0.0,  /* Infinity */              /* psi(1,0) */
      M_PI*M_PI/6.0,                    /* psi(1,1) */
      0.644934066848226436472415,       /* ...      */
      0.394934066848226436472415,
      0.2838229557371153253613041,
      0.2213229557371153253613041,
      0.1813229557371153253613041,
      0.1535451779593375475835263,
      0.1331370146940314251345467,
      0.1175120146940314251345467,
      0.1051663356816857461222010,
      0.0951663356816857461222010,
      0.0869018728717683907503002,
      0.0799574284273239463058557,
      0.0740402686640103368384001,
      0.0689382278476838062261552,
      0.0644937834032393617817108,
      0.0605875334032393617817108,
      0.0571273257907826143768665,
      0.0540409060376961946237801,
      0.0512708229352031198315363,
      0.0487708229352031198315363,
      0.0465032492390579951149830,
      0.0444371335365786562720078,
      0.0425467743683366902984728,
      0.0408106632572255791873617,
      0.0392106632572255791873617,
      0.0377313733163971768204978,
      0.0363596312039143235969038,
      0.0350841209998326909438426,
      0.0338950603577399442137594,
      0.0327839492466288331026483,
      0.0317433665203020901265817,
      0.03076680402030209012658168,
      0.02984853037475571730748159,
      0.02898347847164153045627052,
      0.02816715194102928555831133,
      0.02739554700275768062003973,
      0.02666508681283803124093089,
      0.02597256603721476254286995,
      0.02531510384129102815759710,
      0.02469010384129102815759710,
      0.02409521984367056414807896,
      0.02352832641963428296894063,
      0.02298749353699501850166102,
      0.02247096461137518379091722,
      0.02197713745088135663042339,
      0.02150454765882086513703965,
      0.02105185413233829383780923,
      0.02061782635456051606003145,
      0.02020133322669712580597065,
      0.01980133322669712580597065,
      0.01941686571420193164987683,
      0.01904704322899483105816086,
      0.01869104465298913508094477,
      0.01834810912486842177504628,
      0.01801753061247172756017024,
      0.01769865306145131939690494,
      0.01739086605006319997554452,
      0.01709360088954001329302371,
      0.01680632711763538818529605,
      0.01652854933985761040751827,
      0.01625980437882562975715546,
      0.01599965869724394401313881,
      0.01574770606433893015574400,
      0.01550356543933893015574400,
      0.01526687904880638577704578,
      0.01503731063741979257227076,
      0.01481454387422086185273411,
      0.01459828089844231513993134,
      0.01438824099085987447620523,
      0.01418415935820681325171544,
      0.01398578601958352422176106,
      0.01379288478501562298719316,
      0.01360523231738567365335942,
      0.01342261726990576130858221,
      0.01324483949212798353080444,
      0.01307170929822216635628920,
      0.01290304679189732236910755,
      0.01273868124291638877278934,
      0.01257845051066194236996928,
      0.01242220051066194236996928,
      0.01226978472038606978956995,
      0.01212106372098095378719041,
      0.01197590477193174490346273,
      0.01183418141592267460867815,
      0.01169577311142440471248438,
      0.01156056489076458859566448,
      0.01142844704164317229232189,
      0.01129931481023821361463594,
      0.01117306812421372175754719,
      0.01104961133409026496742374,
      0.01092885297157366069257770,
      0.01081070552355853781923177,
      0.01069508522063334415522437,
      0.01058191183901270133041676,
      0.01047110851491297833872701,
      0.01036260157046853389428257,
      0.01025632035036012704977199,  /* ...        */
      0.01015219706839427948625679,  /* psi(1,99)  */
      0.01005016666333357139524567   /* psi(1,100) */
    };

    double psi_int(const int n)
    {

      if(n <= 0) {
        THROW_DEFAULT_EXCEPTION( FDomainError ); //DOMAIN_ERROR(result);
      }
      else if(n <= PSI_TABLE_NMAX) {
        return psi_table[n];
        //result->val = psi_table[n];
        //result->err = GSL_DBL_EPSILON * fabs(result->val);
        //return GSL_SUCCESS;
      }
      else {
        /* Abramowitz+Stegun 6.3.18 */
        const double c2 = -1.0/12.0;
        const double c3 =  1.0/120.0;
        const double c4 = -1.0/252.0;
        const double c5 =  1.0/240.0;
        const double ni2 = (1.0/n)*(1.0/n);
        const double ser = ni2 * (c2 + ni2 * (c3 + ni2 * (c4 + ni2*c5)));
        //result->val  = log(n) - 0.5/n + ser;
        //result->err  = GSL_DBL_EPSILON * (fabs(log(n)) + fabs(0.5/n) + fabs(ser));
        //result->err += GSL_DBL_EPSILON * fabs(result->val);
        //return GSL_SUCCESS;
        return log((double)n) - 0.5/n +ser;
      }
    }

    double psi_1_int(const int n)
    {
      /* CHECK_POINTER(result) */
      if(n <= 0) {
        THROW_DEFAULT_EXCEPTION( FDomainError ); //DOMAIN_ERROR(result);
      }
      else if(n <= PSI_1_TABLE_NMAX) {
        return psi_1_table[n];
        //result->err = GSL_DBL_EPSILON * result->val;
        //return GSL_SUCCESS;
      }
      else {
        /* Abramowitz+Stegun 6.4.12
         * double-precision for n > 100
         */
        const double c0 = -1.0/30.0;
        const double c1 =  1.0/42.0;
        const double c2 = -1.0/30.0;
        const double ni2 = (1.0/n)*(1.0/n);
        const double ser =  ni2*ni2 * (c0 + ni2*(c1 + c2*ni2));
        return (1.0 + 0.5/n + 1.0/(6.0*n*n) + ser) / n;
        //result->err = GSL_DBL_EPSILON * result->val;
        //return GSL_SUCCESS;
      }
    }




    /* Chebyshev fits from SLATEC code for psi(x)

       Series for PSI        on the interval  0.         to  1.00000D+00
       with weighted error   2.03E-17
       log weighted error  16.69
       significant figures required  16.39
       decimal places required  17.37

       Series for APSI       on the interval  0.         to  2.50000D-01
       with weighted error   5.54E-17
       log weighted error  16.26
       significant figures required  14.42
       decimal places required  16.86

*/

    static double psics_data[23] = {
      -.038057080835217922,
      .491415393029387130,
      -.056815747821244730,
      .008357821225914313,
      -.001333232857994342,
      .000220313287069308,
      -.000037040238178456,
      .000006283793654854,
      -.000001071263908506,
      .000000183128394654,
      -.000000031353509361,
      .000000005372808776,
      -.000000000921168141,
      .000000000157981265,
      -.000000000027098646,
      .000000000004648722,
      -.000000000000797527,
      .000000000000136827,
      -.000000000000023475,
      .000000000000004027,
      -.000000000000000691,
      .000000000000000118,
      -.000000000000000020
    };
    static cheb_series psi_cs = {
      psics_data,
      22,
      -1, 1,
      17
    };


    static double apsics_data[16] = {    
      -.0204749044678185,
      -.0101801271534859,
      .0000559718725387,
      -.0000012917176570,
      .0000000572858606,
      -.0000000038213539,
      .0000000003397434,
      -.0000000000374838,
      .0000000000048990,
      -.0000000000007344,
      .0000000000001233,
      -.0000000000000228,
      .0000000000000045,
      -.0000000000000009,
      .0000000000000002,
      -.0000000000000000 
    };    

    static cheb_series apsi_cs = {
      apsics_data,
      15,
      -1, 1,
      9
    };


    /* digamma for x both positive and negative; we do both
     * cases here because of the way we use even/odd parts
     * of the function
     */
    double
      psi_x(const double x)
      {
        const double y = fabs(x);

        if(x == 0.0 || x == -1.0 || x == -2.0) {
          THROW_DEFAULT_EXCEPTION( FDomainError ); //DOMAIN_ERROR(result);
        }
        else if(y >= 2.0) {
          const double t = 8.0/(y*y)-1.0;
          //gsl_sf_result result_c;
          //cheb_eval_e(&apsi_cs, t, &result_c);
          double result_c = cheb_eval(&apsi_cs, t );
          if(x < 0.0) {
            const double s = sin(M_PI*x);
            const double c = cos(M_PI*x);
            if(fabs(s) < 2.0*F_SQRT_DBL_MIN) {
              THROW_DEFAULT_EXCEPTION( FDomainError ); //DOMAIN_ERROR(result);
            }
            else {
              return log(y) - 0.5/x + result_c- M_PI * c/s;
              //result->err  = M_PI*fabs(x)*F_DBL_EPSILON/(s*s);
              //result->err += result_c.err;
              //result->err += F_DBL_EPSILON * fabs(result->val);
              //return F_SUCCESS;
            }
          }
          else {
            return log(y) - 0.5/x + result_c;
            //result->err  = result_c.err;
            //result->err += F_DBL_EPSILON * fabs(result->val);
            //return F_SUCCESS;
          }
        }
        else { /* -2 < x < 2 */
          //gsl_sf_result result_c;

          if(x < -1.0) { /* x = -2 + v */
            const double v  = x + 2.0;
            const double t1 = 1.0/x;
            const double t2 = 1.0/(x+1.0);
            const double t3 = 1.0/v;
            //cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);
            double result_c = cheb_eval( &psi_cs, 2.0*v-1.0);
            return  -(t1 + t2 + t3) + result_c;
            //result->err  = F_DBL_EPSILON * (fabs(t1) + fabs(x/(t2*t2)) + fabs(x/(t3*t3)));
            //result->err += result_c.err;
            //result->err += F_DBL_EPSILON * fabs(result->val);
            //return F_SUCCESS;
          }
          else if(x < 0.0) { /* x = -1 + v */
            const double v  = x + 1.0;
            const double t1 = 1.0/x;
            const double t2 = 1.0/v;
            //cheb_eval_e(&psi_cs, 2.0*v-1.0, &result_c);
            double result_c = cheb_eval( &psi_cs, 2.0*v-1.0);
            return -(t1 + t2) + result_c;
            //result->err  = F_DBL_EPSILON * (fabs(t1) + fabs(x/(t2*t2)));
            //result->err += result_c.err;
            //result->err += F_DBL_EPSILON * fabs(result->val);
            //return F_SUCCESS;
          }
          else if(x < 1.0) { /* x = v */
            const double t1 = 1.0/x;
            //cheb_eval_e(&psi_cs, 2.0*x-1.0, &result_c);
            double result_c = cheb_eval(&psi_cs, 2.0*x-1.0);
            return -t1 + result_c;
            //result->err  = F_DBL_EPSILON * t1;
            //result->err += result_c.err;
            //result->err += F_DBL_EPSILON * fabs(result->val);
            //return F_SUCCESS;
          }
          else { /* x = 1 + v */
            const double v = x - 1.0;
            return cheb_eval(&psi_cs, 2.0*v-1.0);
          }
        }
      }





    /* generic polygamma; assumes n >= 0 and x > 0
    */
    static double
      psi_n_xg0(const int n, const double x)
      {
        if(n == 0) {
          return psi_x(x);
        }
        else {
          /* Abramowitz + Stegun 6.4.10 */
          //gsl_sf_result ln_nf;
          //gsl_sf_result hzeta;
          //int stat_hz = gsl_sf_hzeta_e(n+1.0, x, &hzeta);
          //int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
          //int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
          //    hzeta.val, hzeta.err,
          //    result);
          //if(GSL_IS_EVEN(n)) result->val = -result->val;
          //return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
          double hz = FMath::zeta::hzeta(n+1.0, x);
          double ln_nf = lnfact((unsigned int) n );
          double result = exp_mult(ln_nf, hz);
          if(!IS_ODD(n)) return -result;
          return result;
        } 
      } 


    double psi_1(const double x)
    {
      /* CHECK_POINTER(result) */

      if(x == 0.0 || x == -1.0 || x == -2.0) {
        THROW_DEFAULT_EXCEPTION( FDomainError );
      }
      else if(x > 0.0)
      {
        return psi_n_xg0(1, x);
      }
      else if(x > -5.0)
      {
        /* Abramowitz + Stegun 6.4.6 */
        int M = (int)-floor(x);
        double fx = x + M;
        double sum = 0.0;
        int m;

        if(fx == 0.0)
          THROW_DEFAULT_EXCEPTION( FDomainError );

        for(m = 0; m < M; ++m)
          sum += 1.0/((x+m)*(x+m));

        {
          return sum + psi_n_xg0(1, fx);
          //result->val += sum;
          //result->err += M * GSL_DBL_EPSILON * sum;
          //return stat_psi;
        }
      }
      else
      {
        /* Abramowitz + Stegun 6.4.7 */
        const double sin_px = sin(M_PI * x);
        const double d = M_PI*M_PI/(sin_px*sin_px);
        double r = psi_n_xg0(1, 1.0-x);
        //result->val = d - r.val;
        //result->err = r.err + 2.0*GSL_DBL_EPSILON*d;
        //return stat_psi;
        return d-r;
      }
    }



    double psi_n(const int n, const double x )
    {
      if(n == 0)
      {
        return psi_x(x);
      }
      else if(n == 1)
      {
        return psi_1(x);
      }
      else if(n < 0 || x <= 0.0) {
        THROW_DEFAULT_EXCEPTION( FDomainError ); //DOMAIN_ERROR(result);
      }
      else {
        //THROW_DEFAULT_EXCEPTION( FNotImplementedException );
#if 0
        gsl_sf_result ln_nf;
        gsl_sf_result hzeta;
        int stat_hz = gsl_sf_hzeta_e(n+1.0, x, &hzeta);
        int stat_nf = gsl_sf_lnfact_e((unsigned int) n, &ln_nf);
        int stat_e  = gsl_sf_exp_mult_err_e(ln_nf.val, ln_nf.err,
            hzeta.val, hzeta.err,
            result);
        if(IS_EVEN(n)) result->val = -result->val;
        return GSL_ERROR_SELECT_3(stat_e, stat_nf, stat_hz);
#endif
       double hz = FMath::zeta::hzeta(n+1.0, x );
       double ln_nf = lnfact((unsigned int)n);
       double result = exp_mult(ln_nf, hz);
       //std::cout << "hz=" << hz << "  ln_nf=" << ln_nf << " res=" << result << std::endl;
       if( !IS_ODD(n)) result = -result;
       return result;
      }
      eassert(false); // never reached
    }


  } // namespace psi

  //} // ns psi


#define F_GAMMA_XMAX 171.0

double pow_int(double x, int n)
{
  double value = 1.0;
  int count = 0;

  if(n<0)
  {
    n = -n;

    if(x == 0.0)
    {
      double u = 1.0/x;
      return( n%2)? u: (u*u); // correct sign of infinity
    }
    x = 1.0/x;
  }

  // repeat squaring method
  // return 0.0^0 = 1.0 so continuous in x
  do
  {
    if(IS_ODD(n)) value *= x;
    n >>= 1;
    x *= x;
    ++count;
  }while(n);
  return value;
}


namespace gamma{

#define F_SF_GAMMA_XMAX (171.0)
  
  double lngamma_sgn(double& sgn, const double x);
  static double gamma_xgthalf(const double x);

  /*
  double lgamma(const double x)
  {
    if(x<0) THROW_EXCEPTION( FException, "x must be >= 0 for lgamma function");
    return ((x + 0.5) * log(x + 5.5) - x - 5.5 
        + log((2.50662827510730 + 190.9551718930764 / (x + 1)
            - 216.8366818437280 / (x + 2) + 60.19441764023333 
            / (x + 3) - 3.08751323928546 / (x + 4) + 0.00302963870525
            / (x + 5) - 0.00001352385959072596 / (x + 6)) / x));
  }
  */
  
  double gamma(const double x)
  {
    double s = (2.50662827510730 + 190.9551718930764 / (x + 1)
        - 216.8366818437280 / (x + 2) + 60.19441764023333
        / (x + 3) - 3.08751323928546 / (x + 4) + 0.00302963870525
        / (x + 5) - 0.00001352385959072596 / (x + 6)) / x;
    if (s < 0) 
      return (-exp((x + 0.5) * log(x + 5.5) - x - 5.5 + log(-s)));
    else 
      return exp((x + 0.5) * log(x + 5.5) - x - 5.5 + log(s));
  }


  /** series for gammastar(x)
   * doubel-precision for x > 10.0
   */
  double gammastar_ser( const double x)
  {
    /*
    Use the Stirling series for the correction to Log(Gamma(x)),
    which is better behaved and easier to compute than the
    regular Stirling series for Gamma(x)
    */
    const double y = 1.0/(x*x);
     
    const double c0 = 1.0/12.0;
    const double c1 =-1.0/360.0;
    const double c2 = 1.0/1260.0;
    const double c3 =-1.0/1680.0;
    const double c4 = 1.0/1188.0;
    const double c5 =-691.0/360360.0;
    const double c6 = 1.0/156.0;
    const double c7 =-3617.0/122400.0;
    const double ser = c0+y*(c1+y*(c2+y*(c3+y*(c4+y*(c5+y*(c6+y*c7))))));
    return exp(ser/x);
  }

/* Chebyshev expansion for log(gamma(x)/gamma(8))
 * 5 < x < 10
 * -1 < t < 1
 */   
static double gamma_5_10_data[24] = {
 -1.5285594096661578881275075214,
  4.8259152300595906319768555035,
  0.2277712320977614992970601978,
 -0.0138867665685617873604917300,
  0.0012704876495201082588139723,
 -0.0001393841240254993658962470,
  0.0000169709242992322702260663,
 -2.2108528820210580075775889168e-06,
  3.0196602854202309805163918716e-07,
 -4.2705675000079118380587357358e-08,
  6.2026423818051402794663551945e-09,
 -9.1993973208880910416311405656e-10,
  1.3875551258028145778301211638e-10,
 -2.1218861491906788718519522978e-11,
  3.2821736040381439555133562600e-12,
 -5.1260001009953791220611135264e-13,
  8.0713532554874636696982146610e-14,
 -1.2798522376569209083811628061e-14,
  2.0417711600852502310258808643e-15,
 -3.2745239502992355776882614137e-16,
  5.2759418422036579482120897453e-17,
 -8.5354147151695233960425725513e-18,
  1.3858639703888078291599886143e-18,
 -2.2574398807738626571560124396e-19
};  
static const cheb_series gamma_5_10_cs = {
  gamma_5_10_data,
  23,
  -1, 1,
  11
};



  
  double gammainv(const double x)
  {
    if(x <= 0.0 && x == floor(x))
    {
      // negative integer
      return 0.0;
    }
    else if(x< 0.5)
    {
      double sgn;
      double ret;
      try{
        ret = FMath::gamma::lngamma_sgn(sgn, x);
      } catch( FDomainError )
      {
        ret = 0.0;
      }
      return ret;
    }
    else
    {
      double g;
      try{
        g = FMath::gamma::gamma_xgthalf(x);
      }CATCH_N_RETHROW(FException);
      return 1.0/g; // hope there is no underflow
    }
    //THROW_DEFAULT_EXCEPTION( FNotImplementedException );
  }

  inline static double lngamma_1_pade(const double eps)
  {
    /* use (2,2) pade for log[gamma[1+eps]]/eps
     * plus a correction series.
     */
    const double n1 = -1.0017419282349508699871138440;
    const double n2 =  1.7364839209922879823280541733;
    const double d1 =  1.2433006018858751556055436011;
    const double d2 =  5.0456274100274010152489597514;
    const double num = (eps + n1) * (eps + n2);
    const double den = (eps + d1) * (eps + d2);
    const double pade = 2.0816265188662692474880210318 * num / den;
    const double c0 =  0.004785324257581753;
    const double c1 = -0.01192457083645441;
    const double c2 =  0.01931961413960498;
    const double c3 = -0.02594027398725020;
    const double c4 =  0.03141928755021455;
    const double eps5 = eps*eps*eps*eps*eps;
    const double corr = eps5 * (c0 + eps*(c1 + eps*(c2 + eps*(c3 + c4*eps))));
    return eps * (pade + corr);
    //result->err = 2.0 * F_DBL_EPSILON * fabs(result->val);
    //return F_SUCCESS;
  }

  inline static double lngamma_2_pade(const double eps)
  {
    /* Use (2,2) Pade for Log[Gamma[2+eps]]/eps
     *    * plus a correction series.
     *       */
    const double n1 = 1.000895834786669227164446568;
    const double n2 = 4.209376735287755081642901277;
    const double d1 = 2.618851904903217274682578255;
    const double d2 = 10.85766559900983515322922936;
    const double num = (eps + n1) * (eps + n2);
    const double den = (eps + d1) * (eps + d2);
    const double pade = 2.85337998765781918463568869 * num/den;
    const double c0 =  0.0001139406357036744;
    const double c1 = -0.0001365435269792533;
    const double c2 =  0.0001067287169183665;
    const double c3 = -0.0000693271800931282;
    const double c4 =  0.0000407220927867950;
    const double eps5 = eps*eps*eps*eps*eps;
    const double corr = eps5 * (c0 + eps*(c1 + eps*(c2 + eps*(c3 + c4*eps))));
    return eps * (pade + corr);
    //result->err = 2.0 * F_DBL_EPSILON * fabs(result->val);
    //return F_SUCCESS;
  }

  /**
   * x = eps near zero
   * gives double-precision for |eps| < 0.02
   *
   * returns eps
   */
  static double lngamma_sgn_0(double &sgn, const double eps)
  {
    // calculate series for g(eps) = Gamma(eps) eps-1/(1+eps)-eps/2
    const double c1  = -0.07721566490153286061;
    const double c2  = -0.01094400467202744461;
    const double c3  =  0.09252092391911371098;
    const double c4  = -0.01827191316559981266;
    const double c5  =  0.01800493109685479790;
    const double c6  = -0.00685088537872380685;
    const double c7  =  0.00399823955756846603;
    const double c8  = -0.00189430621687107802;
    const double c9  =  0.00097473237804513221;
    const double c10 = -0.00048434392722255893;
    const double g6  = c6+eps*(c7+eps*(c8 + eps*(c9 + eps*c10)));
    const double g   = eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*g6)))));

    /* calculate Gamma(eps) eps, a positive quantity */
    const double gee = g + 1.0/(1.0+eps) + 0.5*eps;

    sgn = F_SIGN(eps);

    return log(gee/fabs(eps));
    //lng->err = 4.0 * F_DBL_EPSILON * fabs(lng->val);
    //return GSL_SUCCESS;
  }   


  //namespace detail
  //{


  /* coefficients for gamma=7, kmax=8  Lanczos method */
  static double lanczos_7_c[9] = {
    0.99999999999980993227684700473478,
    676.520368121885098567009190444019,
    -1259.13921672240287047156078755283,
    771.3234287776530788486528258894,
    -176.61502916214059906584551354,
    12.507343278686904814458936853,
    -0.13857109526572011689554707,
    9.984369578019570859563e-6,
    1.50563273514931155834e-7
  };

#define LogRootTwoPi_  0.9189385332046727418

  /* Lanczos method for real x > 0;
   * gamma=7, truncated at 1/(z+8) 
   * [J. SIAM Numer. Anal, Ser. B, 1 (1964) 86]
   */
  static
    double
    lngamma_lanczos(double x)
    {
      int k;
      double Ag;
      double term1, term2;

      x -= 1.0; /* Lanczos writes z! instead of Gamma(z) */

      Ag = lanczos_7_c[0];
      for(k=1; k<=8; k++) { Ag += lanczos_7_c[k]/(x+k); }

      /* (x+0.5)*log(x+7.5) - (x+7.5) + LogRootTwoPi_ + log(Ag(x)) */
      term1 = (x+0.5)*log((x+7.5)/M_E);
      term2 = LogRootTwoPi_ + log(Ag);
      return term1 + (term2 - 7.0);
      //result->err  = 2.0 *   F_DBL_EPSILON * (fabs(term1) + fabs(term2) + 7.0);
      //result->err +=   F_DBL_EPSILON * fabs(result->val);

      //return GSL_SUCCESS;
    }


  /* x near a negative integer
   * Calculates sign as well as log(|gamma(x)|).
   * x = -N + eps
   * assumes N >= 1
   */
  static
    double
    lngamma_sgn_sing(double &sgn, const int N, const double eps)
    {
      if(eps == 0.0) {
        //lng->val = 0.0;
        //lng->err = 0.0;
        //*sgn = 0.0;
        //GSL_ERROR ("error", GSL_EDOM);
        THROW_DEFAULT_EXCEPTION( FDomainError );
      }
      else if(N == 1) {
        /* calculate series for
         * g = eps gamma(-1+eps) + 1 + eps/2 (1+3eps)/(1-eps^2)
         * double-precision for |eps| < 0.02
         */
        const double c0 =  0.07721566490153286061;
        const double c1 =  0.08815966957356030521;
        const double c2 = -0.00436125434555340577;
        const double c3 =  0.01391065882004640689;
        const double c4 = -0.00409427227680839100;
        const double c5 =  0.00275661310191541584;
        const double c6 = -0.00124162645565305019;
        const double c7 =  0.00065267976121802783;
        const double c8 = -0.00032205261682710437;
        const double c9 =  0.00016229131039545456;
        const double g5 = c5 + eps*(c6 + eps*(c7 + eps*(c8 + eps*c9)));
        const double g  = eps*(c0 + eps*(c1 + eps*(c2 + eps*(c3 + eps*(c4 + eps*g5)))));

        /* calculate eps gamma(-1+eps), a negative quantity */
        const double gam_e = g - 1.0 - 0.5*eps*(1.0+3.0*eps)/(1.0 - eps*eps);

        //lng->val = log(fabs(gam_e)/fabs(eps));
        //lng->err = 2.0 * GSL_DBL_EPSILON * fabs(lng->val);
        sgn = ( eps > 0.0 ? -1.0 : 1.0 );
        return log(fabs(gam_e)/fabs(eps));
        //return GSL_SUCCESS;
      }
      else {
        double g;

        /* series for sin(Pi(N+1-eps))/(Pi eps) modulo the sign
         * double-precision for |eps| < 0.02
         */
        const double cs1 = -1.6449340668482264365;
        const double cs2 =  0.8117424252833536436;
        const double cs3 = -0.1907518241220842137;
        const double cs4 =  0.0261478478176548005;
        const double cs5 = -0.0023460810354558236;
        const double e2  = eps*eps;
        const double sin_ser = 1.0 + e2*(cs1+e2*(cs2+e2*(cs3+e2*(cs4+e2*cs5))));

        /* calculate series for ln(gamma(1+N-eps))
         * double-precision for |eps| < 0.02
         */
        double aeps = fabs(eps);
        double c1, c2, c3, c4, c5, c6, c7;
        double lng_ser;
        //gsl_sf_result c0;
        //gsl_sf_result psi_0;
        //gsl_sf_result psi_1;
        //gsl_sf_result psi_2;
        //gsl_sf_result psi_3;
        //gsl_sf_result psi_4;
        //gsl_sf_result psi_5;
        //gsl_sf_result psi_6;
        //psi_2.val = 0.0;
        //psi_3.val = 0.0;
        //psi_4.val = 0.0;
        //psi_5.val = 0.0;
        //psi_6.val = 0.0;
        double c0 = lnfact(N);
        double psi_0 = FMath::psi::psi_int(N+1);
        double psi_1 = FMath::psi::psi_1_int(N+1);
        double psi_2(0.0), psi_3(0.0), psi_4(0.0), psi_5(0.0), psi_6(0.0);
        if(aeps > 0.00001) psi_2 = FMath::psi::psi_n(2, N+1.0 );
        if(aeps > 0.0002)  psi_3 = FMath::psi::psi_n(3, N+1.0 );
        if(aeps > 0.001)   psi_4 = FMath::psi::psi_n(4, N+1.0 );
        if(aeps > 0.005)   psi_5 = FMath::psi::psi_n(5, N+1.0 );
        if(aeps > 0.01)    psi_6 = FMath::psi::psi_n(6, N+1.0 );
        c1 = psi_0;
        c2 = psi_1/2.0;
        c3 = psi_2/6.0;
        c4 = psi_3/24.0;
        c5 = psi_4/120.0;
        c6 = psi_5/720.0;
        c7 = psi_6/5040.0;
        lng_ser = c0-eps*(c1-eps*(c2-eps*(c3-eps*(c4-eps*(c5-eps*(c6-eps*c7))))));

        /* calculate
         * g = ln(|eps gamma(-N+eps)|)
         *   = -ln(gamma(1+N-eps)) + ln(|eps Pi/sin(Pi(N+1+eps))|)
         */
        g = -lng_ser - log(sin_ser);

        //lng->val = g - log(fabs(eps));
        //lng->err = c0.err + 2.0 * GSL_DBL_EPSILON * (fabs(g) + fabs(lng->val));

        sgn = ( IS_ODD(N) ? -1.0 : 1.0 ) * ( eps > 0.0 ? 1.0 : -1.0 );

        return g - log(fabs(eps));
        //return GSL_SUCCESS;
      }
    }




  double lngamma_sgn(double &sgn, const double x)
  {
    if(fabs(x-1.0) < 0.01)
    {
      sgn = 1.0;
      return lngamma_1_pade( x-1.0 );
    }
    else if(fabs(x-2.0) < 0.01)
    {
      sgn = 1.0;
      return lngamma_2_pade(x-2.0);
    }
    else if( x >= 0.5)
    {
      sgn = 1.0;
      return lngamma_lanczos(x);
    }
    else if ( x == 0.0 )
    {
      sgn = 0.0;
      THROW_DEFAULT_EXCEPTION( FDomainError );
    }
    else if( fabs(x) < 0.02)
    {
      return lngamma_sgn_0(sgn, x);
    }
    else if( x > -0.5/(F_DBL_EPSILON*M_PI))
    {
      double z = 1.0-x;
      double s = sin(M_PI*x);
      double as = fabs(s);
      if( s==0.0)
      {
        sgn = 0.0;
        THROW_DEFAULT_EXCEPTION( FDomainError );
      }
      else if(as < M_PI*0.015)
      { // x near negative integer, -N
        if( x < INT_MIN +2.0)
        {
          sgn = 0.0;
          THROW_EXCEPTION( FException, "error" );
        }
        else
        {
          int N = -(int)(x - 0.5);
          double eps = x+N;
          return lngamma_sgn_sing(sgn,N, eps);
        }
      }
      else
      {
        double lg_z = lngamma_lanczos(z);
        sgn = (s > 0.0 ? 1.0: -1.0);
        return M_LNPI -(log(as) +lg_z);
            }
            }
            else
            {
            THROW_EXCEPTION( FException, "eround error" );
            }
            }



            double lngamma(const double x)
            {
            /* CHECK_POINTER(result) */

            if(fabs(x - 1.0) < 0.01) {
            /* Note that we must amplify the errors
             * from the Pade evaluations because of
             * the way we must pass the argument, i.e.
             * writing (1-x) is a loss of precision
             * when x is near 1.
             */
              return lngamma_1_pade(x - 1.0);
              //result->err *= 1.0/(F_DBL_EPSILON + fabs(x - 1.0));
              //return stat;
            }
            else if(fabs(x - 2.0) < 0.01) {
              return lngamma_2_pade(x - 2.0);
              //result->err *= 1.0/(F_DBL_EPSILON + fabs(x - 2.0));
              //return stat;
            }
            else if(x >= 0.5) {
              return lngamma_lanczos(x);
            }
            else if(x == 0.0) {
              THROW_DEFAULT_EXCEPTION( FDomainError ); //DOMAIN_ERROR(result);
            }
            else if(fabs(x) < 0.02) {
              double sgn;
              return lngamma_sgn_0(sgn, x );
            }
            else if(x > -0.5/(F_DBL_EPSILON*M_PI)) {
              /* Try to extract a fractional
               * part from x.
               */
              double z  = 1.0 - x;
              double s  = sin(M_PI*z);
              double as = fabs(s);
              if(s == 0.0) {
                THROW_DEFAULT_EXCEPTION( FDomainError ); //DOMAIN_ERROR(result);
              }
              else if(as < M_PI*0.015) {
                /* x is near a negative integer, -N */
                if(x < INT_MIN + 2.0) {
                  //result->val = 0.0;
                  //result->err = 0.0;
                  //F_ERROR ("error", F_EROUND);
                  THROW_EXCEPTION( FException, "??? error ???" );
                }
                else {
                  int N = -(int)(x - 0.5);
                  double eps = x + N;
                  double sgn;
                  return lngamma_sgn_sing(sgn, N, eps);
                }
              }
              else {
                //gsl_sf_result lg_z;
                double lg_z = lngamma_lanczos(z);
                return M_LNPI - (log(as) + lg_z);
                //result->err = 2.0 * F_DBL_EPSILON * fabs(result->val) + lg_z.err;
                //return F_SUCCESS;
              }
            }
            else {
              /* |x| was too large to extract any fractional part */
              //result->val = 0.0;
              //result->err = 0.0;
              //F_ERROR ("error", F_EROUND);
              THROW_EXCEPTION( FException, "|x| too large to extract fractional part" );
            }
            //THROW_DEFAULT_EXCEPTION( FNotImplementedException );
            }

        /* gamma(x) for x >= 1/2
         * assumes x >= 1/2 
         */   
        static double gamma_xgthalf(const double x) 
        {
          if(x == 0.5) {
            return 1.77245385090551602729817;
            //result->err = F_DBL_EPSILON * result->val;
            //return F_SUCCESS;
          } else if (x <= (FACT_TABLE_MAX + 1.0) && x == floor(x)) {
            int n = (int) floor (x);
            return fact_table[n - 1].f;
            //result->err = F_DBL_EPSILON * result->val;
            //return F_SUCCESS;
          }    
          else if(fabs(x - 1.0) < 0.01) {
            /* Use series for Gamma[1+eps] - 1/(1+eps).
            */
            const double eps = x - 1.0;
            const double c1 =  0.4227843350984671394;
            const double c2 = -0.01094400467202744461;
            const double c3 =  0.09252092391911371098;
            const double c4 = -0.018271913165599812664;
            const double c5 =  0.018004931096854797895;
            const double c6 = -0.006850885378723806846;
            const double c7 =  0.003998239557568466030;
            return 1.0/x + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*c7))))));
            //result->err = F_DBL_EPSILON;
            //return F_SUCCESS;
          }
          else if(fabs(x - 2.0) < 0.01) {
            /* Use series for Gamma[1 + eps].
            */
            const double eps = x - 2.0;
            const double c1 =  0.4227843350984671394;
            const double c2 =  0.4118403304264396948;
            const double c3 =  0.08157691924708626638;
            const double c4 =  0.07424901075351389832;
            const double c5 = -0.00026698206874501476832;
            const double c6 =  0.011154045718130991049;
            const double c7 = -0.002852645821155340816;
            const double c8 =  0.0021039333406973880085;
            return 1.0 + eps*(c1+eps*(c2+eps*(c3+eps*(c4+eps*(c5+eps*(c6+eps*(c7+eps*c8)))))));
            //result->err = F_DBL_EPSILON;
            //return F_SUCCESS;
          }
          else if(x < 5.0) {
            /* Exponentiating the logarithm is fine, as
             * long as the exponential is not so large
             * that it greatly amplifies the error.
             */
            //gsl_sf_result lg;
            double lg = lngamma_lanczos(x);
            return exp(lg);
            //result->err = result->val * (lg.err + 2.0 * F_DBL_EPSILON);
            //return F_SUCCESS;
          }
          else if(x < 10.0) {
            /* This is a sticky area. The logarithm
             * is too large and the gammastar series
             * is not good.
             */
            const double gamma_8 = 5040.0;
            const double t = (2.0*x - 15.0)/5.0;
            //gsl_sf_result c;
            //cheb_eval_e(&gamma_5_10_cs, t, &c);
            //result->val  = exp(c.val) * gamma_8;
            //result->err  = result->val * c.err;
            //result->err += 2.0 * F_DBL_EPSILON * result->val;
            //return F_SUCCESS;
            return exp( cheb_eval( &gamma_5_10_cs, t ))*gamma_8;
                }
                else if(x < F_SF_GAMMA_XMAX) {
                /* We do not want to exponentiate the logarithm
                 * if x is large because of the inevitable
                 * inflation of the error. So we carefully
                 * use pow() and exp() with exact quantities.
                 */
                double p = pow(x, 0.5*x);
                double e = exp(-x);
                double q = (p * e) * p;
                double pre = M_SQRT2 * M_SQRTPI * q/sqrt(x);
                //gsl_sf_result gstar;
                //int stat_gs = gammastar_ser(x, &gstar);
                //result->val = pre * gstar.val;
                //result->err = (x + 2.5) * F_DBL_EPSILON * result->val;
                //return stat_gs;
                double gstar = gammastar_ser(x);
                return pre* gstar;
                }
                else {
                THROW_EXCEPTION( FException, "Overflow" ); //OVERFLOW_ERROR(result);
                }
        }



} // namespace gamma


namespace poch{

  static const double bern[21] = {
    0.0   /* no element 0 */,
    +0.833333333333333333333333333333333e-01,
    -0.138888888888888888888888888888888e-02,
    +0.330687830687830687830687830687830e-04,
    -0.826719576719576719576719576719576e-06,
    +0.208767569878680989792100903212014e-07,
    -0.528419013868749318484768220217955e-09,
    +0.133825365306846788328269809751291e-10,
    -0.338968029632258286683019539124944e-12,
    +0.858606205627784456413590545042562e-14,
    -0.217486869855806187304151642386591e-15,
    +0.550900282836022951520265260890225e-17,
    -0.139544646858125233407076862640635e-18,
    +0.353470703962946747169322997780379e-20,
    -0.895351742703754685040261131811274e-22,
    +0.226795245233768306031095073886816e-23,
    -0.574472439520264523834847971943400e-24,
    +0.145517247561486490186626486727132e-26,
    -0.368599494066531017818178247990866e-28,
    +0.933673425709504467203255515278562e-30,
    -0.236502241570062993455963519636983e-31
  };



  double expm1(double x)
  {
    const double cut = 0.002;

    if(x < F_LOG_DBL_MIN) {
      return -1.0;
    }
    else if(x < -cut) {
      return  exp(x) - 1.0;
    }
    else if(x < cut) {
      return x * (1.0 + 0.5*x*(1.0 + x/3.0*(1.0 + 0.25*x*(1.0 + 0.2*x))));
    }
    else if(x < F_LOG_DBL_MAX) {
      return exp(x) - 1.0;
    }
    else {
      THROW_EXCEPTION( FException, "Oerflow in expm1" );
      return 0;
    }
  }
 


  /* ((a)_x - 1)/x in the "small x" region where
   * cancellation must be controlled.
   *
   * Based on SLATEC DPOCH1().
   */
  /*
     C When ABS(X) is so small that substantial cancellation will occur if
     C the straightforward formula is used, we use an expansion due
     C to Fields and discussed by Y. L. Luke, The Special Functions and Their
     C Approximations, Vol. 1, Academic Press, 1969, page 34.
     C
     C The ratio POCH(A,X) = GAMMA(A+X)/GAMMA(A) is written by Luke as
     C        (A+(X-1)/2)**X * polynomial in (A+(X-1)/2)**(-2) .
     C In order to maintain significance in POCH1, we write for positive a
     C        (A+(X-1)/2)**X = EXP(X*LOG(A+(X-1)/2)) = EXP(Q)
     C                       = 1.0 + Q*EXPREL(Q) .
     C Likewise the polynomial is written
     C        POLY = 1.0 + X*POLY1(A,X) .
     C Thus,
     C        POCH1(A,X) = (POCH(A,X) - 1) / X
     C                   = EXPREL(Q)*(Q/X + Q*POLY1(A,X)) + POLY1(A,X)
     C
     */
  static
    double
    pochrel_smallx(const double a, const double x)
    {
      /*
         SQTBIG = 1.0D0/SQRT(24.0D0*D1MACH(1))
         ALNEPS = LOG(D1MACH(3))
         */
      const double SQTBIG = 1.0/(2.0*M_SQRT2*M_SQRT3*  F_SQRT_DBL_MIN);
      const double ALNEPS = F_LOG_DBL_EPSILON - M_LN2;

      if(x == 0.0) {
        return FMath::psi::psi_x(a);
      }
      else {
        const double bp   = (  (a < -0.5) ? 1.0-a-x : a );
        const int    incr = (int)( (bp < 10.0) ? 11.0-bp : 0 );
        const double b    = bp + incr;
        double dpoch1;
        //gsl_sf_result dexprl;
        //int stat_dexprl;
        double dexprl;
        int i;
        double var    = b + 0.5*(x-1.0);
        double alnvar = log(var);
        double q = x*alnvar;

        double poly1 = 0.0;

        if(var < SQTBIG) {
          const int nterms = (int)(-0.5*ALNEPS/alnvar + 1.0);
          const double var2 = (1.0/var)/var;
          const double rho  = 0.5 * (x + 1.0);
          double term = var2;
          double gbern[24];
          int k, j;

          gbern[1] = 1.0;
          gbern[2] = -rho/12.0;
          poly1 = gbern[2] * term;

          if(nterms > 20) {
            /* NTERMS IS TOO BIG, MAYBE D1MACH(3) IS BAD */
            /* nterms = 20; */
            //result->val = 0.0;
            //result->err = 0.0;
            //F_ERROR ("error", F_ESANITY);
            THROW_EXCEPTION(FException, "sanity error, maybe d1mach(3) is bad" );
          }

          for(k=2; k<=nterms; k++) {
            double gbk = 0.0;
            for(j=1; j<=k; j++) {
              gbk += bern[k-j+1]*gbern[j];
            }
            gbern[k+1] = -rho*gbk/k;

            term  *= (2*k-2-x)*(2*k-1-x)*var2;
            poly1 += gbern[k+1]*term;
          }
        }

        dexprl = expm1(q);
        /*
           stat_dexprl = gsl_sf_expm1_e(q, &dexprl);
           if(stat_dexprl != F_SUCCESS) {
           result->val = 0.0;
           result->err = 0.0;
           return stat_dexprl;
           }
           */
        //    dexprl.val = dexprl.val/q;
        dexprl /= q;
        poly1 *= (x - 1.0);
        dpoch1 = dexprl * (alnvar + q * poly1) + poly1;
        for(i=incr-1; i >= 0; i--) {
          /*
             C WE HAVE DPOCH1(B,X), BUT BP IS SMALL, SO WE USE BACKWARDS RECURSION
             C TO OBTAIN DPOCH1(BP,X).
             */
          double binv = 1.0/(bp+i);
          dpoch1 = (dpoch1 - binv) / (1.0 + x*binv);
        }

        if(bp == a) {
          return dpoch1;
          //result->val = dpoch1;
          //result->err = 2.0 * F_DBL_EPSILON * (fabs(incr) + 1.0) * fabs(result->val);
          //return F_SUCCESS;
        }
        else {
          /*
             C WE HAVE DPOCH1(BP,X), BUT A IS LT -0.5.  WE THEREFORE USE A
             C REFLECTION FORMULA TO OBTAIN DPOCH1(A,X).
             */
          double sinpxx = sin(M_PI*x)/x;
          double sinpx2 = sin(0.5*M_PI*x);
          double t1 = sinpxx/tan(M_PI*b);
          double t2 = 2.0*sinpx2*(sinpx2/x);
          double trig  = t1 - t2;
          //result->val  = dpoch1 * (1.0 + x*trig) + trig;
          //result->err  = (fabs(dpoch1*x) + 1.0) * F_DBL_EPSILON * (fabs(t1) + fabs(t2));
          //result->err += 2.0 * F_DBL_EPSILON * (fabs(incr) + 1.0) * fabs(result->val);
          //return F_SUCCESS;
          return dpoch1 * (1.0+x*trig)+trig;
        }
      }
    }



  // assumes a>0 and a+x > 0
  double lnpoch_pos( const double a, const double x)
  {
    double absx = fabs(x);
    if(absx > 0.1*a || absx*log(std::max(a,2.0)) > 0.1)
    {
      if(a < F_GAMMA_XMAX && a+x < F_GAMMA_XMAX)
      {
        // if it can be done by calculating the gamma functions
        // do it directly because this will be more accurate than
        // computing the substraction of the logs
        double g1 = FMath::gamma::gammainv(a);
        double g2 = FMath::gamma::gammainv(a+x);
        return -log(g2/g1);
      }
      else
      {
        // must do substraction
        double lg1 = FMath::gamma::lngamma(a);
        double lg2 = FMath::gamma::lngamma(a+x);
        return lg2 - lg1;
      }
    }
    else if( absx < 0.1*a && a > 15.0)
    {
      // be careful about the implied subtraction.
      // note that both a+x and a must be large
      // here since a is not small and x is not relatively
      // large. So we calculate using Stirling for Log[gamma(z)]
      // log[gamma(a+x)/gamma(a)] = x(log[a]-1)+(x+a-1/2)log[1+x/a]
      // +(1/(1+eps) -1)   /(12a)
      // -(1/(1+eps))^3 -1)/(360 a^3)
      // + ...
      const double eps = x/a;
      const double den = 1.0 + eps;
      const double d3 = den*den*den;
      const double d5 = d3*den*den;
      const double d7 = d5*den*den;
      const double c1 = -eps/den;
      const double c3 = -eps*(3.0+eps*(3.0+eps))/d3;
      const double c5 = -eps*(5.0+eps*(10.0+eps*(10.0+eps*(5.0+eps))))/d5;
      const double c7 = -eps*(7.0+eps*(21.0+eps*(35.0+eps*(35.0+eps*(21.0+eps*(7.0+eps))))))/d7;
      const double p8 = pow_int(1.0+eps,8);
      const double c8 = 1.0/p8             - 1.0;  /* these need not   */
      const double c9 = 1.0/(p8*(1.0+eps)) - 1.0;  /* be very accurate */
      const double a4 = a*a*a*a;
      const double a6 = a4*a*a;
      const double ser_1 = c1 + c3/(30.0*a*a) + c5/(105.0*a4) + c7/(140.0*a6);
      const double ser_2 = c8/(99.0*a6*a*a) - 691.0/360360.0 * c9/(a6*a4);
      const double ser = (ser_1 + ser_2)/ (12.0*a);

      double term1 = x * log(a/M_E);
      double term2;
      double ln_1peps = func_log_1plusx(eps );  /* log(1 + x/a) */
      term2 = (x + a - 0.5) * ln_1peps;

      return   term1 + term2 + ser;
      //result->err  = F_DBL_EPSILON*fabs(term1);
      //result->err += fabs((x + a - 0.5)*ln_1peps.err);
      //result->err += fabs(ln_1peps.val) * F_DBL_EPSILON * (fabs(x) + fabs(a) + 0.5);
      //result->err += 2.0 * F_DBL_EPSILON * fabs(result->val);
      //return F_SUCCESS;
    }
    else {
      //eassert(false);
      
      //gsl_sf_result poch_rel;
      double poch_rel = pochrel_smallx(a, x);
      double eps = x*poch_rel;
      return func_log_1plusx(eps);
      //int stat_e = gsl_sf_log_1plusx_e(eps, result);
      //result->err  = 2.0 * fabs(x * poch_rel.err / (1.0 + eps));
      //result->err += 2.0 * F_DBL_EPSILON * fabs(result->val);
      //return F_ERROR_SELECT_2(stat_e, stat_p);
    }


  }


  double lnpoch(const double a, const double x)
  {
    if(a<=0.0 || a+x <= 0.0)
      THROW_DEFAULT_EXCEPTION( FDomainError );

    if(x == 0.0 )
      return 1.0;
    else
      return lnpoch_pos(a, x);
  }

} // namespace poch


double legendre_Pmm(const int m, const double x)
{
  if (m==0)
    return 1.0;
  double pmm = 1.0;
  double root_factor = sqrt(1.0-x)*sqrt(1.0+x);
  double fact_coeff = 1.0;
  for(int i=1; i<=m; ++i)
  {
    pmm *= -fact_coeff * root_factor;
    fact_coeff += 2.0;
  }
  return pmm;
}

#ifdef NEW
double legendre(const unsigned int l, const unsigned int m, const double x)
{
  data.clear();
  if(m<0 || l<m || x <-1.0 || x > 1.0) 
    THROW_ECEPTION( FException, "invalid domain" );

  double p_mm = legendre_Pmm(m,x);
  double p_mmp1 = x*(2*m+1)*pmm;
  if(l==m)
  {
    return p_mm;
  }
  else if( l== m+1 )
  {
    retrun p_mmp1;
  }
  else
  {
    double p_ellm2 = p_mm;
    double p_ellm1 = p_mmp1;
    double p_ell = 0.0;
    int ell;

    for(ell=m+2; ell<= l; ++ell)
    {
      p_ell = (x*(2*ell-1)*p_ellm1 -(ell+m-1)*p_ellm2)/(ell-m);
      p_ellm2 = p_ellm1;
      p_ellm1 = p_ell;
    }
    return p_ell;
  } 
}
#else
// static helper
double legendre(const unsigned int l, const unsigned int m, const double x)
{
  double fact, pll(0), pmm, pmmp1, x2;
  unsigned int i,ll;

  //  if ( m< 0|| m>l || fabs(x) > 1.0)
  //      cerr << "bad arguments:" << m << " " << l << " " << x << endl;
  pmm=1.0;
  if( m> 0){
    x2=sqrt((1.0-x)*(1.0+x));
    fact=1.0;
    for(i=1;i<=m;++i)
    {
      pmm *= -fact*x2;
      fact+=2.0;
    }
  }
  if(l==m)
    return pmm;
  else
  {
    pmmp1=x*(2*m+1)*pmm;
    if(l==(m+1))
      return pmmp1;
    else
    {
      for(ll=m+2;ll<=l;++ll)
      {
        pll=(x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
        pmm=pmmp1;
        pmmp1=pll;
      }
      return pll;
    }
  }
}
#endif

// gsl like implementation
void legendrePlmArray(std::vector<double> &results, const int lmax, const int m, const double x)
{
  if(m<0 || lmax < m || x < -1.0 || x > 1.0)
    THROW_EXCEPTION( FException, "domain error" );

  double p_mm = legendre_Pmm(m,x);
  double p_mmp1 = x*(2.0*m+1.0)*p_mm;

  if(lmax ==  m)
  {
    results[0] = p_mm;
    return;
  }
  else if (lmax == m+1)
  {
    results[0] = p_mm;
    results[1] = p_mmp1;
  }
  else
  {
    double p_ellm2 = p_mm;
    double p_ellm1 = p_mmp1;
    double p_ell = 0.0;
    int ell;

    results[0] = p_mm;
    results[1] = p_mmp1;

    for(ell=m+2; ell <= lmax; ++ell)
    {
      p_ell = (x*(2.0*ell-1.0)*p_ellm1 - (ell+m-1)*p_ellm2)/(ell-m);
      p_ellm2 = p_ellm1;
      p_ellm1= p_ell;
      eassert(ell-m < (int)results.size());
      results[ell-m] = p_ell;
    }
  }
  return;
}

void legendrePlArray(std::vector<double> &result, const int lmax, const double x)
{
  if(lmax <0 || x < -1.0 || x > 1.0)
    THROW_EXCEPTION( FException, "domain error" );

  if(lmax == 0)
  {
    result[0] = 1.0;
    return; // ok
  }
  else if(lmax==1)
  {
    result[0] = 1.0;
    result[1] = x;
    return; // ok
  }
  else
  {
    // recurrence l P_l = (2l-1) z P_{l-1} - (l-1) P_{l-2}
    double p_ellm2 = 1.0; // P_0(x)
    double p_ellm1 = x;
    double p_ell = p_ellm1;
    int ell;

    result[0] = 1.0;
    result[1] = x;

    for(ell=2; ell <= lmax; ++ell)
    {
      p_ell = (x*(2*ell-1)*p_ellm1 -(ell-1)*p_ellm2)/ell;
      p_ellm2 = p_ellm1;
      p_ellm1 = p_ell;
      result[ell] = p_ell;
    }
    return; // ok
  }
}

void legendrePlDerivArray(std::vector<double>& result, std::vector<double>& result_deriv, const int lmax, const double x)
{
  legendrePlArray(result, lmax, x);

  if(lmax >= 0) result_deriv[0] = 0.0;
  if(lmax >= 1) result_deriv[1] = 1.0;

  if(true) // result of legendrePlArray test here?
  {
    int ell;
    if(fabs(x-1.0)*(lmax+1.0)*(lmax+1.0) < SQRT_EPSILON)
    {
      // x near 1
      for(ell=2; ell <= lmax; ++ell)
      {
        const double pre = 0.5 * ell*(ell+1.0);
        result_deriv[ell] = pre*(1.0-0.25*(1.0-x)*(ell+2.0)*(ell-1.0));
      }
    }
    else if(fabs(x+1.0)*(lmax+1.0)*(lmax+1.0) < SQRT_EPSILON)
    {
      // x near -1
      for(ell=2; ell<=lmax; ++ell)
      {
        const double sgn=IS_ODD(ell)? 1.0: -1.0;
        const double pre = sgn*0.5*ell*(ell+1.0);
        result_deriv[ell] = pre*(1.0-0.25*(1.0+x)*(ell+2.0)*(ell-1.0));
      }
    }
    else
    {
      const double diff_a = 1.0+x;
      const double diff_b = 1.0-x;

      for(ell=2; ell<=lmax; ++ell)
      {
        result_deriv[ell] = -ell*(x*result[ell]-result[ell-1])/(diff_a*diff_b);
      }
    }
    return; //ok
  }
  // error;
  eassert(false);
}

void legendrePlmDerivArray(std::vector<double>& result, std::vector<double> &result_deriv, const int lmax, const int m, const double x)
{
  if(m<0 || m >lmax)
    THROW_EXCEPTION( FException, "domain error" );

  if(m==0)
  {
    legendrePlDerivArray(result, result_deriv, lmax, x);
    return;
  }
  else
  {
    legendrePlmArray(result, lmax, m, x);

    // if legendre ok
    int ell;
    if(m ==1 && 1.0-fabs(x) < F_DBL_EPSILON)
    {
      /* This divergence is real and comes from the cusp like
       * behavior for m=1, for example P(1,1) = -sqrt(1-x^2)
       */
      THROW_EXCEPTION( FException, "divergence near |x| = 1.0 since m = 1, overflow error" );
    }
    else if(m==2 && (1.0 -fabs(x) < F_DBL_EPSILON))
    {
      if(fabs(x-1.0) < EPSILON)
      {
        for(ell=m; ell <= lmax; ++ell)
          result_deriv[ell-m] = -0.25*x*(ell-1.0)*ell*(ell+1.0)*(ell+2.0);
      }
      else if(fabs(x+1.0) < F_DBL_EPSILON)
      {
        for(ell=m; ell<=lmax; ++ell)
        {
          const double sgn = IS_ODD(ell) ? 1.0 : -1.0;
          result_deriv[ell-m] = -0.25*sgn *x*(ell-1.0)*ell*(ell+1.0)*(ell+2.0);
        }
      }
      return; // oK
    }
    else
    {
      /* m>2 is easier to deal with since the endpoints always vanish */
      if(1.0 - fabs(x) < F_DBL_EPSILON)
      {
        for(ell=m; ell<=lmax; ++ell)
          result_deriv[ell-m] = 0.0;
        return; // ok
      }
      else
      {
        const double diff_a = 1.0+x;
        const double diff_b = 1.0-x;
        result_deriv[0] = -m*x/(diff_a*diff_b)*result[0];
        if(lmax-m>=1) result_deriv[1] = (2.0*m+1.0)*(x*result_deriv[0]+result[0]);
        for(ell=m+2;ell<=lmax; ++ell)
        {
          result_deriv[ell-m]= -(ell*x*result[ell-m]-(ell+m)*result[ell-1-m])/(diff_a*diff_b);
        }
        return; //ok
      }
    }
  }
  // error if we reach here!
  eassert(false);
}



//#if 0 // lnpoch needs to be implemented to use this!

// These functions following from here on are used for spherical harmonics calculation and are already normalized by thone funny factors.

double legendreSHPlm(const int l, int m, const double x)
{
  if(m<0||l<m || x< -1.0 || x > 1.0)
    THROW_DEFAULT_EXCEPTION(FDomainError    );

  if(m==0)
  {
    double result =  legendre_Pl(l,x);
    double pre = sqrt((2.0*l+1.0)/(4.0*M_PI));
    return pre*result;
  }
  else if( x== 1.0 || x == -1.0)
  {
    return 0.0; // ok
  }
  else
  {
    // m> 0 and |x| < 1
    //
    // starting recursion with
    // Y_m^m(x) = sqrt((2m+1)/(4pi m) gamma(m+1/2)/gamma(m))(-1)^m (1-x^2)^(m/2) /pi^(1/4)
    //
    double lncirc;
    double lnpoch_;
    double lnpre_val;
    double sr;
    const double sgn = IS_ODD(m) ? -1.0: 1.0;
    const double y_mmp1_factor = x*sqrt(2.0*m+3.0);
    double y_mm; // y_mm_err;
    double y_mmp1;
    lncirc = func_log_1plusx(-x*x);
    lnpoch_ = FMath::poch::lnpoch(m, 0.5);
    lnpre_val = -0.25*M_LNPI +0.5 *(lnpoch_+m*lncirc);
    double ex_pre = exp(lnpre_val);
    sr = sqrt((2.0+1.0/m)/(4.0*M_PI));
    y_mm = sgn*sr*ex_pre;
    y_mmp1 = y_mmp1_factor*y_mm;

    if(l==m)
    {
      return y_mm;
    }
    else if( l== m+1)
    {
      return y_mmp1;
    }
    else
    {
      double y_ell = 0.0;
      int ell;
      for(ell=m+2; ell <=l; ++ell)
      {
        const double rat1 = (double)(ell-m)/(double)(ell+m);
        const double rat2 = (ell-m-1.0)/(ell+m-1.0);
        const double factor1 = sqrt(rat1*(2*ell+1)*(2*ell-1));
        const double factor2 = sqrt(rat1*rat2*(2*ell+1)/(2*ell-3));
        y_ell = (x*y_mmp1*factor1 -(ell+m-1)*y_mm*factor2)/(ell-m);
        y_mm = y_mmp1;
        y_mmp1 = y_ell;
      }

      return y_ell; //ok
    }
  }
}

void legendreSHPlmArray(std::vector<double>& result, const int lmax, int m, const double x)
{
  if(m<0 || lmax <m || x< -1.0 || x > 1.0)
    THROW_DEFAULT_EXCEPTION( FDomainError     );

  if( m> 0 && (x == 1.0 || x == -1.0))
  {
    int ell;
    for(ell=m; ell<=lmax; ++ell) result[ell-m] = 0.0;
    return; // ok
  }
  else
  {
    double y_mm;
    double y_mmp1;

    if( m== 0)
    {
      y_mm = 0.5/M_SQRTPI;
      y_mmp1 = x*M_SQRT3*y_mm;
    }
    else
    {
      // |x| < 1
      double lncirc;
      double lnpoch;
      double lnpre;
      const double sgn = IS_ODD(m) ? -1.0: 1.0;
      lncirc = func_log_1plusx(-x*x);
      lnpoch = FMath::poch::lnpoch(m, 0.5);
      lnpre= -0.25*M_LNPI+0.5*(lnpoch+m*lncirc);
      y_mm = sqrt((2.0+1.0/m)/(4.0*M_PI))*sgn*exp(lnpre);
      y_mmp1= x*sqrt(2.0*m+3.0)*y_mm;
    }

    if(lmax == m)
    {
      eassert(result.size() >= 1);
      result[0] = y_mm;
      return; // ok;
    }
    else if(lmax == m+1)
    {
      eassert(result.size() >= 2);
      result[0] = y_mm;
      result[1] = y_mmp1;
      return; // ok
    }
    else
    {
      double y_ell;
      int ell;
      result[0] = y_mm;
      result[1] = y_mmp1;

      // compute Y_l^m, l>m+1 upward recurson on l
      for(ell=m+2; ell<=lmax; ++ell)
      {
        const double rat1 = (double)(ell-m)/(double)(ell+m);
        const double rat2 = (ell-m-1.0)/(ell+m-1.0);
        const double factor1 = sqrt(rat1*(2*ell+1)*(2*ell-1));
        const double factor2 = sqrt(rat1*rat2*(2*ell+1)/(2*ell-3));
        y_ell=(x*y_mmp1*factor1 - (ell+m-1)*y_mm*factor2)/(ell-m);
        y_mm = y_mmp1;
        y_mmp1=y_ell;
        eassert(ell-m < (int)result.size());
        result[ell-m] = y_ell;
      }
    }
    return; // ok
  }
}

void legendreSHPlmDerivArray(std::vector<double> &result, std::vector<double>& result_deriv, const int lmax, const int m, const double x)
{
  if(m<0 || lmax <m || x< -1.0 || x > 1.0)
    THROW_DEFAULT_EXCEPTION( FDomainError     );
  if(m==0)
  {
    legendrePlDerivArray(result, result_deriv, lmax, x);
    int ell;
    for(ell=0; ell<=lmax; ++ell)
    {
      const double prefactor = sqrt((2.0*ell+1.0)/(4.0*M_PI));
      result[ell] *= prefactor;
      result_deriv[ell] *= prefactor;
    }
    return; //ok
  }
  else if(m==1)
  {
    legendrePlmDerivArray(result, result_deriv, lmax, m, x);
    int ell;
    for(ell=1; ell<=lmax; ++ell)
    {
      const double prefactor = sqrt((2.0*ell+1.0)/(ell+1.0)/(4.0*M_PI*ell));
      result[ell-1] *= prefactor;
      result_deriv[ell-1] *= prefactor;
    }
    return; //ok
  }
  else
  {
    legendreSHPlmArray(result, lmax, m, x);
    if(true) // check result
    {
      int ell;
      if(1.0-fabs(x) < EPSILON)
      {
        for(ell=m; ell<=lmax; ell++) result_deriv[ell-m] = 0.0;
        return; // ok
      }
      else
      {
        const double diff_a = 1.0+x;
        const double diff_b = 1.0-x;
        result_deriv[0] = -m*x/(diff_a*diff_b)*result[0];
        if(lmax-m >=1) result_deriv[1] = sqrt(2.0*m+3.0)*(x*result_deriv[0]+result[0]);
        for(ell=m+2; ell<=lmax; ++ell)
        {
          const double c1 = sqrt(((2.0*ell+1.0)/(2.0*ell-1.0))*((double)(ell-m)/(double)(ell+m)));
          result_deriv[ell-m] = -(ell*x*result[ell-m]-c1*(ell+m)*result[ell-1-m])/(diff_a*diff_b);
        }
        return; // ok
      }
    }
    eassert(false);
  }
}


double legendre_Plm(const int l, const int m, const double x )
{
  // if m and l are both large, overflows may occur.
  const double dif = l-m;
  const double sum = l+m;
  const double t_d = ( dif == 0.0 ? 0.0 : 0.5*dif *(log(dif)-1.0));
  const double t_s = ( dif == 0.0 ? 0.0 : 0.5*sum *(log(sum)-1.0));

  const double exp_check = 0.5 * log(2.0*l+1.0) + t_d - t_s;

  if(m <0|| l<m || x<-1.0|| x>1.0)
  {
    THROW_DEFAULT_EXCEPTION( FDomainError );
  }
  else if( exp_check < F_LOG_DBL_MIN +10.0)
  {
    THROW_DEFAULT_EXCEPTION( FOverflowError );
  }

  double p_mm = legendre_Pmm(m, x);
  double p_mmp1 = x*(2*m+1)*p_mm;

  if(l==m)
  {
    return p_mm;
  }
  else if(l==m+1)
  {
    return p_mmp1;
  }
  else
  {
    // recurrence (l-m) P(l,m) = (2l-1) z P(l-1,m) -(l+m-1) P(l-2,m);
    // start with P(m,m), P(m+1,m)
    double pellm2 = p_mm;
    double pellm1 = p_mmp1;

    double pell = 0.0;
    int ell;
    for( ell=m+2; ell <= l; ++ell)
    {
      pell = (x*(2*ell-1)*pellm1 - (ell+m-1)*pellm2)/(ell-m);
      pellm2=pellm1;
      pellm1=pell;
    }
    return pell;
  }
}

namespace bessel
{

  namespace{ // anonymous

/* chebyshev expansions for amplitude and phase
   functions used in bessel evaluations

   These are the same for J0,Y0 and for J1,Y1, so
   they sit outside those functions.
*/

static double bm0_data[21] = {
   0.09284961637381644,
  -0.00142987707403484,
   0.00002830579271257,
  -0.00000143300611424,
   0.00000012028628046,
  -0.00000001397113013,
   0.00000000204076188,
  -0.00000000035399669,
   0.00000000007024759,
  -0.00000000001554107,
   0.00000000000376226,
  -0.00000000000098282,
   0.00000000000027408,
  -0.00000000000008091,
   0.00000000000002511,
  -0.00000000000000814,
   0.00000000000000275,
  -0.00000000000000096,
   0.00000000000000034,
  -0.00000000000000012,
   0.00000000000000004
};
const cheb_series bessel_amp_phase_bm0_cs = {
  bm0_data,
  20,
  -1, 1,
  10
};

static double bth0_data[24] = {
  -0.24639163774300119,
   0.001737098307508963,
  -0.000062183633402968,
   0.000004368050165742,
  -0.000000456093019869,
   0.000000062197400101,
  -0.000000010300442889,
   0.000000001979526776,
  -0.000000000428198396,
   0.000000000102035840,
  -0.000000000026363898,
   0.000000000007297935,
  -0.000000000002144188,
   0.000000000000663693,
  -0.000000000000215126,
   0.000000000000072659,
  -0.000000000000025465,
   0.000000000000009229,
  -0.000000000000003448,
   0.000000000000001325,
  -0.000000000000000522,
   0.000000000000000210,
  -0.000000000000000087,
   0.000000000000000036
};
const cheb_series bessel_amp_phase_bth0_cs = {
  bth0_data,
  23,
  -1, 1,
  12
};


static double bm1_data[21] = {
   0.1047362510931285,
   0.00442443893702345,
  -0.00005661639504035,
   0.00000231349417339,
  -0.00000017377182007,
   0.00000001893209930,
  -0.00000000265416023,
   0.00000000044740209,
  -0.00000000008691795,
   0.00000000001891492,
  -0.00000000000451884,
   0.00000000000116765,
  -0.00000000000032265,
   0.00000000000009450,
  -0.00000000000002913,
   0.00000000000000939,
  -0.00000000000000315,
   0.00000000000000109,
  -0.00000000000000039,
   0.00000000000000014,
  -0.00000000000000005,
};
const cheb_series bessel_amp_phase_bm1_cs = {
  bm1_data,
  20,
  -1, 1,
  10
};



static double bth1_data[24] = {
   0.74060141026313850,
  -0.004571755659637690,
   0.000119818510964326,
  -0.000006964561891648,
   0.000000655495621447,
  -0.000000084066228945,
   0.000000013376886564,
  -0.000000002499565654,
   0.000000000529495100,
  -0.000000000124135944,
   0.000000000031656485,
  -0.000000000008668640,
   0.000000000002523758,
  -0.000000000000775085,
   0.000000000000249527,
  -0.000000000000083773,
   0.000000000000029205,
  -0.000000000000010534,
   0.000000000000003919,
  -0.000000000000001500,
   0.000000000000000589,
  -0.000000000000000237,
   0.000000000000000097,
  -0.000000000000000040,
};
const cheb_series bessel_amp_phase_bth1_cs = {
  bth1_data,
  23,
  -1, 1,
  12
};



  }
  
  // based on SLATEC besj0, 1977 version, w. fullerton
  // chebyshev expanasions for bessel functions
  // series for besselJ0 on interval 0.0 to 1.6e+01
  // weighted error: 7.47e-18
  // log weighted error: 17.13
  // significant figures required 16.98
  // decimal places required 17.68
  namespace { // anonymous
#ifdef NOTNR
    static double bj0_data[13] =
    {   0.100254161968939137,
      -0.665223007764405132,
      0.248983703498281314,
      -0.0332527231700357697,
      0.0023114179304694015,
      -0.0000991127741995080,
      0.0000028916708643998,
      -0.0000000612108586630,
      0.0000000009838650793,
      -0.0000000000124235515,
      0.0000000000001265433,
      -0.0000000000000010619,
      0.0000000000000000074
    };

    static cheb_series bj0_cs = {
      bj0_data,
      12,
      -1, 1,
      9
    };
#endif




    /**
     * used to calculate the oscillating bessel function
     *  cos(y - pi/4 + eps )
     */
    double bessel_cos_pi4(double y, double eps)
    {
      const double sy = sin(y);
      const double cy = cos(y);
      const double s = sy + cy;
      const double d = sy - cy;
      double seps;
      double ceps;
      if(fabs(eps) < F_ROOT5_DBL_EPSILON) {
        const double e2 = eps*eps;
        seps = eps * (1.0 - e2/6.0 * (1.0 - e2/20.0));
        ceps = 1.0 - e2/2.0 * (1.0 - e2/12.0);
      }
      else {
        seps = sin(eps);
        ceps = cos(eps);
      }
      return (ceps * s - seps * d)/ M_SQRT2;
    }

    /**
     * used to calculate the oscillating bessel function
     *  sin(y - pi/4 + eps )
     */
    double gsl_sf_bessel_sin_pi4_e(double y, double eps)
    { 
      const double sy = sin(y);
      const double cy = cos(y);
      const double s = sy + cy;
      const double d = sy - cy;
      double seps; 
      double ceps; 
      if(fabs(eps) < F_ROOT5_DBL_EPSILON) {
        const double e2 = eps*eps;
        seps = eps * (1.0 - e2/6.0 * (1.0 - e2/20.0));
        ceps = 1.0 - e2/2.0 * (1.0 - e2/12.0);
      }
      else { 
        seps = sin(eps);
        ceps = cos(eps);
      }
      return (ceps * d + seps * s)/ M_SQRT2;
    }

  } // bessel::anonymous namespace

  double bessel_I0(const double x)
  {
    // NR p 237 code bessi0
    // the results become ugly for larger values of x
    // (100 and up seem to become ugly)
    double ax, ans;
    double y;

    if((ax=fabs(x)) < 3.75)
    {
      y=x/3.75;
      y*=y;
      ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
              +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
    }
    else
    {
      y=3.75/ax;
      ans=(exp(ax)/sqrt(ax))*(0.3984228+y*(0.1328592e-1
            +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
            +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
                  +y*0.392377e-2))))))));
    }
    return ans;
  }

double bessel_I1(const double x)
{
  double ax, ans;
  double y;

  if((ax=fabs(x)) < 3.75)
  {
    y=x/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
              +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3))))));
  }
  else
  {
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
          -y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
          +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
    ans*=(exp(ax)/sqrt(ax));
  }
  return x<0.0 ? -ans : ans;
}


  double bessel_J0(const double x)
  {
#ifdef NOTNR
    /*
    double y = std::abs(x);
    if(y<2.0*F_SQRT_DBL_EPSILON)
    {
      return 1.0;
    }
    else if( y <= 4.0 )
    {
      //THROW_DEFAULT_EXCEPTION(FNotImplementedException);
      return cheb_eval( &bj0_cs, 0.125*y*y-1.0);
    }
    else
    {
      // gsl code
      THROW_DEFAULT_EXCEPTION(FNotImplementedException);
      const double z = 32.0/(y*y)-1.0;
      double ca = cheb_eval( &bessel_amp_phase_bm0_cs, z);
      double ct = cheb_eval( &bessel_amp_phase_bth0_cs, z);
      double cp = bessel_cos_pi4( y, ct/y);
      const double sqrty = sqrt(y);
      const double ampl =(0.75+ca)/sqrty;
      return ampl*cp;
      */
    }
#endif
      // function from numerical receips page 232
      float ax, z;
      double xx, y, ans, ans1, ans2;
      if((ax=fabs(x)) < 8.0)
      {
        y=x*x;
        ans1=57568490574.0+y*(-13362590354.0+y*(651619640.7
              +y*(-11214424.18+y*(77392.33017+y*(-184.9052456)))));
        ans2=57568490411.0+y*(1029532985.0+y*(9494680.718
              +y*(59272.64853+y*(267.8532712+y*1.0))));
        ans=ans1/ans2;
      }
      else
      {
        z=8.0/ax;
        y=z*z;
        xx=ax-0.785398164;
        ans1=1.0+y*(-0.1098628627e-2+y*(0.2734510407e-4
              +y*(-0.2073370639e-5+y*0.2093887211e-6)));
        ans2=-0.1562499995e-1+y*(0.1430488765e-3
            +y*(-0.6911147651e-5+y*(0.7621095161e-6
            -y*0.934945152e-7)));
        ans=sqrt(0.636619772/ax)*(cos(xx)*ans1-z*sin(xx)*ans2);
      }
      return ans;
 
  }
} // bessel functions namespace

} // namespace FMath
