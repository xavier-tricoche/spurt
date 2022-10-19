#ifndef FMathValues_HH
#define FMathValues_HH

#include <math.h>

#ifndef M_E
#define M_E        2.71828182845904523536028747135      /* e */
#endif

#ifndef M_LOG2E
#define M_LOG2E    1.44269504088896340735992468100      /* log_2 (e) */
#endif

#ifndef M_LOG10E
#define M_LOG10E   0.43429448190325182765112891892      /* log_10 (e) */
#endif

#ifndef M_SQRT2
#define M_SQRT2    1.41421356237309504880168872421      /* sqrt(2) */
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2  0.70710678118654752440084436210      /* sqrt(1/2) */
#endif


#ifndef M_SQRT3
#define M_SQRT3    1.73205080756887729352744634151      /* sqrt(3) */
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif

#ifndef M_PI_2
#define M_PI_2     1.57079632679489661923132169164      /* pi/2 */
#endif

#ifndef M_PI_4
#define M_PI_4     0.78539816339744830966156608458      /* pi/4 */
#endif

#ifndef M_SQRTPI
#define M_SQRTPI   1.77245385090551602729816748334      /* sqrt(pi) */
#endif

#ifndef M_2_SQRTPI
#define M_2_SQRTPI 1.12837916709551257389615890312      /* 2/sqrt(pi) */
#endif

#ifndef M_1_PI
#define M_1_PI     0.31830988618379067153776752675      /* 1/pi */
#endif

#ifndef M_2_PI
#define M_2_PI     0.63661977236758134307553505349      /* 2/pi */
#endif

#ifndef M_LN10
#define M_LN10     2.30258509299404568401799145468      /* ln(10) */
#endif

#ifndef M_LN2
#define M_LN2      0.69314718055994530941723212146      /* ln(2) */
#endif

#ifndef M_LNPI
#define M_LNPI     1.14472988584940017414342735135      /* ln(pi) */
#endif

#ifndef M_EULER
#define M_EULER    0.57721566490153286060651209008      /* Euler constant */
#endif

/* other needlessly compulsive abstractions */

#define F_IS_ODD(n)  ((n) & 1)
#define F_IS_EVEN(n) (!(F_IS_ODD(n)))
#define F_SIGN(x)    ((x) >= 0.0 ? 1 : -1)
#define F_SIGN2(x,y)  ((x) >= 0.0 ? y : -y)

/* Return nonzero if x is a real number, i.e. non NaN or infinite. */
#define F_IS_REAL(x) (gsl_finite(x))

/* Define MAX and MIN macros/functions if they don't exist. */

/* plain old macros for general use */
#define F_MAX(a,b) ((a) > (b) ? (a) : (b))
#define F_MIN(a,b) ((a) < (b) ? (a) : (b))

/* function versions of the above, in case they are needed */
double f_max (double a, double b);
double f_min (double a, double b);

/* inline-friendly strongly typed versions */

extern inline int F_MAX_INT (int a, int b);
extern inline int F_MIN_INT (int a, int b);
extern inline double F_MAX_DBL (double a, double b);
extern inline double F_MIN_DBL (double a, double b);
extern inline long double F_MAX_LDBL (long double a, long double b);
extern inline long double F_MIN_LDBL (long double a, long double b);

extern inline int
F_MAX_INT (int a, int b)
{
  return F_MAX (a, b);
}

extern inline int
F_MIN_INT (int a, int b)
{
  return F_MIN (a, b);
}

extern inline double
F_MAX_DBL (double a, double b)
{
  return F_MAX (a, b);
}

extern inline double
F_MIN_DBL (double a, double b)
{
  return F_MIN (a, b);
}

extern inline long double
F_MAX_LDBL (long double a, long double b)
{
  return F_MAX (a, b);
}

extern inline long double
F_MIN_LDBL (long double a, long double b)
{
  return F_MIN (a, b);
}


// reasonable values for epsilons in most cases
//

#define F_DBL_EPSILON        2.2204460492503131e-16
#define F_SQRT_DBL_EPSILON   1.4901161193847656e-08
#define F_ROOT3_DBL_EPSILON  6.0554544523933429e-06
#define F_ROOT4_DBL_EPSILON  1.2207031250000000e-04
#define F_ROOT5_DBL_EPSILON  7.4009597974140505e-04
#define F_ROOT6_DBL_EPSILON  2.4607833005759251e-03
#define F_LOG_DBL_EPSILON   (-3.6043653389117154e+01)

#define F_DBL_MIN        2.2250738585072014e-308
#define F_SQRT_DBL_MIN   1.4916681462400413e-154
#define F_ROOT3_DBL_MIN  2.8126442852362996e-103
#define F_ROOT4_DBL_MIN  1.2213386697554620e-77
#define F_ROOT5_DBL_MIN  2.9476022969691763e-62
#define F_ROOT6_DBL_MIN  5.3034368905798218e-52
#define F_LOG_DBL_MIN   (-7.0839641853226408e+02)

#define F_DBL_MAX        1.7976931348623157e+308
#define F_SQRT_DBL_MAX   1.3407807929942596e+154
#define F_ROOT3_DBL_MAX  5.6438030941222897e+102
#define F_ROOT4_DBL_MAX  1.1579208923731620e+77
#define F_ROOT5_DBL_MAX  4.4765466227572707e+61
#define F_ROOT6_DBL_MAX  2.3756689782295612e+51
#define F_LOG_DBL_MAX    7.0978271289338397e+02

#define F_FLT_EPSILON        1.1920928955078125e-07
#define F_SQRT_FLT_EPSILON   3.4526698300124393e-04
#define F_ROOT3_FLT_EPSILON  4.9215666011518501e-03
#define F_ROOT4_FLT_EPSILON  1.8581361171917516e-02
#define F_ROOT5_FLT_EPSILON  4.1234622211652937e-02
#define F_ROOT6_FLT_EPSILON  7.0153878019335827e-02
#define F_LOG_FLT_EPSILON   (-1.5942385152878742e+01)

#define F_FLT_MIN        1.1754943508222875e-38
#define F_SQRT_FLT_MIN   1.0842021724855044e-19
#define F_ROOT3_FLT_MIN  2.2737367544323241e-13
#define F_ROOT4_FLT_MIN  3.2927225399135965e-10
#define F_ROOT5_FLT_MIN  2.5944428542140822e-08
#define F_ROOT6_FLT_MIN  4.7683715820312542e-07
#define F_LOG_FLT_MIN   (-8.7336544750553102e+01)

#define F_FLT_MAX        3.4028234663852886e+38
#define F_SQRT_FLT_MAX   1.8446743523953730e+19
#define F_ROOT3_FLT_MAX  6.9814635196223242e+12
#define F_ROOT4_FLT_MAX  4.2949672319999986e+09
#define F_ROOT5_FLT_MAX  5.0859007855960041e+07
#define F_ROOT6_FLT_MAX  2.6422459233807749e+06
#define F_LOG_FLT_MAX    8.8722839052068352e+01

#define F_SFLT_EPSILON        4.8828125000000000e-04
#define F_SQRT_SFLT_EPSILON   2.2097086912079612e-02
#define F_ROOT3_SFLT_EPSILON  7.8745065618429588e-02
#define F_ROOT4_SFLT_EPSILON  1.4865088937534013e-01
#define F_ROOT5_SFLT_EPSILON  2.1763764082403100e-01
#define F_ROOT6_SFLT_EPSILON  2.8061551207734325e-01
#define F_LOG_SFLT_EPSILON   (-7.6246189861593985e+00)



#endif /* __F_MATH_H__ */
