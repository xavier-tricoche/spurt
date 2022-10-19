#ifndef FWignerSymbols_hh
#define FWignerSymbols_hh

#include <stdlib.h>
#include <cmath>

#include "eassert.hh"
#include "FException.hh"
#include <iostream>
namespace FMath
{
namespace FWignerSymbols
{

  /** log(n!) */
  double lnfact(int n);
 
  /** checks for triangular condition */
  bool isTriangle( const int ja, const int jb, const int jc );
 
  /** Wigner Rotation */
  double wignerRot( const int two_j, const int two_m1, const int two_m2, const double theta );

  /** 
   * natural logarithm of the delta function (Zare Eq. A2)
   */
  double lndelta( const int two_j1, const int two_j2, const int two_j3);

  /** Wigner 3j-Symbol */
  double ThreeJ( const int j1, const int j2, const int j3, const int m1, const int m2, const int m3);

  /** Wigner 6j-Symbol */
  double SixJ(const int j1, const int j2, const int j3, const int j4, const int j5, const int j6 );

  /** Wigner 9j-Symbol */
  double NineJ(const int j1, const int j2, const int j3,
      const int j4, const int j5, const int j6,
      const int j7, const int j8, const int j9);
  
} // namespace FWignerSymbols
} // namespace FMath

#endif

