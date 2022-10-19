#ifndef FFRT_hh
#define FFRT_hh

#include "FArray.hh"
#include "FSphericalHarmonics.hh"
#include "FSpecialFunctions.hh"

/**
 * Compute the Funk-Radon-Transform in terms of Spherical Harmonics
 * as described by Maxime Descoteaux (INRIA)
 */
class FFunkRadonTransform
{
  public:
  FFunkRadonTransform(int order)
    : transform_(order+1)
    {
      for(unsigned int l=0; l<transform_.size(); ++l)
      {
        transform_[l] = 2.*M_PI*FMath::legendre_Pl( l, 0)/FMath::legendre_Pl(l,1);
      }
    }

  template<class SH>
  void transform(SH &nsh, const SH& sh)
  {
//    nsh.resizeBands(sh.bands());

    try{
    nsh.resize( sh.getOrder());

    typename SH::iterator it;
    while(it.getOrder() <= sh.getOrder())
    {
      nsh(it.getOrder(), it.getMode()) = transform_(it.getOrder())*sh( it.getOrder(), it.getMode());
      it++;
    }
#if 0 
    for(int l=0; l< sh.getOrder(); l+=2)
    {
      for(int m=-l; m<=l; ++m)
        nsh(l,m) = transform_(l)*sh(l,m);
    }
#endif
    }CATCH_N_RETHROW(FException);
  }

  template<class SH>
  FMatrix getMatrix(const SH& sh) const
  {
    FMatrix m( sh.sz(sh.getOrder()), sh.sz(sh.getOrder()));
    int i=0;
    typename SH::iterator it;
    while ( it.getOrder()<= sh.getOrder() )
    {
      m( i,i ) = transform_( it.getOrder() );
      it++; ++i;
    }
    return m;
  }

  private:
  FArray transform_;
};

#endif
