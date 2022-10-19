#ifndef FTemplateSphericalHarmonics_HH
#define FTemplateSphericalHarmonics_HH

#include "stdAliases.hh"
#include "FArray.hh"
#include "FTensor.hh"

#include <complex>
#include <iostream>

#include "FMatrix.hh"

//#include <iostream>
#include <iosfwd>
//#include <algorithm>
#include <iterator>
#include "eassert.hh"
#include "FMathValues.hh"

using namespace std;

template<class T> inline bool isReal(const T&){ return false; }
inline bool isReal(const double&){ return true; }

/** 
 * A spherical harmonic base class.
 * 
 * Data is kept here, further information is stored in struct SHTraits.
 */
template<class SHTraits>
class FSHBase : public SHTraits
{
#ifndef NODEBUG
#define CHECK
#else
#define CHECK \
  eassert( (this->order==-1 && this->data.size()==0) || (SHTraits::sz(this->order) == this->data.size()) );
#endif
 
  public:
    typedef FSHBase<SHTraits> this_type;
    typedef typename SHTraits::value_type value_type;// type of the SH coefficients
    typedef typename SHTraits::func_type  func_type; // type of the surface function
    typedef SHTraits traits;
    
    FSHBase(){order = -1;}
    FSHBase( const FSHBase& rhs ) : data( rhs.data ), order(rhs.order) {CHECK;}
    FSHBase( const int order ) : data( SHTraits::sz(order) ), order(order) {CHECK;}
    FSHBase( const FArray& data ) : data( data ), order( orderFromSize( data.size() ) ){CHECK;}

    static int orderFromSize( int size )
    {
      int i=0;
      while ( true )
      {
        if ( SHTraits::sz( i ) == size)
        {
          return i;
        }
        ++i;
      }
    }
    
    int getOrder() const { CHECK; return order; }
    unsigned int nbComponents() const { CHECK; return data.size(); }
    unsigned int memSize()  const { CHECK; return sizeof(this) + data.size()*sizeof(value_type); }

    void resize(unsigned int order) { CHECK; this->order=order; data.resize( SHTraits::sz(order), value_type(0)); CHECK; }
    void clear() { for(unsigned int i=0; i< data.size(); ++i) data[i]=value_type(0); }
    
    value_type  operator()(const unsigned int l, const int m) const { CHECK; return data[ SHTraits::idx(l,m) ]; }
    value_type& operator()(const unsigned int l, const int m)       { CHECK; return data[ SHTraits::idx(l,m) ]; }

    func_type evaluate(const double theta, const double phi) const {CHECK; return SHTraits::evaluate( *this, theta, phi );}

    // not yet implemented for all types
    func_type derivate_theta(const double theta, const double phi) const {CHECK; return SHTraits::derivate_theta( *this, theta, phi );}
    func_type derivate_phi(const double theta, const double phi) const {CHECK; return SHTraits::derivate_phi( *this, theta, phi );}

    template< class SHTraits_ >
    friend std::ostream& operator<<(std::ostream& os, const FSHBase< SHTraits_ >& sh);

    template< class OTHER >
    void set( const OTHER &sh )
    {
      copySH( *this, sh );
    }
    
    //template<class A, class B> friend void copyFromSH(A&, const B&);


    /**
     * scale spherical harmonic by a constant factor
     */
    this_type& operator*=(const double scale)
    {
      typename std::vector<value_type>::iterator it;
      it = data.begin();
      for(;it != data.end(); ++it)
      {
        (*it) *= scale;
      }
      return *this;
    }

    /**
     * add two spherical harmonic, increase size of *this if
     * rhs is larger than we are
     */
    this_type& operator+=(const this_type& rhs)
    {
      if(rhs.getOrder() > this->getOrder())
      {
        std::vector<value_type> data2;
        data.swap(data2);
        this.resize( rhs.getOrder() );
         for(typename std::vector<value_type>::iterator it = data.begin(), i1 = data2.begin(), i2=rhs.data.begin();
            i2 != rhs.data.end(); ++it,++i1, ++i2)
        {
          (*it) = *i1 + *i2;
        }
      }
      else
      {
        eassert( rhs.getOrder() <= this.getOrder() );
        for(typename std::vector<value_type>::iterator it = data.begin(), i2=rhs.data.begin();
            i2 != rhs.data.end(); ++it,++i2)
        {
          (*it) += *i2;
        }
      }
      return *this;
    }
      
   void setFromArray( const FArray& a )
   {
     eassert( a.size() == data.size());
     for(unsigned int i=0; i< a.size(); ++i)
     {
       data[i] = a[i];
     }
   }

   void getAsArray( FArray& a ) const
   {
     eassert( isReal( value_type() ) );
     a.resize( data.size() );
     for(unsigned int i=0; i< data.size(); ++i)
     {
       a[i] = data[i];
     }
   }
    
  private:
    std::vector<value_type> data;
    int order;
};


template< class SHTraits>
std::ostream& operator<<(std::ostream& os, const FSHBase< SHTraits> & sh)
{
  os << "[ ";
  std::copy( 
      sh.data.begin(), 
      sh.data.end(), 
      std::ostream_iterator< typename FSHBase< SHTraits >::value_type >( os, " " )
      );
  os << "]";
  return os;
}


struct FRealSymmetricSHTraits
{
  typedef double value_type;
  typedef double func_type;
  static unsigned int idx(unsigned int l, int m) { eassert(l%2==0); eassert( -(int)l <=m && m <=(int)l); return (l*l+l)/2+m; }
  static unsigned int sz(unsigned int order)     { eassert(order%2==0); return idx(order, order)+1; }
  static value_type evaluate(const FSHBase<FRealSymmetricSHTraits>& sh, const double theta, const double phi);

  struct iterator
  {
    iterator(): l(0), m(0)
    {}

   iterator& operator++(int)
    {
      if(m==l)
      {
        l+=2; m=-l;
      }
      else
        ++m;
      return *this;
    }

   int getOrder() const { return l; }
   int getMode()  const { return m; }

    private:
   int l;
   int m;
  };
};

/**
 * \todo test!
 * Real valued spherical harmonics
 */
struct FRealSHTraits
{
  typedef double value_type;
  typedef double func_type;
  static unsigned int idx(unsigned int l, int m) { eassert( -(int)l<=m && m <=(int)l ); return l*l+m+l; }
  static unsigned int sz(unsigned int order)     { return FRealSHTraits::idx(order+1, -order-1); }
  static value_type evaluate(const FSHBase<FRealSHTraits>& sh, const double theta, const double phi);

  struct iterator
  {
    iterator(): l(0), m(0)
    {}

   iterator& operator++(int)
    {
      if(m==l)
      {
        ++l; m=-l;
      }
      else
        ++m;
      return *this;
    }

   unsigned int getIndex() const { return FRealSHTraits::idx(l,m); }

   int getOrder() const { return l; }
   int getMode()  const { return m; }
    private:
   int l;
   int m;
  };

};


/**
 * Spherical harmonics as defined by Maxime Descoteaux (INRIA) 
 */
struct FDescoteauxSHTraits
{
  typedef FDescoteauxSHTraits self_type;
  typedef FSHBase<self_type> base_type;
  typedef double value_type;
  typedef double func_type;
  static unsigned int idx(unsigned int l, int m) { eassert( l % 2 == 0); eassert( -(int)l<=m && m <=(int)l ); return (l*l+l+2)/2+m-1; }
  static unsigned int sz(unsigned int order)     { if ( order %2 == 1 ) return 0; return FDescoteauxSHTraits::idx(order, order)+1; }
  static value_type evaluate(const base_type& sh, const double theta, const double phi);

  struct iterator 
  {
    iterator() : l(0), m(0)
    {}

   iterator& operator++(int)
    {
      if(m==l)
      {
        l+=2; m=-l;
      }
      else
        ++m;
      return *this;
    }

   unsigned int getIndex() const { return FDescoteauxSHTraits::idx(l,m); }

   int getOrder() const { return l; }
   int getMode()  const { return m; }

    private:
   int l;
   int m;
  };

  static value_type derivate_theta(const base_type& base, double theta, double phi);
  static value_type derivate_phi(const base_type& base, const double theta, const double phi);
};

/** comples spherical harmonics as defined on MathWorld */
struct FComplexSHTraits
{
  typedef std::complex<double> value_type;
  typedef std::complex<double> func_type;

  static unsigned int idx( unsigned int l, int m) { eassert( -(int)l<=m && m<=(int)l); return l*l+m+l; }
  static unsigned int sz( unsigned int order ) { return FComplexSHTraits::idx(order+1, -order-1); }
  static value_type evaluate( const FSHBase<FComplexSHTraits>& sh, const double theta, const double phi );


  struct iterator 
  {
    iterator() : l(0), m(0)
    {}

   iterator& operator++(int)
    {
      if(m==l)
      {
        l+=1; m=-l;
      }
      else
        ++m;
      return *this;
    }

   unsigned int getIndex() const { return FComplexSHTraits::idx(l,m); }

   int getOrder() const { return l; }
   int getMode()  const { return m; }
    private:
   int l;
   int m;
  };


};

//! some typedefs for cleaner notation
typedef FSHBase<FRealSymmetricSHTraits>  FRealSymmetricSphericalHarmonic;
typedef FSHBase<FRealSHTraits>           FRealSphericalHarmonic;
typedef FSHBase<FDescoteauxSHTraits>     FDescoteauxSphericalHarmonic;
typedef FSHBase<FComplexSHTraits>        FComplexSphericalHarmonic;


/**
 * compute the matrix to evaluate a spherical harmonic of order order at all positions thetaPhi
 */
template<class FSHType>
inline void computeForwardMatrix( FMatrix& B, const std::vector<std::pair<double, double> > &thetaPhi, int order )
{
  B.resize( thetaPhi.size(), FSHType::sz(order) );
  typename FSHType::iterator it;
  for(unsigned int i=0; i< FSHType::sz(order); ++i)
  {
    FSHType sh; sh.resize(order);
    sh( it.getOrder(), it.getMode() ) = 1.0;
    for(unsigned int j=0; j< thetaPhi.size(); ++j)
    {
      B(j,i) = sh.evaluate( thetaPhi[j].first, thetaPhi[j].second );
    }
    it++;
  }
}


/**
 * compute the matrix to project input data at directions thetaPhi to a spherical harmonic of order order
 */
template<class FSHType>
inline void computeMatrix( FMatrix &BI, const std::vector<std::pair<double, double> > &thetaPhi, int order )
{
  eassert( thetaPhi.size() >= FSHType::sz(order) );

  FMatrix M;
  computeForwardMatrix<FSHType>( M, thetaPhi, order );

#ifdef VERBOSE
  std::cout << M << std::endl;
#endif

  // invert Matrix B to BI ( pseudo inverse )
  FMatrix U(M); // copy matrix, it will be modified
  FMatrix v(FSHType::sz(order), FSHType::sz(order));
  FVector w(FSHType::sz(order));
  U.svdcmp( w, v );
  FMatrix d(w.getDimension()); // build diagonal matrix from vector
  for(unsigned int c=0; c<w.getDimension();c++)
  {
    if(w[c]==0.0) d(c,c) = 0.0; // never happens in my calculation
    else d(c,c) = 1./w[c];//1.0 / w[c];
  }

  BI = v*d*U.transpose();
#ifdef VERBOSE
  std::cout << "Matching Matrix: " << std::endl;
  std::cout << BI << std::endl;
#endif

};








inline
void copySH( FRealSphericalHarmonic &to, const FRealSphericalHarmonic & from )
{
  // untested
  /*
  to.data = from.data;
  to.order = from.order;
  */
  FRealSphericalHarmonic::iterator it;
  to.resize( from.getOrder() );
  while(it.getOrder() < from.getOrder() )
  {
    to( it.getOrder(), it.getMode()) = from( it.getOrder(), it.getMode());
    it++;
  }
}

inline
void copySH( FRealSphericalHarmonic &to, const FRealSymmetricSphericalHarmonic & from )
{
  // untested
  FRealSymmetricSphericalHarmonic::iterator it;
  to.resize( from.getOrder() );
  to.clear();
  while(it.getOrder() < from.getOrder() )
  {
    to( it.getOrder(), it.getMode()) = from( it.getOrder(), it.getMode());
    it++;
  }
}

inline
void copySH( FComplexSphericalHarmonic &to, const FDescoteauxSphericalHarmonic & from )
{
  // untested
  to.resize( from.getOrder() );
  to.clear();

  for(int l=0; l<= from.getOrder(); l+=2)
  {
    for(int m=-l; m<=l; ++m)
    {
      if(m < 0)
      {
        if( (-m)%2 != 0)
          to(l,m) += -M_SQRT1_2*from(l,m);
        else
          to(l,m) += M_SQRT1_2*from(l,m);
        to(l,-m) += M_SQRT1_2*from(l,m);
      }
      else if( m> 0)
      {
        if( (m)%2 != 0)
          to(l,m) += std::complex<double>(0., M_SQRT1_2)*from(l,m);
        else
          to(l,m) += std::complex<double>(0., -M_SQRT1_2)*from(l,m);
        to(l,-m) += std::complex<double>(0., M_SQRT1_2)*from(l,m);
      }
      else
      {
        to(l,m) += 1.*from(l,m);
      }
      to(l,m) +=1.*from(l,m);
    }
  }
}


inline
void copySH( FComplexSphericalHarmonic &to, const FRealSphericalHarmonic & from )
{
  // untested
  FRealSphericalHarmonic::iterator it;
  while( it.getOrder() <= from.getOrder() )
  {
    const unsigned int l = it.getOrder();
    const int m = it.getMode();
    if(m < 0 )
    {
      to(l,m) = std::complex<double>( from(l,-m)/2., -from(l,m));
    }
    else if( m > 0 )
    {
      to(l,m) = std::complex<double>( from( l, -m), from(l,-m)/2. );
    }
    else // m == 0
    {
      to(l, 0 ) =  from( l, 0 ); // real only
    }
    it++;
  }
}

inline
void copySH( FRealSphericalHarmonic & to, const FComplexSphericalHarmonic & from )
{
  // untested
  std::cerr << "Warning: copying complex to real spherical harmonic, information may be lost" << std::endl;

  to.resize( from.getOrder() );
  FRealSphericalHarmonic::iterator it;
  while( it.getOrder() <= from.getOrder() )
  {
    const unsigned int l = it.getOrder();
    const int m = it.getMode();
    if(m < 0 )
    {
      to(l,m) = -from(l,m).imag() + from(l,-m).imag();
    }
    else if( m > 0 )
    {
      to(l,m) = from(l,m).real() + from(l,-m).real();
    }
    else // m == 0
    {
      to(l, 0 ) =  from( l, 0 ).real(); // real only
    }
    it++;
  }
}

inline
void copySH( FRealSymmetricSphericalHarmonic &to, const FRealSphericalHarmonic & from )
{
  // untested
  std::cerr << "Warning: copying real  to real symmetric spherical harmonic, information may be lost" << std::endl;
  
  FRealSymmetricSphericalHarmonic::iterator it;
  while(it.getOrder() <= from.getOrder() )
  {
    const unsigned int l = it.getOrder();
    const int m = it.getMode();
    to(l,m) = from(l,m);
    it++;
  }
}

inline
void copySH( FComplexSphericalHarmonic &to, const FRealSymmetricSphericalHarmonic &from )
{
  // untested
  FRealSphericalHarmonic rsh;
  rsh.set( from );
  to.set( rsh );
}

#endif
