#include "FSphericalHarmonics.hh"

#include <math.h>

#include <iostream>
#include "FException.hh"

#include "FMatrix.hh"
#include "eassert.hh"

#include "FTensor.hh"
#include "FSpecialFunctions.hh"
#define VERBOSE
//using namespace std;
using namespace FMath;

#if 0
//
// these are rotations of the first 2 bands around the
// x axis by 90 degree.
//
// It is possible to construct arbitrary rotations of 
// these first bands with this matrices and those of
// getZRotationMatrix(...)
//
double Xm90[] =
{
  1,  0, 0, 0,  0, 0, 0, 0, 0,
  
  0,  0, 1, 0,  0, 0, 0, 0, 0,
  0, -1, 0, 0,  0, 0, 0, 0, 0,
  0,  0, 0, 1,  0, 0, 0, 0, 0,
  
  0,  0, 0, 0,  0, 0, 0, 1, 0,
  0,  0, 0, 0,  0,-1, 0, 0, 0,
  0,  0, 0, 0,  0, 0,-0.5,0,-sqrt(3)/2,
  0,  0, 0, 0, -1, 0, 0, 0, 0,
  0,  0, 0, 0,  0, 0,-sqrt(3)/2,0,0.5
    };

double Xp90[] =
{
  1,  0, 0, 0,  0, 0, 0, 0, 0,  
  
  0,  0,-1, 0,  0, 0, 0, 0, 0,
  0,  1, 0, 0,  0, 0, 0, 0, 0,
  0,  0, 0, 1,  0, 0, 0, 0, 0,
  
  0,  0, 0, 0,  0, 0, 0,-1, 0,
  0,  0, 0, 0,  0,-1, 0, 0, 0,
  0,  0, 0, 0,  0, 0,-0.5,0,-sqrt(3)/2,
  0,  0, 0, 0,  1, 0, 0, 0, 0,
  0,  0, 0, 0,  0, 0,-sqrt(3)/2,0,0.5
};

#endif



// rotate a point (theta,phi) 90 degree around x-axis
void rotX90(double &theta2, double&phi2, double theta, double phi)
{
  const double st= sin(theta);
  phi2   = atan2(-cos(theta), st*cos(phi));
  phi2   = fmod( phi2+2*M_PI, 2*M_PI);
  theta2 = acos (st*sin(phi));
}

// rotate a point (theta,phi) back to original coordinate system
void invRotX90(double &theta2, double&phi2, double theta, double phi)
{
  const double st= sin(theta);
  const double cp= cos(phi);
  phi2   = atan2(cos(theta), cp*st);
  theta2 = acos (-cp*st);
}


// to simplify factorial calculation I only calculate the fraction
// (l-m)!/(l+m)! = 1./[(l-m+1)*..*(l+m)]
// static helper
double FSphericalHarmonics::calcdenomfactorial( unsigned int from, unsigned int to)
{
    double value = 1.0;
    for(unsigned int i=from; i<=to; ++i)
    {
        value *= i;
    }
    return value;
}


// static
double FSphericalHarmonics::funcK(unsigned int l, int m)
{
    if(m<0) m = -m;
    double num = (2*l+1);//*factorial(l-m);
    double denom = 4.*M_PI*calcdenomfactorial(l-m+1, l+m);//*factorial(l+m);
    return sqrt(num/denom);
}

// static
std::complex<double> FSphericalHarmonics::complexY(unsigned int l, int m, double theta, double phi)
{
    if( m > 0 )
    {
        double factor = M_SQRT2 * FSphericalHarmonics::funcK(l,m) * FMath::legendre(l, m, cos(theta));
        return std::complex<double>( factor * cos(m*phi), factor * sin(m*phi) );
    }
    else if( m < 0 )
    {
        // Y_l,-m(\theta,\phi) = (-1)^m Y_l,^m(\theta, \phi) ??
        double factor = M_SQRT2 * FSphericalHarmonics::funcK(l,-m) * FMath::legendre(l, -m, cos(theta));
//        double sign = ((-m)%2 == 0) ? 1. : -1;
        return std::complex<double>( factor * cos(m*phi), -factor * sin(-m*phi) );
    }
    else
    {
        return FSphericalHarmonics::funcK(l,0) * FMath::legendre(l,0, cos(theta)); // only real valued result.
    }
}

#if 0
// static
double FSphericalHarmonics::imagY(unsigned int l, int m, double theta, double phi)
{
    std::complex<double> val = FSphericalHarmonics::complexY(l,m,theta,phi);
    return val.imag();
}

// static
double FSphericalHarmonics::realY(unsigned int l, int m, double theta, double phi)
{
    std::complex<double> val = FSphericalHarmonics::complexY(l,m,theta,phi);
    return val.real();
}

//
double FSphericalHarmonics::squaredY(unsigned int l, int m, double theta, double phi)
{
    std::complex<double> val = FSphericalHarmonics::complexY(l,m,theta,phi);
    return val.real()*val.real()+val.imag()*val.imag();
}
#endif
// static
double FSphericalHarmonics::therealY(unsigned int l, int m, double theta, double phi)
{
    double ret=0;
    if( m > 0 )
    {
        ret = M_SQRT2 * FSphericalHarmonics::funcK(l, m) * cos( m * phi ) * FMath::legendre(l, m, cos(theta));
    }
    else if (m < 0)
    {
        ret = M_SQRT2 * FSphericalHarmonics::funcK(l,-m) * sin( m * phi ) * FMath::legendre(l,-m, cos(theta)); // sin has a -, FIXME?
    }
    else // m == 0
    {
        ret = FSphericalHarmonics::funcK(l, 0) * FMath::legendre(l,0, cos(theta));
    }
    return ret;
}
//#if 0
#ifndef NEW_SH
// static
double FSphericalHarmonics::derivativeY(int l, int m, double theta, double phi, bool dir_is_theta)
{
  double ret = 0;
  if(!dir_is_theta)
  {
    if( m > 0 )
    {
      ret = M_SQRT2 * FSphericalHarmonics::funcK(l, m) * (-m)* sin( m * phi ) * FMath::legendre(l, m, cos(theta));
    }
    else if( m < 0 )
    {
      ret = M_SQRT2 * FSphericalHarmonics::funcK(l,-m) * m * cos( m* phi ) * FMath::legendre(l,-m, cos(theta));
    }
    else
    {
      ret = 0;
    }
  }
  else
  {
    if( m > 0 )
    {
      ret = M_SQRT2 * FSphericalHarmonics::funcK(l, m) * cos( m * phi ) * sin(theta) * FMath::legendreDerivate( l, m, cos(theta) );
    }
    else if( m < 0 )
    {
      ret = M_SQRT2 * FSphericalHarmonics::funcK(l, m) * sin( m * phi ) * sin(theta) * FMath::legendreDerivate( l,-m, cos(theta) );
    }
    else
    {
      ret = FSphericalHarmonics::funcK(l, 0) * sin(theta) * FMath::legendreDerivate( l, 0, cos(theta) );
    }
  }
  return ret;
}

void FSphericalHarmonics::derivativeY(FVector& deriv, int l, int m, double theta, double phi)
{
  deriv.resize(2);
  deriv[0] = derivativeY(l,m,theta, phi, true);
  deriv[1] = derivativeY(l,m,theta, phi, false);
}

// use_theta = true: derivation along theta direction
double FSphericalHarmonics::derivate( double theta, double phi, bool use_theta ) const
{
  if(use_theta)
  {
  double resultt(0);
  for(unsigned int l=0; l< nbBands; ++l)
  {
    for(int m=-l; m<=(int)l; ++m)
    {
      resultt += (*this)(l,m) * derivativeY(l,m,theta, phi, true);
    }
  }
  return -resultt;
  }

  double resultp(0);
  for(unsigned int l=0; l< nbBands; ++l)
  {
    for(int m=-l; m<=(int)l; ++m)
    {
      resultp += (*this)(l,m) * derivativeY(l,m,theta, phi, false);
    }
  }
  return resultp;
}
#endif // NEW_SH
//#endif // if 0
void FSphericalHarmonics::derivative(FVector& deriv, double theta, double phi) const
{
#ifdef NEW_SH
  deriv.resize(2);
  double ctheta = cos(theta);
  double stheta = sin(theta);

  // check how close to the z-axis we are
  if(1.0 - fabs(ctheta) >= 10*F_DBL_EPSILON)
  {
    // ok, we are fine, go on
    std::vector<double> vals(nbBands);
    std::vector<double> deri(nbBands);
    double rphi  =0.0;
    double rtheta=0.0;

    FMath::legendreSHPlmDerivArray( vals, deri, nbBands-1, 0, ctheta);
    for(unsigned int l=0; l<nbBands; ++l)
    {
      rtheta += (*this)(l,0)*deri.at(l);
      rphi   += (*this)(l,0)*vals.at(l);
    }
    for(unsigned int m=1; m<nbBands; ++m)
    {
      const double smphi = sin(m*phi);
      const double cmphi = cos(m*phi);
      FMath::legendreSHPlmDerivArray( vals, deri,nbBands-1, m, ctheta);
      const double f= IS_ODD(m) ? 1.0 : -1.0;
      for(unsigned int l=0; l<nbBands-m; ++l)
      {
        rtheta += (*this)(l+m,m)*deri[l] *cmphi*stheta;
        rtheta += (*this)(l+m,-m)*f*deri[l] *smphi*stheta;
        rphi   += (*this)(l+m,m)*vals[l]*(-m)*smphi;
        rphi   += (*this)(l+m,-m)*f*vals[l]*(m)*cmphi;
      }
    }
    deriv(1) = rphi;
    deriv(0) = rtheta;
  }
  else
  {
    // oops, not so fine
    std::cerr << "FIXME: implement rotated derivatives!!!" << std::endl;

    double thetaO = theta, phiO = phi;
    rotX90( theta, phi, thetaO, phiO);
   
    if(!shrot)
    {
      FMatrix rotation;
      FMatrix smallrot(3,3);
      smallrot(0,1) =-1.0;
      smallrot(1,0) = 1.0;
      smallrot(2,2) = 1.0;
      getRotationMatrix( rotation, smallrot, nbBands ); 
      shrot = new FSphericalHarmonics( rotation * (*this) );
    }

    std::cout << "Rotated sh: " << (*shrot) << std::endl;
    
    // ok, now we should be fine
    std::vector<double> vals(nbBands);
    std::vector<double> deri(nbBands);
    double rphi  =0.0;
    double rtheta=0.0;

    FMath::legendreSHPlmDerivArray( vals, deri, nbBands-1, 0, ctheta);
    for(unsigned int l=0; l<nbBands; ++l)
    {
      rtheta += (*shrot)(l,0)*deri.at(l);
      rphi   += (*shrot)(l,0)*vals.at(l);
    }
    for(unsigned int m=1; m<nbBands; ++m)
    {
      const double smphi = sin(m*phi);
      const double cmphi = cos(m*phi);
      FMath::legendreSHPlmDerivArray( vals, deri,nbBands-1, m, ctheta);
      const double f= IS_ODD(m) ? 1.0 : -1.0;
      for(unsigned int l=0; l<nbBands-m; ++l)
      {
        rtheta += (*shrot)(l+m,m)*deri[l] *cmphi*stheta;
        rtheta += (*shrot)(l+m,-m)*f*deri[l] *smphi*stheta;
        rphi   += (*shrot)(l+m,m)*vals[l]*(-m)*smphi;
        rphi   += (*shrot)(l+m,-m)*f*vals[l]*(m)*cmphi;
      }
    }

    deriv(1) = rphi;
    deriv(0) = rtheta;
    
    return;
  }

#else
  deriv.resize(2);
  double resultt(0);
  double resultp(0);
  for(unsigned int l=0; l< nbBands; ++l)
  {
    for(int m=-l; m<=(int)l; ++m)
    {
      resultt += (*this)(l,m) * derivativeY(l,m,theta, phi, true);
      resultp += (*this)(l,m) * derivativeY(l,m,theta, phi, false);
    }
  }
  deriv[0] = resultt;
  deriv[1] = resultp;
#endif
}

double FSphericalHarmonics::evaluate(double theta, double phi) const
{
#ifdef NEW_SH
  const double ctheta = cos(theta);
  double result = 0.0;
  std::vector<double> vals(nbBands);
  FMath::legendreSHPlmArray( vals, nbBands-1, 0, ctheta);
  eassert( vals.size() == nbBands );
  for(unsigned int l=0; l<nbBands; ++l)
    result += (*this)(l,0)*vals[l];
  
  for(unsigned int m=1; m<nbBands; ++m)
  {
    const double smphi = sin(m*phi);
    const double cmphi = cos(m*phi);
    FMath::legendreSHPlmArray( vals, nbBands-1, m, ctheta);
    const double f= IS_ODD(m) ? 1.0 : -1.0;
    for(unsigned int l=0; l<nbBands-m; ++l)
    {
      result += (*this)(l+m,m)*vals[l] *cmphi;
      result += (*this)(l+m,-m)*f*vals[l] *smphi;
    }
  }
  return result;
#else
    double result(0);
    for(unsigned int l = 0; l< nbBands; ++l)
    {
        for(int m=-l; m<=(int)l; ++m)
        {
          double f = 1;
//          if((abs(m))%2 == 0) f=-1;
            result += f*(*this)(l,m) * therealY(l,m,theta, phi);
        }
    }
    return result;
#endif
}


const double& FSphericalHarmonics::operator()(unsigned int l, int m) const
{
    return comp[index(l,m)];
}


double& FSphericalHarmonics::operator()(unsigned int l, int m)
{
    return comp[index(l,m)];
}


FSphericalHarmonics::FSphericalHarmonics()
{
    nbBands = 0; // an invalid spherical harmonic which evaluates to constant zero
    shrot = 0;
}

FSphericalHarmonics::FSphericalHarmonics(unsigned int l) : /*FArray( l*(l+2)+1 ), */ nbBands(l+1), shrot(0)
{
    comp.resize( l*(l+2)+1 );
    //FArray::resize( l*(l+2)+1 );
//    std::cout << "Spherical harmonic with size: " << size() << " created."<< std::endl;
}

// static helper, const membor only needed for debugging
unsigned int FSphericalHarmonics::index(unsigned int l, int m)// const
{
    eassert(((int)l >= m) && ((int)l >= -m));
    // on level l we have : 2*l+1 values
    // the middle value of level l is: l;
    // up to level l including level l we have : l*(l+1) + l + 1  ==  (l+1)*(l+1) values
    // up to level l-1 we have (l)*(l) values, thus the first index of level l is l*l
    // so element (l,m) is at position:
    unsigned int index = (l)*(l)+( m + l);               //l*(l+1)+m;
#ifdef VERBOSE
//    if(index >= size()){
//        std::cerr << "ERROR: l=" << l << "m=" << m << "   => " << index << " but size=" << size() << std::endl;
//    }
#endif
    return index;
}

// iterative, thus a bit slow.
void FSphericalHarmonics::revindex(int &l, int &m, int size)// const
{
    l=0;
    while(true)
    {
        if( size <= 2*l+1) 
        {
            m=size-1;
            return;
        }
        size -= 2*l+1;
        l++;
    }   
}


void FSphericalHarmonics::precomputeMatch(FMatrix &BI, const std::vector<FVector> &vecs) const
{
    FSphericalHarmonics::precomputeMatch(BI,vecs, this->nbBands-1);
}

void FSphericalHarmonics::precomputeMatch(FMatrix &BI, const std::vector<FVector> &vecs, int order)
{
    int size = FSphericalHarmonics::index(order+1,-order-1);
#ifndef NODEBUG
    if(vecs.size() < (unsigned int)size)
      THROW_EXCEPTION( FInvalidDimensionException, "Number of vectors must be larger than size of FSphericalHarmonics" );
#endif
    // create the matrix
    FMatrix B(vecs.size(), size);
//    FVector b(vecs.size());

#ifdef NEW_SH
    {
      int nbBands = order+1;
      double theta, phi, r;
      unsigned int j; // tmp index
      std::vector<double> leg(size);  // buffer for legendre polynomials
     for(unsigned int i=0; i< vecs.size(); ++i)
      {
        if(vecs[i].getDimension() != 3) THROW_EXCEPTION( FInvalidDimensionException, "Dimension must be 3" );
        toThetaPhi(theta, phi, r, vecs[i]);
        double ct = cos(theta);

        FMath::legendreSHPlmArray( leg, nbBands-1, 0, ct);
        for(int l=0; l< nbBands; ++l)
        {
          unsigned int j=index(l,0);
          B(i,j) = leg[l]; // fixed
        }

        
        for(int m=1; m<nbBands; ++m)
        {
          double smp = sin(m*phi);
          double cmp = cos(m*phi);
          const double f= IS_ODD(m) ? 1.0 : -1.0;
          FMath::legendreSHPlmArray( leg, nbBands-1, m, ct);
          for(int l=m; l< nbBands; ++l)
          {
            j = index(l, m);
            B(i,j) = leg[l-m]*cmp;
            j = index(l,-m);
            B(i,j) = f*leg[l-m]*smp;
          }
        }
      }
    }
#else
    {
        double theta, phi, r;

        // create the matrix
        // If many values are calculated, this matrix should be calculated only once
        // for a given set of sample points. 
        for(unsigned int i=0; i< vecs.size(); ++i)
        {
            if(vecs[i].getDimension() != 3) THROW_EXCEPTION( FInvalidDimensionException, "Dimension must be 3" );
            toThetaPhi(theta, phi, r, vecs[i]);
            for(int l=0; l< order+1; ++l)
            {
                for(int m=-l; m<=(int)l; ++m)
                {
                    unsigned int j = index(l,m);
                    B(i,j) = therealY(l,m,theta,phi);
                }
            }
        }
    }
#endif
//#ifndef NODEBUG
//    std::cout << B << std::endl;
    FMatrix B_(B);
    FMatrix B_T(B_.transposed());
    FMatrix BB = B_T*B;
    FVector evals;
    std::vector< FArray > evecs;
    BB.getEigenSystem( evals, evecs );
    std::cout << "ES: " << evals << std::endl;
//#endif                    

    // invert Matrix B to BI ( pseudo inverse )
    FMatrix U(B); // copy matrix, it will be modified
    FMatrix v(size, size);
    FVector w(size);
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
}

void FSphericalHarmonics::precomputeSymmetricMatch(FMatrix &BI, const std::vector<FVector> &vecs) const
{
    FSphericalHarmonics::precomputeSymmetricMatch(BI,vecs, this->nbBands-1);
}

/* static */
void FSphericalHarmonics::precomputeSymmetricMatch(FMatrix &BI, const std::vector<FVector> &vecs, int order)
{
  int redorder = order/2;
  int size = 2*redorder*(redorder+1)+(redorder+1);
  int fullsize = FSphericalHarmonics::index(order+1,-order-1);
#ifndef NODEBUG
  std::cout << "size = " << size <<  "   vecs.size() = " << vecs.size() << std::endl;
  std::cout << "fullsize = " << fullsize << std::endl;
  std::cout << "order = " << order << " redorder = " << redorder << std::endl;
#endif
#ifndef NODEBUG
    if(vecs.size() < (unsigned int)size)
      THROW_EXCEPTION( FInvalidDimensionException, "Number of vectors must be larger than size of FSphericalHarmonics" );
#endif
  
    FMatrix B(vecs.size(), size);
#ifdef NEW_SH
    {
      int nbBands = order+1;
      //std::cout << "nbbands = " << nbBands << std::endl;
      double theta, phi, r;
      unsigned int j; // tmp index
      std::vector<double> leg(size);  // buffer for legendre polynomials
      for(unsigned int i=0; i< vecs.size(); ++i)
      {
        if(vecs[i].getDimension() != 3) THROW_EXCEPTION( FInvalidDimensionException, "Dimension must be 3" );
        toThetaPhi(theta, phi, r, vecs[i]);
        double ct = cos(theta);

        FMath::legendreSHPlmArray( leg, nbBands-1, 0, ct);
        for(int l=0; l<= order; l=l+2)
        {
          int ell = l/2;
          int j = (2*ell*ell+ell);
          //std::cout << "l=" << l << " j=" << j << std::endl;
          B(i,j) = leg[l]; // fixed
        }


        for(int m=1; m<=order; ++m)
        {
          double smp = sin(m*phi);
          double cmp = cos(m*phi);
          const double f= IS_ODD(m) ? 1.0 : -1.0;
          FMath::legendreSHPlmArray( leg, nbBands-1, m, ct);
          for(int l=m; l<= order; l++)
          {
            if(l%2 != 0) continue;
            int ell=(l)/2;
            //if(ell==0) continue;
            j = 2*ell*ell+ell+m;
            //std::cout << "+l=" << l << " ell=" << ell <<  " m=" << m << " j=" << j << " ell=" << ell << std::endl;
            B(i,j) = leg[2*ell]*cmp;
            j = 2*ell*ell+ell-m;
            //std::cout << "-l=" << l << " m=" << m << " j=" << j << " ell=" << ell << std::endl;
            B(i,j) = f*leg[2*ell]*smp;
          }
        }
      }
    }
#else
    {
      double theta, phi, r;
      for(unsigned int i=0; i< vecs.size(); ++i)
      {
        if(vecs[i].getDimension() != 3) THROW_EXCEPTION( FInvalidDimensionException, "Dimension must be 3" );
        toThetaPhi(theta, phi, r, vecs[i]);
        unsigned int j = 0;
        int l=0;
        for(l=0; l<= order; ++l)
        {
          if(l%2 == 1) continue;
          for(int m=-l; m<=(int)l; ++m)
          {
            B(i,j) = therealY(l,m,theta,phi);
            ++j;
          }
        }
      }
    }
#endif

#ifndef NODEBUG
    std::cout << B << std::endl;
    FMatrix B_(B);
    FMatrix B_T(B_.transposed());
    FMatrix BB = B_T*B;
    FVector evals;
    std::vector< FArray > evecs;
    BB.getEigenSystem( evals, evecs );
    std::cout << "ES: " << evals << std::endl;
#endif                    

    FMatrix BT(B.transposed());
    FMatrix mBI = BT*B;
    try{
      mBI.invert();
    }catch(FMatrixSingularException &e )
    {
      e.addTraceMessage("Problem not solvable, check gradient directions or gradient values");
      RETHROW_EXCEPTION(e);
    }
    mBI*=BT;

    FMatrix FullBI(fullsize, mBI.getDimensionX());

    int si = 0;
    for(unsigned int l=0; l<=(unsigned)order; ++l)
    {
      if(l%2 == 1) continue;
      for(int m=-l; m<=(int)l;++m)
      {
        unsigned int i= index(l,m);
        for(unsigned int j=0; j< vecs.size(); ++j)
          FullBI( i, j) = mBI(si,j);
        ++si;
      }
    }
#ifndef NODEBUG
#ifdef VERBOSE
    std::cout << mBI << std::endl;
    std::cout << FullBI << std::endl;
#endif
#endif

    BI = FullBI;
}


/*
void FSphericalHarmonics::precomputeComplexMatch(FMatrix &BI, const std::vector<FVector> &vecs, std::complex<double> values, int order)
{
  int size = FSphericalHarmonics::index(order+1,-order-1);
#ifndef NODEBUG
  eassert( vecs.size() > (unsigned int) size);
#endif
  FMatrix B(vecs.size()*2, size);
//  FVector b(vecs.size()*2);
  {
    double theta, phi, r;
    for(unsigned int i=0; i < vecs.size(); ++i)
    {
      eassert(vecs[i].getDimension() == 3);
      for(int l=0; l < order+1;++l)
      {
        for(int m=-l; m=(int)l; ++m)
        {
          toThetaPhi( theta, phi, r, vecs[i]);
          unsiged int j = index(l,m);
          std::complex<double> base= complexY(l,m,theta,phi);
          B(i*2  ,j) = base.real();
          B(i*2+1,j) = base.imag();
        }
      }
    }
  }
  FMatrix U(B);
  FMatrix v(size,size);
  FVector w(size);
  U.svdcmp( w, v );
  FMatrix d(w.getDimension());
  for(unsigned int c=0; c<w.getDimension();++c)
  {
    if(w[c]=0.0) d(c,c) = 0.0;
    else d(c,c) =1./w[c];
  }
  BI = v*d*U.transpose();
}

*/        

// rather inefficient but working.
// SVD would be a bit better afaik.
void FSphericalHarmonics::match(const std::vector<FVector> &vecs, const FArray& b)
{
    //    std::cout << *this << std::endl;


    //    std::cout << "Matrix:"
    //        << mat << std::endl;
    //    std::cout << b << std::endl;
#ifdef NOSVD
    {
        FMatrix matT = mat.transposed();
        FVector rhs = matT*b;
        //        std::cout << "r:" << rhs << std::endl;
        FMatrix A = matT*mat;
        //        std::cout << "matT*mat" << A << std::endl;
        try{ 
            A.invert();
            //            std::cout << "A.invert()" << A << std::endl;
        }
        catch(FMatrixSingularException &e)
        {
            solved = false;
            return false;
        }
        comp = A*rhs;
    }
#else
    FMatrix BI;
    precomputeMatch(BI,vecs);
    comp = BI*b;
#endif

//        this->threshold( 1e-9 );
}

void FSphericalHarmonics::matchSymmetric( const std::vector<FVector> &vecs, const FArray &b)
{
  FMatrix BI;
  precomputeSymmetricMatch( BI, vecs);
  comp = BI*b;
}

void FSphericalHarmonics::match(const FMatrix &BI, const FArray& vals)
{
    comp = BI*vals;
}
/*
void FSphericalHarmonics::symmetricMatch(const FMatrix &BI, const FArray& vals)
{
    comp = BI*vals;
}
*/

std::ostream& operator <<(std::ostream& o,const FSphericalHarmonics& sh)
{
    o << "[";
    for(int l=0; l< (int)sh.nbBands; ++l)
    {
        o << "[";
        for(int m=-l; m<=l; ++m)
        {
            o << sh(l,m) << "\t ";
        }
        o << "]";
        if( l < (int)sh.nbBands-1 ) o << std::endl;
    }
    o << "]" << std::endl;
    return o;
}



///////////////////////////////////////////////////////////////////////
/*
   FRealSphericalHarmonics::FRealSphericalHarmonics(unsigned int l)
   : FSphericalHarmonics()
   {
//    resize( (l+1)*(l+2)/2 );

nbBands = l+1;
#ifdef VERBOSE
std::cout << "Real Spherical harmonic with size: " << size() << " created."<< std::endl;
#endif
}



unsigned int FRealSphericalHarmonics::index(unsigned int l, int m) const
{
//    std::cout << "index real" << std::endl;
// l+1 values in level l
// (l+1)*(l+2)/2 values up to and including level l
// l*(l+1)/2 is the first index in level l
// m is always >= 0
#ifndef NDEBUG
eassert( m>= 0);
#endif
return l*(l+1)/2+m; 
}

double FRealSphericalHarmonics::evaluate(double theta, double phi) const 
{
//    std::cout << "evaluating real" << std::endl;
double result = 0;
for(unsigned int l = 0; l< nbBands; ++l)
for(int m=0; m<=(int)l; ++m)
{
result += (*this)(l,m) * realY(l,m,theta, phi);
}
return result;
}


bool FRealSphericalHarmonics::match(const std::vector<FVector> &vecs, const std::vector<double> vals)
{
//    std::cout << "matching real" << std::endl;
bool solved=true;
FMatrix mat(vals.size(), this->size());
FVector b(vals.size());

double theta, phi, r;

// create the matrix
// If many values are calculated, this matrix should be calculated only once
// for a given set of sample points. 
for(unsigned int i=0; i< vals.size(); ++i)
{
eassert(vecs[i].getDimension() == 3); // we are in 3D!
for(unsigned int l=0; l< nbBands; ++l)
{
for(int m=0; m<=(int)l; ++m)
{
toThetaPhi(theta, phi, r, vecs[i]);
unsigned int j = index(l,m);
//                std::cout << "index: "<< i << " " << j<< std::endl;
mat(i,j) = realY(l,m, theta, phi);
}
}
}



for(unsigned int i=0; i< vals.size(); ++i)
b[i] = vals[i];

//    std::cout << "Matrix:"
//              << mat << std::endl;
//    std::cout << b << std::endl;

{
    FMatrix matT = mat.transposed();
    FVector r = matT*b;
    //    std::cout << "r:" << r << std::endl;
    FMatrix A = matT*mat;
    //    std::cout << "matT*mat" << A << std::endl;
    try{ 
        A.invert();
        //        std::cout << "A.invert()" << A << std::endl;
    }
    catch(FMatrixSingularException &e)
    {
        solved = false;
        return false;
    }

    (*static_cast<FArray*>(this)) = A*r;
    //    cout << "Solution: " << (*this) << endl;
}

return solved;
}

///////////////////////////////////////////////////////////////////////

    FSymmetricSphericalHarmonics::FSymmetricSphericalHarmonics(unsigned int l)
: FSphericalHarmonics()
{
    nbBands = l+1;
    //    int valid_lines= l/2;
    comp.resize( (((l/2)*2)+1)*(l/2+1));
#ifdef VERBOSE
    std::cout << "Symmetric Spherical harmonic with size: " << size() << " created."<< std::endl;
#endif
}



unsigned int FSymmetricSphericalHarmonics::index(unsigned int l, int m) const
{
    //    std::cout << "index real" << std::endl;
    // on level l we have : 2*l+1 values
    // the middle value of level l is: l;
    // up to level l including level l we have : (((l/2)*2)*(l/2+1)+l/2+1 values
#ifndef NDEBUG
    eassert( l % 2 == 0);
#endif
    int line = l/2;
    int pline= line-1;
    return 2*pline*(pline+1)+pline+1+l+m;
}

double FSymmetricSphericalHarmonics::evaluate(double theta, double phi) const 
{
    //    std::cout << "evaluating real" << std::endl;
    double result = 0;
    for(unsigned int l = 0; l< nbBands; l+=2)
        for(int m=0; m<=(int)l; ++m)
        {
            result += (*this)(l,m) * realY(l,m,theta, phi);
        }
    return result;
}


bool FSymmetricSphericalHarmonics::match(const std::vector<FVector> &vecs, const std::vector<double> vals)
{
    std::cout << "matching symmetric" << std::endl;
    bool solved=true;
    FMatrix mat(vals.size(), this->size());
    FVector b(vals.size());

    double theta, phi, r;

    // create the matrix
    // If many values are calculated, this matrix should be calculated only once
    // for a given set of sample points. 
    for(unsigned int i=0; i< vals.size(); ++i)
    {
        eassert(vecs[i].getDimension() == 3); // we are in 3D!
        for(unsigned int l=0; l< nbBands; l+=2)
        {
            for(int m=-l; m<=(int)l; ++m)
            {
                toThetaPhi(theta, phi, r, vecs[i]);
                unsigned int j = index(l,m);
                std::cout << "index: "<< i << " " << j<< std::endl;
                mat(i,j) = imagY(l,m, theta, phi);
            }
        }
    }



    for(unsigned int i=0; i< vals.size(); ++i)
        b[i] = vals[i];

    std::cout << "Matrix:"
        << mat << std::endl;
    std::cout << b << std::endl;

    {
        FMatrix matT = mat.transposed();
        FVector r = matT*b;
        //    std::cout << "r:" << r << std::endl;
        FMatrix A = matT*mat;
        std::cout << "matT*mat" << A << std::endl;
        try{ 
            A.invert();
            std::cout << "A.invert()" << A << std::endl;
        }
        catch(FMatrixSingularException &e)
        {
            solved = false;
            return false;
        }

        (*static_cast<FArray*>(this)) = A*r;
        //    cout << "Solution: " << (*this) << endl;
    }
    return solved;
}
*/

FTensor FSphericalHarmonics::toTensor() const
{
  FTensor tensor;
  if(nbBands > 1)
    tensor.resizeTensor(3, this->nbBands-1);
  else
    tensor.resizeTensor(1, 0); // stupid FAnToM convention,
                      // dimension should simply be ignored, but isn't everywhere
  
  switch( this->nbBands )
  {
    case 1: // Tensor of order 0 (scalar)
      {
        tensor() = comp[0] / sqrt(4.*M_PI);
      }break;
    case 3: // Tensor of order 2
      {
        // warning: this only handles symmetric part!
#ifndef NODEBUG
        double delta=0;
        for(unsigned int i=1; i< 4; ++i)
        {
          delta += fabs(comp[i]);
        }
        if(delta > 1e-5)
        {
          THROW_EXCEPTION(FException, "Spherical harmonic contains antisymmetric parts");
        }
#endif
        static double f1 = 1./sqrt(4.*M_PI);
        static double f2 = sqrt(5./16/M_PI);
        static double f3 = sqrt(15./8/M_PI);

        tensor(0,0) = (*this)[0]*f1 - f2* (*this)[6] + f3*(*this)[8];
        tensor(1,1) = (*this)[0]*f1 - f2* (*this)[6] - f3*(*this)[8];
        tensor(2,2) = (*this)[0]*f1 + 2.*f2 * (*this)[6];
        tensor(0,1) = tensor(1,0) = -f3*(*this)[4];
        tensor(0,2) = tensor(2,0) = -f3*(*this)[7];
        tensor(1,2) = tensor(2,1) = +f3*(*this)[5];
      }break;
    default:
      THROW_DEFAULT_EXCEPTION(FNotImplementedException);
      break;
  }
  return tensor;
}

const FArray& FSphericalHarmonics::toArray() const
{
  return comp;
}

void FSphericalHarmonics::getEigenSystem(FVector& vals, FVector v[3]) const
{
#ifndef NODEBUG
  if(this->nbBands < 2) // we need 3 bands to get the tensor
    THROW_DEFAULT_EXCEPTION(FInvalidDimensionException);
#endif
  eassert( this->size() >= 9);
  
  FArray data(comp);
  data.resize(9,true);// resize but keep content
  FSphericalHarmonics sh(data);
  FTensor tensor=sh.toTensor();
  tensor.getEigenSystem(vals, v );
}

FSphericalHarmonics::FSphericalHarmonics(const FTensor& tensor)
{
  setFromTensor( tensor, true, true );
}

void FSphericalHarmonics::setFromTensor(const FTensor& t, bool resize, bool clear)
{
#ifndef NODEBUG
  if(t.getDimension() != 3) THROW_DEFAULT_EXCEPTION(FInvalidDimensionException);
#endif
  
  unsigned int l = t.getOrder();
  if(resize)
  {
    comp.resize( l*(l+2)+1 );
    nbBands = l+1; 
  }

#ifndef NODEBUG
  if(nbBands-1 < t.getOrder()) THROW_DEFAULT_EXCEPTION(FInvalidDimensionException);
#endif

  if(clear)
  {
    for(unsigned int i=0; i< size(); ++i)
    {
      comp[i] = 0;
    }
  }

  switch( t.getOrder() )
  {
    case 0:
      {
    comp[0] = sqrt(4.*M_PI)*t();
      }break;
    case 1:
      {
        THROW_DEFAULT_EXCEPTION(FNotImplementedException);
        /*
        comp[0] = t(0);
        comp[1] = t(1);
        comp[2] = t(2);
        */
      }break;
    case 2:
      {
        comp[0] = sqrt(4.*M_PI)/3.*( t(0,0) + t(1,1) + t(2,2));
        comp[8] = sqrt(2.*M_PI/15.)*(t(0,0)-t(1,1)); // Re 2,2
        comp[4] = sqrt(2.*M_PI/15.)*(-2.*t(0,1));    // Im 2,2
        comp[7] = sqrt(8.*M_PI/15.)*(-t(0,2));       // Re 2,1
        comp[5] = sqrt(8.*M_PI/15.)*(t(1,2));        // Im 2,1
        comp[6] = -sqrt(4.*M_PI/45.)*(t(0,0)+t(1,1)-2*t(2,2)); // 2,2
      }break;
    default:
      THROW_DEFAULT_EXCEPTION(FNotImplementedException);
      break;
  }
}

// working, tested
void FSphericalHarmonics::getComplexCoefficients( std::vector< std::complex< double > >& data ) const
{
  data.resize(this->size());
  int i=0;
  for(unsigned int l=0; l< nbBands; ++l)
  {
    for(int m=-l; m<=(int)l; ++m)
    {
      if(m<0)
        data[i] = std::complex< double >( comp[i-2*m]/2., -comp[i]/2.    );
      if(m>0)
        data[i] = std::complex< double >( comp[i]/2.    ,  comp[i-2*m]/2.);
      if(m==0)
        data[i] = comp[i]; // only real 
      ++i;
    }
  }
}

// working, tested
void FSphericalHarmonics::setFromComplexCoefficients( std::vector< std::complex< double > > cdata )
{
  unsigned int i=0;
  int l=0;
  do
  {
    for(int m=-l; m<=(int)l; ++m)
    {
      if(m<0)
        comp[i] = ( -cdata[i].imag()+ cdata[i-2*m].imag() );
      if(m==0)
        comp[i] = ( cdata[i].real() );
      if(m>0)
        comp[i] = ( cdata[i].real()+cdata[i-2*m].real() );
      ++i;
    }
    ++l;
  }while(i<cdata.size());
}

// just a guess, _NOT_TESTED_ !!!!!!!!!!
void FSphericalHarmonics::setComplexFromComplexCoefficients( std::vector< std::complex< double > > cdata )
{
  unsigned int i=0;
  int l=0;
  do{
    for(int m=-l; m<=(int)l; ++m)
    {
      if(m<0)
        comp[i] = ( cdata[i].imag() + cdata[i-2*m].imag() );
      if(m==0)
        comp[i] = ( cdata[i].imag() );
      if(m>0)
        comp[i] = ( cdata[i].real() - cdata[i-2*m].real() );
      ++i;
    }
    ++l;
  }while(i<cdata.size());
}

// tested, working
/* static */
std::complex< double > FSphericalHarmonics::evaluateComplex( const std::vector< std::complex< double > > &cdata, double theta, double phi )
{
  int L;
  int M;
  revindex(L,M, cdata.size());

  unsigned int nbBands = (unsigned int)L+1;
  std::complex < double > result;
  int i=0;
  for(unsigned int l = 0; l< nbBands; ++l)
    for(int m=-l; m<=(int)l; ++m)
    {
      result += cdata[i] * FSphericalHarmonics::complexY(l,m,theta, phi);
      ++i;
    }
  return result;
}

void FSphericalHarmonics::resizeBands(unsigned int nbBands, bool clear)
{
#ifndef NODEBUG
  if(nbBands < 1) THROW_DEFAULT_EXCEPTION(FInvalidDimensionException);
#endif
  this->nbBands = nbBands;
  comp.resize( (nbBands-1)*(nbBands+1)+1, !clear );
}


// static
void FSphericalHarmonics::getZRotationMatrix( FMatrix& r, double theta, int bands )
{

  //std::cout << "pass" << std::endl;
  int size = (bands-1)*(bands+1)+1;
  //std::cout << size << std::endl;
  r.resize(size,size, false);
  int pos=0;
  int sz =1;
  for(int l=0; l<bands; ++l)
  {
    //std::cout << "r= " << r << std::endl;
    if(l == 0 )
    {
      r(0,0) = 1.0;
    }
    else
    {
      //std::cout << pos << " " << size << " " << sz <<std::endl;
      FMatrix s;
      s.resize(sz,sz);
      r.getSubMatrix( s, pos,pos );
      pos += sz;
      r.setSubMatrix( pos+1, pos+1, s);
    
    
    sz = 2*l+1;

    double ct = cos(theta*l);
    double st = sin(theta*l);
    
    r(pos,pos) = ct;
    r(pos+0, pos+sz-1) = st;
    r(pos+sz-1, pos) = -st;
    r(pos+sz-1, pos+sz-1) = ct;
    } 
  }
}



namespace FMath
{

// tangent vector along increasing theta direction
void toTangentTheta(FVector &v, double theta, double phi)
{
  v.resize(3);
  v[0] = -cos(phi)*cos(theta);
  v[1] = -sin(phi)*cos(theta);
  v[2] = sin(theta);
}

// tangent along increasing phi direction
// (which is independant of theta parameter)
void toTangentPhi(FVector &v, double /*theta*/, double phi)
{
  v.resize(3);
  v[0] = -sin(phi);
  v[1] = cos(phi);
  v[2] = 0.;
}

void toVector(FVector &v, double theta, double phi, double r=1.0)
{
    v.resize(3);
    v[0] = sin(theta)*cos(phi)*r;
    v[1] = sin(theta)*sin(phi)*r;
    v[2] = cos(theta)*r;
}


void toThetaPhi(double &theta, double &phi, double &r, const FVector& v)
{
//    cout << v << endl;
    if(v.getDimension() != 3) THROW_EXCEPTION( FInvalidDimensionException, "Dimension must be 3." );
    r = v.norm();
    eassert( r > 1E-10);
    phi = atan2(v[1],v[0]);
    phi = fmod( phi+2*M_PI, 2*M_PI);
    theta=acos(v[2]/r);
//    cout << phi << " " << theta << " " << r << endl;
}

// Rotation matrix computation as defined in
// "Rotation Matrices for Real Spherical Harmonics"
// by Joseph Ivanic and Klaus Ruedenberg
// in "Journal of Physical Chemistry 1996, Vol. 100
// pages 6342-6347
// and corrections in
// "Journal of Physical Chemistry A", Vol. 102,
// No 45 (1998, 1999) pages 9099f
// as shown, too, in spherical-harmonic-lighting.pdf
// "Spherical Harmonic Lighting, the gritty details"


//static FMatrix ROT;
/**
 * compute values using the knowledge of M^(l-1) stored in Ml_1
 */
double sub_P(const int i, const int l, const int a, const int b, const FMatrix& Ml_1, const FMatrix &ROT )
{
  try{
    //std::cout << "sub_P( " << i << ", " << l << ", " << a << ", " << b << ")" << std::endl;
    int off = l-1;

    // quite a lot argument checks
    // Do not try to call this function with anything where sub_{u|v|w} == 0, this would cause
    // unwanted parameters here... funny though his is nowhere documented in the paper.
    eassert(ROT.getDimensionX() == 3);
    eassert(ROT.getDimensionY() == 3);
    eassert( -1 <= i && i <= 1 );
    eassert( Ml_1.getDimensionX() == Ml_1.getDimensionY() );
    eassert( (int)Ml_1.getDimensionX() >=3 );
    eassert( (int)Ml_1.getDimensionY() >=3 );
    eassert( (int)Ml_1.getDimensionX() >a+off );
    eassert( (int)Ml_1.getDimensionX() >=b+off );
    eassert( a+off >= 0);
    eassert( (int)b+off >= -1);

    if( abs(b) < l )
    {
      return ROT(i+1,1) * Ml_1(a+(off),b+(off)); // i,0
    }
    if( b == l )
    {
      return ROT(i+1,1+1)* Ml_1(a+(off),l-1+off) - ROT(i+1,-1+1)* Ml_1(a+off,-l+1+off);
    }

    return ROT(i+1,1+1)* Ml_1(a+(off),-l+1+(off)) + ROT(i+1,-1+1)*Ml_1(a+off,(int)l-1+(off));

  }CATCH_N_RETHROW( FException );
}

inline double sub_U( const int l, const int m, const int n, const FMatrix& last, const FMatrix& rot )
{
  return sub_P(0,l,m,n, last,rot);
  //  if(m==0) return sub_P(0,l,0,n,last);
  //  else if(m>0) return sub_P(0,l,m,n,last);
  //  else if(m<0) return sub_P(0,l,m,n,last);
}

inline double sub_V(const int l, const int m, const int n, const FMatrix &last, const FMatrix &rot )
{
  if(m == 0)
    return sub_P(1,l,1,n,last,rot) + sub_P(-1,l,-1,n,last,rot);
  else if(m>0)
  {
    double sqr = (m==1) ? M_SQRT2 : 1.0;
    double ret = sqr*sub_P(1,l,m-1,n,last,rot);
    if(m != 1) return ret - sub_P(-1,l,-m+1,n,last,rot);
    return ret;
  }
  else
  {
    double sqr = (m==-1) ? M_SQRT2 : 1.0; // this seems to be wrong in the papers?
    double ret = sqr*sub_P(-1,l,-m-1,n,last,rot);
    //    double ml = (m == -l) ? 2.0 : 1.0;
    if(m != -1) return ret + sub_P(1,l,m+1,n,last,rot);
    return ret;
  }
}

inline double sub_W(const int l, const int m, const int n, const FMatrix& last, const FMatrix& rot )
{
  if(m>0)
    return sub_P(1,l,m+1,n,last,rot)+sub_P(-1,l,-m-1,n, last,rot);
  if(m<0)
    return sub_P(1,l,m-1,n,last,rot)-sub_P(-1,l,-m+1,n, last,rot);

  // m==0 never happens because sub_w( m==0 ) == 0, thus
  // sub_W should never be evaluated!!!
  eassert(false);
  return 0.0;
}

double sub_u(const int l, const int m, const int n )
{
  double a1 = l+m;
  double a2 = l-m;
  if(abs(n)==l)
  {
    return sqrt((a1*a2)/((2.*l)*(2.*l-1)));
  }
  double a3 = l+n;
  double a4 = l-n;
  return sqrt((a1*a2)/(a3*a4));
}

double sub_v(const int l, const int m, const int n)
{
  double a1 = (m==0) ? 2 :1;
  double a2 = l+abs(m)-1;
  double a3 = l+abs(m);
  double a6 = (m==0) ? -1: 1;

  if(abs(n) == l )
  {
    return 0.5*sqrt((a1*a2*a3)/(2.*l)/(2.*l-1))*a6;
  }
  double a4 = l+n;
  double a5 = l-n;

  return 0.5*sqrt((a1*a2*a3)/(a4*a5))*a6;
}

double sub_w(const int l, const int m, const int n)
{
  //we shouldn't consider calling this for m==0 because it returns zero
  //  eassert( m != 0);
  if(m==0) return 0.0; // this is the (1-delta(m,0)) term

  double a1 = l-abs(m)-1;
  double a2 = l-abs(m);

  if(abs(n)==l)
  {
    return -0.5*sqrt((a1*a2)/((2.*l)*(2.*l-1)));
  }
  double a3 = l+n;
  double a4 = l-n;
  return - 0.5*sqrt((a1*a2)/a3/a4);
}




void getMatrices(std::vector<FMatrix> &mats, const FMatrix& rot, const int l)
{
  eassert( l>= 0);
  FMatrix mat0(1,1);
  mat0(0,0) = 1.0;
  mats.push_back( mat0 );

//  FMatrix mat1(3,3);
  mats.push_back( rot );

  for(int i=2; i<=l; ++i)
  {
    int size= 2*i+1;
    FMatrix mat(size,size);
    /*const*/ FMatrix &m_=mats.back();
//    std::cout << "last sizes: " << m_.getDimensionX() << " " << m_.getDimensionY() << std::endl;

    // std::cout << "back dimension: " << m_.getDimensionX() << std::endl;
    eassert( (int)m_.getDimensionX() == 2*(i-1)+1);
    eassert( (int)m_.getDimensionY() == 2*(i-1)+1);

    for(int x=0;x<size; ++x)
    {
      for(int y=0; y<size;++y)
      {
        int m=x-i;
        int n=y-i;
//      std::cout << "m,n: " << m << " " << n << std::endl;
        // eq. (8.1)
        double u = sub_u(i,m,n);
        double v = sub_v(i,m,n);
        double w = sub_w(i,m,n);
//        if( w==0 ) std::cout << "wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww=0" << std::endl;
//        if( u==0 ) std::cout << "uuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuuu=0" << std::endl;
//        if( v==0 ) std::cout << "vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv=0" << std::endl;

//        std::cout << u << " " << v << " " << w << std::endl;
        if( u != 0)
        {
          double U = sub_U(i,m,n,m_,rot);
          mat(x,y) += u * U;
//          std::cout << "U=" << U << std::endl;
        }
        if ( v != 0)
        {
//          std::cout << "sub_V" << std::endl;
          double V = sub_V(i,m,n,m_,rot);
//          std::cout << "-sub_V" << std::endl;
          mat(x,y) += v * V;
//          std::cout << "V=" << V << std::endl;
        }
        if( m != 0 && w != 0)
        {
//          std::cout << " w = " << w << std::endl;
          double W = sub_W(i,m,n,m_,rot);
//          std::cout << "W=" << W << std::endl;
          mat(x,y) += w * W; 
        }
      }
    }
    mats.push_back( mat );
  }
}

void getRotationMatrix( FMatrix& mat, const FMatrix& rotation, const int L)
{
  std::vector<FMatrix> mats;
  getMatrices( mats, rotation, L );
  eassert((int)mats.size() == L+1);
  int size = L*(L+2) +1;
  std::cout << size << std::endl;
  mat.resize( size, size, false );
  int off =0;
  for(int l=0;l<=L; ++l)
  {
    //std::cout << l << "/" << L << " off=" << off << " size=" << mats[l].getDimensionX() << " " << mats[l].getDimensionY() << std::endl;
    mat.setSubMatrix(off, off, mats[l]);
    off += mats[l].getDimensionX();
  }
}

void getROTMatrix( FMatrix &rot, double alpha, double beta, double gamma )
{
  rot(0,0) = cos(alpha)*cos(gamma)-sin(alpha)*sin(gamma)*cos(beta);
  rot(0,1) = sin(alpha)*sin(beta);
  rot(0,2) = cos(alpha)*sin(gamma)+sin(alpha)*cos(gamma)*cos(beta);
  rot(1,0) = sin(gamma)*sin(beta);
  rot(1,1) = cos(beta);
  rot(1,2) = -cos(gamma)*sin(beta);
  rot(2,0) = -cos(alpha)*sin(gamma)*cos(beta)-sin(alpha)*cos(gamma);
  rot(2,1) = cos(alpha)*sin(beta);
  rot(2,2) = cos(alpha)*cos(gamma)*cos(beta)-sin(alpha)*sin(gamma);
}

} // namespace FMath

FSphericalHarmonics::~FSphericalHarmonics()
{
  if(shrot) delete shrot;
  shrot = 0;
};


