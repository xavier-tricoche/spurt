//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTensor.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:38:16 $
// Author:    $Author: garth $
// Version:   $Revision: 1.14 $
//
//--------------------------------------------------------------------------- 

#include "FTensor.hh"
#include <iostream>
#include "binio.hh"

#include <cmath>

//---------------------------------------------------------------------------

FTensor::~FTensor()
{
}

//---------------------------------------------------------------------------

std::ostream& operator<< (std::ostream& os, const FTensor& t)
{
  os <<"[ ";

  if( t.comp == 0 )
      os << " empty ";
  else
  for (unsigned char i=0; i<(unsigned char)FTensor::pow (t.dimension,
							 t.order); i++)
    os << t.comp[i] << " ";

  os <<" ]";
  return os;
}

//---------------------------------------------------------------------------

std::istream& binread( std::istream& in, FTensor &t )
{
  unsigned int dimension;
  binread_raw( in, &dimension );

  unsigned int order;
  binread_raw( in, &order );

  unsigned int size;
  binread_raw( in, &size );

  std::vector<double> values;
  values.resize( size );

  for( unsigned int i = 0; i < size; ++i )
    {
      binread_raw( in, &values[i] );
    }
  
  FTensor tmp( dimension, order, values );
  t = tmp;

  return in;
}

std::ostream& binwrite( std::ostream& out, const FTensor& t )
{
  unsigned int dimension = t.getDimension();
  unsigned int order = t.getOrder();
  std::vector<double> values;
  t.getValues( values );
  unsigned int size = values.size();
  
  binwrite_raw( out, &dimension );
  binwrite_raw( out, &order );
  binwrite_raw( out, &size );

  for( unsigned int i = 0; i < size; ++i )
  {
      binwrite_raw( out, &values[i] );
  }

  return out;

}

//---------------------------------------------------------------------------

void FTensor::getEigenSystem2(FVector& vals, FVector v[3])
{
#ifndef NODEBUG
  if (order != 2)
    THROW_EXCEPTION( FInvalidDimensionException, "ERROR: invalid order (has to be 2)!");
  
  if (dimension != 3)
    THROW_EXCEPTION( FInvalidDimensionException, "ERROR: invalid dimension (has to be 3)!");
#endif
#define T ( *this )

  vals.resize( 3 );
  double I1 = T( 0, 0 ) + T( 1, 1 ) + T( 2, 2 );
  double I2 = T( 0, 0 ) * T( 1, 1 ) 
            + T( 1, 1 ) * T( 2, 2 )
            + T( 2, 2 ) * T( 0, 0 )
            - ( T( 0,1 )*T( 0,1 )+T( 0,2 )*T( 0,2 )+T( 1,2 )*T( 1,2 ) );
  double I3 = T( 0,0 )*T( 1,1 )*T( 2,2 )
            + 2.*T( 0,1 ) *T( 1,2 )*T( 2,0 )
            -( T( 0,0 )*T( 1,2 ) + T( 1,1 )*T( 2,0 )+T( 2,2 )*T( 0,1 ) );

  double V = ( I1/3. )*( I1/3. ) - I2/3.;
  double S = ( I1/3. )*( I1/3. )*( I1/3. )-I1*I2/6.+I3/2.;
  double phi= acos( S/V*sqrt( 1./V ) )/3.;

  vals[ 0 ] = I1/3.+2.*sqrt( V )*cos( phi );
  vals[ 1 ] = I1/3.-2.*sqrt( V )*cos( M_PI/3.+phi );
  vals[ 2 ] = I1/3.-2.*sqrt( V )*cos( M_PI/3.-phi );
  for ( int i=0; i<2; ++i )
  {
    double A = T( 0,0 )- vals[ i ];
    double B = T( 1,1 )- vals[ i ];
    double C = T( 2,2 )- vals[ i ];

    v[ i ] = FVector( 3 );
    v[ i ][ 0 ] = ( T( 0,1 )*T( 1,2 )-B*T( 0,2 ) ) * ( T( 0,2 )*T( 1,2 )-C*T( 0,1 ) ); // FIXME!!!
    v[ i ][ 1 ] = ( T( 0,2 )*T( 1,2 )-C*T( 0,2 ) ) * ( T( 0,2 )*T( 0,1 )-A*T( 1,2 ) );
  }
  v[ 2 ] = FVector( 3 );
  v[ 2 ] = v[ 0 ].crossProductConst( v[ 1 ] );
#undef T
}

//===========================================================================

#ifdef OUTLINE
#include "FTensor.icc"
#endif

//---------------------------------------------------------------------------
