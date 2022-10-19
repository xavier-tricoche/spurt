//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTensorSet.cc,v $
// Language:  C++
// Date:      $Date: 2004/05/19 12:08:03 $
// Author:    $Author: hlawit $
// Version:   $Revision: 1.8 $
//
//--------------------------------------------------------------------------- 

#include "FTensorSet.hh"
#include "math.h"

#include "eassert.hh"
#include <cassert>

namespace {
inline positive power(positive d,positive p)
{
  positive ret=1;
  for(;p!=0;--p)
    ret*=d;
  return ret;
}
}

//---------------------------------------------------------------------------

using namespace std;

//---------------------------------------------------------------------------

FTensorSet::~FTensorSet()
{
}

//---------------------------------------------------------------------------

const FString& FTensorSet::getClassName() const
{
    static FString ugly = "FTensorSet";
    return ugly;
}

//---------------------------------------------------------------------------
FTensorSet::FTensorSet(positive dimension, positive order, 
		       boost::shared_ptr< FanyArray<double> >  ldata,
		       const std::string& aName
		       )
  :name(aName),
   data(ldata),
   tensorOrder(order),
   tensorDimension(dimension),
   compSize(power(dimension,order))
{
    // sanity check memorial

    bufvector.resize(compSize);

    if ( data->size() % compSize )
	THROW_EXCEPTION(FException,
			" given array size was not "
			"a multiple of the tensor size" );
}

//---------------------------------------------------------------------------
FTensorSet::FTensorSet( positive dimension, positive order, 
			vector<double>& ldata,
			const std::string& aName )
  : name(aName),    
    data(new FanyVector<double> (ldata) ), // this does an internal swap, so ldata is deleted
    tensorOrder(order),
    tensorDimension(dimension),
    compSize(power(dimension,order))
{
    // sanity check memorial

    bufvector.resize(compSize);

    if ( data->size() % compSize )
	THROW_EXCEPTION(FException,
			" given array size was not "
			"a multiple of the tensor size" );


}

//---------------------------------------------------------------------------

FTensorSet::FTensorSet( positive dimension, positive order, 
			const std::vector<FTensor>& tensors,
			const std::string& aName )
  : name(aName)
  , tensorOrder( order )
  , tensorDimension( dimension )
  , compSize( power(dimension,order))
{
    // sanity check memorial

    bufvector.resize(compSize);

    data.reset(new FanyVector<double>( compSize*tensors.size()) );

    FAssignableAnyArray<double>::iterator d = boost::shared_dynamic_cast< FanyVector<double> >(data)->begin();

    for( vector<FTensor>::const_iterator i=tensors.begin(); 
	 i!=tensors.end(); ++i )
    {
	const FArray &a = (*i);

	assert( i->getDimension() == dimension );
	assert( i->getOrder() == order );
	assert( a.size() == compSize );

	for( unsigned int j=0; j<compSize; j++, d++ )
	    *d = a[j];
    }
}

//---------------------------------------------------------------------------

FTensorSet::FTensorSet( const FTensorSet &ts, const std::string& aName )
  : FObject(*this),
    name(aName),
    data( ts.data ), // this does an internal swap, so ldata is deleted
    tensorOrder(ts.tensorOrder),
    tensorDimension(ts.tensorDimension),
    compSize(ts.compSize)
{
}

FTensorSet::FTensorSet( positive dimension, positive order, const FTensorSet &ts, const std::string& aName )
  : FObject(*this),
    name(aName),    
    data(ts.data ),
    tensorOrder(order),
    tensorDimension(dimension),
    compSize(power(dimension,order))
{
  // share data with changed attributes
  // at least the size of the tensors should be the same,
  // otherwise the grid changed as well
  assert( compSize == ts.compSize );
}

//---------------------------------------------------------------------------

positive FTensorSet::getDimension() const
{
    return tensorDimension;
}

//---------------------------------------------------------------------------

positive FTensorSet::getOrder() const
{
    return tensorOrder;
}

//---------------------------------------------------------------------------

positive FTensorSet::getNbTensorComponents() const
{
    return compSize;
}

//---------------------------------------------------------------------------

const string& FTensorSet::getName() const
{
    return name;
}

//---------------------------------------------------------------------------

positive FTensorSet::getNbTensors() const
{
   eassert( compSize != 0 );
   return ((positive) data->size()) / compSize;
}

//---------------------------------------------------------------------------

void FTensorSet::getTensor(FTensor& result, const FIndex& tensorId) const
{
    positive id = getNbTensorComponents()*tensorId.getIndex();

#ifndef NODEBUG
    if (id >= (positive) data->size()) 
    {
      THROW_DEFAULT_EXCEPTION( FIndexOutOfBoundsException );
    }
#endif

    result.resizeTensor(tensorDimension,tensorOrder);

    // this version uses the data array bufvector that prevents this
    // part of the code to be reentrant
    // In addition to that, it is a copy operation more than we
    // should need
    data->get_range( id, compSize ,&bufvector[0] );
    result.setValues(bufvector);
#if 0
    // I'd prefer it this way, but it seems that we do not ensure that
    // tensors lie in the same block when storing the data e.g. in an
    // mmapped array, so we can not guarantee this to work in all
    // cases
    double *first = &(*data)[id];
    result.setValues( first );
#endif
}

//---------------------------------------------------------------------------

positive FTensorSet::memSize() const
{
    positive memsize = data->memSize() + sizeof(FTensorSet);
    memsize += bufvector.size()*sizeof(double);
    return memsize;
}


//===========================================================================
