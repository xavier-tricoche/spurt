//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTensorSet.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:12 $
// Author:    $Author: garth $
// Version:   $Revision: 1.11 $
//
//--------------------------------------------------------------------------- 

#ifndef __FTensorSet_hh
#define __FTensorSet_hh

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>
#include "stdAliases.hh"
#include "FObject.hh"
#include "FString.hh"
#include "FException.hh"
#include "FTensor.hh"
#include "FIndex.hh"
#include "FRefTensor.hh"
#include "FanyArray.hh"

//===========================================================================


/**
 * FTensorSet is a set of mathematical tensors. The data is stored in
 * a simple array of doubles. If one asks for a tensor the data is
 * linked to the returnTensor which is then given back to the user.
 *
 *   \ingroup DataSet
 */
class FTensorSet : public FObject
{
public:

  /** 
   *\par Description:
   *Constructor: provides a tensors set of dimension \b dim and order 
   *\b order. A name can also be given, namely \b name. Default value
   *for this parameter is "noName". 
   *\pre
   *order <= 2, dimension = {2,3}
   *\post
   * data is empty
   *\exception
   *FInvalidDimensionException
   *\param data:
   * components of the tensors
   *\param
   *dimension: dimension of the tensors.
   *\param
   *order: order of the tensors.
   *\param
   *aName: name of the tensors set.
   */
  FTensorSet(positive dimension, positive order, 
	     vector<double> & ldata,
	     const std::string& aName = "noName");

  FTensorSet(positive dimension, positive order, 
	     boost::shared_ptr< FanyArray<double> >  ldata,
	     const std::string& aName = "noName" );

  FTensorSet( positive dimension, positive order, 
	      const std::vector<FTensor>& tensors,
	      const std::string& aName = "noName" );


  // copy constructors that share an internal pointer to the data array
  FTensorSet( const FTensorSet &ts, const std::string& name = "noName" );
  FTensorSet( positive dimension, positive order, const FTensorSet &ts, 
      const std::string& name = "noName" );
  
  /** 
   *\par Description:
   *Destructor.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  ~FTensorSet();

  /** 
   *\par Description:
   *Gets the class name.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  const FString& getClassName() const; 

  /** 
   *\par Description:
   *Gives the dimension of the tensors set.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  positive getDimension() const;

  /** 
   *\par Description:
   *Gives the order of the tensors set.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  positive getOrder() const;

   /** 
   *\par Description:
   *Return the number of tensor components.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  positive getNbTensorComponents() const;

  /** 
   *\par Description:
   *Gives the number of allocated tensors in the tensors set.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  positive getNbTensors() const;

  /** 
   *\par Description:
   *Gives the name of the tensors set.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  const std::string& getName() const; 

  /** 
   *\par Description:
   *Returns the tensor with index \b tensorId.
   *\pre
   *tensorId.getIndex() < getNbTensors()
   *\post
   *none
   *\exception
   *FIndexOutOfBoundsException
   *\param
   *result: returned tensor value.
   *\param
   *tensorId: index of the requested tensor.
   */
  void getTensor(FTensor& result, const FIndex& tensorId) const;

  /**
   * Access to the raw data.
   * This should only be used in code that is bound close to the
   * internal data as too much knowledge of the internal data
   * structure (at least more than good softwar engineering allows)
   * has to be known here.
   *
   * Use in code that needs speed, only and try to mark it clearly
   * as such.
   */
  boost::shared_ptr<const FanyArray<double> > getRawData() const
  {
    return data;
  }
  
  virtual positive memSize() const;

private:

  // Name of the FTensorSet.
  std::string name;

  // Tensor data
  boost::shared_ptr<FanyArray<double> > data;

  positive tensorOrder, tensorDimension;
  positive compSize;

  //vector for the buffering of one tensor entry,
  //used in getTensor
  mutable vector<double> bufvector;
};

#endif // __FTensorSet_hh
