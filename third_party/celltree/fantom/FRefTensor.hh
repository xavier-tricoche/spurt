//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FRefTensor.hh,v $
// Language:  C++
// Date:      $Date: 2003/03/27 08:29:52 $
// Author:    $Author: garth $
// Version:   $Revision: 1.5 $
//
//--------------------------------------------------------------------------- 

#ifndef __FRefTensor_hh
#define __FRefTensor_hh

#include "FException.hh"
#include "FString.hh"
#include "FTensor.hh"
#include "FIndex.hh"
#include <vector>


/**
 *FRefTensor is an exact copy of a FTensor
 *(not the components, but the pointer on them is copied)
 */
class FRefTensor:public FTensor
{
public:
  /**
   *Default constructor
   */
  FRefTensor();

  /**
   *Copy constructor
   *doesn't copy the components of the Tensor, 
   *but only the pointer on them
   */
  FRefTensor(FTensor&);

  /**
   *Constructor that initializes the Fields directly
   *\param order
   *Order of the referenced tensor  
   *\param dimension 
   *Dimension of the referenced tensor
   *\param Comp
   *pointer on the double array
   *holding the Components of the referenced Tensor
   */
  FRefTensor(int Dimension,int Order,double*Comp);
  
  /**
   *sets the pointer on a new set in FTensorSet
   *\param Comp: pointer on a valid array of doubles(or NIL)
   *\pre 
   *Comp has at least Dimension^Order entries
   *\post
   *none
   *\exception
   *none
   */
  void setPointerOnComponents(double*Comp);


  /**
   *sets the order,dimension and pointer of the tensor   
   *\pre 
   *Comp has at least pow(Dimension,Order) entries
   *\post
   *none
   *\exception
   *none
   *
   *\param Dimension:
   * new dimension of the tensor
   *
   *\param Order:
   * new order of the tensor
   *
   *\param Comp:
   * pointer on array with at least pow(Dimension,Order) entries
   *
   *\param Size:
   * pow(Dimension,Order)
   */
  void setDimensionOrderPointerAndSize(int Dimension, int Order,double*Comp,int Size=0);


  /** 
   *\par Description:
   *sets the values of the referenced tensor to the ones in T.
   *\pre
   *none
   *\post
   *none
   *\param
   *T: Tensor from where values are copied
   */
  FTensor& operator = (const FTensor&T);

  /**
   *Destructor
   */
  ~FRefTensor();
};

#ifndef OUTLINE
#include "FRefTensor.icc"
#endif

#endif // __FRefTensor_hh
