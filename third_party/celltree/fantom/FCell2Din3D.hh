//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCell2Din3D.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:01 $
// Author:    $Author: garth $
// Version:   $Revision: 1.14 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCell2Din3D_hh
#define __FCell2Din3D_hh

#include "FObject.hh"
#include "stdAliases.hh"
#include "FException.hh"
#include "FString.hh"
#include "FIndex.hh"
#include "FPosition.hh"
#include "FPositionSet.hh"
#include "FTensor.hh"
#include "FTensorSet.hh"
#include "FCell.hh"

#include <list>

#include <utility>


//===========================================================================

/**
 *The FCell2Din3D class provides an abstract 2Din3D geometric cell lying in a 3D space.
 */
class FCell2Din3D : public FCell
{
public:

  /** 
   *\par Description:
   *Constructor
   */
  FCell2Din3D (const geoInfo&g,
	       FIndex*const i,
	       FRefArray*const p,
	       double *const pd,
	       FRefTensor*const t) : FCell(g,i,p,pd,t){ };
  
  /** 
   *\par Description:
   *Destructor
   */
  ~FCell2Din3D() { };
  
  /**
   *\par Description:
   * returns the cell dimension
   */
  virtual positive getDimension() const;
  
};

#endif // __FCell2Din3D_hh

