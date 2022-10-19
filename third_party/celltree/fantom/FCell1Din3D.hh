//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCell1Din3D.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:01 $
// Author:    $Author: garth $
// Version:   $Revision: 1.12 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCell1Din3D_hh
#define __FCell1Din3D_hh

#include "FObject.hh"
#include "FIndex.hh"
#include "FCell.hh"
#include "stdAliases.hh"

#include <list>
#include <vector>
#include <utility>


//===========================================================================

/**
 *The FCell1Din3D class provides an abstract 1D geometric cell lying in a 1D space.
 */
class FCell1Din3D : public FCell
{
public:
  /** 
   *\par Description:
   *Constructor
   */
  FCell1Din3D (const geoInfo& g, FIndex* const i, FRefArray* const p, 
	   double * const pd, FRefTensor* const t) 
    : FCell(g,i,p,pd,t)
  {  }
  
  /** 
   *\par Description:
   *Destructor
   */
  virtual ~FCell1Din3D ();

};

#endif // __FCell1Din3D_hh

