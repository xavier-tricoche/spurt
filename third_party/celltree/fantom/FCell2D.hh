//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCell2D.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:01 $
// Author:    $Author: garth $
// Version:   $Revision: 1.12 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCell2D_hh
#define __FCell2D_hh

#include "FObject.hh"
#include "FIndex.hh"
#include "FCell.hh"
#include "FAMSingularPoint.hh"
#include "stdAliases.hh"

#include "FBoundingBox.hh"

#include <list>
#include <vector>
#include <utility>


//class defined in FCell2D.cc
//where edgeToGo stores its data
class FCellEdgeInfo;

//===========================================================================

/**
 *The FCell2D class provides an abstract 2D geometric cell lying in a 2D space.
 */
class FCell2D : public FCell
{
public:
  /** 
   *\par Description:
   *Constructor
   */
  FCell2D (const geoInfo& g, FIndex* const i, FRefArray* const p, 
	   double * const pd, FRefTensor* const t) 
    : FCell(g,i,p,pd,t)
  { edgeInfo=0; }
  

  /**
   * help function for cell locator:
   * get edge where line through start
   * in direction dir exits the cell
   *\pre
   * start is inside cell or the ray described by start,dir
   * intersects cell
   * \param start
   * start position in actual cell
   * \param dir
   * direction in which line goes
   *\retval length
   * parameter for line : start + length * dir = cutpoint with face
   *\return
   * Id of next edge where line leaves the cell
   */
  FIndex edgeToGo(double&length,
		  const FArray& start, 
		  const FArray& dir) const;

  /** 
   *\par Description:
   *Destructor
   */
  virtual ~FCell2D ();

protected:

  mutable FCellEdgeInfo * edgeInfo;

};

#endif // __FCell2D_hh

