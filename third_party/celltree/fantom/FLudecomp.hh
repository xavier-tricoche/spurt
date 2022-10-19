//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FLudecomp.hh,v $
// Language:  C++
// Date:      $Date: 2003/02/06 10:22:32 $
// Author:    $Author: garth $
// Version:   $Revision: 1.6 $
//
//--------------------------------------------------------------------------- 

#ifndef __FLudecomp_hh
#define __FLudecomp_hh

#include "FMatrix.hh"

//===========================================================================

/** 
 * The FLudecomp class provides an implementation of Ludecomp algorithm.
 */

class FLudecomp : public FMatrix
{
private:
  /** Default constructor. (shouldn´t be used...)
   */
  FLudecomp();

public:
  /** Copy constructor.
   * \\\par Description:
   *   Constructs a matrix with the same dimension and components as
   *   the given argument.
   * \param
   *   {\tt matrix} matrix to copy.
   */
  FLudecomp(const FMatrix& matrix);

  /** Destructor
   */
  ~FLudecomp();

  /** Solver.
   * \\\par Description:
   *   Solves the previously given Matrix for the FVector in the argument
   * \param
   *   {\tt d} FVector to solve the matrix for (holds the result afterwards).
   */
  void solve(FVector& d) const ;

private:
  
  positive* indx;   // used to store the row permutation
  double d;      // ...becomes 1 if equal number of Row interchanges, -1 else
};

//===========================================================================

#ifndef OUTLINE
#include "FLudecomp.icc"
#endif

#endif // __FLudecomp_hh
