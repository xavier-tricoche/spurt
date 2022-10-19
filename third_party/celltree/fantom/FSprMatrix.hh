//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSprMatrix.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:54:49 $
// Author:    $Author: garth $
// Version:   $Revision: 1.12 $
//
//--------------------------------------------------------------------------- 


//---------------------------------------------------------------------
// Sorry if some things aren't as effective as they could be
// but for reasons of time effort and error minimizing this is 
// taken as found in Numerical Recipes Chapter 2.7
//
//
// --------------------------------------------------------------------
// IMPORTANT NOTE :::::::::::::::::::::::::::::::::::::::::::::::::::::
// --------------------------------------------------------------------
//
// this class is not foolproof, so read the interface description !!!
//
//---------------------------------------------------------------------

#ifndef __FSprMatrix_hh
#define __FSprMatrix_hh

#include <map>
#include <vector>
#include <iosfwd>

#include "FException.hh"

class FArray;
typedef FArray FVector;
class FMatrix;


//===========================================================================

/**
 * The FSprMatrix class provides an implementation of an MxM Square Sparse
 * Matrix. Since it is incompatible with the methods operating on regular
 * matrices, an iterative solver for linear systems is included, using the
 * special nature of the matrix representation.
 * {\bf IMPORTANT NOTE}
 * this class is not foolproof, so read the interface description !!!
 * (especially the efficient fill)
 * --------------------------------------------------------------------
 */

class FSprMatrix
{
private:

  /** Default constructor.
   */
  FSprMatrix();

public:

  /** Copy constructor.
   * \param m
   *  matrix to copy.
   */
  FSprMatrix(const FSprMatrix& m);

  /** Copy constructor (uses Matrix as input).
   * \param m
   *   matrix to copy.
   * \param tresh
   *   components whose absolute value is less than this treshold are
   *   are considered zero.
   */
  FSprMatrix(const FMatrix& m, double thresh);

  /**
   * \par Description:
   *   To-Use-Constructor (use THIS one to build a sparse Matrix).
   * \param m
   *   Dimension of the square Matrix.
   * \param fill_estimate
   *   By giving this value between 0 and 1, you can speed up the
   *   construction of the matrix, since memory is allocated in advance
   *   according to your estimate. The better your estimate, the fewer
   *   reallocations will be necessary.
   */
  FSprMatrix(unsigned int m, double fill_estimate);

  /** Destructor.
   * \pre
   * Matrix exists.
   * \post
   * Matrix no longer exists.
   */
  ~FSprMatrix();

  /**
   * \par Description:
   *   Resets the matrix by deallocating the used space and redefining
   *   dimension.
   * \param m
   *   Dimension of the square Matrix.
   */
  void _reset_ (unsigned int m);
  
  /**
   * \par Description:
   *   Scalar multiplication (binary).
   * \param lambda
   *   Scalar to multiply the matrix with.
   */
  FSprMatrix operator*(const double lambda);

  /**
   * \par Description:
   *   Scalar multiplication (unary).
   * \param lambda
   *   Scalar to multiply.
   */
  FSprMatrix& operator*=(const double lambda);

  /**
   * \par Description:
   *   Matrix_by_Vector Multiplication 
   * \pre
   *   Vector dimension has to be equal to matrix-dimension.
   * \param v
   *   Right hand side operand.
   * \exception
   *   FInvalidDimensionException
   */
  FArray operator*(const FVector &v);

  /**
   * \par Description:
   *   transpose_of_Matrix by Vector Multiplication 
   * \pre
   *   vector dimension has to be equal to matrix-dimension.
   * \param v
   *   Right hand side operand.
   * \exception
   *   FInvalidDimensionException
   */
  FArray tmult(FVector &v);

  /**
   * \par Description:
   *   Sets the minimum absolute value entries have to be to be
   *   considered non-zero. The absolute value of every entry that
   *   is added into the matrix will first be compared to this
   *   treshold and ignored if less.
   *   !!!! does not affect the existing object !!!!
   * \param t
   *   Treshold value.
   */
  void setTreshold(double t);

  /**
   * \par Description:
   *   Transposition of a matrix
   * \pre
   *   Matrix is valid.
   * \post
   *   Matrix is transposed, the initial matrix does not exist any more.
   * \return
   *   Reference to transposed self
   */
  FSprMatrix& transpose(void);

  /**
   * \par Description:
   *   Substitute for Index Operator
   * \return
   *   Element m,n points to
   */
  double getElement(unsigned int m, unsigned int n) const;

#if 0
  /**
   * \par Description:
   *   Substitute for Index Operator.
   *   Vector has to be filled top-left --> right-down
   *   otherwise throws exception
   *  !!!!!!!! use only this OR only cueElement() !!!!!!!!!!
   * \param m,n
   *   Indiex of element to access.
   * \param val
   *   value to insert
   * \exception FException
   *   Thrown if the required fill order is violated.
   * \return
   *   The given value is returned, no other sense.
   *
   *
   *
   *   This is currently not working!
   */
  double setElement(unsigned int m, unsigned int n, double val);
#endif

  /**
   * \par Description:
   *   Efficiently fill matrix by collecting the entries to put them
   *   all at once. Once you begin using cueElement(), you must insert
   *   all matrix elements this way and finish by calling commitCue().
   *   You may insert the elements in any order. Until commitCue() is
   *   called, no other operations on the matrix are allowed.
   *   (exceptions are thrown)
   * \param
   *   {\tt m,n} index of element to add, {\tt val} value to insert
   */
  void cueElement (unsigned int m, unsigned int n, double val);

  /**
   * \par Description:
   *   Generates sparse form of data previously collected by cueElement(...)
   * \pre
   *   The matrix elements have been defined via cueElement().
   * \post
   *   The sparse form of the matrix has been generated.
  */
  void commitCue ( void );

  /**
   * \par Description:
   *   Solves the matrix for the given vector v using the biconjugate
   *   gradient method.
   * \param v
   *   Vector to solve this matrix for.
   * \param itol
   *   Value between 1 and 4, defines the convergence criteria used.
   *   Sorry, consult the source to find out which effect the value have.
   * \param tol
   * \param itmax
   *   The maximum number of iterations that should be taken to reach
   *   a good result.
   * \param iter
   * \param err
   * \return
   *   copy of FArray holding the result
   */
  FArray solve_bcg( FVector & v, unsigned int itol, double tol, unsigned int itmax, unsigned int & iter, double & err);

  /**
   * \par Description:
   *    returns the allocated space (number of array-entries)
   */
  unsigned int space ( void ) const;

  /**
   * \par Description:
   *    returns the dimension of the Sparce matrix (just to be compatible
   *    with the other matrix classes...
   */
  unsigned int getDimension(void);

  /**
   * \par Description:
   *    No need to explain that... no aligning is done and be careful
   *    when trying it on 10000x10000 matrices... *g*
   */
  friend std::ostream& operator<< (std::ostream& os, const FSprMatrix& matrix);
  
private:
  
  void asolve(FVector& b,FVector& x); // used by solve_bcg
  double snrm(FVector& sx, unsigned int itol); // used by solve_bcg

  std::map<std::pair<unsigned int, unsigned int>, double>* _map;  //! used for efficient fillment of the matrix

  double EPS; //! the treshold for nonzero-determination
  unsigned int dimension; //! the dimension of the matrix
  std::vector<unsigned int> ija; //! the array holding indices of nonzero matrix entries
  std::vector<double> sa; //! the entries themselves
  
  unsigned int fill_m, fill_n; //! used to determine if the matrix is filled in the right way
  
};

#endif // __FSprMatrix_hh
