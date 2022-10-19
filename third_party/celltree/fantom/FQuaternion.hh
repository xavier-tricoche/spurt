//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FQuaternion.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:54:48 $
// Author:    $Author: garth $
// Version:   $Revision: 1.8 $
//
//--------------------------------------------------------------------------- 

#ifndef __FQuaternion_hh
#define __FQuaternion_hh

#include "FObject.hh"
#include "FMatrix.hh"

//===========================================================================

/** 
 * The FQuaternion class provides an implementation of an 
 * 4-dimensional Quaternion
 */
class FQuaternion : public FObject
{
public:

  /** Constructs Quaternion {\tt n}.
   * \post
   *   Components being initialized to zero.
   */
  FQuaternion();

  /** Constructor
   * \post
   *   Components being initialized to value.
   */
  FQuaternion (const double &value);

  /** Copy constructor.
   * \\\par Description:
   *   Constructs a Quaternion with the same values as the given argument.
   * \param
   *   {\tt vector} vector to copy.
   */
  FQuaternion (const FVector&);

  /** Copy constructor.
   * \\\par Description:
   *   Constructs a Quaternion with the same values as the given argument.
   * \param
   *   {\tt quaternion} Quaternion to copy.
   */
  FQuaternion (const FQuaternion&);

  /** Destructor.
   */
  ~FQuaternion();

  /** getClassName
   */
  const FString& getClassName () const;

  /// Undocumented.
  const double& operator[](unsigned int c) const;
  /// Undocumented.
  double& operator[] (unsigned int c);

  friend std::ostream& operator<<(std::ostream& os, const FQuaternion& a);

  /** Quaternion multiplication.
   * \\\par Description:
   *   Multiplication of this quaternion with qR.
   * \param
   * Right quaternion qR.
   * \param
   * Destination.
   */
  void mult(const FQuaternion &qR, FQuaternion& dest);

  /**
   * \brief
   * Conjugate quaternion.
   * \pre
   * None.
   * \post
   * None.
   * \param
   * Quaternion q to store conjugate of this in.
   * \exception
   * None.
   */
  void conj(FQuaternion& q) const;

  /**
   * \brief
   * Construct rotation matrix from (possibly non-unit) quaternion.
   * Assumes matrix is used to multiply column vector on the left:
   * vNew = mat x vOld.toMatri Works correctly for right-handed coordinate 
   * system and right-handed rotations.
   * \pre
   * None.
   * \post
   * None.
   * \exception
   * None.
   * \param
   * Matrix m to store resultant rotation matrix in.
   */
  void toMatrix(FMatrix& m);
  
  /**
   * \brief
   * Constant unit == (0,0,0,1) quaternion.
   * \pre
   * None.
   * \post
   * None.
   * \return
   * Unit quaternion.
   */
  static const FQuaternion qOne ();

private:
  
  double comp[4];
};

//===========================================================================
#ifndef OUTLINE
#include "FQuaternion.icc"
#endif
//=========================================================================== 

#endif // __FQuaternion_hh
