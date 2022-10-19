//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FObject.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:24:52 $
// Author:    $Author: garth $
// Version:   $Revision: 1.8 $
//
//--------------------------------------------------------------------------- 

#ifndef __FObject_hh
#define __FObject_hh

#include "FString.hh"

using namespace std;

//===========================================================================

/** 
 * The FObject class is the abstract base class of all other FAnToM classes. 
 */

class FObject
{
protected:
 
  /** default constructor
   */
  FObject();


  /** copy constructor
   */
  FObject(const FObject &object);

public:

  /** destructor
   */
  virtual ~FObject();

  /** Returns the class name as string.
   *  \return class name   
   */
  virtual const FString& getClassName() const;

  virtual FObject* Clone() const;
};

//===========================================================================
#ifndef OUTLINE
#include "FObject.icc"
#endif
//===========================================================================

#endif // __FObject_hh
 
