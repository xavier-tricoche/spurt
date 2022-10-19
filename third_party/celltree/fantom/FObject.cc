//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FObject.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:24:52 $
// Author:    $Author: garth $
// Version:   $Revision: 1.5 $
//
//--------------------------------------------------------------------------- 

#include "FObject.hh"
#include "FString.hh"

#ifdef OUTLINE
#include "FObject.icc"
#endif

//---------------------------------------------------------------------------

FObject::~FObject()
{
#ifdef FOBJECT_DEBUG
  FObjectRegistry::registry.unregisterObject(this);
#endif
}

//---------------------------------------------------------------------------

const FString& FObject::getClassName() const
{
  static const FString className("FObject");

  return className;
}

//===========================================================================
