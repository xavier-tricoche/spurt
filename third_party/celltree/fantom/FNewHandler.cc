//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FNewHandler.cc,v $
// Language:  C++
// Date:      $Date: 2003/03/27 08:41:31 $
// Author:    $Author: garth $
// Version:   $Revision: 1.4 $
//
//--------------------------------------------------------------------------- 

#include "FNewHandler.hh"
#include "FException.hh"

//---------------------------------------------------------------------------

void FNewHandler::new_handler()
{
  throw FOutOfMemoryException();
}

//===========================================================================
