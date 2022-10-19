//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: eassert.cc,v $
// Language:  C++
// Date:      $Date: 2003/09/16 07:06:31 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#include "eassert.hh"

#include <cstring>
#include <cstdlib>

const char *__make_assertion_msg( const char *expr,
				  const char *file,
				  unsigned int line )
{
    static char msg[512];
    
    snprintf( msg, 512, "\"%s\" at %s, line %u",
	      expr, file, line );

    return msg;
}

