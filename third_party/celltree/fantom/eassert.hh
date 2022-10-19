//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: eassert.hh,v $
// Language:  C++
// Date:      $Date: 2003/09/16 07:06:32 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __eassert_hh
#define __eassert_hh

#include <stdexcept>
#include <cstdio>

//--------------------------------------------------------------------------- 

const char* __make_assertion_msg( const char* e, const char* f, 
				  unsigned int l );


class assertion_failure: public std::logic_error
{
public:
    assertion_failure( const char* e, const char* f, unsigned int l ) :
	std::logic_error( __make_assertion_msg( e, f, l ) )
    {
    }

    ~assertion_failure() throw() {}
};

//--------------------------------------------------------------------------- 

#ifdef WIN32
#define eassert(expr) \
    ((void) ((expr) ? 0 : throw assertion_failure( #expr, __FILE__, __LINE__ )))
#else

#define eassert(expr)  \
    ((void) ((expr) ? 0 : throw assertion_failure( __STRING(expr), \
                                                   __FILE__, __LINE__ )))
#endif
#endif //__eassert_hh
