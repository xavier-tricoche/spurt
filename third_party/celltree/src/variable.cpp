#include <netcdf.h>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <cstdio>

#include "variable.hpp"

inline static void nccheck( int result )
{
    if( result != NC_NOERR )
        throw std::runtime_error( nc_strerror(result) );
}
    
// -------------------------------------------------------------------------
