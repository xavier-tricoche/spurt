/*************************************************************************
* celltree: boundary interval hierarchy tree-based solution for point 
*           location and interpolation in unstructured meshes. This
*           technique was first proposed in:
*           "Fast, Memory-Efficient Cell location in Unstructured Grids
*            for Visualization"
*                by Christop Garth and Kenneth I Joy, VisWeek 2011.
* Author: Christoph Garth
* Copyright (C) 2011 Christoph Garth
* All rights reserved
**************************************************************************/

#ifndef __netcdf_util_hpp
#define __netcdf_util_hpp

// #define NETCDF_CPP_INTERFACE

#ifndef NETCDF_CPP_INTERFACE
#include <netcdf.h>
#else
#include <netcdfcpp.h>
#endif


namespace netcdf_util {
    
struct netcdf_exception: public std::exception
{
    int error;
    
    netcdf_exception( int e ) : error(e)
    {
    }
    
    const char* what() const throw()
    {
        return nc_strerror(error);
    }
};

static void nccheck( int result )
{
    if( result != NC_NOERR )
        throw netcdf_exception( result );
}

} // netcdf_util

#endif
