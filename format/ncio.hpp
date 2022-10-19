#ifndef __ncio_hpp
#define __ncio_hpp

#include <netcdf.h>
#include <cstdio>
#include <vector>
#include <stdexcept>
#include <cstdarg>

#include <netcdf.h>

#include  <math/fixed_vector.hpp>

template<typename T>
struct netcdf_traits
{
};

template<> struct netcdf_traits<float>
{
    enum { components = 1, nctype = NC_FLOAT };
};

template<> struct netcdf_traits<double>
{
    enum { components = 1, nctype = NC_DOUBLE };
};

template<> struct netcdf_traits<int>
{
    enum { components = 1, nctype = NC_INT };
};

template<> struct netcdf_traits<unsigned int>
{
    enum { components = 1, nctype = NC_INT };
};

template<typename T, size_t N> 
struct netcdf_traits< nvis::fixed_vector<T,N> >
{
    enum { components = N, nctype = netcdf_traits<T>::nctype };
};

template<typename T, unsigned int N> 
struct netcdf_traits< nvis::fixed_vector<T,N> >
{
    enum { components = N, nctype = netcdf_traits<T>::nctype };
};

// template<typename T, unsigned int N, unsigned int M> 
// struct netcdf_traits< nvis::fixed_matrix<T,N,M> >
// {
//     enum { components = N*M, nctype = netcdf_traits<T>::nctype };
// };

// -------------------------------------------------------------------------

namespace ncio {

inline void check( int result )
{
    if( result != NC_NOERR )
        throw std::runtime_error( nc_strerror(result) );
}

inline void message( int result, const char* path, const char* msg, ... )
{
    va_list argp;
        
    fprintf( stderr, "NetCDF %s (%s): ", result == NC_NOERR ? "warning" : "error", path );
    
    va_start( argp, msg );
    vfprintf( stderr, msg, argp );
    va_end( argp );
    
    fprintf( stderr, "\n" );
    
    check( result );
}

inline const char* netcdf_typename( nc_type type )
{
    switch( type )
    {
    case NC_BYTE:   return "NC_BYTE";
    case NC_CHAR:   return "NC_CHAR";
    case NC_SHORT:  return "NC_SHORT";
    case NC_INT:    return "NC_INT";
    case NC_FLOAT:  return "NC_FLOAT";
    case NC_DOUBLE: return "NC_DOUBLE";
    default:        return "UNKNOWN";
    }
}

// -------------------------------------------------------------------------

inline void create( const char* path )
{ 
    int ncid;
    
    check( nc_create( path, 0, &ncid ) );
    check( nc_close( ncid ) );
}

template<typename T>
void put( const char* path, const std::vector<T>& data, const char* name )
{
    int  ncid;
    bool define = false;

    try
    {
        int result = nc_open( path, NC_WRITE, &ncid );
        
        if( result == ENOENT )
        {
            check( nc_create( path, NC_NOCLOBBER, &ncid ) );
            define = true;
        }
        else
            check( result );
        
        printf( "NetCDF file %s opened for writing, ncid = %d\n", path, ncid );
        
        int varid;
        result = nc_inq_varid( ncid, name, &varid );
        
        if( result == NC_NOERR )
        {
            int ndims;
            check( nc_inq_varndims( ncid, varid, &ndims ) );
            
            int dimid[ndims];
            check( nc_inq_vardimid( ncid, varid, dimid ) );
            
            size_t nelem = 1, length;
            
            for( unsigned int i=0; i<ndims; ++i )
            {
                check( nc_inq_dimlen( ncid, dimid[i], &length ) );
                nelem *= length;
            }
            
            if( nelem != netcdf_traits<T>::components * data.size() )
                message( NC_EVARSIZE, path, "writing %d elements to existing variable \"%s\" of total size %d", 
                         data.size(), name, nelem );
                
            nc_type type;
            check( nc_inq_vartype( ncid, varid, &type ) );
            
            if( type != (nc_type)netcdf_traits<T>::nctype )
                message( NC_NOERR, path, "writing %s to existing variable \"%s\" of type %s",
                         netcdf_typename((nc_type)netcdf_traits<T>::nctype),
                         name, netcdf_typename(type) );
        }
        else if( result == NC_ENOTVAR )
        {
            char dimname[2][NC_MAX_NAME];
            int ndims = 1, dimid[2];
            size_t dimlen[2];
            
            snprintf( dimname[0], NC_MAX_NAME, "%s_len", name );
            dimlen[0] = data.size();
            
            if( netcdf_traits<T>::components > 1 )
            {
                ++ndims;
                
                snprintf( dimname[1], NC_MAX_NAME, "%s_dim", name );
                dimlen[1] = netcdf_traits<T>::components;
            }
            
            for( unsigned int i=0; i<ndims; ++i )
            {
                result = nc_inq_dimid( ncid, name, &dimid[i] );
                
                if( result == NC_NOERR )
                {
                    size_t length;
                    check( nc_inq_dimlen( ncid, dimid[i], &length ) );
                    
                    if( length != dimlen[i] )
                        message( NC_EDIMSIZE, path, "creating dimension \"%s\" of size %ld conflicts" \
                                                    " with existing dimension of size %ld",
                                                    dimname[i], dimlen[i], length );
                }
                else if( result == NC_EBADDIM )
                    dimid[i] = -1;
                else 
                    check( result );
            }

            if( !define )
                check( nc_redef( ncid ) );
            
            for( unsigned int i=0; i<ndims; ++i )
            {
                if( dimid[i] > 0 )
                    continue;
                    
                check( nc_def_dim( ncid, dimname[i], dimlen[i], &dimid[i] ) );
            }            
            
            check( nc_def_var( ncid, name, (nc_type)netcdf_traits<T>::nctype, ndims, dimid, &varid ) );
            check( nc_enddef( ncid ) );
        }
        else
            check( result );

        nc_type type = (nc_type)netcdf_traits<T>::nctype;
        
        switch( type )
        {
        case NC_BYTE:   check( nc_put_var_uchar( ncid, varid, (unsigned char*)&*(data.begin() ) ) ); break;
        case NC_CHAR:   check( nc_put_var_schar( ncid, varid, (signed char*)&*(data.begin() ) ) ); break;
        case NC_SHORT:  check( nc_put_var_short( ncid, varid, (short*)&*(data.begin() ) ) ); break;
        case NC_INT:    check( nc_put_var_int( ncid, varid, (int*)&*(data.begin() ) ) );break;
        case NC_FLOAT:  check( nc_put_var_float( ncid, varid, (float*)&*(data.begin() ) ) ); break;
        case NC_DOUBLE: check( nc_put_var_double( ncid, varid, (double*)&*(data.begin() ) ) ); break;
        }
        
        check( nc_close( ncid ) );
    }    
    catch( std::runtime_error& e )
    {
        nc_abort( ncid );
        throw;
    }
}

template<typename T>
void get( const char* path, std::vector<T>& data, const char* name )
{
    int ncid;
    
    try
    {
        check( nc_open( path, NC_NOWRITE, &ncid ) );

        int varid;
        check( nc_inq_varid( ncid, name, &varid ) );
        
        int ndims;
        check( nc_inq_varndims( ncid, varid, &ndims) );
        
        int dimids[ndims];
        check( nc_inq_vardimid( ncid, varid, dimids ) );
        
        size_t dimlen[ndims], nelem = 1;
        
        for( unsigned int i=0; i<ndims; ++i )
        {
            check( nc_inq_dimlen( ncid, dimids[i], &dimlen[i] ) );
            nelem *= dimlen[i];
        }

        if( ndims > 1 )
        {
            if( netcdf_traits<T>::components > 1 )
            {
                if( netcdf_traits<T>::components != dimlen[ndims-1] )
                    message( NC_EVARSIZE, path, "reading %d-component data from variable \"%s\", " \
                                                "but last dimension has size %ld\n",
                                                netcdf_traits<T>::components, name, dimlen[ndims-1] );
            }
        }
        else
        {
            if( nelem % netcdf_traits<T>::components )
                message( NC_EVARSIZE, path, "reading %d-component data from variable \"%s\", " \
                                            "total size %ld is not a multiple of %d",
                                            netcdf_traits<T>::components, name, 
                                            nelem, netcdf_traits<T>::components );
        }
        
        nelem /= netcdf_traits<T>::components;
        
        nc_type vartype;
        check( nc_inq_vartype( ncid, varid, &vartype ) );

        if( vartype != (nc_type)netcdf_traits<T>::nctype )
            message( NC_NOERR, path, "reading %s from variable \"%s\" of type %s",
                       netcdf_typename((nc_type)netcdf_traits<T>::nctype), 
                       name, netcdf_typename(vartype) );
    
        data.resize( nelem );
        
        int type = netcdf_traits<T>::nctype; 
        
        switch( type )
        {
        case NC_BYTE:   check( nc_get_var_uchar( ncid, varid, (unsigned char*)&*(data.begin() ) ) ); break;
        case NC_CHAR:   check( nc_get_var_schar( ncid, varid, (signed char*)&*(data.begin() ) ) ); break;
        case NC_SHORT:  check( nc_get_var_short( ncid, varid, (short*)&*(data.begin() ) ) ); break;
        case NC_INT:    check( nc_get_var_int( ncid, varid, (int*)&*(data.begin() ) ) );break;
        case NC_FLOAT:  check( nc_get_var_float( ncid, varid, (float*)&*(data.begin() ) ) ); break;
        case NC_DOUBLE: check( nc_get_var_double( ncid, varid, (double*)&*(data.begin() ) ) ); break;
        }
        
        check( nc_close( ncid ) );
    }
    catch( std::runtime_error& )
    {
        nc_close( ncid );
        throw;
    }
};

} // namespace ncio

#endif // __ncio_hpp
