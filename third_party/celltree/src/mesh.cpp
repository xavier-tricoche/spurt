#include <netcdf.h>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include <cstdio>

#include "mesh.hpp"

inline static void nccheck( int result )
{
    if( result != NC_NOERR )
        throw std::runtime_error( nc_strerror(result) );
}
    
// -------------------------------------------------------------------------

#if 0
void readDLR( const char* file, mesh& m )
{
    int ncid;
    nccheck( nc_open( file, NC_NOWRITE, &ncid ) );

    // find out number of points, cells, and indices

    int pdimid;
    nccheck( nc_inq_dimid( ncid, "no_of_points", &pdimid ) );

    size_t npoints;
    nccheck( nc_inq_dimlen( ncid, pdimid, &npoints ) );

    size_t ncells = 0, nindices = 0;

    int cdimid;

    if( NC_NOERR == nc_inq_dimid( ncid, "no_of_tetraeders", &cdimid ) )
    {
        size_t dimlen;
        nccheck( nc_inq_dimlen( ncid, cdimid, &dimlen ) );
    
        ncells   += dimlen;
        nindices += dimlen * 4;
    }
    
    if( NC_NOERR == nc_inq_dimid( ncid, "no_of_hexaeders", &cdimid ) )
    {
        size_t dimlen;
        nccheck( nc_inq_dimlen( ncid, cdimid, &dimlen ) );
    
        ncells   += dimlen;
        nindices += dimlen * 8;
    }
    
    if( NC_NOERR == nc_inq_dimid( ncid, "no_of_prisms", &cdimid ) )
    {
        size_t dimlen;
        nccheck( nc_inq_dimlen( ncid, cdimid, &dimlen ) );
    
        ncells   += dimlen;
        nindices += dimlen * 6;
    }

    if( NC_NOERR == nc_inq_dimid( ncid, "no_of_pyramids", &cdimid ) )
    {
        size_t dimlen;
        nccheck( nc_inq_dimlen( ncid, cdimid, &dimlen ) );
    
        ncells   += dimlen;
        nindices += dimlen * 5;
    }
    
    std::cout << npoints << " points, " << ncells << " cells\n";

    printf( "mesh memory size = %.2fMB\n", 
            (npoints*sizeof(float)*3 + nindices*sizeof(uint32_t) + 
             ncells*(sizeof(uint32_t)+sizeof(unsigned char))) / 1048576.0 );

    // ---

    m.points.resize( npoints );
    
    {
        float* tmp = new float[npoints];
        int varid;

        nccheck( nc_inq_varid( ncid, "points_xc", &varid ) );
        nccheck( nc_get_var_float( ncid, varid, tmp ) );
        
        for( unsigned int i=0; i<npoints; ++i )
            m.points[i].coord[0] = tmp[i];

        nccheck( nc_inq_varid( ncid, "points_yc", &varid ) );
        nccheck( nc_get_var_float( ncid, varid, tmp ) );

        for( unsigned int i=0; i<npoints; ++i )
            m.points[i].coord[1] = tmp[i];

        nccheck( nc_inq_varid( ncid, "points_zc", &varid ) );
        nccheck( nc_get_var_float( ncid, varid, tmp ) );
        
        for( unsigned int i=0; i<npoints; ++i )
            m.points[i].coord[2] = tmp[i];
            
        delete[] tmp;
    }
    
    m.cells.resize( ncells );
    m.indices.resize( nindices );
    
    std::vector<mesh::cell>::iterator ci = m.cells.begin();
    
    int varid, dimid[2];
    size_t dimsz[2], nind, start = 0;
    
    if( NC_NOERR == nc_inq_varid( ncid, "points_of_hexaeders", &varid ) )
    {
        nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
        nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
        nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );
    
        nccheck( nc_get_var_int( ncid, varid, (int*)&(m.indices[start]) ) );
    
        for( unsigned int i=0; i<dimsz[0]; ++i, start += 8 )
            *(ci++) = mesh::cell( HEXAHEDRON, start );
    }
    
    if( NC_NOERR == nc_inq_varid( ncid, "points_of_tetraeders", &varid ) )
    {
        nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
        nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
        nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );
    
        nccheck( nc_get_var_int( ncid, varid, (int*)&(m.indices[start]) ) );
    
        for( unsigned int i=0; i<dimsz[0]; ++i, start += 4 )
            *(ci++) = mesh::cell( TETRAHEDRON, start );
    }
    
    if( NC_NOERR == nc_inq_varid( ncid, "points_of_prisms", &varid ) )
    {
        nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
        nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
        nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );
    
        nccheck( nc_get_var_int( ncid, varid, (int*)&(m.indices[start]) ) );
    
        for( unsigned int i=0; i<dimsz[0]; ++i, start += 6 )
            *(ci++) = mesh::cell( PRISM, start );
    }
    
    if( NC_NOERR == nc_inq_varid( ncid, "points_of_pyramids", &varid ) )
    {
        nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
        nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
        nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );

        nccheck( nc_get_var_int( ncid, varid, (int*)&(m.indices[start]) ) );

        for( unsigned int i=0; i<dimsz[0]; ++i, start += 5 )
            *(ci++) = mesh::cell( PYRAMID, start );
    }
     
    assert( ci == m.cells.end() );
}
#endif

// -------------------------------------------------------------------------
