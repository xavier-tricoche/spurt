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

#ifndef __datasetDLR_hpp
#define __datasetDLR_hpp

#include <iostream>
#include <fstream>

//#include <libgen.h>
#include <string>
#include <cassert>
#include <cstring>
#include <stdexcept>
#include <algorithm>

#include "dataset.hpp"
#include "netcdf_util.hpp"

namespace CellTree {
    
    
/* datasetDLR: import CFD mesh and associated scalar and vector attributes
 *             in NetCDF format used e.g., by Tau solver at Deutsches Zentrum
 *             fuer Luft- und Raumfahrt, Goettigen, Germany
 */ 
class datasetDLR: public dataset
{
public:
    
    datasetDLR(const std::string& filename) {
        load_from_file(filename);
    }
    
    datasetDLR(const std::string& gridname, const std::string& valuename) {
        timestep ts;
        ts.file = valuename;
        ts.time = 0.0;
        m_timesteps.push_back(ts);
    }
    
    void load_from_file( const std::string& filename )
    {
        // check that we have the right format
        std::string suffix = filename.substr( filename.rfind( '.' )+1, std::string::npos );
        if (suffix != "dlr" && suffix != "cfg") 
            throw std::runtime_error("dataset descriptor does not correspond to a DLR dataset");
        
        std::ifstream in( filename.c_str() );

        if( !in.good() )
            throw std::runtime_error( "cannot read DLR dataset descriptor " + std::string(filename) );

        char* tmp = strdup( filename.c_str() );
        free( tmp );
        
        m_basepath += '/';
    
        // start parsing and loading
        in >> m_gridfile;
        
        timestep ts;
        
        while( in.good() )
        {
            ts.file.clear();
            ts.time = 0.0;
            
            in >> ts.file >> ts.time;

            if( in.fail() )
                break;

            m_timesteps.push_back( ts );
        }

        in.close();
        
        if( m_timesteps.empty() )
            throw std::runtime_error( "could not read timesteps" );
        
        std::sort( m_timesteps.begin(), m_timesteps.end(), sort_by_time() );
    }
    
    virtual mesh* read_mesh() const
    {
        double* points = 0;
        unsigned int* indices = 0;
        mesh::cell* cells = 0;
    
        using namespace netcdf_util;
        
        try
        {
#ifndef NETCDF_CPP_INTERFACE
            int ncid;
            nccheck( nc_open( m_gridfile.c_str(), NC_NOWRITE, &ncid ) );
#else
            NcFile data_file(m_gridfile.c_str());
#endif

            size_t npoints;
#ifndef NETCDF_CPP_INTERFACE
            int pdimid;
            nccheck( nc_inq_dimid( ncid, "no_of_points", &pdimid ) );
            nccheck( nc_inq_dimlen( ncid, pdimid, &npoints ) );
#else
            npoints = data_file.get_dim("no_of_points")->size();
#endif

            size_t ncells = 0, nindices = 0;
#ifndef NETCDF_CPP_INTERFACE            
            int cdimid;
            if( NC_NOERR == nc_inq_dimid( ncid, "no_of_hexaeders", &cdimid ) )
            {
                size_t dimlen;
                nccheck( nc_inq_dimlen( ncid, cdimid, &dimlen ) );
    
                ncells   += dimlen;
                nindices += dimlen * 8;
            }
#else
            ncells += data_file.get_dim("no_of_hexaeders")->size();
            nindices += 8*ncells;
#endif

#ifndef NETCDF_CPP_INTERFACE               
            if( NC_NOERR == nc_inq_dimid( ncid, "no_of_tetraeders", &cdimid ) )
            {
                size_t dimlen;
                nccheck( nc_inq_dimlen( ncid, cdimid, &dimlen ) );
    
                ncells   += dimlen;
                nindices += dimlen * 4;
            }
#else
            ncells += data_file.get_dim("no_of_tetraeders")->size();
            nindices += 4*ncells;
#endif
    
#ifndef NETCDF_CPP_INTERFACE               
            if( NC_NOERR == nc_inq_dimid( ncid, "no_of_prisms", &cdimid ) )
            {
                size_t dimlen;
                nccheck( nc_inq_dimlen( ncid, cdimid, &dimlen ) );
    
                ncells   += dimlen;
                nindices += dimlen * 6;
            }
#else
            ncells += data_file.get_dim("no_of_prisms")->size();
            nindices += 6*ncells;
#endif

#ifndef NETCDF_CPP_INTERFACE               
            if( NC_NOERR == nc_inq_dimid( ncid, "no_of_pyramids", &cdimid ) )
            {
                size_t dimlen;
                nccheck( nc_inq_dimlen( ncid, cdimid, &dimlen ) );
    
                ncells   += dimlen;
                nindices += dimlen * 5;
            }
#else
            ncells += data_file.get_dim("no_of_pyramids")->size();
            nindices += 6*ncells;
#endif
    
            // ---

            points = new double[npoints*3];
            double* tmp = new double[npoints];
    
            try
            {
#ifndef NETCDF_CPP_INTERFACE
                int varid;            
                nccheck( nc_inq_varid( ncid, "points_xc", &varid ) );
                nccheck( nc_get_var_double( ncid, varid, tmp ) );
#else
                data_file.get_var("points_xc")->get(tmp);
#endif
        
                for( unsigned int i=0; i<npoints; ++i )
                    points[3*i+0] = tmp[i];

#ifndef NETCDF_CPP_INTERFACE
                nccheck( nc_inq_varid( ncid, "points_yc", &varid ) );
                nccheck( nc_get_var_double(ncid, varid, tmp));
#else
                data_file.get_var("points_yc")->get(tmp);
#endif

                for( unsigned int i=0; i<npoints; ++i )
                    points[3*i+1] = tmp[i];

#ifndef NETCDF_CPP_INTERFACE
                nccheck( nc_inq_varid( ncid, "points_zc", &varid ) );
                nccheck( nc_get_var_double( ncid, varid, tmp ) );
#else
                data_file.get_var("points_zc")->get(tmp);
#endif
        
                for( unsigned int i=0; i<npoints; ++i )
                    points[3*i+2] = tmp[i];
            
                delete[] tmp;
            }
            catch( ... )
            {
                delete[] tmp;
                throw;
            }
            
            // ---
            
            cells = new mesh::cell[ncells];
            indices = new unsigned int[nindices];
    
            mesh::cell* ci = cells;
    
            int varid, dimid[2];
            size_t dimsz[2], nind, start = 0;
    

#ifndef NETCDF_CPP_INTERFACE
            if( NC_NOERR == nc_inq_varid( ncid, "points_of_hexaeders", &varid ) )
            {
                nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
                nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
                nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );
    
                nccheck( nc_get_var_int( ncid, varid, (int*)(indices+start) ) );
#else
            if (true) {
                NcVar* hexa_var = data_file.get_var("points_of_hexaeders");
                dimsz[0] = hexa_var->get_dim(0)->size();
                // dimsz[1] = hexa_var->get_dim(1);
                hexa_var->get((int*)(indices+start));
#endif
    
                for( unsigned int i=0; i<dimsz[0]; ++i, start += 8 )
                    *(ci++) = mesh::cell( HEXAHEDRON, start );
            }

#ifndef NETCDF_CPP_INTERFACE    
            if( NC_NOERR == nc_inq_varid( ncid, "points_of_tetraeders", &varid ) )
            {
                nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
                nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
                nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );
    
                nccheck( nc_get_var_int( ncid, varid, (int*)(indices+start) ) );
#else
            if (true) {
                NcVar* tet_var = data_file.get_var("points_of_tetraeders");
                dimsz[0] = tet_var->get_dim(0)->size();
                // dimsz[1] = tet_var->get_dim(1);
                tet_var->get((int*)(indices+start));
#endif
                for( unsigned int i=0; i<dimsz[0]; ++i, start += 4 )
                    *(ci++) = mesh::cell( TETRAHEDRON, start );
            }
    
#ifndef NETCDF_CPP_INTERFACE    
            if( NC_NOERR == nc_inq_varid( ncid, "points_of_prisms", &varid ) )
            {
                nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
                nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
                nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );
    
                nccheck( nc_get_var_int( ncid, varid, (int*)(indices+start) ) );
#else
            if (true) {
                NcVar* pri_var = data_file.get_var("points_of_prisms");
                dimsz[0] = pri_var->get_dim(0)->size();
                // dimsz[1] = pri_var->get_dim(1);
                pri_var->get((int*)(indices+start));
#endif                
    
                for( unsigned int i=0; i<dimsz[0]; ++i, start += 6 )
                    *(ci++) = mesh::cell( PRISM, start );
            }

#ifndef NETCDF_CPP_INTERFACE        
            if( NC_NOERR == nc_inq_varid( ncid, "points_of_pyramids", &varid ) )
            {
                nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
                nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
                nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );

                nccheck( nc_get_var_int( ncid, varid, (int*)(indices+start) ) );
#else
            if (true) {
                NcVar* py_var = data_file.get_var("points_of_pyramids");
                dimsz[0] = py_var->get_dim(0)->size();
                // dimsz[1] = py_var->get_dim(1);
                py_var->get((int*)(indices+start));
#endif
                for( unsigned int i=0; i<dimsz[0]; ++i, start += 5 )
                    *(ci++) = mesh::cell( PYRAMID, start );
            }

            return new mesh( npoints, ncells, nindices, points, cells, indices );
        }
        catch( ... )
        {
            delete[] points;
            delete[] cells;
            delete[] indices;
            
            throw;
        }
    };

    virtual variable* read_scalar_variable( unsigned int timestep, const std::string& name ) const
    {
        using namespace netcdf_util;
        
        timestep = 0;
        std::string file = m_timesteps[timestep].file;

#ifndef NETCDF_CPP_INTERFACE    
        int ncid;
        nccheck(nc_open(file.c_str(), NC_NOWRITE, &ncid));
#else
        NcFile data_file(file.c_str());
#endif

        double *data = 0, *tmp = 0;
        unsigned int size = 0;

        try {            
            int ndims;

#ifndef NETCDF_CPP_INTERFACE    
            int varid;
            nccheck(nc_inq_varid(ncid, name.c_str(), &varid));
            nccheck(nc_inq_varndims(ncid, varid, &ndims));
#else
            NcVar* var = data_file.get_var(name.c_str());
            ndims = var->num_dims();
#endif

            if (ndims > 1)
                throw std::runtime_error("multi-dimensional NetCDF variables not supported");
            size_t vsize;
#ifndef NETCDF_CPP_INTERFACE    
            int dimid;
            nccheck(nc_inq_vardimid(ncid, varid, &dimid));

            nccheck(nc_inq_dimlen(ncid, dimid, &vsize));
#else
            vsize = var->get_dim(0)->size();
#endif

            size = vsize;
            data = new double[size];
            tmp = new double[size];

#ifndef NETCDF_CPP_INTERFACE    
            nccheck(nc_get_var_double(ncid, varid, tmp));
#else
            var->get(tmp);
#endif

            for (unsigned int i = 0; i < size; ++i)
                data[i] = tmp[i];

            delete[] tmp;
        }
        catch (...)
        {
            delete[] data;
            delete[] tmp;

            throw;
        }
        return new variable(1, size, data);
    }
    
    virtual variable* read_vector_variable( unsigned int timestep, const std::string& name ) const
    {
        using namespace netcdf_util;
        
        timestep = 0;

        std::string file = m_timesteps[timestep].file;
        
#ifndef NETCDF_CPP_INTERFACE    
        int ncid;
        nccheck( nc_open( file.c_str(), NC_NOWRITE, &ncid ) );
#else
        NcFile data_file(file.c_str());
#endif

        double *data = 0, *tmp = 0;
        unsigned int size = 0;

        try
        {
            const char* prefix[3] = { "x_", "y_", "z_" };

            for( unsigned int d=0; d<3; ++d )
            {
                int ndims;
                
#ifndef NETCDF_CPP_INTERFACE    
                int varid;
                nccheck( nc_inq_varid( ncid, (prefix[d]+name).c_str(), &varid ) );
                nccheck( nc_inq_varndims( ncid, varid, &ndims ) );
#else
                NcVar* var = data_file.get_var((prefix[d]+name).c_str());
                ndims = var->num_dims();
#endif
        
                if( ndims > 1 )
                    throw std::runtime_error( "multi-dimensional NetCDF variables not supported" );

                size_t vsize;
#ifndef NETCDF_CPP_INTERFACE        
                int dimid;    
                nccheck( nc_inq_vardimid( ncid, varid, &dimid ) );
                nccheck( nc_inq_dimlen( ncid, dimid, &vsize ) );
#else
                vsize = var->get_dim(0)->size();
#endif
                
                if( d == 0 )    
                {
                    size = vsize;
                    data = new double[3*size];
                    tmp  = new double[size];
                }
                else if( size != vsize )
                    throw std::runtime_error( "invalid variable size" );

#ifndef NETCDF_CPP_INTERFACE            
                nccheck( nc_get_var_double( ncid, varid, tmp ) );
#else
                var->get(tmp);
#endif

                for( unsigned int i=0; i<size; ++i )
                    data[3*i+d] = tmp[i];
            }
            
            delete[] tmp;
        }
        catch( ... )
        {
            delete[] data;
            delete[] tmp;
            
            throw;
        }
        
        return new variable( 3, size, data );
    }
    
    virtual double get_timestep( unsigned int num ) const 
    {
        return m_timesteps[num].time;
    }
    
    virtual unsigned int get_num_timesteps() const
    {
        return m_timesteps.size();
    }
        
protected:

    struct timestep {
        double      time;
        std::string file;
    };

    struct sort_by_time {
        bool operator()( const timestep& t0, const timestep& t1 ) const {
            return t0.time < t1.time;
        }
    };

    std::vector<timestep> m_timesteps;
    std::string           m_gridfile;
    std::string           m_basepath;
};

} // namespace CellTree


#endif // __datasetDLR_hpp
