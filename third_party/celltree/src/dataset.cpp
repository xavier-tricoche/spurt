#include <iostream>
#include <fstream>

#include <netcdf.h>
#include <libgen.h>
#include <string>
#include <cassert>
#include <cstring>
#include <stdexcept>
#include <algorithm>

#include "dataset.hpp"

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

// -------------------------------------------------------------------------

class datasetVTK: public celltree::dataset 
{
public:
    datasetVTK( const std::string& filename) 
    {
    }
    
    virtual celltree::mesh* read_mesh() const
    {
        throw std::runtime_error("Not implemented");        
    }
    
    virtual celltree::variable* read_scalar_variable( unsigned int timestep, const std::string& name ) const
    {
        throw std::runtime_error("Not implemented");
    }
    
    virtual celltree::variable* read_vector_variable( unsigned int timestep, const std::string& name ) const
    {
        throw std::runtime_error("Not implemented");
    }
    
    virtual double get_timestep( unsigned int num ) const 
    {
        throw std::runtime_error("Not implemented");
    }
    
    virtual unsigned int get_num_timesteps() const
    {
        throw std::runtime_error("Not implemented");
    }
    
};

// -------------------------------------------------------------------------

class datasetDLR: public celltree::dataset
{
public:
    
    datasetDLR( const std::string& filename )
    {
        std::ifstream in( filename.c_str() );

        if( !in.good() )
            throw std::runtime_error( "cannot read DLR dataset descriptor " + std::string(filename) );

        char* tmp = strdup( filename.c_str() );
        m_basepath = dirname( (char*)tmp );
        free( tmp );
        
        m_basepath += '/';
    
        // start parsing and loading
        in >> m_gridfile;

        if( m_gridfile[0] != '/' )
            m_gridfile = m_basepath + m_gridfile;
        
        timestep ts;
        
        while( in.good() )
        {
            ts.file.clear();
            ts.time = 0.0;
            
            in >> ts.file >> ts.time;

            if( in.fail() )
                break;

            if( ts.file[0] != '/' )
                ts.file = m_basepath + ts.file;

            m_timesteps.push_back( ts );
        }

        in.close();
        
        if( m_timesteps.empty() )
            throw std::runtime_error( "could not read timesteps" );
        
        std::sort( m_timesteps.begin(), m_timesteps.end(), sort_by_time() );
    }
    
    virtual celltree::mesh* read_mesh() const
    {
        float* points = 0;
        unsigned int* indices = 0;
        celltree::mesh::cell* cells = 0;
        
        try
        {
            int ncid;
            nccheck( nc_open( m_gridfile.c_str(), NC_NOWRITE, &ncid ) );

            int pdimid;
            nccheck( nc_inq_dimid( ncid, "no_of_points", &pdimid ) );

            size_t npoints;
            nccheck( nc_inq_dimlen( ncid, pdimid, &npoints ) );

            size_t ncells = 0, nindices = 0;

            int cdimid;

            if( NC_NOERR == nc_inq_dimid( ncid, "no_of_hexaeders", &cdimid ) )
            {
                size_t dimlen;
                nccheck( nc_inq_dimlen( ncid, cdimid, &dimlen ) );
    
                ncells   += dimlen;
                nindices += dimlen * 8;
            }
    
            if( NC_NOERR == nc_inq_dimid( ncid, "no_of_tetraeders", &cdimid ) )
            {
                size_t dimlen;
                nccheck( nc_inq_dimlen( ncid, cdimid, &dimlen ) );
    
                ncells   += dimlen;
                nindices += dimlen * 4;
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
    
            // ---

            points = new float[npoints*3];
            float* tmp = new float[npoints];
    
            try
            {
                int varid;

                nccheck( nc_inq_varid( ncid, "points_xc", &varid ) );
                nccheck( nc_get_var_float( ncid, varid, tmp ) );
        
                for( unsigned int i=0; i<npoints; ++i )
                    points[3*i+0] = tmp[i];

                nccheck( nc_inq_varid( ncid, "points_yc", &varid ) );
                nccheck( nc_get_var_float( ncid, varid, tmp ) );

                for( unsigned int i=0; i<npoints; ++i )
                    points[3*i+1] = tmp[i];

                nccheck( nc_inq_varid( ncid, "points_zc", &varid ) );
                nccheck( nc_get_var_float( ncid, varid, tmp ) );
        
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
            
            cells = new celltree::mesh::cell[ncells];
            indices = new unsigned int[nindices];
    
            celltree::mesh::cell* ci = cells;
    
            int varid, dimid[2];
            size_t dimsz[2], nind, start = 0;
    
            if( NC_NOERR == nc_inq_varid( ncid, "points_of_hexaeders", &varid ) )
            {
                nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
                nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
                nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );
    
                nccheck( nc_get_var_int( ncid, varid, (int*)(indices+start) ) );
    
                for( unsigned int i=0; i<dimsz[0]; ++i, start += 8 )
                    *(ci++) = celltree::mesh::cell( celltree::HEXAHEDRON, start );
            }
    
            if( NC_NOERR == nc_inq_varid( ncid, "points_of_tetraeders", &varid ) )
            {
                nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
                nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
                nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );
    
                nccheck( nc_get_var_int( ncid, varid, (int*)(indices+start) ) );
    
                for( unsigned int i=0; i<dimsz[0]; ++i, start += 4 )
                    *(ci++) = celltree::mesh::cell( celltree::TETRAHEDRON, start );
            }
    
            if( NC_NOERR == nc_inq_varid( ncid, "points_of_prisms", &varid ) )
            {
                nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
                nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
                nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );
    
                nccheck( nc_get_var_int( ncid, varid, (int*)(indices+start) ) );
    
                for( unsigned int i=0; i<dimsz[0]; ++i, start += 6 )
                    *(ci++) = celltree::mesh::cell( celltree::PRISM, start );
            }
    
            if( NC_NOERR == nc_inq_varid( ncid, "points_of_pyramids", &varid ) )
            {
                nccheck( nc_inq_vardimid( ncid, varid, dimid ) );
                nccheck( nc_inq_dimlen( ncid, dimid[0], &dimsz[0] ) );
                nccheck( nc_inq_dimlen( ncid, dimid[1], &dimsz[1] ) );

                nccheck( nc_get_var_int( ncid, varid, (int*)(indices+start) ) );

                for( unsigned int i=0; i<dimsz[0]; ++i, start += 5 )
                    *(ci++) = celltree::mesh::cell( celltree::PYRAMID, start );
            }

            return new celltree::mesh( npoints, ncells, nindices, points, cells, indices );
        }
        catch( ... )
        {
            delete[] points;
            delete[] cells;
            delete[] indices;
            
            throw;
        }
    };

    virtual celltree::variable* read_scalar_variable( unsigned int timestep, const std::string& name ) const
    {
        std::string file = m_timesteps[timestep].file;
        
        int ncid;
        nccheck( nc_open( file.c_str(), NC_NOWRITE, &ncid ) );

        float *data = 0;
        unsigned int size = 0;

        try
        {
            int varid;
            nccheck( nc_inq_varid( ncid, name.c_str(), &varid ) );
        
            int ndims;
            nccheck( nc_inq_varndims( ncid, varid, &ndims ) );
        
            if( ndims > 1 )
                throw std::runtime_error( "multi-dimensional NetCDF variables not supported" );
        
            int dimid;    
            nccheck( nc_inq_vardimid( ncid, varid, &dimid ) );
          
            size_t vsize;
            nccheck( nc_inq_dimlen( ncid, dimid, &vsize ) );
            
            size = vsize;
            data = new float[size];
    
            nccheck( nc_get_var_float( ncid, varid, data ) );
        }
        catch( ... )
        {
            delete[] data;
            
            throw;
        }
        
        return new celltree::variable( 1, size, data );
    }
    
    virtual celltree::variable* read_vector_variable( unsigned int timestep, const std::string& name ) const
    {
        // timestep = 0;

        std::string file = m_timesteps[timestep].file;
        
        int ncid;
        nccheck( nc_open( file.c_str(), NC_NOWRITE, &ncid ) );

        float *data = 0, *tmp = 0;
        unsigned int size = 0;

        try
        {
            const char* prefix[3] = { "x_", "y_", "z_" };

            for( unsigned int d=0; d<3; ++d )
            {
                int varid;
                nccheck( nc_inq_varid( ncid, (prefix[d]+name).c_str(), &varid ) );
        
                int ndims;
                nccheck( nc_inq_varndims( ncid, varid, &ndims ) );
        
                if( ndims > 1 )
                    throw std::runtime_error( "multi-dimensional NetCDF variables not supported" );
        
                int dimid;    
                nccheck( nc_inq_vardimid( ncid, varid, &dimid ) );
            
                size_t vsize;
                nccheck( nc_inq_dimlen( ncid, dimid, &vsize ) );
                
                if( d == 0 )    
                {
                    size = vsize;
                    data = new float[3*size];
                    tmp  = new float[size];
                }
                else if( size != vsize )
                    throw std::runtime_error( "invalid variable size" );
    
                nccheck( nc_get_var_float( ncid, varid, tmp ) );

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
        
        return new celltree::variable( 3, size, data );
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
        float       time;
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

// -------------------------------------------------------------------------

#include <format/avtNek5000FileFormat.hpp>

class datasetNek5000: public celltree::dataset
{
public:
    
    datasetNek5000( const std::string& metafile ) : fmt( metafile.c_str() )
    {
    }
        
    virtual celltree::mesh* read_mesh() const
    {  
        std::vector<float> coords;
        std::vector<unsigned int> ids;
        fmt.GetMesh(coords, ids);
        
        int dim = fmt.GetDimension();
        assert(dim == 3);
        
        int ncells = ids.size() / 8;
        int nindices = ids.size();
        int npts = coords.size()/3;
        
        celltree::mesh::cell* cells = new celltree::mesh::cell[ncells];
        unsigned int* indices = new unsigned int[nindices];
        float* points = new float[npts*3];
        
        std::copy(coords.begin(), coords.end(), points);
        std::copy(ids.begin(), ids.end(), indices);
        

        celltree::mesh::cell* ci = cells;
        for (size_t i=0; i<ncells; ++i) {
            *(ci++) = celltree::mesh::cell(celltree::HEXAHEDRON, 8*i);
        }
        
        return new celltree::mesh( npts, ncells, nindices, points, cells, indices );
    };

    virtual celltree::variable* read_scalar_variable( unsigned int timestep, const std::string& name ) const
    {
        throw std::runtime_error( "not implemented" );
    }
    
    virtual celltree::variable* read_vector_variable( unsigned int timestep, const std::string& name ) const
    {   
        if( name != "velocity" )
            throw std::runtime_error( "unknown variable " + name );
        std::vector<float> values;
        
        fmt.GetVectorVar(timestep,  values);
        
        float* vals = new float[values.size()];
        return new celltree::variable(3, values.size()/3, vals);
    }
    
    virtual double get_timestep( unsigned int num ) const 
    {
        std::vector<double> t;
        fmt.GetTimes( t );
        
        return t[num];
    }
    
    virtual unsigned int get_num_timesteps() const
    {
        return fmt.GetNTimesteps();
    }
        
protected:

    mutable avtNek5000FileFormat fmt;
};

// -------------------------------------------------------------------------

celltree::dataset* celltree::dataset::create( const std::string& filename )
{
    std::string suffix = filename.substr( filename.rfind( '.' )+1, std::string::npos );
    
    if( suffix == "dlr" || suffix == "cfg" )
        return new datasetDLR( filename );
    else if( suffix == "nek3d" )
        return new datasetNek5000( filename );
    else if ( suffix == "vtk" || suffix == "VTK")
        return new datasetVTK( filename );
    else
        throw std::runtime_error( "could not identify a dataset reader for \"" + suffix + "\" suffix" );
}
