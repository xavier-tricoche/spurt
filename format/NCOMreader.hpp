#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <math/fixed_vector.hpp>
#include <stdexcept>
#include <sstream>

#include <cassert>
#include <list>

#include <stdlib.h>
#include <string>

#include <fstream>
#include <unistd.h>
#include <netcdf.h>

namespace {
inline void check_nc( int status ) {
    if( status != NC_NOERR )
        throw std::runtime_error(nc_strerror(status));
}
}

namespace spurt {

namespace NCOMreader {

void load_dim_value(size_t& nlon, size_t& nlat, int file) {
    int ndimid;
    
    try {
    check_nc(nc_inq_dimid(file, "lat", &ndimid));
    check_nc(nc_inq_dimlen(file, ndimid, &nlat));

    check_nc(nc_inq_dimid(file, "lon", &ndimid));
    check_nc(nc_inq_dimlen(file, ndimid, &nlon));
    }
    catch( std::runtime_error& e) {
        std::cout << "exception caught in load_dim_value: "
            << e.what() << '\n';
    }
}

void load_coordinates(std::vector<double>& lat, std::vector<double>& lon,
                      int file) {
                          
    size_t nlat, nlon;
    load_dim_value(nlon, nlat, file);
    
    int id;
    lat.resize(nlat);
    lon.resize(nlon);
    size_t count = nlat, start = 0;
    ptrdiff_t stride = 1, imap = 1;
    
    try {
    check_nc(nc_inq_varid(file, "lat", &id));
    check_nc(nc_get_varm_double(file, id, &start, &count,
                                &stride, &imap, &lat[0]));
                    
    count = nlon;
    start = 0;            
    check_nc(nc_inq_varid(file, "lon", &id));
    check_nc(nc_get_varm_double(file, id, &start, &count,
                                &stride, &imap, &lon[0]));
    }
    catch( std::runtime_error& e) {
        std::cout << "exception caughtin load_coordinates: "
            << e.what() << '\n';
    }
}

void var_info(int file, const std::string& var_name) {
    int  var_id;                        /* variable ID */
    nc_type var_type;                   /* variable type */
    int  var_ndims;                     /* number of dims */
    int  var_dimids[NC_MAX_VAR_DIMS];   /* dimension IDs */
    int  var_natts;                     /* number of attributes */
        
    check_nc(nc_inq_varid (file, var_name.c_str(), &var_id));
    check_nc(nc_inq_var (file, var_id, 0, &var_type, &var_ndims, var_dimids,
                         &var_natts));
                         
    std::cout << var_name << " has id " << var_id << '\n';
    std::cout << var_name << " has type " << var_type << '\n';
    std::cout << var_name << " has " << var_ndims << " dimensions\n";
    for (int i=0; i<var_ndims; ++i) {
        std::cout << var_name << "[" << i << "] has id " << var_dimids[i] << '\n';
    }
    std::cout << var_name << " has " << var_natts << " attributes\n";
}

void load_velocity(std::vector<nvis::vec2>& uv, 
                   size_t nlat, size_t nlon, int file) {
    int id;
    uv.resize(nlat*nlon);
    size_t start[4] = {0, 0, 0, 0};
    size_t count[4] = {1, 1, nlat, nlon};
    
    std::vector<double> u(nlat*nlon);
    std::vector<double> v(nlat*nlon);
    
    try {
        check_nc(nc_inq_varid(file, "water_u", &id));
        check_nc(nc_get_vara_double(file, id, start, count,
                                    reinterpret_cast<double *>(&u[0])));
        check_nc(nc_inq_varid(file, "water_v", &id));
        check_nc(nc_get_vara_double(file, id, start, count,
                                    reinterpret_cast<double *>(&v[0])));
         
        // ugly but smarter strided / mapped version of nc_get_var is
        // a pain to use                        
        for (int i=0; i<u.size(); ++i) {
            uv[i][0]=u[i];
            uv[i][1]=v[i];
        }
    }
    catch( std::runtime_error& e) {
        std::cout << "exception caught in load_velocity: " << e.what()
            << '\n';
    }
}

void load_variable(std::vector<double>& v, const std::string& name, size_t nlat, size_t nlon, int file) {
    int id;
    v.resize(nlat*nlon);
        
    try {
        check_nc(nc_inq_varid(file, name.c_str(), &id));
        check_nc(nc_get_var_double(file, id,
                                   reinterpret_cast<double *>(&v[0])));
    }
    catch( std::runtime_error& e) {
        std::cout << "exception caught in load_variable: " << e.what()
            << '\n';
    }
    
}

double load_time(int file) {
    int id;
    size_t start = 0;
    double t;
    try {
        check_nc( nc_inq_varid( file, "time", &id ) );
        check_nc( nc_get_var1_double( file, id, &start, &t ) );
    }
    catch( std::runtime_error& e ) {
        std::cout << "exception caught in load_velocity: " << e.what()
            << '\n';
    }
    return t;
}

void load_dataset(const std::string& filename, 
                  std::vector<double>& lat, std::vector<double>& lon,
                  std::vector<nvis::vec2>& vel, double& t) {
    int file;
    check_nc( nc_open( filename.c_str(), NC_NOWRITE, &file ) );
    load_coordinates(lat, lon, file);
    load_velocity(vel, lat.size(), lon.size(), file);
    t=load_time(file);
}

void load_dataset(const std::string& filename, size_t& nlat, size_t& nlon,
                  std::vector<nvis::vec2>& vel, double& t) {
    int file;
    check_nc( nc_open( filename.c_str(), NC_NOWRITE, &file ) );
    load_dim_value(nlon, nlat, file);
    load_velocity(vel, nlat, nlon, file);
    t=load_time(file);
}

void load_dataset(const std::string& filename, std::vector<double>& lat, std::vector<double>& lon,
                  std::vector<nvis::vec2>& vel, const std::vector<std::string>& var_names, 
                  std::vector< std::vector<double> >& variables, double& t) {
    int file;
    check_nc( nc_open( filename.c_str(), NC_NOWRITE, &file ) );
    load_coordinates(lat, lon, file);
    load_velocity(vel, lat.size(), lon.size(), file);
    t=load_time(file);
    variables.resize(var_names.size());
    for (int i=0; i<var_names.size(); ++i) {
        load_variable(variables[i], var_names[i], lat.size(), lon.size(), file);
    }
    
}

} // NCOMreader

} // spurt
