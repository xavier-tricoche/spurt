#ifndef __LOGICAL_TO_PHYSICAL_HPP__
#define __LOGICAL_TO_PHYSICAL_HPP__

#include <vector>
#include <fstream>
#include <math/fixed_vector.hpp>
#include <hdf5.h>
#include <poincare/metric.hpp>

#define CHECK_POINTS

namespace xavier {

class logical2physical {
    int index(int i, int j) const {
        return i + j*size[0];
    }
    
public:
    logical2physical(const std::string& h5file) {
        hid_t file_id = H5Fopen(h5file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        
        hid_t dataset_id = H5Dopen1(file_id, "cartGrid");
        
        if (dataset_id < 0) {
            throw std::runtime_error("couldn't get \"cartGrid\" dataset");
        }
        
        hid_t dataspace_id = H5Dget_space(dataset_id);
        
        if (dataspace_id < 0) {
            throw std::runtime_error("couldn't get \"cartGrid\" dataspace");
        }
        
        
        int ndims = H5Sget_simple_extent_ndims(dataspace_id);
        
        if (ndims != 4) {
            throw std::runtime_error("dataspace \"cartGrid\" is not of dimension 4");
        }
        
        hsize_t dims[4];
        H5Sget_simple_extent_dims(dataspace_id, dims, NULL);
        H5Sclose(dataspace_id);
        
        std::reverse(dims, dims + 4);
        std::copy(dims + 1, dims + 4, size.begin());
        std::reverse(dims, dims + 4);
        
        std::cout << "grid dimensions: " << size[0] << 'x' << size[1] << 'x' << size[2] << '\n';
        hsize_t npoints = dims[0] * dims[1] * dims[2];
        
        std::cout << "npoints = " << npoints << '\n';
        
        points = new nvis::fvec3[npoints];
        H5Dread(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, points);
        H5Dclose(dataset_id);
        
        box.reset();
        for (int i=0 ; i<size[0]*size[1] ; ++i) {
            nvis::fvec3 p = points[i];
            box.add(nvis::vec2(p[0], p[2]));
        }
        
#ifdef CHECK_POINTS
        std::fstream pts("points.txt", std::ios::out);
        for (int i=0 ; i<size[0]*size[1] ; ++i) {
            pts << points[i][0] << " " << points[i][2] << " 0\n";
        }
        pts.close();
#endif
        
        metric.bounds() = nvis::bbox2(nvis::vec2(0,0), nvis::vec2(size[1]-1, size[0]-1));
        metric.periodic(0) = true;
        metric.periodic(1) = false;
    }
    
    nvis::vec2 operator()(const nvis::vec2& x) const {
        nvis::vec2 y = metric.modulo(x);
        nvis::vec2 z(y[1], y[0]);
        int i = floor(z[0]);
        int j = floor(z[1]);
        nvis::vec2 loc = z - nvis::vec2(i,j);
        double v = loc[0];
        double u = loc[1];
        nvis::fvec3 tmp = (1-u)*(1-v)*points[index(i,j)] + u*(1-v)*points[index(i+1,j)] +
                          u*v*points[index(i+1,j+1)] + (1-u)*v*points[index(i,j+1)];
        return nvis::vec2(tmp[0], tmp[2]);
    }
    
    const nvis::bbox2& bounds() const {
        return box;
    }
    
protected:
    nvis::fvec3*    points;
    nvis::uvec3     size;
    map_metric      metric;
    nvis::bbox2     box;
};


}

#endif