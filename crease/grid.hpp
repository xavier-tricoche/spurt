#ifndef __GRID__HPP__
#define __GRID__HPP__

#include <teem/nrrd.h>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <image/nrrd_wrapper.hpp>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>
#include <string>

namespace spurt {
typedef std::pair< nvis::vec3, nvis::vec3 > Edge;
typedef std::pair< unsigned int, unsigned int > EdgeId;

// unique identification of a voxel face
struct FaceId {
    FaceId() : is(4) {}
    FaceId(unsigned int i0, unsigned int i1,
           unsigned int i2, unsigned int i3);
    FaceId(const FaceId& fid);
    
    int operator<(const FaceId& fid) const;
    int operator==(const FaceId& fid) const;
    FaceId& operator=(const FaceId& fid);
    
    std::vector< unsigned int > is;
};

std::ostream& operator<<(std::ostream& os, const FaceId& fid);

// sampling grid
struct Grid {
    Grid(const Nrrd* nrrd, unsigned int upsampling = 1);
    
    unsigned int id(unsigned int i, unsigned int j, unsigned int k) const;
    nvis::vec3 operator()(unsigned int i, unsigned int j, unsigned int k) const;
    nvis::vec3 operator()(const nvis::vec3& p) const;
    nvis::vec3 globalc(const nvis::vec3& p) const;
    
    void voxel(std::vector< nvis::vec3 >& v,
               unsigned int i, unsigned int j, unsigned int k) const;
               
    nvis::bbox3 bounds;
    nvis::ivec3 size;
    nvis::vec3 d;
};

struct Slice {
    Slice(const Nrrd* nrrd, unsigned int dim, unsigned int pos,
          unsigned int upsampling = 1);
          
    unsigned int id(unsigned int i, unsigned int j) const;
    nvis::vec3 operator()(unsigned int i, unsigned int j) const;
    
    nvis::bbox3 bounds;
    nvis::ivec3 size;
    nvis::vec3 d;
    unsigned int dim;
    double z;
};

namespace {
const unsigned int idx[][3] = { {1, 2, 0}, {0, 2, 1}, {0, 1, 2} };
};


// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

inline
FaceId::FaceId(unsigned int i0, unsigned int i1,
               unsigned int i2, unsigned int i3)
    : is(4)
{
    is[0] = i0;
    is[1] = i1;
    is[2] = i2;
    is[3] = i3;
    std::sort(is.begin(), is.end());
}

inline
FaceId::FaceId(const FaceId& fid)
    : is(4)
{
    for (unsigned int i = 0 ; i < 4 ; i++) {
        is[i] = fid.is[i];
    }
}

inline
FaceId& FaceId::operator=(const FaceId& fid)
{
    for (unsigned int i = 0 ; i < 4 ; i++) {
        is[i] = fid.is[i];
    }
    return *this;
}

inline
int FaceId::operator<(const FaceId& fid) const
{
    return (is[0] < fid.is[0] ||
            (is[0] == fid.is[0] &&
             (is[1] < fid.is[1] ||
              (is[1] == fid.is[1] &&
               (is[2] < fid.is[2] ||
                (is[2] == fid.is[2] &&
                 is[3] < fid.is[3]))))));
}

inline
int FaceId::operator==(const FaceId& fid) const
{
    return (is[0] == fid.is[0] &&
            is[1] == fid.is[1] &&
            is[2] == fid.is[2] &&
            is[3] == fid.is[3]);
}

inline
std::ostream& operator<<(std::ostream& os, const spurt::FaceId& fid)
{
    os << "[ " << fid.is[0] << ", " << fid.is[1] << ", "
       << fid.is[2] << ", " << fid.is[3] << "]";
    return os;
}

inline
Grid::Grid(const Nrrd* nrrd, unsigned int upsample)
{
    assert(nrrd->dim >= 3);
    bounds = spurt::nrrd_utils::get_bounds<3>(nrrd);
    
    for (unsigned int i = 0 ; i < 3 ; i++) {
        size[i] = nrrd->axis[nrrd->dim-3+i].size * upsample;
        d[i] = (bounds.max()[i] - bounds.min()[i]) / (size[i] - 1);
    }
}

inline
unsigned int
Grid::
id(unsigned int i, unsigned int j, unsigned int k) const
{
    return i + size[0]*(j + size[1]*k);
}

inline
nvis::vec3
Grid::
operator()(unsigned int i, unsigned int j, unsigned int k) const
{
    return bounds.min() + nvis::vec3(i,j,k)*d;
}

inline
nvis::vec3
Grid::
operator()(const nvis::vec3& p) const
{
    return bounds.min() + p*d;
}

inline
void Grid::voxel(std::vector< nvis::vec3 >& v,
                 unsigned int i, unsigned int j, unsigned int k) const
{
    v.resize(8);
    v[0] = (*this)(i, j, k);
    v[1] = (*this)(i + 1, j, k);
    v[2] = (*this)(i + 1, j + 1, k);
    v[3] = (*this)(i, j + 1, k);
    v[4] = (*this)(i, j, k + 1);
    v[5] = (*this)(i + 1, j, k + 1);
    v[6] = (*this)(i + 1, j + 1, k + 1);
    v[7] = (*this)(i, j + 1, k + 1);
}

inline
Slice::Slice(const Nrrd* nrrd, unsigned int dim, unsigned int pos,
             unsigned int upsample)
{
    this->dim = dim;
    
    NrrdAxisInfo* info;
    assert(nrrd->dim >= 3);
    assert(dim < nrrd->dim);
    assert(pos < upsample*nrrd->axis[nrrd->dim-3+dim].size);
    bounds = spurt::nrrd_utils::get_bounds<3>(nrrd);
    
    for (unsigned int j = 0 ; j < 3 ; j++) {
        unsigned int i = idx[dim][j];
        const NrrdAxisInfo& axis = nrrd->axis[nrrd->dim-3+i];
        
        size[j] = axis.size * upsample;
        d[j] = (bounds.max()[j] - bounds.min()[j]) / (size[j] - 1);
    }
    
    z = bounds.min()[2] + pos * d[2];
}

inline
unsigned int
Slice::
id(unsigned int i, unsigned int j) const
{
    return i + size[0]*j;
}

inline
nvis::vec3
Slice::
operator()(unsigned int i, unsigned int j) const
{
    double u = bounds.min()[0] + i * d[0];
    double v = bounds.min()[1] + j * d[1];
    nvis::vec3 p;
    p[idx[dim][0]] = u;
    p[idx[dim][1]] = v;
    p[idx[dim][2]] = z;
    return p;
}
};

#endif








