#ifndef __CONVERT_HPP__
#define __CONVERT_HPP__

#include <teem/nrrd.h>
#include <vector>
#include <iostream>

namespace xavier
{
    void vector2nrrd( Nrrd* nrrd, const std::vector< double >& data,
        unsigned int nx, unsigned int ny,
        double dx=1, double dy=1 );

    void vector2nrrd( Nrrd* nrrd, const std::vector< double >& data,
        unsigned int nx, unsigned int ny, unsigned int nz,
        double dx=1, double dy=1, double dz=1 );
};

inline void xavier::
vector2nrrd( Nrrd* nrrd, const std::vector< double >& data,
    unsigned int nx, unsigned int ny,
    double dx, double dy )
{
    if ( nx*ny != data.size() )
    {
        std::cout << "INVALID SIZE in vector2nrrd!!" << std::endl;
        throw "incompatible sizes in vector2nrrd";
    }
    double maxx = dx*( double )( nx-1 );
    double maxy = dy*( double )( ny-1 );
    double dz = std::max( dx, dy );

    double *_copy = ( double * )calloc( 2*nx*ny, sizeof( double ) );
    for ( unsigned int i=0 ; i<data.size() ; i++ )
    {
        _copy[i] = _copy[i+nx*ny] = data[i];
    }

    nrrdWrap_va( nrrd, _copy, nrrdTypeDouble, 3, nx, ny, 2 );  
    nrrdAxisInfoSet_va( nrrd, nrrdAxisInfoKind, nrrdKindSpace, 
        nrrdKindSpace, nrrdKindSpace );
    nrrdAxisInfoSet_va( nrrd, nrrdAxisInfoMin, 0, 0, -dz/2 );
    nrrdAxisInfoSet_va( nrrd, nrrdAxisInfoMax, maxx, maxy, dz/2 );
    nrrdAxisInfoSet_va( nrrd, nrrdAxisInfoSpacing, dx, dy, dz );
    nrrdAxisInfoSet_va( nrrd, nrrdAxisInfoCenter, nrrdCenterNode, 
        nrrdCenterNode, nrrdCenterNode );
}

inline void xavier::
vector2nrrd( Nrrd* nrrd, const std::vector< double >& data,
    unsigned int nx, unsigned int ny, unsigned int nz,
    double dx, double dy, double dz )
{
    if ( nx*ny*nz != data.size() )
    {
        throw "incompatible sizes in vector2nrrd";
    }

    double *_copy = new double[nx*ny*nz];
    std::copy( data.begin(), data.end(), _copy );

    nrrdWrap_va( nrrd, _copy, nrrdTypeDouble, 3, nx, ny, nz );  
    nrrdAxisInfoSet_va( nrrd, nrrdAxisInfoKind, nrrdKindSpace, 
        nrrdKindSpace, nrrdKindSpace );
    nrrdAxisInfoSet_va( nrrd, nrrdAxisInfoSpacing, dx, dy, dz );
    nrrdAxisInfoSet_va( nrrd, nrrdAxisInfoMin, 0, 0, 0 );
}

#endif
