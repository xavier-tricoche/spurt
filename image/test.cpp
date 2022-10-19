#include <vector>
#include <iostream>
#include "convert.hpp"
#include "probe.hpp"

using namespace spurt;
using namespace gage_interface;

int main( int argc, char* argv[] )
{
    unsigned int w, h;
    float dx, dy;

    w = atoi( argv[1] );
    h = atoi( argv[2] );

    dx = atof( argv[3] );
    dy = atof( argv[4] );

    std::cout << "w = " << w << ", h = " << h << ", dx = " << dx 
        << ", dy = " << dy << std::endl;

    std::vector< double > data( w*h );
    for ( unsigned int i=0 ; i<w ; i++ )
        for ( unsigned int j=0 ; j<h ; j++ )
    {
        data[i+j*w] = ( double )i/( double )( w-1 );
    }

    double minx, maxx, miny, maxy;
    minx = miny = 0;
    maxx = dx*( w-1 );
    maxy = dy*( h-1 );

    Nrrd *nrrd = nrrdNew();
    vector2nrrd( nrrd, data, w, h, dx, dy );

    scalar_wrapper image( nrrd );

    srand48( time(0) );
    vec2 x;
    double val, err;

    std::cout << "checking interpolation" << std::endl;
    unsigned int N = atoi( argv[5] );
    for ( unsigned int n=0 ; n<N ; n++ )
    {
        double u,v;
        u = drand48();
        v = drand48();

        x[0] = u*(w-1);
        x[1] = v*(h-1);

        image.value( x, val );
        err += fabs( ( val-u )/val );
    }

    err /= ( double )N;
    std::cout << "mean error = " << err << std::endl;

    std::cout << "checking boundaries" << std::endl;
    double _minx, _maxx, _miny, _maxy;
    _minx = _miny = 1000000;
    _maxx = _maxy = -_minx;
    for ( unsigned int i=0 ; i<10*N ; i++ )
    {
        double u = -1+drand48()*(w+1);
        double v = -1+drand48()*(h+1);

        vec2 x( u, v );
        if ( image.value( x, val ) )
        {
            if ( u<_minx ) _minx = u;
            else if ( u>_maxx ) _maxx = u;

            if ( v<_miny ) _miny = v;
            else if ( v>_maxy ) _maxy = v;
        }
    }

    std::cout << "minx = " << _minx << ", maxx = " << _maxx << std::endl;
    std::cout << "miny = " << _miny << ", maxy = " << _maxy << std::endl;

    return 0;
}
