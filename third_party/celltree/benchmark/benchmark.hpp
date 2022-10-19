#ifndef __benchmark_hpp
#define __benchmark_hpp

#include <vector>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <stdexcept>

#include <wall_timer.hpp>
#include <progress.hpp>

#include "dopri5.hpp"

// all benchmarks time out after 10 minutes
const double benchmark_timeout = 600.0;

static wall_timer benchmark_timer;

// -------------------------------------------------------------------------

template<typename Intp>
void benchmark_begin( Intp& intp )
{
    benchmark_timer.restart();
}

// -------------------------------------------------------------------------

template<typename Intp>
void benchmark_end( const char* name, Intp& intp, unsigned int neval ) 
{
    double time = benchmark_timer.elapsed();

    if( time < benchmark_timeout )
    {
        double rate = (double)neval / time;
    
        fprintf( stdout, "%20s | %6.2fs | %11.2f | %10u\n", 
                 name, time, rate, neval );
    }
    else
    {
        fprintf( stdout, "%20s | %7s | %11s | %10s\n", 
                 name, "--", "--", "--" );
    }
};

// -------------------------------------------------------------------------

static
void benchmark_header()
{
    fprintf( stdout, "benchmarks: %8s | %7s | %11s | %10s \n",
             "name", "time", "rate", "neval" );
}

// -------------------------------------------------------------------------

template<typename Intp>    
void points( Intp& intp, const std::string& prefix )
{
    typedef typename Intp::coord_type coord_type;
    typedef typename Intp::value_type value_type;
    
    // read points from file
    std::ifstream in( (prefix+".rpts").c_str() );

    if( !in.good() )
	    throw std::runtime_error( "could not read point file " + 
	                              prefix + ".rpts" );

    std::vector<coord_type> pts;
    std::istream_iterator<coord_type> iin(in), ieos;

    std::copy( iin, ieos, std::back_inserter(pts) );
    in.close();

    const unsigned int npts = pts.size() / 3;

    // interpolate at these locations
    benchmark_begin( intp );
    progress pr( "points", npts );

    for( unsigned int i=0; i<npts; ++i, ++pr )
    {
        value_type res[3];
        intp( 0.0f, &pts[3*i], res );
        
        if( pr.elapsed() > benchmark_timeout )
            break;
    }

    pr.finished();
    benchmark_end( "points", intp, npts );
}

// -------------------------------------------------------------------------

template<typename Intp>
void plane( Intp& intp, const std::string& prefix )
{
    typedef typename Intp::coord_type coord_type;
    typedef typename Intp::value_type value_type;
    
    // read in plane description and resolution from file
    std::ifstream in( (prefix + ".plane").c_str() );
    
    unsigned int nx, ny;    
    in >> nx >> ny; 

    coord_type orig[3], vec0[3], vec1[3];
    
    in >> orig[0] >> orig[1] >> orig[2] 
       >> vec0[0] >> vec0[1] >> vec0[2] 
       >> vec1[0] >> vec1[1] >> vec1[2];

    if( in.fail() )
        throw std::runtime_error( "could not read plane information from " +
                                  prefix + ".plane" );
    
    in.close();

    // perform regular sampling over the extent of the plane
    benchmark_begin( intp );
    progress pr( "plane", nx*ny );

    for( int y=0; y<ny; ++y )
    {
        for( int x=0; x<nx; ++x, ++pr )
        {
            coord_type pos[3] = {
                orig[0] + (x+0.5)/nx * vec0[0] + (y+0.5)/ny * vec1[0],
                orig[1] + (x+0.5)/nx * vec0[1] + (y+0.5)/ny * vec1[1],
                orig[2] + (x+0.5)/nx * vec0[2] + (y+0.5)/ny * vec1[2]
            };
            
            value_type res[3];
          
            intp( 0.0f, pos, res );

            if( pr.elapsed() > benchmark_timeout )
                break;
        }

        if( pr.elapsed() > benchmark_timeout )
            break;
    }
    
    pr.finished();
    benchmark_end( "plane", intp, nx*ny );
}

// -------------------------------------------------------------------------

template<typename Intp>
void volume( Intp& intp, const std::string& prefix )
{
    typedef typename Intp::coord_type coord_type;
    typedef typename Intp::value_type value_type;

    // read volume description from file
    std::ifstream in( (prefix + ".vol").c_str() );
    
    unsigned int nx, ny, nz;    
    coord_type min[3], max[3];
    
    in >> nx >> ny >> nz 
       >> min[0] >> max[0] 
       >> min[1] >> max[1] 
       >> min[2] >> max[2];

    if( in.fail() )
        throw std::runtime_error( "could not read plane information from " + 
                                  prefix + ".vol" );
    
    in.close();

    // perform regular sampling over the extent of the volume
    benchmark_begin( intp );
    progress pr( "volume", nx*ny*nz );

    for( int z=0; z<nz; ++z )
    {
        for( int y=0; y<ny; ++y )
        {
            for( int x=0; x<nx; ++x, ++pr )
            {
                coord_type pos[3] = {         
                    min[0] + (max[0]-min[0])*(x+0.5)/nx,
                    min[1] + (max[1]-min[1])*(y+0.5)/ny,
                    min[2] + (max[2]-min[2])*(z+0.5)/nz,
                };
                
                value_type res[3];

                intp( 0.0f, pos, res );

                if( pr.elapsed() > benchmark_timeout )
                    break;
            }

            if( pr.elapsed() > benchmark_timeout )
                break;
        }

        if( pr.elapsed() > benchmark_timeout )
            break;
    }
    
    pr.finished();
    benchmark_end( "volume", intp, nx*ny*nz );
}

// -------------------------------------------------------------------------

template<typename Intp>
void streamlines( Intp& intp, const std::string& prefix )
{
    typedef typename Intp::coord_type coord_type;
    typedef typename Intp::value_type value_type;
    
    // read streamline parameters from file
    unsigned int N, K;
    float h, a[3], b[3];

    std::ifstream in( (prefix+".slseed").c_str() );

    in >> N >> K >> h
       >> a[0] >> a[1] >> a[2]
       >> b[0] >> b[1] >> b[2];

    if( in.fail() )
        throw std::runtime_error( "couldn't read streamline seed line from " +
                                  prefix + ".slseed" );

    in.close();

    // ---

    unsigned int neval = 0;

    benchmark_begin( intp );
    progress pr( "streamlines", N );

    for( int n=0; n<N; ++n, ++pr )
    {
        dopri5<coord_type,3> intg;

        intg.y[0]  = a[0] + (n+0.5)/N * (b[0]-a[0]);
        intg.y[1]  = a[1] + (n+0.5)/N * (b[1]-a[1]);
        intg.y[2]  = a[2] + (n+0.5)/N * (b[2]-a[2]);

        intg.t_max = K*h;
        intg.t     = 0.0;
        intg.h     = 0.0;

        while( true )
        {
            if( intg.step( intp ) != dopri5<coord_type,3>::OK )
                break;

            if( pr.elapsed() > benchmark_timeout )
                break;
        }

        neval += intg.n_eval;

        if( pr.elapsed() > benchmark_timeout )
            break;
    }

    pr.finished();
    benchmark_end( "streamline", intp, neval );
}

// -------------------------------------------------------------------------

template<typename Intp>
void run_benchmarks( Intp& intp, const std::string& prefix )
{
    benchmark_header();
    
    points( intp, prefix );
    plane( intp, prefix );
    volume( intp, prefix );
    streamlines( intp, prefix );
    
    fprintf( stdout, "\n" );
}

#endif // __benchmark_hpp