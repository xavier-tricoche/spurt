#include <cstdio>
#include <dataset.hpp>
#include <celltree.hpp>
#include <celltree_builder.hpp>
#include <interpolator.hpp>

#include "timer.hpp"
#include "cuda_interpolator.hpp"

// -------------------------------------------------------------------------

int main( int argc, char** argv )
{
    if( argc < 2 )
    {
        printf( "usage: %s <dataset>\n", argv[0] );
        exit( EXIT_FAILURE );
    }
    
    try 
    {
        std::string datafile = argv[1];
        std::string prefix = datafile.substr( 0, datafile.rfind( '.' ) );

        // load dataset and create celltree
        dataset* ds = dataset::create( datafile );
        
        const mesh*     m = ds->read_mesh();
        const variable* v = ds->read_vector_variable( 0, "velocity" );
              
        celltree ct;
        celltree_builder builder;
        builder.build( ct, *m );
        
        // read point list from file
        std::ifstream in( (prefix+".rpts").c_str() );

        if( !in.good() )
	        throw std::runtime_error( "could not read point file " + 
	                                  prefix + ".rpts" );
        
        std::vector<float3> points;

        while( true )
        {
            float3 p;
            in >> p.x >> p.y >> p.z;
            
            if( !in.good() )
                break;
                
            points.push_back( p );
        }

        in.close();
        
        // interpolate points on device
        thrust::host_vector<float3> result_gpu( points.size() );
        
        {
            cuda_interpolator ci( m, v, &ct );
            progress pr( "interpolating on GPU", points.size() );

            event_timer timer;

            thrust::device_vector<float3> dpoints = points;
            thrust::device_vector<float3> dresult( points.size() );
        
            timer.tic();
        
            thrust::transform( dpoints.begin(), dpoints.end(), 
                               dresult.begin(), ci.get_interpolator() );
                        
            timer.toc();
                               
            result_gpu = dresult;
            
            fprintf( stdout, "interpolating on GPU: %7.2fms (%.2fMeval/s)\n",
                     timer.elapsed(), 
                     (double)points.size() / timer.elapsed() * 1e-3 );
        }
        
        // interpolate points on host
        std::vector<float3> result_cpu( points.size() );
        
        {
            interpolator intp( m, v, ct );
            progress pr( "interpolating on CPU", points.size() );
            
            for( unsigned int i=0; i<points.size(); ++i, ++pr )
            {
                if( !intp( 0.0, &points[i].x, &result_cpu[i].x ) )
                    result_cpu[i] = make_float3( 0, 0, 0 );
            }
            
            fprintf( stdout, "interpolating on CPU: %7.2fms (%.2fMeval/s)\n",
                     pr.elapsed() * 1e3,
                     (double)points.size() / pr.elapsed() * 1e-6 );
        }
        
        // compare results
        unsigned int errorHist[8] = { 0, 0, 0, 0, 0, 0, 0, 0 };
        double       errorSum = 0.0;

        for( unsigned int i=0; i<points.size(); ++i )
        {
            if( length( result_cpu[i] ) > 0 )
            {
                float dist = length( result_cpu[i] - result_gpu[i] ) /
                             length( result_cpu[i] );

                int log10dist = floor( -log10( dist ) );
            
                if( log10dist < 0 ) 
                    log10dist = 0;
                if( log10dist > 7 )
                    log10dist = 7;
            
                ++errorHist[log10dist];
            
                errorSum += dist;
            }
        }

        fprintf( stdout, "comparison:\n" );
        
        for( unsigned int i=0; i<8; ++i )
            fprintf( stdout, "    rel. difference < 1e-%u: %u\n", 
                     i, errorHist[i] );

        fprintf( stdout, "    avg. difference       : %.2e\n",
                 errorSum / points.size() );
    }
    catch( std::exception& e )
    {
        printf( "exception: %s\n", e.what() );
        exit( EXIT_FAILURE );
    }
    
    exit( EXIT_SUCCESS );    
}