#include <cstdio>

#include <dataset.hpp>
#include <interpolator.hpp>
#include <celltree_builder.hpp>

#include <wall_timer.hpp>

#include "memstat.hpp"
#include "benchmark.hpp"

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
        
        dataset* ds = dataset::create( datafile );
        
        const mesh*     m = ds->read_mesh();
        const variable* v = ds->read_vector_variable( 0, "velocity" );

        std::size_t datasize = 
            mesh_traits<mesh>::memsize( *m ) + 
            variable_traits<variable>::memsize( *v );

        fprintf( stdout, "\n"
                         "dataset: %u points\n"
                         "         %u cells\n"
                         "         %ld bytes\n\n",
                 m->npoints, m->ncells, datasize );
  
        wall_timer timer;
        memstat_reset();
      
        celltree ct;
        
        {
            celltree_builder builder;
        
            builder.m_buckets  = 5;
            builder.m_leafsize = 8;

            builder.build( ct, *m );
        }

        fprintf( stdout, "interp:  %.2fs build time\n"
                         "         %ld bytes allocated (%d%% overhead)\n"
                         "         %ld bytes allocated max.\n\n",
                         timer.elapsed(), 
                         memstat_cur(), 
                         (int)rint( 100.0 * memstat_cur() / datasize ),
                         memstat_max() );
        
        interpolator intp( m, v, ct );

        // ---
        
        run_benchmarks( intp, prefix );
    }
    catch( std::exception& e )
    {
        printf( "exception: %s\n", e.what() );
        exit( EXIT_FAILURE );
    }
    
    exit( EXIT_SUCCESS );    
}
