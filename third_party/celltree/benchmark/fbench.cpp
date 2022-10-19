#include "FDefaultInterpolator.hpp"
#include "FCellTreeInterpolator.hpp"

#include <FTensorFieldReader.hh>
#include <FTensorField.hh>
#include <FGrid3DArbitrary.hh>
#include <FCellDefinitions.hh>
#include <FPositionSet.hh>
#include <cstdio>

#include <dataset.hpp>
#include <interpolator.hpp>
#include <celltree_builder.hpp>

#include <wall_timer.hpp>
#include <mesh_traits.hpp>

#include "memstat.hpp"
#include "benchmark.hpp"

#include "fantom_mesh_traits.hpp"

template<typename Interpolator>
void run_fantom_benchmarks( shared_ptr<FTensorField> tf, 
                            const std::string& prefix )
{
    typedef mesh_traits< shared_ptr<FTensorField> > mtraits;

    wall_timer timer;
    memstat_reset();

    Interpolator intp( tf );

    fprintf( stdout, "interp:  %.2fs build time\n"
                     "         %ld bytes allocated (%d%% overhead)\n"
                     "         %ld bytes allocated max.\n\n",
                     timer.elapsed(), 
                     memstat_cur(), 
                     (int)rint( 100.0 * memstat_cur() / mtraits::memsize(tf) ),
                     memstat_max() );

    run_benchmarks( intp, prefix );    
}

// --------------------------------------------------------------------------

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
     
        typedef mesh_traits< shared_ptr<FTensorField> > mtraits;
     
        shared_ptr<FTensorField> tf = FTensorFieldReader::loadDLR( argv[1] );
   
        std::size_t datasize = mtraits::memsize( tf );
   
        fprintf( stdout, "\n"
                         "dataset: %ld points\n"
                         "         %ld cells\n"
                         "         %ld bytes\n\n",
                 mtraits::points_size( tf ),
                 mtraits::cells_size( tf ), 
                 datasize );
  
        fprintf( stdout, "FCellTreeInterpolator:\n" );
        run_fantom_benchmarks<FCellTreeInterpolator>( tf, prefix );

        fprintf( stdout, "FDefaultInterpolator:\n" );
        run_fantom_benchmarks<FDefaultInterpolator>( tf, prefix );
    }
    catch( std::exception& e )
    {
        printf( "exception: %s\n", e.what() );
        exit( EXIT_FAILURE );
    }
    
    exit( EXIT_SUCCESS );    
}
