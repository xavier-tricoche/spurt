#include <vtkDataSet.h>
#include <vtkCellLocator.h>
#include <vtkModifiedBSPTree.h>
#include <vtkUnstructuredGrid.h>

#include "vtkCellTree.hpp"
#include "vtkDLRReader.hpp"
#include "vtkInterpolator.hpp"
#include "vtk_mesh_traits.hpp"

#include <wall_timer.hpp>
#include "memstat.hpp"
#include "benchmark.hpp"

// --------------------------------------------------------------------------

template<typename Locator>
void run_vtk_benchmarks( vtkDataSet* ds, const std::string& prefix )
{
    typedef mesh_traits<vtkDataSet*> mtraits;

    wall_timer timer;
    memstat_reset();
    
    Locator* loc = Locator::New();
    loc->SetDataSet( ds );
    loc->LazyEvaluationOff();
    loc->CacheCellBoundsOff();
    loc->BuildLocator();
            
    fprintf( stdout, "interp:  %.2fs build time\n"
                     "         %ld bytes allocated (%d%% overhead)\n"
                     "         %ld bytes allocated max.\n\n",
                     timer.elapsed(), 
                     memstat_cur(), 
                     (int)rint( 100.0 * memstat_cur() / mtraits::memsize( ds ) ),
                     memstat_max() );

    vtkInterpolator* intp = new vtkInterpolator( ds, loc );
      
    run_benchmarks( *intp, prefix );      
    
    delete intp;
    loc->Delete();
    
    printf( " --- debug: %ld bytes still allocated\n", memstat_cur() );
}

// --------------------------------------------------------------------------

int main( int argc, char** argv )
{
    typedef mesh_traits<vtkDataSet*> mtraits;
     
    if( argc < 2 )
    {
        printf( "usage: %s <dataset>\n", argv[0] );
        exit( EXIT_FAILURE );
    }
    
    try 
    {
        std::string datafile = argv[1];
        std::string prefix = datafile.substr( 0, datafile.rfind( '.' ) );
     
        vtkDLRReader* dr = vtkDLRReader::New();
        dr->SetFileName( argv[1]);
        dr->Update();
        
        vtkDataSet* ds = dr->GetOutput();
        ds->Register( 0 );
        dr->Delete();
        
        fprintf( stdout, "\n"
                         "dataset: %ld points\n"
                         "         %ld cells\n"
                         "         %ld bytes\n\n",
                 mtraits::points_size( ds ),
                 mtraits::cells_size( ds ), 
                 mtraits::memsize( ds ) );
 
         fprintf( stdout, "vtkCellTree:\n" );
         run_vtk_benchmarks<vtkCellTree>( ds, prefix );

         fprintf( stdout, "vtkCellLocator:\n" );
         run_vtk_benchmarks<vtkCellLocator>( ds, prefix );

         fprintf( stdout, "vtkModifiedBSPTree:\n" );
         run_vtk_benchmarks<vtkModifiedBSPTree>( ds, prefix );
    }
    catch( std::exception& e )
    {
        printf( "exception: %s\n", e.what() );
        exit( EXIT_FAILURE );
    }
    
    exit( EXIT_SUCCESS );    
}
