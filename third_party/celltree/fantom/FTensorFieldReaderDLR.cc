#include "FTensorSet.hh"
#include "FPositionSet3DArbitrary.hh"
#include "FGrid3DArbitrary.hh"
#include "FGrid2Din3DArbitrary.hh"
#include "FCellDefinitions3DUnstructured.hh"
#include "FCellDefinitions2Din3DUnstructured.hh"
#include "FCellDefinitions3DTriangulation.hh"
#include "FTensorField.hh"

#include "FException.hh"
#include "FStatusReporter.hh"

#include "FTensorFieldReader.hh"
#include "FDataSet.hh"

#include <libgen.h>
#include <stdexcept>
#include <fstream>
#include <netcdf.h>
#include <cassert>
#include <iostream>

using namespace std;

namespace {
  const char *NPOINTS = "no_of_points";
  const char *POINTSX = "points_xc";
  const char *POINTSY = "points_yc";
  const char *POINTSZ = "points_zc";

  struct TCellConv
  { 
    const char *const name;
    FCell::CellType type;
  };

  static const TCellConv cellCore[] = {
    { "points_of_hexaeders", FCell::ARBITRARY_HEX },
    { "points_of_tetraeders", FCell::TETRAHEDRON },
    { "points_of_prisms", FCell::PRISM },
    { "points_of_pyramids", FCell::PYRAM },
    { 0, FCell::UNDEFINED } //delimiter
  };

  static const TCellConv cellBdry[] = {
    { "points_of_surfacetriangles", FCell::TRIANGLE_3D },
    { "points_of_surfacequadrilaterals", FCell::QUADRILATERAL_3D },
    { 0, FCell::UNDEFINED } //delimiter
  };

  typedef shared_ptr<FTensorSet>   tsPtr;
  typedef shared_ptr<FTensorField> tfPtr;

  //-----------------------------------------------------------------------------

  inline void check_nc( int status )
  {
    if( status != NC_NOERR )
      throw FException( nc_strerror(status) );
  }

} // anonymous namespace
//-----------------------------------------------------------------------------

namespace FTensorFieldReader {

shared_ptr<FPositionSet> loadPositionSetDLR( const std::string& file )
{
    int id, nPointsDimId, ncid;
    size_t nPoints;

    check_nc( nc_open( file.c_str(), NC_NOWRITE, &ncid ) );
    check_nc( nc_inq_dimid( ncid, NPOINTS, &nPointsDimId ) );
    check_nc( nc_inq_dimlen( ncid, nPointsDimId, &nPoints ) );

    vector<double> points( nPoints*3 );

    size_t count = nPoints,start = 0;
    ptrdiff_t stride = 1, imap = 3;    
    
    theStatusRep->say( "loading position set" );
    
    check_nc( nc_inq_varid( ncid, POINTSX, &id ) );

    check_nc( nc_get_varm_double( ncid, id, &start, &count,
				  &stride, &imap, &points[0] ) );

    check_nc( nc_inq_varid( ncid, POINTSY, &id ) );

    check_nc( nc_get_varm_double( ncid, id, &start, &count,
				  &stride, &imap, &points[1] ) );

    check_nc( nc_inq_varid( ncid, POINTSZ, &id ) );
    
    check_nc( nc_get_varm_double( ncid, id, &start, &count,
				  &stride, &imap, &points[2] ) );
    
    nc_close( ncid );

    return shared_ptr<FPositionSet>( new FPositionSet3DArbitrary(points) );
}

//-----------------------------------------------------------------------------

namespace {
  struct dimLen
  { 
    size_t nCells, nPointsPerCell; 
    int varId; 
  };
}

shared_ptr<FCellDefinitions> loadCellDefsDLR( const std::string& file )
{
    int allCells = 0, allCellsSize = 0, nPointsDimId, ncid;
    size_t nPoints;
    bool onlyTriangles=true;

    vector<dimLen> dimLens;

    const TCellConv *cellConv = cellCore;

    check_nc( nc_open( file.c_str(), NC_NOWRITE, &ncid ) );
    check_nc( nc_inq_dimid( ncid, NPOINTS, &nPointsDimId ) );
    check_nc( nc_inq_dimlen( ncid, nPointsDimId, &nPoints ) );

    theStatusRep->say( "loading cell definitions (dims)" );

    for( int i=0; cellConv[i].name!=0; i++ )
    {
        dimLen dl;
        dl.nCells = 0;
        
        if( nc_inq_varid( ncid, cellConv[i].name, 
			  &dl.varId ) == NC_NOERR )
        {
            if(cellConv[i].name!="points_of_surfacetriangles")
                onlyTriangles=false;
	      
            int dimIds[2];
	  
            check_nc( nc_inq_vardimid( ncid, dl.varId, dimIds ) );
            check_nc( nc_inq_dimlen( ncid, dimIds[0], &dl.nCells) );
            check_nc( nc_inq_dimlen( ncid, dimIds[1], &dl.nPointsPerCell ) );
	  
            allCells += dl.nCells;
            allCellsSize += dl.nPointsPerCell * dl.nCells;
        }
        
        dimLens.push_back( dl );
    }

    if( !allCells )
        return shared_ptr<FCellDefinitions>();
    
    vector<FIndex> allVertices( allCellsSize );
    typedef vector< pair<FCell::CellType, unsigned int> > typeVector;

    typeVector allTypes( allCells+1 );
    typeVector::iterator k = allTypes.begin(), kend;

    unsigned int actInd = 0;
  
    for( int i=0; cellConv[i].name!=0; i++ )
    {
	if( dimLens[i].nCells )
	{
	    size_t start[] = { 0, 0 };
	    size_t count[] = { dimLens[i].nCells, dimLens[i].nPointsPerCell };
    
	    assert( sizeof(long) == sizeof(FIndex) );

	    char name[256];
	    nc_inq_varname( ncid, dimLens[i].varId, name );

	    theStatusRep->say( "loading cell definitions ("+string(name)+")" );

	    check_nc( nc_get_vara_long( ncid, dimLens[i].varId,
					start, count, 
					(long*)&allVertices[actInd] ) );
	    
	    FCell::CellType t = cellConv[i].type;
	    kend = k + dimLens[i].nCells;        

	    for( ; k!=kend; k++,actInd+=dimLens[i].nPointsPerCell )
	    {
		k->first = t;
		k->second = actInd;
	    }
	}
    }

    nc_close( ncid );

    k->first=FCell::UNDEFINED;
    k->second=actInd;

    // if( boundary )
    //   if(onlyTriangles)
	// return shared_ptr<FCellDefinitions>
	//   ( new FCellDefinitions3DTriangulation( nPoints, 
	// 					 "triangular boundary", 
	// 					 allVertices,true) );
    //   else
	// return shared_ptr<FCellDefinitions>
	//   ( new FCellDefinitions2Din3DUnstructured( nPoints, "boundary",
	// 					    allTypes, allVertices) );
    // else
      return shared_ptr<FCellDefinitions>
	( new FCellDefinitions3DUnstructured( nPoints, "core", 
					      allTypes,allVertices ) );
}

//-----------------------------------------------------------------------------

tsPtr loadTensorSetDLR( const std::string& file )
{
    int nvars, nPointsDimId, ncid;
    size_t nPoints;

    check_nc( nc_open( file.c_str(), NC_NOWRITE, &ncid ) );
    check_nc( nc_inq_nvars( ncid, &nvars ) );

    check_nc( nc_inq_dimid( ncid, NPOINTS, &nPointsDimId ) );
    check_nc( nc_inq_dimlen( ncid, nPointsDimId, &nPoints ) );


    for( int i=0; i<nvars; i++ )
    {
	char name[256];
	check_nc( nc_inq_varname( ncid, i, name ) );
		
	int ndims;
	
	check_nc( nc_inq_varndims( ncid, i, &ndims ) );

	assert( ndims<1000 );
	assert( ndims>0 );

	int *dimid = new int[ndims];

	check_nc( nc_inq_vardimid( ncid, i, &(dimid[0]) ) );

	if( dimid[0] == nPointsDimId )
	{
	    if( ndims == 1 )
	    {
            if( strncmp( "x_", name, 2 ) == 0 )
            {
                theStatusRep->say( "loading tensor sets ("+string(name+2)+")" );
		    
                vector<double> data( nPoints*3 );
                size_t start = 0, count = nPoints;
                ptrdiff_t stride = 1, imap = 3;
		    
                check_nc( nc_get_varm_double( ncid, i, &start, &count,
                                              &stride, &imap, &data[0] ) );
		    
                name[0]='y';
                check_nc( nc_inq_varid( ncid, name, &i ) );
		    
                check_nc( nc_get_varm_double( ncid, i, &start, &count,
                                              &stride,&imap,&data[1] ) );
		    
                name[0]='z';
                check_nc( nc_inq_varid( ncid, name, &i ) );
		    
                check_nc( nc_get_varm_double( ncid, i, &start, &count,
                                              &stride, &imap, &data[2] ) );
		    
                // heuristic: any vector field must be the velocity field
                // it is to be the first field added to the data set
                nc_close( ncid );
                return tsPtr( new FTensorSet( 3, 1, data, name+2 ));
            }
            else
            {
                continue;
                // theStatusRep->say( "loading tensor sets ("+string(name)+")" );
		    
                // vector<double> data( nPoints );
		    
                // size_t start = 0, count = nPoints;
		    
                // check_nc( nc_get_vara_double( ncid, i, &start, 
                //                               &count, &data[0] ) ); 
		    
                // tsets.push_back( shared_ptr<FTensorSet>
                //                  ( new FTensorSet( 1, 0, data, name ) ) );
            }
	    }
	    else
	    {
            continue;
            // theStatusRep->say( "loading tensor sets ("+string(name)+")" );

            // size_t tsize=1,tdim;
            // size_t *start = new size_t[ndims];
            // size_t *count = new size_t[ndims];

		
            // check_nc( nc_inq_dimlen( ncid, dimid[1], &tdim ) );
		
            // for( int j=1; j<ndims; ++j )
            // {
            //     size_t tmp;
            //     check_nc( nc_inq_dimlen( ncid, dimid[j], &tmp ) );
            //     assert(tmp==tdim);
            //     tsize *= tdim;

            //     start[j] = 0;
            //     count[j] = tdim;
            // }

            // start[0] = 0;
            // count[0] = nPoints;

            // cout << "tsize: " << tsize << endl;

            // vector<double> data( nPoints*tsize );

            // check_nc( nc_get_vara_double( ncid, i, &start[0], &count[0] ,&data[0] ) );
            // tsets.push_back( shared_ptr<FTensorSet>
            //                  ( new FTensorSet( tdim, ndims-1, data, name ) ) );

            // delete [] start;
            // delete [] count;
	    }
	
	}

    delete [] dimid;

    }
    
    return tsPtr();
}

//-----------------------------------------------------------------------------

tfPtr loadDLR( const string& file )
{
    std::string gridfile, pvalfile;

    std::ifstream in( file.c_str() );
    in >> gridfile >> pvalfile;

    if( !in )
        throw std::runtime_error( "Unable to read DLR file" );
    
    in.close();
    
    std::string path;
    
    {
        char* tmp = strdup( file.c_str() );
        path = dirname( tmp );
        free( tmp );
    }
    
    // gridfile = path + '/' + gridfile;
    // pvalfile = path + '/' + pvalfile;
    
    shared_ptr<FPositionSet> ps = 
        FTensorFieldReader::loadPositionSetDLR( gridfile );

    shared_ptr<FCellDefinitions> cd =
        FTensorFieldReader::loadCellDefsDLR( gridfile );

    shared_ptr<FTensorSet> ts =
        FTensorFieldReader::loadTensorSetDLR( pvalfile );

    shared_ptr<FGrid> grid = FGrid::constructSuitableGrid( "inner grid", cd, ps );

    return tfPtr( new FTensorField( "field", ts, grid ) );
}

}
