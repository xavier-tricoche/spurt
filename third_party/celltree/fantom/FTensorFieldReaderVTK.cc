//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTensorFieldReaderVTK.cc,v $
// Language:  C++
// Date:      $Date: 2004/07/01 13:42:24 $
// Author:    $Author: garth $
// Version:   $Revision: 1.9 $
//
//---------------------------------------------------------------------------

#ifndef NO_DATASET_VTK

// Cell definitions types
#include "FCellDefinitions2DRectilinear.hh"
#include "FCellDefinitions2DStructured.hh"
#include "FCellDefinitions2DRectilinear.hh"
#include "FCellDefinitions2DUnstructured.hh"
#include "FCellDefinitions3DRectilinear.hh"
#include "FCellDefinitions3DStructured.hh"
#include "FCellDefinitions3DRectilinear.hh"
#include "FCellDefinitions3DUnstructured.hh"
#include "FCellDefinitions3DTriangulation.hh"
#include "FCellDefinitions2Din3DUnstructured.hh"
#include "FCellDefinitions2DPoints.hh"
#include "FCellDefinitions3DPoints.hh"
#include "FCellDefinitions3DLines.hh"
#include "FCellDefinitions2DLines.hh"

// Cell types
#include "FQuadrilateralCell2D.hh"
#include "FHexahedronCell.hh"

// PositionSet types
#include "FPositionSet.hh"
#include "FPositionSet2DRectilinear.hh"
#include "FPositionSet3DRectilinear.hh"
#include "FPositionSet2DArbitrary.hh"
#include "FPositionSet3DArbitrary.hh"
#include "FPositionSet2DCurvilinear.hh"
#include "FPositionSet3DCurvilinear.hh"

#include "FPositionDistributionChecker.hh"
#include "FTensorField.hh"
#include "FGrid.hh"

#include "FDummyArray.hh"

#include "FTensorFieldReader.hh"
#include "FStatusReporter.hh"
#include "FException.hh"

#include <cmath>
#include <iostream>
#include <fstream>
#include <cassert>
#ifdef __APPLE__
// basename(...)
#include <libgen.h>
#endif
#ifdef WIN32
// basename(...)
#include "util.hh"
#endif

#include "FEndianess.hh"

#include "eassert.hh"
using namespace std;

#define VTK_VERTEX           1
#define VTK_POLY_VERTEX      2
#define VTK_LINE             3
#define VTK_POLY_LINE        4
#define VTK_TRIANGLE         5
#define VTK_TRIANGLE_STRIP   6
#define VTK_POLYGON          7
#define VTK_PIXEL            8
#define VTK_QUAD             9
#define VTK_TETRA           10
#define VTK_VOXEL           11
#define VTK_HEXAHEDRON      12
#define VTK_WEDGE           13
#define VTK_PYRAMID         14

#define VTK_FANTOM_PRISM    20
#define VTK_FANTOM_PYRAMID  30

//---------------------------------------------------------------------------


#include "parseutils.hh"

using boost::spirit::ch_p;
using boost::spirit::real_p;
using boost::spirit::int_p;
using boost::spirit::uint_p;
using boost::spirit::str_p;
using boost::spirit::assign;
using boost::spirit::append;
using boost::spirit::rule;
using boost::spirit::lexeme_d;
using boost::spirit::phrase_scanner_t;
using boost::spirit::graph_p;

// -------------------------------------------------------------------------
/**
 * anonymous namespace for FTensorFieldReaderVTK
 */
namespace {
template<typename T> std::string to_string( const T& t )
{
    std::ostringstream o;
    o << t;

    return o.str();
}

template<typename T> std::string type_to_string( const T&)
{
  return "unknown";
}

template<> std::string type_to_string( const float& )
{
  return "float";
}

template<> std::string type_to_string( const double& )
{
  return "double";
}

template<> std::string type_to_string( const unsigned int& )
{
  return "unsigned int";
}

template<> std::string type_to_string( const int& )
{
  return "int";
}

// -------------------------------------------------------------------------

/**
 * a parse helper that can parse binary data, too
 */
struct vtk_parse_helper : public parse_helper
{
  public:
    vtk_parse_helper( const std::string&path, bool bigEndian ) : parse_helper( path )
      , bigEndian(bigEndian)
        , eh(bigEndian
            ? FEndianess::EndianessHelper::BigEndian
            : FEndianess::EndianessHelper::LittleEndian )
      {
        binary = false; // init, will be changed later
      }

  private:
    /**
     * \todo binary reading is not completely implemented. some of the data is only parsed
     * in big endian format as we don't handle dataIsBigEndian bit through correctly at the
     * moment.
     *
     * Rectilinear grids and grids described by ORIGIN thingy should at least work in binary, too
     * when stored in big endian format
     */
    template< typename T, typename TTarget >
      void readBinary( std::vector<TTarget> &data, const unsigned int size )
      {
        data.clear();

        T buffer;
        //eassert(sizeof(buffer) == 4);
        char *access = (char*)&buffer;
        //ph.goto_sol();
        const char*pos=  this->get_pos();

        while(data.size() < size )
        {
          if(pos+sizeof(T) > this->end())
            THROW_EXCEPTION( FException, "Not enough binary data in file." );
          for(unsigned int i=0; i< sizeof(T); ++i)
          {
            access[i]=pos[i];
          }
          if(data.size() < 3) std::cout << type_to_string( T()) << "  ** " << buffer;
          eh.toHost(access, sizeof(T));
          if(data.size() < 3) std::cout << "  ** " << buffer << std::endl;
          data.push_back( (TTarget)buffer );
          pos += sizeof(T);
          this->goto_pos( pos );
        }
      }

  public:
    template<typename type>
    void parseArray( std::vector<double> &data, const unsigned int size )
    {
      if(binary)
      {
        readBinary<type, double>(data, size );
      }
      else
      {
        while( data.size() < size )
          this->parse( *real_p[append(data)] );
      }
    }

    template<typename type>
    void parseArray( std::vector<unsigned int> &data, const unsigned int size )
    {
      if(binary)
      {
        readBinary<type, unsigned int>( data, size );
      }
      else
      {
        while(data.size() < size )
          this->parse( *uint_p[append(data)]);
      }
    }

  public:
    bool bigEndian;
    bool binary;
    FEndianess::EndianessHelper eh;
};

struct header_parser
{
    enum _grid_type
    {
	unknown,
	structured,
	unstructured,
	rectilinear,
    structured_points,
    polydata
    } grid_type;

    enum _storage_type
    {
      ASCII, BINARY
    } storage_type;

    header_parser( parse_helper& ph )
    {
	grid_type = unknown;

	ph.goto_sol();
	ph.goto_sol();

	// preamble
    ph.parse(
          str_p( "BINARY")[set(storage_type, BINARY)]
        | str_p( "ASCII" )[set(storage_type, ASCII)]);

	// grid type
	ph.parse(
		 ( str_p( "DATASET" )|
		   str_p( "dataset" )
		   )
		 >>
		 ( str_p( "STRUCTURED_GRID" )
		   [set(grid_type,structured)] |
		   str_p( "UNSTRUCTURED_GRID" )
		   [set(grid_type,unstructured)] |
		   str_p( "POLYDATA" )
		   [set(grid_type,polydata)] |
		   str_p( "RECTILINEAR_GRID" )
		   [set(grid_type,rectilinear)] |
           str_p( "STRUCTURED_POINTS" )
           [set(grid_type,structured_points)]) );
    }
};

// -------------------------------------------------------------------------

// Parse so called "structured points", i.e. an axis aligned rectilinear
// grid with fixed spacing and arbitrary origin.
//
// DATASET STRUCTURED_POINTS
// DIMENSIONS dimx dimy dimz
// ORIGIN ox oy oz
// SPACING sx sy sz
struct structured_points_grid_parser
{
  unsigned int size_x, size_y, size_z;
  double ox, oy, oz, sx, sy, sz;

  std::vector<double> px, py, pz;

  structured_points_grid_parser( parse_helper& ph )
  {
    // grid dimension
    ph.parse( str_p( "DIMENSIONS" ) >>
        uint_p[assign(size_x)] >>
        uint_p[assign(size_y)] >>
        uint_p[assign(size_z)]);

     // aspect ratio
    ph.try_parse( (str_p( "SPACING" ) | str_p( "ASPECT_RATIO" )) >>
        real_p[assign(sx)] >>
        real_p[assign(sy)] >>
        real_p[assign(sz)] );

    // these lines are in some files we get, commented out.
    // as this class does not handle comments, I have to do it this way
    while(ph.try_parse( ch_p('#')))
    {
    ph.goto_sol();
    }

    // origin
    ph.parse( str_p( "ORIGIN" ) >>
        real_p[assign(ox)] >>
        real_p[assign(oy)] >>
        real_p[assign(oz)] );
    std::cout << "got origin " << ox << " " << oy << " " << oz << std::endl;
    // spacings
    ph.try_parse( (str_p( "SPACING" ) | str_p( "ASPECT_RATIO" )) >>
        real_p[assign(sx)] >>
        real_p[assign(sy)] >>
        real_p[assign(sz)] );

    ph.goto_sol();

    // some sanity checks
    if( size_x <= 0 || size_y <= 0 || size_z <= 0)
     throw parse_failure( "STRUCTURED_POINTS: DIMENSIONS must be larger than zero." );
    if( sx <= 0. || sy <= 0. || sz <= 0.)
      throw parse_failure( "STRUCTURED_POINTS: SPACINGs must be larger than zero." );


    // for easier handling later, initialize some point structures
    // for the grid to construct
    for(unsigned int i=0; i< size_x; ++i)
      px.push_back(ox+ sx*i);
    for(unsigned int i=0; i< size_y; ++i)
      py.push_back(oy+ sy*i);
    for(unsigned int i=0; i< size_z; ++i)
      pz.push_back(oz+ sz*i);
  }
};


// -------------------------------------------------------------------------

struct rectilinear_grid_parser
{
    unsigned int size_x, size_y, size_z;
    std::vector<double> xc, yc, zc;

    rectilinear_grid_parser( vtk_parse_helper& ph )
    {
	// grid dimension
	ph.parse( str_p( "DIMENSIONS" ) >>
		  uint_p[assign(size_x)] >>
		  uint_p[assign(size_y)] >>
		  uint_p[assign(size_z)] );

    std::cout << "1" << std::endl;
	// point data
	unsigned int np;
    int type=0;
	ph.parse( (
		   str_p( "X_COORDINATES" )|
		   str_p( "x_coordinate" )
		   )
		 >>
            uint_p[assign(np)] >>
           ( str_p( "float" )[set(type,0)] | str_p( "double" )[set(type, 1)] ) );
      std::cout << "Np=" << np << std::endl;

      if( np != size_x )
          throw parse_failure( "VTK parse failure reading X_COORDINATES" );

      std::cout << "1" << std::endl;
      xc.reserve( np );

      ph.goto_sol();
      if(type == 0)
        ph.parseArray<float>( xc, np );
      else
        ph.parseArray<double>( xc, np );
      eassert( xc.size() == np );
      ph.goto_sol();

      for(unsigned int i=0; i< 10 && i< np; ++i)
        std::cout << "-++ " << xc[i] << std::endl;

      ph.parse( (
		 str_p( "Y_COORDINATES" )|
		 str_p( "y_coordinate" )
		 )
		>>
            uint_p[assign(np)] >>
            ( str_p( "float" )[set(type, 0)] | str_p( "double" )[set(type, 1)] ) );

      std::cout << "2" << std::endl;
      if( np != size_y )
          throw parse_failure( "VTK parse failure reading Y_COORDINATES" );

      yc.reserve( np );
      ph.goto_sol();
      if(type == 0)
        ph.parseArray<float>( yc, np );
      else
        ph.parseArray<double>( yc, np );
      eassert( yc.size() == np );
      ph.goto_sol();

      std::cout << "3" << std::endl;

      ph.parse( (
		 str_p( "Z_COORDINATES" )|
		 str_p( "z_coordinate" )
		 )
		 >>
           uint_p[assign(np)] >>
           ( str_p( "float" )[set(type, 0)] | str_p( "double" )[set(type, 1)] ) );

      if( np != size_z )
          throw parse_failure( "VTK parse failure reading Z_COORDINATES" );

      std::cout << "4" << std::endl;
      zc.reserve( np );
      ph.goto_sol();
      if(type == 0)
        ph.parseArray<float>( zc, np );
      else
        ph.parseArray<double>( zc, np );
      eassert( zc.size() == np );
      ph.goto_sol();

      std::cout << "done grid" << std::endl;
      }
  };

  // -------------------------------------------------------------------------

  struct structured_grid_parser
  {
      unsigned int size_x, size_y, size_z;
      std::vector<double> coords;

      structured_grid_parser( vtk_parse_helper& ph )
      {
      // grid dimension
      ph.parse( str_p( "DIMENSIONS" ) >>
           uint_p[assign(size_x)] >>
           uint_p[assign(size_y)] >>
           uint_p[assign(size_z)] );

#ifndef NODEBUG
      cout << "dim ok " << size_x << ' ' << size_y << ' ' << size_z << '\n';
#endif
      // point data
      unsigned int np;
      int type;
      ph.parse( str_p( "POINTS" ) >>
           uint_p[assign(np)] >>
           ( str_p( "float" )[set(type, 0)] | str_p( "double" )[set(type,1)] ) );

#ifndef NODEBUG
      cout << "points ok " << np << '\n';
#endif

      coords.reserve( 3*np );

      if(type == 0)
        ph.parseArray<float>( coords, np*3);
      else
        ph.parseArray<double>( coords, np*3);
      eassert( coords.size() == np*3);

      ph.goto_sol();
      }
  };

  // -------------------------------------------------------------------------

/**
 * Struct for parsing unstructtured grid for VTK ASCII
 */

struct unstructured_grid_parser
{
  unsigned int nbPos, nbCells;
  std::vector<double>       coords;
  std::vector<unsigned int> cells;
  std::vector<unsigned int> celltypes;

    /**
     * Funtion for parsing unstructtured grid for VTK ASCII
     */
    unstructured_grid_parser( vtk_parse_helper& ph )
      {
        int type=0;
      ph.parse( str_p( "POINTS" ) >>
            uint_p[assign(nbPos)] >>
            ( str_p( "float" )[set(type, 0 )] | str_p( "double" )[set(type, 1)] ) );

      ph.goto_sol();
      /**
       * \todo parse binary here
       */
      coords.reserve( 3*nbPos );

      if(type == 0 )
        ph.parseArray<float>( coords, 3*nbPos );
      else
        ph.parseArray<double>( coords, 3*nbPos );
      eassert( coords.size() == 3*nbPos );
      ph.goto_sol();

      unsigned int nbInts;
      ph.parse( (str_p( "CELLS" ) | str_p( "POLYGONS" ) ) >>
            uint_p[assign(nbCells)] >>
            uint_p[assign(nbInts)] );

      std::cout << "CELLS " << nbCells << " " << nbInts << std::endl;
      ph.goto_sol();
      cells.reserve( nbInts );

      ph.parseArray<unsigned int>( cells, nbInts );
      eassert( cells.size() == nbInts );

      ph.goto_sol();
      ph.parse( str_p( "CELL_TYPES" ) >>
            uint_p[assign(nbInts)] );

      std::cout << "CELL_TYPES " << nbInts << std::endl;
      ph.goto_sol();

      if( nbInts != nbCells )
          throw parse_failure( "VTK parse failure reading CELL_TYPES" );

      celltypes.reserve( nbInts );

      ph.parseArray<unsigned int>( celltypes, nbInts );
      eassert( celltypes.size() == nbInts );

      ph.goto_sol();
      }
  };

  struct polydata_parser
  {
    unsigned int nbPos, nbCells, nbVertices, nbLines, nbStrips;
    std::vector<double> coords;

    std::vector<unsigned int> vertices;
    std::vector<unsigned int> lines;
    std::vector<unsigned int> triangle_strips;

    std::vector<unsigned int> cells;
    std::vector<unsigned int> celltypes;

    polydata_parser( vtk_parse_helper& ph )
    {

      int type;
      /** \todo parse binary here */
      if ( ph.try_parse( str_p( "POINTS" ) >>
          uint_p[assign(nbPos)] >>
          ( str_p( "float" )[set(type, 0)]| str_p( "double" )[set(type,1)] ) )
          )
      {
        coords.reserve( 3*nbPos );

        if ( type == 0)
          ph.parseArray<float>(coords, 3*nbPos );
        else
          ph.parseArray<double>(coords, 3*nbPos );
      }
      else
      {
        throw parse_failure( "cannot parse POINTS in file" );
      }

      unsigned int nbInts;
      ph.goto_sol();
      if ( ph.try_parse( str_p( "VERTICES") >> uint_p[assign(nbVertices)] >> uint_p[assign(nbInts)] ) )
      {
        std::cout << "parsing VERTICES" << std::endl;
        vertices.reserve( nbInts );
        ph.parseArray<unsigned int>( vertices, nbInts );
        ph.goto_sol();
      }
      if ( ph.try_parse( str_p( "LINES") >> uint_p[assign(nbLines)] >> uint_p[assign(nbInts)] ) )
      {
        std::cout << "parsing LINES" << std::endl;
        lines.reserve( nbInts );
        ph.parseArray<unsigned int>( lines, nbInts );
        ph.goto_sol();
      }
      if ( ph.try_parse( str_p( "POLYGONS") >> uint_p[assign(nbLines)] >> uint_p[assign(nbInts)] ) )
      {
        cells.reserve( nbInts );
        ph.parseArray<unsigned int>( cells, nbInts );
        ph.goto_sol();
      }
      if ( ph.try_parse( str_p( "TRIANGLE_STRIPS") >> uint_p[assign(nbLines)] >> uint_p[assign(nbInts)] ) )
      {
        triangle_strips.reserve( nbInts );
        ph.parseArray<unsigned int>( triangle_strips, nbInts );
      }
    }
  };

  // -------------------------------------------------------------------------

  /** This reads our point data.
   * Currently the following VTK types are supported:
   * SCALARS: simple, single valued data
   * VECTORS: 2 or 3 Components, depending on third component of Position Set (see flat)
   * TENSORS: 4 or 9 Components, second order tensor. (depending on third component
   *          of position set (see flat)
   * NORMALS: handled just like vectors, we don't care about normalization as proposed
   *          by VTK documentation
   * TEXTURE_COORDINATES: handled as ARRAY of arbitrary dimension (see below)
   *
   * Additions to VTK File format only used in FAnToM:
   * MULTIVECTORS: Clifford algebra object, 4 or 8 components, depending on flat
   * ARRAYS:  Vector of arbitrary dimension. Can be used to store any of the previous types
   *          of data. Syntax same as TEXTURE_COORDINATES:
   *          ARRAY dim {float|double}
   *
   * valid examples paresed here are:
   * POINT_DATA 1234
   * SCALARS some_funny_name iGnoReD_DataTyPE
   * 0 1 2 3 (...) 1233
   *
   * or
   *
   * POINT_DATA 3
   * TEXTURE_COORDINATES texme_name 2
   * 10 11
   * 20 21
   * 30 31
   *
   */
  struct point_data_parser
  {
      enum _data_type
      {
	scalar,
	vector,
	tensor,
	multivector,
	array
      } data_type;

      int dim;
      int order;

      unsigned int comp;
      unsigned int        nbEntries;
      std::vector<double> data;
      std::string name;

      bool empty;

      enum type_
      {
        t_short,
        t_char,
        t_double,
        t_float
      };
      type_ type ;


      point_data_parser( vtk_parse_helper& ph,
                 bool flat)
        : empty( false )
    {
      if ( !
	ph.try_parse( (
		   str_p( "POINT_DATA" )|
		   str_p( "point_data" )
		   )
		  >>
		  uint_p[assign(nbEntries)] )
    )
      {
        std::cerr << "no POINT_DATA field found, maybe that's intented..." << std::endl;
        empty = true;
        return;
      }
    // for name parsing: turn of whitespace skipping and accept all graphical chars
    rule< phrase_scanner_t > np = lexeme_d[+(graph_p)];

    // currently numComp used in Scalars is skipped and ignored (i.e. assumed to be 1)
    // as the syntax would really become ugly then.
	ph.parse(
	  (( str_p( "SCALARS" ) | str_p( "SCALAR" ))[set(data_type,scalar)]
	     >> (np)[assign(name)]
         >> (str_p("short")[set(type, t_short)]
             |str_p("char")[set(type, t_char)]
             |str_p("double")[set(type, t_double)]
             |str_p("float")[set(type, t_float)]))
	  | (( str_p( "VECTORS" ) | str_p( "VECTOR" ))[set(data_type,vector)]
	     >> (np)[assign(name)]
         >> (str_p("short")[set(type, t_short)]
             |str_p("char")[set(type, t_char)]
             |str_p("double")[set(type, t_double)]
             |str_p("float")[set(type, t_float)])	)
	  | (( str_p( "TENSORS" )  | str_p( "TENSOR" ))[set(data_type,tensor)]
	     >> (np)[assign(name)]
         >> (str_p("short")[set(type, t_short)]
             |str_p("char")[set(type, t_char)]
             |str_p("double")[set(type, t_double)]
             |str_p("float")[set(type, t_float)]))
             | (( str_p( "NORMALS" )  | str_p( "NORMAL" ))[set(data_type,vector)]
	     >> (np)[assign(name)]
         >> (str_p("short")[set(type, t_short)]
             |str_p("char")[set(type, t_char)]
             |str_p("double")[set(type, t_double)]
             |str_p("float")[set(type, t_float)]))
             | (( str_p( "COLOR_SCALARS" ) | str_p( "COLOR_SCALAR" ))[set(data_type,array), set(type, t_char)]
	     >> (np)[assign(name)] >> int_p[assign(comp)])
             | (( str_p( "TEXTURE_COORDINATES") | str_p("TEXTURE_COORDINATE"))[set(data_type,array)]
	     >> (np)[assign(name)] >> int_p[assign(comp)]
         >> (str_p("short")[set(type, t_short)]
             |str_p("char")[set(type, t_char)]
             |str_p("double")[set(type, t_double)]
             |str_p("float")[set(type, t_float)]))
      | (( str_p( "MULTIVECTORS" ) | str_p( "MULTIVECTOR" ))[set(data_type,multivector)]
	     >> (np)[assign(name)]
         >> (str_p("short")[set(type, t_short)]
             |str_p("char")[set(type, t_char)]
             |str_p("double")[set(type, t_double)]
             |str_p("float")[set(type, t_float)]))
             | (( str_p( "ARRAYS" ) | str_p( "ARRAY" ))[set(data_type,array)]
	     >> (np)[assign(name)] >> int_p[assign(comp)]
         >> (str_p("short")[set(type, t_short)]
             |str_p("char")[set(type, t_char)]
             |str_p("double")[set(type, t_double)]
             |str_p("float")[set(type, t_float)]))
          );
	// usually here is a data type field specifying float|double|char, we ignore that
    // hoping that it is in the same line ...(!)
	ph.goto_sol();

	if( ph.try_parse( (
			   str_p("LOOKUP_TABLE")|
			   str_p("lookup_table")
			   )
			  ) )
	    ph.goto_sol();

#ifndef NODEBUG
    std::cout << "name: " << name << std::endl;
#endif
	switch( data_type )
	{
	case scalar:
	    comp = 1; dim = 1; order = 0;
	    break;
	case vector:
	    dim = 3; comp = dim; order = 1;
	    break;
	case tensor:
	    dim = flat ? 2 : 3; comp =  dim * dim; order = 2;
	    break;
	case multivector:
	    dim = flat ? 4 : 8; comp = dim; order = 1;
	    break;
  case array:
      dim = comp; order = 1;// comp already set
      break;
	}

	data.reserve( comp*nbEntries );

#ifndef NODEBUG
	cout << "data ok " << comp << ' ' << nbEntries << '\n';
#endif
    {
      // read data binary
      switch(type)
      {
        case t_char:
          {
            ph.parseArray<char>( data, comp*nbEntries );
            //readBinary<char>( ph, dataIsBigEndian, data, comp*nbEntries );
            eassert( data.size() == comp*nbEntries );
          }break;
        case t_double:
          {
            ph.parseArray<double>( data, comp*nbEntries );
            //readBinary<float>( ph, dataIsBigEndian, data, comp*nbEntries );
            eassert( data.size() == comp*nbEntries );
          }break;
        case t_float:
          {
            // FLOAT (4 Byte)
            ph.parseArray<float>( data, comp*nbEntries );
            //readBinary<float>( ph, dataIsBigEndian, data, comp*nbEntries );
          } break;
        case t_short:
          {
            // SHORT (2 Byte)
            ph.parseArray<short>(data, comp*nbEntries );
          } break;
        default:
          THROW_EXCEPTION( FNotImplementedException, "Data type not implemented for binary, currently only short and float supported");
      }
    }

    if( data.size() < comp*nbEntries)
    {
      switch ( data_type )
      {
        case vector:
          std::cout << 2*nbEntries << std::endl;
          std::cout << comp*nbEntries << std::endl;
          std::cout << data.size() << std::endl;
          if( data.size() == 2*nbEntries )
          {
            std::cout << "dim == 2, OK. fixed." << std::endl;
            //data.resize( 2*nbEntries );
            dim = 2; comp = 2;
          }
          else
            THROW_EXCEPTION(FException, "Dimension=2 of data or number of data values do not match.");
          break;
        default:
          THROW_EXCEPTION(FException, "Dimension of data or number of data values do not match.");
      }
    }

    if(!ph.binary)
	ph.goto_sol();
    };
} ;

// -------------------------------------------------------------------------

} // end anonymous namespace


/* do not use from reader itself, because it needs to skip code by hand */
/* flat forces field to be 2D */
shared_ptr<FTensorSet> FTensorFieldReader::loadVTKTensorSetOnly( const std::string& path, bool flat, bool dataIsBigEndian )
{
    vtk_parse_helper ph( path, dataIsBigEndian );
    header_parser hp( ph );
    ph.binary = (hp.storage_type == header_parser::BINARY);
    const char* pos= ph.get_pos();

    // FIXME: any faster method here than trying out every position?
    while( !ph.try_parse( (
			   str_p( "POINT_DATA" )|
			   str_p( "point_data" )
			   )
			  ))
    {
      ph.goto_pos( ++pos );
    } // go on
    ph.goto_pos( pos );

    theStatusRep->say( "parsing point data" );
    point_data_parser pdp( ph, flat );

    // if 2d vectors, need to discard dummy 3rd component
    if( pdp.data_type == point_data_parser::vector && flat && pdp.dim == 3 )
    {
    	assert( pdp.data.size() % 3 == 0);
    	pdp.dim = 2;
	for( unsigned int i=0, k=0, j=0 ; i<pdp.nbEntries ; i++ )
	{
	    pdp.data[j++] = pdp.data[k++];
	    pdp.data[j++] = pdp.data[k++];
	    k++;
	}

	pdp.data.resize( pdp.dim*pdp.nbEntries );
    }

#ifndef NODEBUG
    cout << "after pdp: " << pdp.dim << ' ' << pdp.nbEntries << '\n';
#endif
     // extract base for grid and field names from the file name
    string base = basename( path.c_str() );
    base += string("(") + pdp.name + ")";

    shared_ptr<FTensorSet> tset( new FTensorSet( pdp.dim, pdp.order, pdp.data, base) );

    return tset;
}

shared_ptr<FTensorField> FTensorFieldReader::loadVTK( const std::string& path,
						      bool triangulate, bool dataIsBigEndian, bool enforce2D  )
{
    // these will be filled in as we go
    // making them shared_ptrs right from the start
    // everything will be freed correctly on error
    shared_ptr<FPositionSet>     pset;
    shared_ptr<FCellDefinitions> cdef;
    shared_ptr<FTensorSet>       tset;

    // flat denotes if the field is 2D
    bool flat = enforce2D;

    // setup parser helper
    vtk_parse_helper ph( path, dataIsBigEndian );

    // header parsing --------------------------------------
    header_parser hp( ph );

    ph.binary = (hp.storage_type == header_parser::BINARY);
#ifndef NODEBUG
    cout << "header ok\n";
    std::cout << "grid type is : " << hp.grid_type << std::endl;
#endif
    // grid parsing ----------------------------------------
    switch( hp.grid_type )
    {
    case header_parser::rectilinear:
    {
	theStatusRep->say( "parsing rectilinear grid" );
	rectilinear_grid_parser rp( ph );

	if( rp.size_z == 1 )
		flat = true;

	if( flat )
	{
	    pset  = shared_ptr<FPositionSet>
		( new FPositionSet2DRectilinear( rp.xc, rp.yc ) );

	    cdef = shared_ptr<FCellDefinitions>
		( new FCellDefinitions2DRectilinear( rp.size_x-1, rp.size_y-1,
						     "2D structured grid",
						     triangulate ) );
	}
	else
	{
	    pset  = shared_ptr<FPositionSet3DRectilinear>
		( new FPositionSet3DRectilinear( rp.xc, rp.yc, rp.zc ) );

	    cdef = shared_ptr<FCellDefinitions3DRectilinear>
		( new FCellDefinitions3DRectilinear( rp.size_x-1,
						     rp.size_y-1,
						     rp.size_z-1,
						     "3D structured grid",
						     triangulate ) );
	}
	break;
    }
    case header_parser::structured_points:
    {
      structured_points_grid_parser sp( ph );
      if( sp.size_z == 1 )
        flat = true;
      if( flat )
      {
        pset  = shared_ptr<FPositionSet>
          ( new FPositionSet2DRectilinear( sp.px, sp.py ) );

        cdef = shared_ptr<FCellDefinitions>
          ( new FCellDefinitions2DRectilinear( sp.size_x-1, sp.size_y-1,
                                               "2D Structured Points",
                                               triangulate ) );
      }
      else
      {


        // ignore flat stuff from above, curently only handle 3D
        pset = shared_ptr<FPositionSet3DRectilinear>
          ( new FPositionSet3DRectilinear( sp.px, sp.py, sp.pz ));
        cdef = shared_ptr<FCellDefinitions3DRectilinear>
          ( new FCellDefinitions3DRectilinear( sp.size_x-1,
                                               sp.size_y-1,
                                               sp.size_z-1,
                                               "3D Structured Points",
                                               triangulate) );
      }
      break;
    }
    case header_parser::structured:
    {
	theStatusRep->say( "parsing structured grid" );
	structured_grid_parser sp( ph );
#ifndef NODEBUG
	cout << "grid ok\n";
#endif

	if( sp.size_z == 1 )
	    flat = true;

	FPositionDistributionChecker checker;
	checker.initializeDimensions( sp.size_x, sp.size_y, sp.size_z );

	unsigned int n = 3 * sp.size_x * sp.size_y *sp.size_z;

	for( unsigned int i=0; i<n; i+=3 )
	    checker.checkStep( &(sp.coords[i]) );

	bool depend[3][3];
	bool rectilinear = true;
	checker.getDependencies( depend );

	if( checker.isStructured() )
	{
	    unsigned int dim = flat ? 2 : 3;

	    for( unsigned int i=0; rectilinear && i<dim; i++ )
		for( unsigned int j=0; rectilinear && j<dim; j++ )
		    if( i!=j && depend[i][j] )
			rectilinear = false;
	}

#ifndef NODEBUG
	cout << "flat: " << flat << '\n';
#endif

	if( checker.isStructured() && rectilinear )
	{
	    vector<double> xc, yc, zc;

	    checker.getCoordinates1D(xc, 0);
	    checker.getCoordinates1D(yc, 1);

	    if( flat )
	    {
		pset = shared_ptr<FPositionSet>(
		    new FPositionSet2DRectilinear( xc, yc ) );

		cdef = shared_ptr<FCellDefinitions>(
		    new FCellDefinitions2DRectilinear( sp.size_x-1,
						       sp.size_y-1,
						       "2D rectilinear grid",
						       triangulate ) );
	    }
	    else
	    {
		checker.getCoordinates1D(zc, 2);

		pset = shared_ptr<FPositionSet>(
		    new FPositionSet3DRectilinear( xc, yc, zc ) );

		cdef = shared_ptr<FCellDefinitions>(
		    new FCellDefinitions3DRectilinear( sp.size_x-1,
						       sp.size_y-1,
						       sp.size_z-1,
						       "3D rectilinear grid",
						       triangulate ) );
	    }
	}
	else
	{
	    if( flat )
	    {
		unsigned int nbPos = sp.size_x * sp.size_y;
		// resort the comp array to eleminate the z-coordinates for 2D
		for( unsigned int i=0, k=0, j=0 ; i<nbPos ; i++ )
		{
		    sp.coords[j++] = sp.coords[k++];
		    sp.coords[j++] = sp.coords[k++];
		    k++;
		}

		sp.coords.resize( 2*nbPos );

		pset = shared_ptr<FPositionSet>(
		    new FPositionSet2DCurvilinear( sp.coords,
						   sp.size_x,
						   sp.size_y ) );

		cdef = shared_ptr<FCellDefinitions>
		    ( new FCellDefinitions2DStructured( sp.size_x-1,
							sp.size_y-1,
							"2D structured grid",
							triangulate ) );
	    }
	    else
	    {
		pset = shared_ptr<FPositionSet>(
		    new FPositionSet3DCurvilinear( sp.coords,
						   sp.size_x,
						   sp.size_y,
						   sp.size_z ) );

		cdef = shared_ptr<FCellDefinitions>
		    ( new FCellDefinitions3DStructured( sp.size_x-1,
							sp.size_y-1,
							sp.size_z-1,
							"3D structured grid",
							triangulate ) );
	    }
	}

	break;
    }
    case header_parser::polydata:
    {
      theStatusRep->say( "parsing polydata" );
      polydata_parser polydp( ph );
      if ( polydp.lines.size() > 0 )
      {
        cdef = shared_ptr<FCellDefinitions>
          ( new FCellDefinitions3DLineStrips( polydp.coords.size()/3, "vtk_lines", polydp.lines ) );
        pset = shared_ptr<FPositionSet>( new FPositionSet3DArbitrary( polydp.coords ) );
      }
      else if ( polydp.coords.size() > 0 )
      {
		cdef = shared_ptr<FCellDefinitions>
		    ( new FCellDefinitions3DPoints( polydp.coords.size()/3, "vtk_points" ) );
        pset = shared_ptr<FPositionSet>(
        new FPositionSet3DArbitrary( polydp.coords ) );
      }
    }
    break;
    case header_parser::unstructured:
    {
	theStatusRep->say( "parsing unstructured grid" );
	unstructured_grid_parser usp( ph );

	flat = true;

	// if there is any point with a non-zero z-component,
	// the pointset is assumed three-dimensional
	for( unsigned int i=2; i<usp.coords.size(); i+=3 )
	    if( usp.coords[i] != 0.0 )
	    {
		flat = false;
		break;
	    }

	vector< pair<FCell::CellType, unsigned int> > types;
	vector< FIndex > verts;

	// transform usp.cells for FCellDefinitions
	// cells is  : nbVerts v1 v2 ... vN nbVerts v1 ...
	// need verts: v1 v2 ... vN v1 ...

	types.resize( usp.nbCells+1 );
	verts.resize( usp.cells.size() - usp.nbCells );

	// r: read index
	// w: write index
	// c: cell index
	unsigned int r=0, w=0;
	bool triangulated=true;
	bool point_cells =false;
	bool line_cells=false;


	// these are to output warnings if unsupported types are found
	bool hasVTKPixels = false;
	bool hasVTKVoxels = false;

	for( unsigned int c=0; c<usp.nbCells; c++ )
	{
	    unsigned int n = usp.cells[r++];

	    FCell::CellType type;


#define MYASSERT(a,b) {if(((a)!= (b))){ std::ostringstream oss; oss << "Wrong number of values in cell type: " << a << " but should be " << b; THROW_EXCEPTION( FException, oss.str()); }}

	    switch( usp.celltypes[c] )
	    {
	      case VTK_TRIANGLE:
		type = flat ? FCell::TRIANGLE_2D : FCell::TRIANGLE_3D;
		MYASSERT( n, 3 );
		break;
	      case VTK_PIXEL:
		type = flat ? FCell::QUADRILATERAL_2D : FCell::QUADRILATERAL_3D;
		triangulated=false;
		FQuadrilateralCell2D::fromToPixelEnum( &usp.cells[r] );
		hasVTKPixels=true;
		MYASSERT( n, 4 );
		break;
	      case VTK_QUAD:
		type = flat ? FCell::QUADRILATERAL_2D : FCell::QUADRILATERAL_3D;
		triangulated=false;
		MYASSERT( n, 4 );
		break;
	      case VTK_TETRA:
		type = FCell::TETRAHEDRON;
		MYASSERT( n, 4 );
		break;
	      case VTK_VOXEL:
		type = FCell::ARBITRARY_HEX;
		hasVTKVoxels = true;
		FHexahedronCell::fromToVoxelEnum( &usp.cells[r] );
		MYASSERT( n, 8 );
		break;
	      case VTK_HEXAHEDRON:
		type = FCell::ARBITRARY_HEX;
		MYASSERT( n, 8 );
		break;
	      case VTK_WEDGE: // VTK_WEDGE VTK Documentation
	      case VTK_FANTOM_PRISM: // FAnToM-Style
		type = FCell::PRISM;
		MYASSERT( n, 6 );
		break;
	      case VTK_PYRAMID: // VTK Documentation
	      case VTK_FANTOM_PYRAMID: // other VTK Docs and FAnToM-Style
		type = FCell::PYRAM;
		MYASSERT( n, 5 );
		break;
		// the following types are rather untested:
	      case VTK_VERTEX:
		type = flat ? FCell::POINT_2D : FCell::POINT_3D;
		point_cells =true;
		MYASSERT( n, 1 );
        break;
          case VTK_LINE:
        type = flat ? FCell::LINE_2D : FCell::LINE_3D;
        line_cells = true;
        MYASSERT( n, 2 );
		break;
	      default:
		throw std::runtime_error(
		  "unknown cell type " + to_string( usp.celltypes[c] ) );
	    }
#undef MYASSERT
	    types[c].first = type;
	    types[c].second = w;
	eassert( w <= verts.size() - n );
        for( unsigned int i=0; i<n; ++i )
          verts[w++] = usp.cells[r++];
	}


    if(hasVTKPixels)
    {
      std::cerr << "WARNING: VTK_PIXEL read, but enumeration problem may occur, please use VTK_QUAD enumeration of points instead." << std::endl;
    }

    if(hasVTKVoxels)
    {
      std::cerr << "WARNING: VTK_VOXEL read, but enumeration problem may occur, please use VTK_HEXAHEDRON enumeration of points instead." << std::endl;
    }

	// last element of types is bogus but must be set for
	// cell definitions
	types.back() =
	    pair<FCell::CellType,unsigned int>( FCell::UNDEFINED, w );

	if( flat )
	{
	    // resort the comp array to eleminate the z-coordinates for 2D
	    for( unsigned int i=0, k=0, j=0 ; i<usp.nbPos ; i++ )
	    {
		usp.coords[j++] = usp.coords[k++];
		usp.coords[j++] = usp.coords[k++];
		k++;
	    }

	    usp.coords.resize( 2*usp.nbPos );

	    pset = shared_ptr<FPositionSet>(
		new FPositionSet2DArbitrary( usp.coords ) );

        if(!point_cells)
        {
	    cdef.reset( new FCellDefinitions2DUnstructured( usp.nbPos,
						    "2D unstructured grid",
						    types, verts ) );
        }
        else
        {
          cdef.reset( new FCellDefinitions2DPoints( usp.nbPos, "2D Points"));
        }
	}
	else
	{
	    pset = shared_ptr<FPositionSet>(
		new FPositionSet3DArbitrary( usp.coords ) );

	    if( types[0].first == FCell::TRIANGLE_3D || types[0].first == FCell::QUADRILATERAL_3D)
	    {
	      if(triangulated)
		cdef.reset( new FCellDefinitions3DTriangulation( usp.nbPos,
							    "2D in 3D unstructured grid",
							     verts, true ) );
	      else
		cdef.reset( new FCellDefinitions2Din3DUnstructured( usp.nbPos,
							    "2D in 3D unstructured grid",
							    types, verts ) );
	    }
	    else
	    {
          if(!point_cells)
          {
		    cdef.reset( new FCellDefinitions3DUnstructured( usp.nbPos,
							"3D unstructured grid",
							types, verts ) );
          }
          else
          {
            cdef.reset( new FCellDefinitions3DPoints( usp.nbPos, "3D Points" ) );
          }
	    }
	}
	break;
    }
    default:
	throw parse_failure( "VTK parse failure: unknown grid type" );
	break;
    }

    bool tensors3D = false;

     // extract base for grid and field names from the file name
    string base = basename( path.c_str() );


    // data parsing ----------------------------------------
    theStatusRep->say( "parsing point data" );
    point_data_parser pdp( ph, flat && !tensors3D );

    if ( !pdp.empty )
    {
    base += string("(") + pdp.name + ")";
    //unsigned int dim(0), order(0);

    /*
    switch( pdp.data_type )
    {
    case point_data_parser::scalar:
	order = 0; dim = 1; break;
    case point_data_parser::vector:
	order = 1; dim = flat ? 2 : 3; break;
    case point_data_parser::tensor:
	order = 2; dim = flat ? 2 : 3; break;
    case point_data_parser::multivector:
	order = 1; dim = flat ? 4 : 8; break;
    case point_data_parser::array:
  order = 1; dim = pdp.comp; break;
    }
    */


    // if 2d vectors requested by flat-flag, need to discard dummy 3rd component
    if( flat && pdp.data_type == point_data_parser::vector
        && pdp.dim == 3 )
    {
      pdp.dim = 2;
	for( unsigned int i=0, k=0, j=0 ; i<pdp.nbEntries ; i++ )
	{
	    pdp.data[j++] = pdp.data[k++];
	    pdp.data[j++] = pdp.data[k++];
	    k++;
	}


	pdp.data.resize( pdp.dim*pdp.nbEntries );
    }

#ifndef NODEBUG
    cout << "after pdp: " << pdp.dim << ' ' << pdp.nbEntries << '\n';
#endif
    tset = shared_ptr<FTensorSet>(
	new FTensorSet( pdp.dim, pdp.order, pdp.data, base) );
    }
    else
    {
      // a dummy scalar tensor field
      shared_ptr<FDummyArray> da(  new FDummyArray( 0.0, pset->getNbPositions() ) );
      tset.reset( new FTensorSet( 1, 0, da) );
    }
#ifndef NODEBUG
    cout << "number of tensors: "<< tset->getNbTensors() << ' ';
    cout << "number ofositions(): "<< pset->getNbPositions()<<endl;
#endif

    // construct the data structure
    shared_ptr<FGrid> grid = FGrid::constructSuitableGrid( base + " grid",
							   cdef, pset );

    // whoo, we're done
    return shared_ptr<FTensorField>( new FTensorField( base, tset, grid ) );
}


#endif
