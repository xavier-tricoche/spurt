//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile$
// Language:  C++
// Date:      $Date$
// Author:    $Author$
// Version:   $Revision$
//
//---------------------------------------------------------------------------

#include "FCellDefinitions3DLines.hh"
#include "FNeighborhoodDataUnstructured.hh"

#include "FLineCell3D.hh"

#include <cassert>
#include "eassert.hh"
using namespace std;

//---------------------------------------------------------------------------

FCellDefinitions3DLines::
FCellDefinitions3DLines( positive nb_pos,
				 const std::string& newname,
				 std::vector<FIndex> & vertices,
				 bool /*buildneighborhooddata*/ )
    : FCellDefinitions( newname )
{
#ifndef NODEBUG
    cout << "nb of vertices in cell def= " << vertices.size() << endl;
#endif

    eassert( vertices.size() );
    eassert( !( vertices.size()%2) );
    eassert( nb_pos );

    nbCells = vertices.size()/2;
    nbPos = nb_pos;

    cell_vertices.swap( vertices );

//    if ( buildneighborhooddata )
//	neighborData
//	    = new FNeighborhoodDataLines( this, cell_vertices );
//
    neighborData = new FNeighborhoodDataUnstructured(this );
}

//---------------------------------------------------------------------------

FCellDefinitions3DLines::~FCellDefinitions3DLines()
{
    // neighborData is deleted by base class
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions3DLines::getClassName() const
{
  static FString name("FCellDefinitions3DLines");
  return name;
}


//---------------------------------------------------------------------------

void FCellDefinitions3DLines::getCellVerticesIndices( const FIndex& cellId,
							      std::vector< FIndex >& vertices ) const
{
    eassert( cellId < nbCells );

    vertices.clear();

    for( positive i=cellId.getIndex()*2 ; i<(cellId.getIndex()+1)*2 ; ++i )
	vertices.push_back( cell_vertices[i] );
}

//---------------------------------------------------------------------------

void FCellDefinitions3DLines::getCellType( const FIndex& cellId,
						   FCell::CellType& cell_type ) const
{
    eassert( cellId < nbCells );

    cell_type = FCell::LINE_3D;
}

//---------------------------------------------------------------------------

shared_ptr<FCell> FCellDefinitions3DLines::getCellTorso(const FIndex& cellId) const
{
    eassert( cellId < nbCells );

    return shared_ptr<FCell>
	( new FLineCell3D
	  ( vector< FIndex >
	    (cell_vertices.begin()+2*cellId.getIndex(),
	     cell_vertices.begin()+2*( cellId.getIndex()+1 ) )
	      )
	    );
}

//---------------------------------------------------------------------------

positive FCellDefinitions3DLines::memSize() const
{
    return
      cell_vertices.size()*sizeof(cell_vertices[0])
      +
      neighborData->memSize()
      +
      sizeof(*this);
}





//---------------------------------------------------------------------------


//---------------------------------------------------------------------------

FCellDefinitions3DLineStrips::
FCellDefinitions3DLineStrips( positive nb_pos,
				 const std::string& newname,
				 std::vector<FIndex> & vertices,
				 bool /*buildneighborhooddata*/ )
    : FCellDefinitions( newname )
{
#ifndef NODEBUG
    cout << "nb of vertices in cell def= " << vertices.size() << endl;
#endif

    eassert( vertices.size() );
    eassert( nb_pos );

    nbPos = nb_pos;

    cell_vertices.swap( vertices );

    nbCells = 0;
    unsigned int i;
    for ( i=0; i< cell_vertices.size(); )
    {
      nbCells += cell_vertices[ i ]-FIndex( 1 );
      starts.push_back( i );
      for ( unsigned int j=i+1; j < i + cell_vertices[ i ]; ++j )
      {
        cells.push_back( j );
        eassert( cell_vertices[ j ] < nbPos );
      }
      eassert( cell_vertices[ i+ cell_vertices[ i ] ] < nbPos );
      i+= cell_vertices[ i ] + FIndex( 1 );
    }
    eassert( i == cell_vertices.size() );
    eassert( cells.size() == nbCells );

//    if ( buildneighborhooddata )
//	neighborData
//	    = new FNeighborhoodDataLines( this, cell_vertices );
//
    neighborData = new FNeighborhoodDataUnstructured(this );
}

FCellDefinitions3DLineStrips::
FCellDefinitions3DLineStrips( positive nb_pos,
				 const std::string& newname,
				 std::vector<unsigned int> & vertices,
				 bool /*buildneighborhooddata*/ )
    : FCellDefinitions( newname )
{
#ifndef NODEBUG
    cout << "nb of vertices in cell def= " << vertices.size() << endl;
#endif

    eassert( vertices.size() );
    eassert( nb_pos );

    nbPos = nb_pos;

    std::copy( vertices.begin(), vertices.end(), std::back_inserter( cell_vertices ) );
    vertices.clear();

    nbCells = 0;
    unsigned int i;
    for ( i=0; i< cell_vertices.size(); )
    {
      nbCells += cell_vertices[ i ]-FIndex( 1 );
      starts.push_back( i );
      for ( unsigned int j=i+1; j < i + cell_vertices[ i ]; ++j )
      {
        cells.push_back( j );
        eassert( cell_vertices[ j ] < nbPos );
      }
      eassert( cell_vertices[ i+ cell_vertices[ i ] ] < nbPos );
      i+= cell_vertices[ i ] + FIndex( 1 );
    }
    eassert( i == cell_vertices.size() );

    eassert( cells.size() == nbCells );

//    if ( buildneighborhooddata )
//	neighborData
//	    = new FNeighborhoodDataLines( this, cell_vertices );

    std::cout << "Cells: " << std::endl;
    for ( unsigned int i=0; i< cells.size(); ++i )
    {
      std::cout << cells[ i ] << std::endl;
    }

    std::cout << "Starts: " << std::endl;
    for ( unsigned int i=0; i< starts.size(); ++i )
    {
      std::cout << starts[ i ] << std::endl;
    }

    std::vector<FIndex> v;
    for ( unsigned int i=0; i< nbCells; ++i )
    {
      std::cout << "cellid=" << i << std::endl;
      getCellVerticesIndices( i, v );
      for ( unsigned int j=0; j< v.size(); ++j )
        std::cout << v[ j ] << " ";
      std::cout << std::endl;
    }
    neighborData = new FNeighborhoodDataUnstructured(this );
}
//---------------------------------------------------------------------------

FCellDefinitions3DLineStrips::~FCellDefinitions3DLineStrips()
{
    // neighborData is deleted by base class
}

//---------------------------------------------------------------------------

void FCellDefinitions3DLineStrips::getCellVerticesIndices( const FIndex& cellId,
							      std::vector< FIndex >& vertices ) const
{
    eassert( cellId < nbCells );

    vertices.resize(2);

	vertices[ 0 ] = cell_vertices[ cells[ cellId.getIndex() ] ];
	vertices[ 1 ] = cell_vertices[ cells[ cellId.getIndex() ] + FIndex( 1 ) ];
}

//---------------------------------------------------------------------------

void FCellDefinitions3DLineStrips::getCellType( const FIndex& cellId,
						   FCell::CellType& cell_type ) const
{
    eassert( cellId < nbCells );

    cell_type = FCell::LINE_3D;
}

//---------------------------------------------------------------------------

shared_ptr<FCell> FCellDefinitions3DLineStrips::getCellTorso(const FIndex& cellId) const
{
    eassert( cellId < nbCells );

    return shared_ptr<FCell>
	( new FLineCell3D( &cell_vertices[ cells[ cellId.getIndex() ] ] ) );
//	  ( vector< FIndex >
//	    (cell_vertices.begin()+cells[ cellId.getIndex() ],
//         cell_vertices.begin()+cells[ cellId.getIndex() ] + 1 )
//	      )
//	    );
}

//---------------------------------------------------------------------------

positive FCellDefinitions3DLineStrips::memSize() const
{
    return
      cell_vertices.size()*sizeof(cell_vertices[0])
      +
      starts.size()*sizeof( starts[ 0 ] )
      +
      cells.size()*sizeof( cells[ 0 ] )
      +
      neighborData->memSize()
      +
      sizeof(*this);
}
