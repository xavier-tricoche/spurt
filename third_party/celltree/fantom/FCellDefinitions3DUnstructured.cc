//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions3DUnstructured.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:03 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#include "FCellDefinitions3DUnstructured.hh"
#include "FNeighborhoodDataUnstructured.hh"

#include "FTetrahedronCell.hh"
#include "FAxisParallelTetCell.hh"
#include "FArbitraryHexahedronCell.hh"
#include "FAxisParallelHexCell.hh"
#include "FPrismCell.hh"
#include "FPyramidCell.hh"
#include "FLineCell3D.hh"

#include "FException.hh"

#include <cassert>
#include <sstream>

using namespace std;

//--------------------------------------------------------------------------- 

FCellDefinitions3DUnstructured::FCellDefinitions3DUnstructured( positive nb_pos, 
								const std::string& name, 
								std::vector< pair<FCell::CellType, unsigned int> >& types,
								vector<FIndex> & vertices )
    : FCellDefinitions( name ), cell_types(), cell_vertices(), 
      nbHexaCells ( 0 ), nbTetraCells( 0 )
{
    nbPos = nb_pos;
    assert( nbPos );

    cell_types.swap( types );
    cell_vertices.swap( vertices );

    // cell_types contains a dummy last cell
    nbCells = cell_types.size()-1;

    neighborData = 0;
	// = new FNeighborhoodDataUnstructured( this, cell_types, cell_vertices );

    // count hexahedrons and tetrahedrons
    vector< pair<FCell::CellType,unsigned int> >::iterator 
	it    = cell_types.begin(),
	itEnd = cell_types.end();

    itEnd--; // last element is no cell

    for( ; it != itEnd; it++ )
    {
      if(it->first == FCell::AXIS_PARALLEL_HEX
	 ||it->first == FCell::AXIS_PARALLEL_TET
	 ||it->first == FCell::AXIS_PARALLEL_TRI_2D
	 ||it->first == FCell::AXIS_PARALLEL_QUAD_2D)
      {
	ostringstream msg;
	msg<<"This celltype (AXIS_PARALLEL_...) is not allowed for unstructured 3D cell definitions!"<<endl;
	msg<<it->first<<endl;
	
	
	THROW_EXCEPTION(FInvalidCellException,msg.str().c_str());
      }
      if( it->first == FCell::TETRAHEDRON )
	nbTetraCells++;
      else
	nbHexaCells++;
    }
}

//--------------------------------------------------------------------------- 

FCellDefinitions3DUnstructured::~FCellDefinitions3DUnstructured()
{
    // neighborData is deleted by base class
}

//--------------------------------------------------------------------------- 

const FString& FCellDefinitions3DUnstructured::getClassName() const
{
  static FString name("FCellDefinitions3DUnstructured");
  return name;
}

//--------------------------------------------------------------------------- 

unsigned int FCellDefinitions3DUnstructured::getNbHexaCells() const
{
    return nbHexaCells;
}

//--------------------------------------------------------------------------- 

unsigned int FCellDefinitions3DUnstructured::getNbTetraCells() const
{
    return nbTetraCells;
}

//--------------------------------------------------------------------------- 

void FCellDefinitions3DUnstructured::getCellVerticesIndices( const FIndex& cellId, 
							     std::vector<FIndex>& vertices ) const
{
    assert( cellId < nbCells );
    
    vertices.clear();
    
    for( positive i=cell_types[cellId].second;
	 i < cell_types[(positive)cellId+1].second ; 
	 // cellIndex + 1 is always valid thanks an additional last element 
	 // in cell_types
	 i++)

	vertices.push_back( cell_vertices[i] );
}

//--------------------------------------------------------------------------- 

void FCellDefinitions3DUnstructured::getCellType( const FIndex& cellId, 
						  FCell::CellType& cell_type ) const
{
    assert( cellId < nbCells );

    cell_type = cell_types[cellId].first;
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions3DUnstructured::memSize() const
{
//     unsigned int tmp;

//     cout << "Size of FCellDefinitions3DUnstructured:" << endl;
//     cout << "cell_types (" << cell_types.capacity() << ") : " 
// 	 << cell_types.capacity() * sizeof( pair< FCell::CellType, unsigned int > )
// 	 << endl;
//     cout << "cell_vertices (" << cell_vertices.capacity() << ") : "
// 	 << cell_vertices.capacity() * sizeof( FIndex )
// 	 << endl;
//     cout << "neighbordata : " << neighborData->size() << endl;
//     tmp = 
// 	cell_types.capacity() * sizeof( pair< FCell::CellType, unsigned int > ) +
// 	cell_vertices.capacity() * sizeof( FIndex ) +
// 	neighborData->size();
    
//     cout << "total is " << tmp << endl << endl;
//     return tmp;

  return 
    cell_vertices.capacity() * sizeof (cell_vertices[0] )
    +
    cell_types.capacity() * sizeof (cell_types[0] )
    +
    (neighborData ? neighborData->memSize() : 0)
    +
    sizeof(*this);
    
}

//--------------------------------------------------------------------------- 

shared_ptr<FCell> FCellDefinitions3DUnstructured::getCellTorso( const FIndex& cellId ) const
{
    assert( cellId < nbCells );

    FCell *cell;

    switch( cell_types[cellId].first )
    {
      // case FCell::AXIS_PARALLEL_TET:
	
// 	cell = new FAxisParallelTetCell
// 	  (  
// 	    vector<FIndex>
// 	    ( cell_vertices.begin() + cell_types[ (positive)cellId   ].second,
// 	      cell_vertices.begin() + cell_types[ (positive)cellId+1 ].second ) 
// 	    );
// 	break;

      case FCell::TETRAHEDRON:
	
	cell = new FTetrahedronCell
	  (  
	    vector<FIndex>
	    ( cell_vertices.begin() + cell_types[ (positive)cellId   ].second,
	      cell_vertices.begin() + cell_types[ (positive)cellId+1 ].second ) 
	    );
	break;
	
      case FCell::PRISM:
	
	// sanity check
	assert( cell_types[(positive)cellId+1].second - cell_types[cellId].second == 6 );
	
	cell =  new FPrismCell ( &(*(cell_vertices.begin() + cell_types[ cellId ] .second)) );
	break;

      case FCell::PYRAM:
      {
	// sanity check
	assert( cell_types[(positive)cellId+1].second - cell_types[cellId].second == 5 );
	
	//vector<FIndex>::const_iterator it
	//  = cell_vertices.begin() + cell_types[cellId].second;
	
	cell = new FPyramidCell
	  ( vector<FIndex> 
	    ( cell_vertices.begin() + cell_types[ (positive)cellId   ].second,
	      cell_vertices.begin() + cell_types[ (positive)cellId+1 ].second ) );
	//	FIndex indices[] = { it[0],it[4],it[1],it[3],it[4],it[2] };
	//cell = new FPrismCell( indices );
	break;
      }
      
      case FCell::ARBITRARY_HEX:
	
	cell = new FArbitraryHexahedronCell
	  (  
	    vector<FIndex> 
	    ( cell_vertices.begin() + cell_types[ (positive)cellId   ].second,
	      cell_vertices.begin() + cell_types[ (positive)cellId+1 ].second ) 
	    );
	break;	
      case FCell::AXIS_PARALLEL_HEX:
	
	cell = new FAxisParallelHexCell
	  (  
	    vector<FIndex> 
	    ( cell_vertices.begin() + cell_types[ (positive)cellId   ].second,
	      cell_vertices.begin() + cell_types[ (positive)cellId+1 ].second ) 
	    );
	break;
      
    // this is a hack useful for CGNS files,
    // we should think of what to do with FLineCell2D/FLineCell3D somwhere else
    //
      case FCell::LINE_3D:
                 cell = new FLineCell3D
                   (
                    vector<FIndex>(
                      cell_vertices.begin() + cell_types[ (positive)cellId ].second,
                      cell_vertices.begin() + cell_types[ (positive)cellId+1 ].second
                      )
                   );
    break;	
      default:
	cout<<"-- for the moment we only support the cell types"<<endl<<"-- for CellDefinitions3DUnstructured::getCellTorso()"
	    <<endl<<"-- FCell::{TETRAHEDRON|PRISM|PYRAM|ARBITRARY_HEX|AXIS_PARALLEL_HEX|LINE_3D}"<<endl;
	cout<<"-- Current cell's type: "<<cell_types[cellId].first<<endl;
	FException e("for the moment we only support the cell types\n\
                      FCell::{TETRAHEDRON|PRISM|PYRAM|ARBITRARY_HEX|AXIS_PARALLEL_HEX|LINE_3D}for CellDefinitions3DUnstructured::getCellTorso()"); 
	throw e;
    }

    return shared_ptr<FCell>( cell );
}

//--------------------------------------------------------------------------- 

const FNeighborhoodData* FCellDefinitions3DUnstructured::getNeighborhoodData() const
{
    if( !neighborData )
        neighborData = new FNeighborhoodDataUnstructured( this, cell_types, cell_vertices );
    
    return neighborData;
}
