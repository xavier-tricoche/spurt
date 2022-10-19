//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:01 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef F_CELL_DEFINITIONS_UNSTRUCTURED_SINGLETYPE_HH
#define F_CELL_DEFINITIONS_UNSTRUCTURED_SINGLETYPE_HH


#include <vector>
#include <string>
#include <algorithm>
#include <boost/shared_ptr.hpp>

#include "FIndex.hh"
#include "FCellDefinitions.hh"
#include "FNeighborhoodDataStructured3D.hh"
#include "FanyMMappedArray.hh"

using namespace boost;
using namespace std;

/**
 * this is a celldefinitions class containing cells of only one type
 *
 * \ingroup TensorField
 */

template<class cellClass,FCell::CellType cellTypeId,unsigned vertsPerCell>
class FCellDefinitionsUnstructuredSingleType:public FCellDefinitions
{
public:

  FCellDefinitionsUnstructuredSingleType( const std::string& newname, 
				     shared_ptr< FanyArray< positive > > vertexIndices,
				     positive noPos)
    :FCellDefinitions(newname),indices(vertexIndices)
  {
    nbCells = vertexIndices->size() / vertsPerCell;
    nbPos = noPos;
  }

  ~FCellDefinitionsUnstructuredSingleType()
  {
    //only to prevent assertion in FCellDefinition
    neighborData = new FNeighborhoodDataStructured3D(1,1,1,false);
  }



  void getCellVerticesIndices( const FIndex& cellId, 
			       std::vector< FIndex >& vertices ) const
  {    
    vertices.resize(vertsPerCell);
  
    positive off=positive(cellId) * vertsPerCell;
  
    for(positive i=0;i<vertsPerCell;i++)
	vertices[i]=(*indices)[i+off];
  }

  inline const vector<FIndex>& getCellVerticesIndices( const FIndex& cellId ) const
  {    
    indbuf.resize(vertsPerCell);
 
    positive off=positive(cellId) * vertsPerCell;

    for(positive i=0;i<vertsPerCell;i++)
        indbuf[i]=(*indices)[i+off];
    
    return indbuf;    
  }

  /**
   * get cell with only vertex indices set
   */
  shared_ptr<FCell> getCellTorso( const FIndex& cellId ) const
  {
    return shared_ptr<FCell> 
      ( new cellClass(  getCellVerticesIndices(cellId)	 ) );
	
  }

  /** 
   *\par Description:
   *Returns the type of the cell with cellId in cellType
   */
  void getCellType( const FIndex& cellId, 
		    FCell::CellType& cellType) const
  { cellType=cellTypeId; }
  
  
  positive memSize() const
  { return indices->size()*sizeof(positive); }

  void getDistribution(vector<positive>& sizes,
		       vector<string> & names) const
  {
    const FanyMMappedArray<positive> * a =
      dynamic_cast<const FanyMMappedArray<positive> * >( indices.get() );
    
    if(!a)
      THROW_EXCEPTION(FException,"not a distributed array" );
    
    a->arr->getBlocksPerFile(sizes);

    for(unsigned i=0;i<sizes.size();i++)
      sizes[i] /= vertsPerCell;

    a->arr->getFileNames(names);
   

  }
  
protected:
  
  shared_ptr< FanyArray<positive> > indices;

  //buffer var for use in getCellVerticesIndices
  mutable vector<FIndex> indbuf;
};















#endif // F_CELL_DEFINITIONS_UNSTRUCTURED_SINGLETYPE_HH


