//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions3DUnstructured.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:03 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellDefinitions3DUnstructured_hh
#define __FCellDefinitions3DUnstructured_hh

#include "FCellDefinitions.hh"
#include "FCell.hh"

/**
 * This is a celldefinitions class containing arbitrary cells
 *
 * \ingroup TensorField
 */
class FCellDefinitions3DUnstructured : public FCellDefinitions
{
public:

  /**
   * This constructor is quite tricky: Both provided vectors will be
   * taken over by the data structure. The second one must contain 
   * the sequential list of all cell vertices (with redundancy). The 
   * first one provides all cell types along with the position of the
   * corresponding vertices in the second array. Note that an additional
   * element is set a the end to a dummy cell type pointing after the
   * last element in the vertices' array (see \ref dummyelement). This is precisely the kind of
   * data structure used in VTK for unstructured grids. To (try to) 
   * understand how this works, have a look at the function 
   * readUnstructuredGeometry() in FDataSetLoadAlgorithm.cc.
   */
  FCellDefinitions3DUnstructured( positive nb_pos, 
				  const std::string& name, 
				  std::vector< pair<FCell::CellType, unsigned int> >& types,
				  vector<FIndex> & vertices );

  ~FCellDefinitions3DUnstructured();

  /** 
   *\par Description:
   *Gets the class name.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  const FString& getClassName() const;

  /**
   * Returns the cell's vertices' indices.
   * \subsection dummyelement 
   * cellIndex + 1 is always valid thanks an additional last element 
   * in cell_types
   */
  void getCellVerticesIndices( const FIndex& cellId, 
			       std::vector<FIndex>& vertices ) const;

  void getCellType( const FIndex& cellId, FCell::CellType& cell_type ) const;

  shared_ptr<FCell> getCellTorso( const FIndex& cellId ) const;

  unsigned int getNbHexaCells() const;
  unsigned int getNbTetraCells() const;
  
  positive memSize() const;

  virtual const FNeighborhoodData* getNeighborhoodData() const;

private:

  void buildNeighborhoodData() const;

  std::vector< pair<FCell::CellType, unsigned int> > cell_types;
  std::vector<FIndex> cell_vertices; 

  unsigned int nbHexaCells, nbTetraCells;
};

#endif
