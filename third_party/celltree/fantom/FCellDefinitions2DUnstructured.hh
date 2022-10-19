//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions2DUnstructured.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:02 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellDefinitions2DUnstructured_hh
#define __FCellDefinitions2DUnstructured_hh

#include "FCellDefinitions.hh"
#include "FCell.hh"

/**
 * This is a celldefinitions class containing arbitrary 2D cells
 *
 * \ingroup TensorField
 */
class FCellDefinitions2DUnstructured : public FCellDefinitions
{
public:

  /**
   * This constructor is quite tricky: Both provided vectors will be
   * taken over by the data structure. The second one must contains 
   * the sequential list of all cell vertices (with redundancy). The 
   * first one provides all cell types along with the position of the
   * corresponding vertices in the second array. Note that an additional
   * element is set in the end to a dummy cell type pointing after the
   * last element in the vertices' array. This is precisely the kind of
   * data structure used in VTK for unstructured grids. To (try to) 
   * understand how this works, have a look at the function 
   * readUnstructuredGeometry() in FDataSetLoadAlgorithm.cc.
   */
  FCellDefinitions2DUnstructured( positive nb_pos, const std::string& name,
				  vector< pair<FCell::CellType,unsigned int> >& types,
				  vector< FIndex > & vertices );

  ~FCellDefinitions2DUnstructured();

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

#ifdef OLD_CODE
  void addNewCell( FCell::CellType myCellType, vector< FIndex > vertices );
#endif
  
  void completeCellDefinitions();

  void getCellVerticesIndices( const FIndex& cellId, 
			       vector< FIndex >& vertices ) const;

  void getCellType( const FIndex& cellId, FCell::CellType& cell_type ) const;
  
  shared_ptr<FCell> getCellTorso( const FIndex& cellId ) const;

  positive memSize() const;
private:

  vector< pair< FCell::CellType, unsigned int > > cell_types;
  vector< FIndex > cell_vertices; 
  
 unsigned int nbQuadrilateralCells, nbTriangleCells;
protected:
 /**
   * Just for inheritance
   **/
  FCellDefinitions2DUnstructured( positive nb_pos, 
				      const std::string& name);
};

#endif
