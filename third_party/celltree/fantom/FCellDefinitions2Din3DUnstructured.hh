//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions2Din3DUnstructured.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:02 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellDefinitions2Din3DUnstructured_hh
#define __FCellDefinitions2Din3DUnstructured_hh

#include "FCellDefinitions.hh"
#include "FCell.hh"
/**
 * This is a celldefinitions class containing 2D cells with 3D vertices
 *
 * \ingroup TensorField
 */
class FCellDefinitions2Din3DUnstructured : public FCellDefinitions
{
public:

  /**
   * This constructor is quite tricky: Both provided vectors will be
   * taken over by the data structure. The second one must contain
   * the sequential list of all cell vertices (with redundancy). The 
   * first one provides all cell types along with the position of the
   * corresponding vertices in the second array. Note that an additional
   * element is set a the end to a dummy cell type pointing after the
   * last element in the vertices' array. This is precisely the kind of
   * data structure used in VTK for unstructured grids. To (try to) 
   * understand how this works, have a look at the function 
   * readUnstructuredGeometry() in FDataSetLoadAlgorithm.cc.
   */
  FCellDefinitions2Din3DUnstructured( positive nb_pos, 
				      const std::string& name,
				      std::vector< pair<FCell::CellType,unsigned int> > &types,
				      std::vector<FIndex>& vertices,
				      bool checkConsistency=false);
  
  virtual ~FCellDefinitions2Din3DUnstructured();

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
  
  virtual void getCellVerticesIndices( const FIndex& cellId, 
			       std::vector<FIndex>& vertices ) const;
  
  virtual void getCellType( const FIndex& cellId, FCell::CellType& cell_type ) const;
  
  virtual shared_ptr<FCell> getCellTorso( const FIndex& cellId ) const;

  virtual positive memSize() const;

private:

  std::vector< pair<FCell::CellType, unsigned int> > cell_types;
  std::vector<FIndex> cell_vertices; 

  unsigned int nbQuadrilateralCells, nbTriangleCells;
protected:
  /**
   * Just for inheritance
   **/
  FCellDefinitions2Din3DUnstructured( positive nb_pos, 
				      const std::string& name);
};

#endif

