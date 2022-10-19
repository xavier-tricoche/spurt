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

#ifndef __FCellDefinitions3DTriangulation_hh
#define __FCellDefinitions3DTriangulation_hh

#include "FCellDefinitions2Din3DUnstructured.hh"
#include "FCell.hh"
/**
 * This is a celldefinitions class containing cells of only triangular type
 * 
 * \ingroup TensorField
 */
class FCellDefinitions3DTriangulation : public FCellDefinitions2Din3DUnstructured
{
public:

  /**
   * As opposed to the constructor of the general 2Din3DUnstructured cell 
   * definitions, this constructor only expects triangles and does not require
   * explicit definition of cell types.
   */
  FCellDefinitions3DTriangulation( positive nb_pos, 
				   const std::string& name,
				   std::vector<FIndex>& vertices,
				   bool buildneighborhooddata=true);
  
  ~FCellDefinitions3DTriangulation();

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

  
  void getCellVerticesIndices( const FIndex& cellId, 
			       std::vector<FIndex>& vertices ) const;
  
  void getCellType( const FIndex& cellId, FCell::CellType& cell_type ) const;
  
  shared_ptr<FCell> getCellTorso( const FIndex& cellId ) const;

  positive memSize() const;

  // only use, when really needed
  const std::vector<FIndex>& getVertices() const { return cell_vertices; }
private:

  std::vector<FIndex> cell_vertices; 
};

#endif

