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

#ifndef __FCellDefinitions2DPoints_hh
#define __FCellDefinitions2DPoints_hh

#include "FCellDefinitions.hh"
#include "FCell.hh"
/**
 * This is a celldefinitions class containing cells of only triangular type
 *
 * \ingroup TensorField
 */
class FCellDefinitions2DPoints : public FCellDefinitions
{
public:

  /**
   * As opposed to the constructor of the general 2Din2DUnstructured cell 
   * definitions, this constructor only expects triangles and does not require
   * explicit definition of cell types.
   */
  FCellDefinitions2DPoints( positive nb_pos, 
				   const std::string& name);
  
  ~FCellDefinitions2DPoints();

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

};

#endif

