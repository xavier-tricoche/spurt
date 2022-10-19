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

#ifndef __FCellDefinitions3DPoints_hh
#define __FCellDefinitions3DPoints_hh

#include "FCellDefinitions2Din3DUnstructured.hh"
#include "FCell.hh"
/**
 * This is a celldefinitions class containing cells of only triangular type
 *
 * \ingroup TensorField
 */
class FCellDefinitions3DPoints : public FCellDefinitions
{
public:

  /**
   * As opposed to the constructor of the general 2Din3DUnstructured cell 
   * definitions, this constructor only expects triangles and does not require
   * explicit definition of cell types.
   */
  FCellDefinitions3DPoints( positive nb_pos, 
				   const std::string& name);
  
  virtual ~FCellDefinitions3DPoints();

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

};

#endif

