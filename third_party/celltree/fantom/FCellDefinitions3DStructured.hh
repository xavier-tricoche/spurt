//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions3DStructured.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:03 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellDefinitions3DStructured_hh
#define __FCellDefinitions3DStructured_hh

#include "FCellDefinitions.hh"

/**
 * This is a celldefinitions class containing 3D cells for a structured grid
 *
 * \ingroup TensorField
 */
class FCellDefinitions3DStructured : public FCellDefinitions
{
public:
  /*
   * This is the constructor for this celldefinition.
   * \param sizeX,sizeY,sizeZ:
   * number Of Cells in x,y or z direction
   */
  FCellDefinitions3DStructured( positive sizeX, positive sizeY,
				positive sizeZ, const std::string& newname,
				bool triang );

  ~FCellDefinitions3DStructured();

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

  shared_ptr<FCell> getCellTorso( const FIndex& cellId ) const;

  void getCellType( const FIndex& cellId, FCell::CellType& cell_type ) const;

  bool getTriangulated() const;

  positive memSize() const;

  positive getDimensionX() const;
  positive getDimensionY() const;
  positive getDimensionZ() const;

  /**
   * Lookup the index of a certain cell adressed by x,y,z offsets
   */
  FIndex getCellIndex(unsigned int x, unsigned int y, unsigned int z) const
  { 
    if(x>= getDimensionX() || y >= getDimensionY() || z >= getDimensionZ()) return FIndex(); // invalid cell
    return x+getDimensionX()*(y+getDimensionY()*z);
  }

protected:

  bool triangulated;

private:

  positive nbX, nbY, nbZ;
};

#endif
