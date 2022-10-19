//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions2DStructured.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:01 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellDefinitions2DStructured_hh
#define __FCellDefinitions2DStructured_hh

#include "FCellDefinitions.hh"

/**
 * This is a celldefinitions class containing 2D cells for a structured grid
 *
 * \ingroup TensorField
 */
class FCellDefinitions2DStructured : public FCellDefinitions
{
public:

  /*
   *\param sizeX,sizeY:
   *number Of Cells in x or y direction
   */
  FCellDefinitions2DStructured( positive sizeX, positive sizeY,
				const std::string& newname, bool triang );

  ~FCellDefinitions2DStructured();

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

  void getCellType( const FIndex& cellId, FCell::CellType& cell_type) const;

  bool getTriangulated() const;

  positive memSize() const;

  positive getDimensionX() const;
  positive getDimensionY() const;

  /**
   * Lookup the index of a certain cell adressed by x,y offsets
   */
  FIndex getCellIndex(unsigned int x, unsigned int y)
  { 
    if(x>= getDimensionX() || y >= getDimensionY()) return FIndex(); // invalid cell
    return x+getDimensionX()*y;
  }

protected:

  bool triangulated;

private:

  positive nbX, nbY;
};

#endif
