//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions3DRectilinear.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:02 $
// Author:    $Author: garth $FAxisParallelQuadCell2D
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellDefinitions3DRectilinear_hh
#define __FCellDefinitions3DRectilinear_hh

#include "FCellDefinitions3DStructured.hh"

/**
 * This is a celldefinitions class containing cells for rectilinear grids
 *
 * \ingroup TensorField
 */
class FCellDefinitions3DRectilinear : public FCellDefinitions3DStructured
{
public:

  /*
   *\param sizeX,sizeY,sizeZ:
   *number Of Cells in x,y or z direction
   */
  FCellDefinitions3DRectilinear( positive sizeX, positive sizeY, 
				 positive sizeZ, const std::string& newname, 
				 bool triang );

  ~FCellDefinitions3DRectilinear();

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

  shared_ptr<FCell> getCellTorso( const FIndex& cellId ) const;

  void getCellType( const FIndex& cellId, FCell::CellType& cell_type ) const;

private:

};

#endif
