//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions2DRectilinear.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:01 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellDefinitions2DRectilinear_hh
#define __FCellDefinitions2DRectilinear_hh

#include "FCellDefinitions2DStructured.hh"

/**
 * This is a celldefinitions class containing 2D cells for a rectilinear grid
 *
 * \ingroup TensorField
 */
class FCellDefinitions2DRectilinear : public FCellDefinitions2DStructured
{
public:

  /*
   *\param sizeX,sizeY:
   *number Of Cells in x or y direction
   */
  FCellDefinitions2DRectilinear( positive sizeX, positive sizeY,
				 const std::string& newname, bool triang );

  ~FCellDefinitions2DRectilinear();

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

  void getCellType( const FIndex& cellId, FCell::CellType& cell_type) const;

};

#endif
