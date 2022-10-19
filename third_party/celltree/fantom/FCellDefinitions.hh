//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:01 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCellDefinitions_hh
#define __FCellDefinitions_hh

#include <vector>
#include <string>
#include <boost/shared_ptr.hpp>

#include "FIndex.hh"

#include "FNeighborhoodData.hh"
#include "FCell.hh"

using namespace boost;

/**
 * FCellDefinitions is the abstract base class for specializations
 * containing the connectivity of points to define cells on a set of
 * positions.
 *
 * \ingroup TensorField
 */
class FCellDefinitions : public FObject
{
public:

  FCellDefinitions( const std::string& newname );

  virtual ~FCellDefinitions();

  /**
   * \par Description: 
   * Get a cell's vertices indices.
   */
  virtual void getCellVerticesIndices( const FIndex& cellId, 
				       std::vector< FIndex >& vertices ) const = 0; 

  /**
   * \par Description:
   * Get a cell's torso, i.e. the FCell itself.
   */
  virtual shared_ptr<FCell> getCellTorso( const FIndex& cellId ) const = 0;

  /** 
   *\par Description:
   *Returns the type of the cell with cellId in cellType
   */
  virtual void getCellType( const FIndex& cellId, 
			    FCell::CellType& cellType) const = 0;

  /** 
   *\par Description:
   *Returns the celldefinition's number of cells
   */
  positive getNbCells() const;
  
  /** 
   *\par Description:
   *Returns celldefinition's number of positions
   */
  positive getNbPositions() const;

  /** 
   *\par Description:
   *Returns celldefinition's neighborhood information
   */
  virtual const FNeighborhoodData* getNeighborhoodData() const;

  /** 
   *\par Description:
   * Returns the celldefinition's size in memory.
   */
  virtual positive memSize() const = 0;

  
  /**
   *\par if cellDefinition is distributed,
   * get distribution of the cells
   */
  virtual void getDistribution(vector<positive>& sizes,
			       vector<string> & names) const;

  const std::string& getName() const;

protected:

  positive nbCells, nbPos;

  // neighborhood information
  mutable FNeighborhoodData *neighborData;

  std::string name;
};

#endif



