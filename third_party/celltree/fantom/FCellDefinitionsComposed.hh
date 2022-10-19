//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions3DUnstructured.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:03 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef F_CELL_DEFINITIONS_COMPOSED
#define F_CELL_DEFINITIONS_COMPOSED



#include "FCellDefinitions.hh"
#include "FCell.hh"

/**
 * FCellDefinitionsComposed
 *
 * \ingroup TensorField
 */

class FCellDefinitionsComposed : public FCellDefinitions
{
public:

  /**
   */
  FCellDefinitionsComposed( positive nb_pos, 
			    const std::string& name, 
			    const vector< shared_ptr< FCellDefinitions > >& lcelldefs
			    );

  ~FCellDefinitionsComposed();

  void getCellVerticesIndices( const FIndex& cellId, 
			       std::vector<FIndex>& vertices ) const;

  void getCellType( const FIndex& cellId, FCell::CellType& cell_type ) const;

  shared_ptr<FCell> getCellTorso( const FIndex& cellId ) const;
  
  positive memSize() const;

  void getDistribution(vector<positive>& sizes,
		       vector<string> & names) const ;
 
  void setNbData(FNeighborhoodData*nbdat);


private:
  
  vector< shared_ptr< FCellDefinitions > > cellDefs;


  vector< positive > startIndices;
  inline positive getPartId(positive cellId) const;
    
};

#endif // F_CELL_DEFINITIONS_COMPOSED
