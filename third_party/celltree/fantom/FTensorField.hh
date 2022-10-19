//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTensorField.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef __FTensorField_hh
#define __FTensorField_hh

//===========================================================================

class FIndex;

#include "FObject.hh"
#include "FHashTable.hh"

#include "stdAliases.hh"
#include "FTensorSet.hh"
#include "FGrid.hh"
#include "FFeature.hh"
#include "FAnalysisModuleData.hh"
#include "FAMSingularPoints.hh"

#include <string>
/**
 * This struct provides shortcuts to most of the important information
 * of a tenso field.
 *
 * \ingroup TensorField
 */
struct FTensorFieldInfo
{
    positive nbPos;
    positive posDim;
    positive nbCells;
    positive nbTensors;
    positive tensorOrder;
    positive tensorDim;
    positive tensorComp;

    bool isCellBased;

    shared_ptr<FCellDefinitions> cellDef;
    shared_ptr<FGrid>            grid;
    shared_ptr<FPositionSet>     posSet;
    shared_ptr<FTensorSet>       tensorSet;
    
    FAnalysisModuleData *analysis;
};


/**
 * The FTensorField class represents a discrete tensor field defined over
 * a bounded domain. 
 * In dataSet terms, it is the combination of a grid (deriving from FGrid) 
 * with a tensor set (an instance of FTensorSet). It provides functions that 
 * rely on the availability of both grid and tensor set (such as 
 * interpolation). Grid and tensor set can be accessed via shared pointers.
 * For quick access to all constituting components of a field, an 
 * FTensorFieldInfo can be obtained via getConvenienceInfo()
 */
class FTensorField : public FObject
{
public:

  /**
   *\par Description:
   *value-constructor: generates a TensorField with given values.
   *\exception
   *FNullPointerAssignmentException: newpTensorSet || newpGrid are NULL
   */
  FTensorField( const std::string& newname, 
	        shared_ptr<FTensorSet> newpTensorSet,
		shared_ptr<FGrid> newpGrid,
		bool isCellBased=0);

  /** 
   *\par Description:
   *Destructor.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  ~FTensorField();

  /** Returns the class name as string
   * \return class name
   */
  const FString& getClassName() const;

  /** 
   *\par Description:
   *Returns name which identifies this tensor field
   */
  const std::string& getName() const;

  /** 
   *\par Description:
   *Sets name which identifies this tensor field
   */
  void setName(const std::string& newname);

  /** 
   *\par Description:
   *Returns pointer to FTensorSet which defines the tensor values
   *according to the defined FGrid
   */
  shared_ptr<FTensorSet> getTensorSet() const;

  /** 
   *\par Description:
   *Returns underlying FGrid identifier
   *\pre
   *FTensorField is not empty
   *\post
   *none
   */
  shared_ptr<FGrid> getGrid() const;

  /** 
   *\par Description:
   *returns cell (along with its tensor information) with given index
   *\pre
   *cellID < numberOfCells
   *\post
   *none
   */
  shared_ptr<FCell> getCell( const FIndex& cellID ) const;

  /** 
   *\par Description:
   *returns cell (along with its tensor information) containing given position
   *\pre
   *position lies in the grid
   *\post
   *none
   */
  shared_ptr<FCell> getCell( const FPosition& position) const;

  void getConvenienceInfo( FTensorFieldInfo& info ) const;

  /** 
   *\par Description:
   *Returns contained FAnalysisModuleData
   */
  FAnalysisModuleData* getAnalysisModuleData() const;
  
  // ----------------- INTERPOLATION & DERIVATIVES ---------------------------
  
  /**
   *\par Description:
   *interpolates grid at position and returns result tensor
   *\exception
   *FException: error during interpolation
   */
  virtual bool interpolate( FTensor& outResult,
			    const FPosition& position ) const;
  
  /**
   *\par Description:
   *calculates derivatives at position and returns result tensor
   *\exception
   *FException: error during derivatives calculation
   */
  virtual bool derivatives( FTensor& outResult,
			    const FPosition& position ) const;


  /**
   *\par Description:
   *computes singularities lying in every cell and completes 
   *AnalysisModuleData
   */
  void computeSingularities() const;

  /*
   * \return
   * approximate size in bytes
   */
  virtual positive memSize() const;
  
  // feature list
  std::list< shared_ptr<const FFeature> > features;

  /**
   * \returns true if field is cellBased
   */
  bool isCellBased() const; 
  
private:


  // cell cache
  mutable shared_ptr<FCell> lastCell;
  mutable FIndex lastCellID;
  
  // pointers to FTensorSet ans FGrid defining the FTensorField
  shared_ptr<FTensorSet> pTensorSet;
  shared_ptr<FGrid> pGrid;

  // analysis module data (we don't need the interface of FAnalysisModule
  // here since we are concerned with a particular grid and tensor set).
  mutable FAnalysisModuleData analysisModuleData;
  std::string name;

  // hash table used for caching
  mutable FHashTable hashtable;
  bool cellBased;
};

#endif // __FTensorField_hh
/*! 
  @ingroup DataSet
  \defgroup TensorField TensorField (The Tensor Field and its components)

  \brief
  This submodule contains the FTensorField class with all its components: cell definition types, tensor sets types, grid types and position set types.

  \par WARNING
  Be aware that this module consists of many more classes that are not listed here because noone has added them yet.
*/
