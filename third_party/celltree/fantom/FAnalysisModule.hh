//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAnalysisModule.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:44 $
// Author:    $Author: garth $
// Version:   $Revision: 1.13 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAnalysisModule_hh
#define __FAnalysisModule_hh

#include "FObject.hh"
#include "FString.hh"
#include "FIndex.hh"
#include "FException.hh"

#include "FAnalysisModuleData.hh"

#include <map>

class FAMLevels;

/**
 * This is the link betwen a given FLevel, FTimestep, FTensorSet and
 * the Analysis module data for that particular configuration
 */
   

class FAnalysisModule : public FObject
{
  //=== Constructors ========================================================
public:

  /**
   * \brief
   *generates an empty Module
   * \pre
   *none
   * \post
   * none
   * \exception
   * none
   * \param 
   * none
   * \return 
   * none
   */
  FAnalysisModule();
  
  /**
   * \brief
   * Initializes a new Module with current pointers
   * \pre
   * none
   * \post
   * getModule() can be used
   * \exception
   * FObjectNotEmptyException: ERROR: FAnalysisModuleData already existed
   * \param 
   *inLevel iterator of the current level
   * \param 
   *inTimestep iterator of the current timestep
   *\param
   *inTensorSet iterator of the current tensorset
   * \return 
   * none
   */
  FAnalysisModule(const FIndex& inLevel, const FIndex& inTimestep,
                  const FIndex& inTensorSet);


  //=== Destructor =========================================================

public:
  /**
   * \brief
   *Destructor
   * \pre
   *none
   * \post
   *all FAnalyseModules belonging to this module are deleted
   */
  ~FAnalysisModule();


  //=== Member Functions ====================================================

public:
    const FString& getClassName() const;

  /**
   * \brief
   *returns a pointer to the current Module data
   * \pre
   *level, timestep, tensorset have to be set first
   * \post
   * none
   * \exception
   * FEmptyObjectException: there is no current Module
   * \param 
   *none
   */
  FAnalysisModuleData* getModuleData() const;

  /**
   * \brief
   *returns a pointer to the module data of the given level, timestep, tensorset combination
   * \pre
   *none
   * \post
   *a new FAnalysisModuleData object might be allocated
   * \exception
   *none 
   * \param 
   *inLevel index of the current level
   * \param 
   *inTimestep index of the current timestep
   *\param
   *inTensorSet index of the current tensorset
   */
  FAnalysisModuleData* getModuleData(const FIndex& inLevel, 
				     const FIndex& inTimestep,
                                     const FIndex& inTensorSet) const;

  /**
   * \brief
   *sets current level
   * \pre
   *none
   * \post
   * none
   * \exception
   * none
   * \param 
   *inLevel the new current level index
   */
  void setLevel(const FIndex& inLevel);

  /**
   * \brief
   *sets current Timestep
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \param 
   *inTStep the new current timeStep index
   */
  void setTimestep(const FIndex& inTStep);

  /**
   * \brief
   *sets current TensorSet
   * \pre
   *none
   * \post
   * none
   * \exception
   * none
   * \param 
   *inTSet the new current TensorSet index
   */
  void setTensorSet(const FIndex& inTSet);
  
private:
  /**
   * \brief
   * checks if a FAnalysisModuleData object for the given combination of
   * level, timestep, tensorset exists.
   * \pre
   *none
   * \post
   *none 
   * \exception
   * none
   * \param 
   *inLevel index of the current level
   * \param 
   *inTimestep index of the current timestep
   *\param
   *inTensorSet index of the current tensorset
   */
  bool find(const FIndex& inLevel, const FIndex& inTimestep,
            const FIndex& inTensorSet) const;

  
  mutable std::map<FIndex, std::map<FIndex, std::map <FIndex, FAnalysisModuleData*> > > theData;

  FIndex currentLevel;
  FIndex currentTimestep;
  FIndex currentTensorSet;

  mutable bool validTimestep, validLevel, validTensorSet;
};

//===========================================================================
#ifndef OUTLINE
#include "FAnalysisModule.icc"
#endif
//=========================================================================== 

#endif // __FAnalysisModule_hh
