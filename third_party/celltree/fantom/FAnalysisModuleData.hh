//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAnalysisModuleData.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:43 $
// Author:    $Author: garth $
// Version:   $Revision: 1.33 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAnalysisModuleData_hh
#define __FAnalysisModuleData_hh

#include <vector>
#include "FString.hh"
#include "FAMSingularities.hh"
#include "FAMSingularEdges.hh"
#include "FAMBifurcations.hh"
#include "FAMSingularPaths.hh"
#include "FAMSeparatrices.hh"
#include "FAMGridBoundary.hh"
#include "FAMInOutFlowSeparatrix.hh"
#include "FIndex.hh"
#include "FAMElement.hh"


#include <boost/shared_ptr.hpp>
#include <list>
#include "stdAliases.hh"

/**
 * The analysis module data class is the interface to the various features
 * offered by the analysis module. Adding and retrieving information on a
 * certain tensor field is done through the methods defined here.
 *
 * The meaning of the "changed" states is not fully cleared and must be
 * discussed again after the redesign of the data set component.
 */

class FAnalysisModuleData : public FObject
{
  //=== Constructors ========================================================
public:

  /**
   * \brief
   *generates an empty Module
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \param 
   * none
   * \return 
   * none
   */
  FAnalysisModuleData();
  

  //=== Destructor =========================================================

  /**
   * \brief
   *Destructor
   * \pre
   *none
   * \post
   *all FAnalyseModules belonging to this module are deleted
   */
  ~FAnalysisModuleData();


  //=== Member Functions ====================================================

public:
  /// Undocumented
  const FString& getClassName() const;
  void addAMObject (boost::shared_ptr<FAMElement> add_object);
  void getAMObject  (boost::shared_ptr<FAMElement> get_object) const;
  boost::shared_ptr<FAMElement> getAMObject  (std::string objName) const;
  const std::list< boost::shared_ptr<FAMElement> >& getAllAMObjects () const;
  std::vector< std::string > getAllClassNames() const;
  void deleteAMObject (std::string className);
  void deleteAllAMObjects ();
  void print () const;
  positive size() const;


private:
  std::list< boost::shared_ptr<FAMElement> > storedObjects;
};

#endif // __FAnalysisModuleData_hh
