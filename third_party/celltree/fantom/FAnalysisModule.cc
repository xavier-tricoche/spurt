//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAnalysisModule.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:43 $
// Author:    $Author: garth $
// Version:   $Revision: 1.11 $
//
//---------------------------------------------------------------------------

#include "FAnalysisModule.hh"

#ifdef OUTLINE
#define inline
#include "FAnalysisModule.icc"
#endif

using namespace std;

//---------------------------------------------------------------------------

FAnalysisModule::~FAnalysisModule()
{
  map<FIndex, map<FIndex, map <FIndex, FAnalysisModuleData*> > >::iterator 
    lIter = theData.begin();
  while (lIter != theData.end()){
    map<FIndex, map <FIndex, FAnalysisModuleData*> >::iterator 
      tStepIter = lIter->second.begin();
    while (tStepIter != lIter->second.end()){
      map <FIndex, FAnalysisModuleData*>::iterator 
	tSetIter= tStepIter->second.begin();
      while (tSetIter != tStepIter->second.end()){
        delete tSetIter->second;
        tSetIter++;
      }
      tStepIter++;
    }
    
    lIter++;
  }  
} 

const FString& FAnalysisModule::getClassName() const
{
  static FString className("AnalysisModule");

  return className;
}


FAnalysisModuleData* FAnalysisModule::getModuleData(const FIndex& inLevel,
                                                    const FIndex& inTimestep,
                                                    const FIndex& inTensorSet) const
{
  if (!find(inLevel, inTimestep, inTensorSet)){
    theData[inLevel][inTimestep][inTensorSet] = new FAnalysisModuleData;
    validTimestep = validLevel = validTensorSet = true;
  }

  return theData[inLevel][inTimestep][inTensorSet];
}


