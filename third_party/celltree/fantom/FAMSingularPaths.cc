//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularPaths.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:42 $
// Author:    $Author: garth $
// Version:   $Revision: 1.6 $
//
//--------------------------------------------------------------------------- 

#include "FAMSingularPaths.hh"
#include "FException.hh"

using namespace std;


FAMSingularPaths::FAMSingularPaths()
{
  empty = true;
}

FAMSingularPaths::~FAMSingularPaths()
{
}

const FString& FAMSingularPaths::getClassName() const
{
  static FString name("FAMSingularPaths");

  return name;
}

void FAMSingularPaths::getSingularPaths(std::vector<FAMSingularPath>& outPath) const
{
  outPath.clear();
  outPath = paths;
}

void FAMSingularPaths::setSingularPaths(const std::vector<FAMSingularPath>& inPath)
{
  paths.clear();
  paths = inPath;
  if (inPath.size())
    empty = false;
}

void FAMSingularPaths::addSingularPath(const FAMSingularPath& inPath)
{
  paths.push_back(inPath);
  empty = false;
}

void FAMSingularPaths::addSingularPaths(const std::vector<FAMSingularPath>& inPath)
{
  for (positive i=0 ; i<inPath.size() ; i++)
    paths.push_back(inPath[i]);

  if (inPath.size())
    empty = false;
}

//---------------------------------------------------------------------------

void FAMSingularPaths::save(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException ); 
}

//---------------------------------------------------------------------------

void FAMSingularPaths::load(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException ); 
}

//---------------------------------------------------------------------------
