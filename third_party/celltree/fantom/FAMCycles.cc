//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMCycles.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:41 $
// Author:    $Author: garth $
// Version:   $Revision: 1.10 $
//
//--------------------------------------------------------------------------- 

#include "FAMCycles.hh"
#include "FException.hh"

using namespace std;


FAMCycles::FAMCycles()
  : FAMElement ()
{
  empty = true;
}

FAMCycles::~FAMCycles()
{
}

const FString& FAMCycles::getClassName() const
{
  static FString name("FAMCycles");

  return name;
}

void FAMCycles::getCycles(std::vector<FAMCycle>& outCyc) const
{
  outCyc.clear();
  outCyc = points;
}

void FAMCycles::setCycles(const std::vector<FAMCycle>& inCyc)
{
  points.clear();
  points = inCyc;
  if (inCyc.size())
    empty = false;

  change (true);
}

void FAMCycles::addCycle(const FAMCycle& inCyc)
{
  points.push_back(inCyc);
  empty = false;

  change (true);
}

void FAMCycles::addCycles(const std::vector<FAMCycle>& inCyc)
{
  for (vector<FAMCycle>::const_iterator iter = inCyc.begin();
       iter != inCyc.end();
       iter++)
    points.push_back(*iter);

  if (inCyc.size())
    empty = false;

  change (true);
}

//---------------------------------------------------------------------------

void FAMCycles::save(const FString& /*fileNameBase*/) {
  THROW_EXCEPTION( FNotImplementedException, "FAMCycles::save not implemented YET... ask Tom :)");
}

//---------------------------------------------------------------------------

void FAMCycles::load(const FString& /*fileNameBase*/) {
  THROW_EXCEPTION( FNotImplementedException, "FAMCycles::load not implemented YET... ask Tom :)");
}

//---------------------------------------------------------------------------
