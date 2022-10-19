//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMBifurcations.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:40 $
// Author:    $Author: garth $
// Version:   $Revision: 1.9 $
//
//--------------------------------------------------------------------------- 

#include "FAMBifurcations.hh"
#include "FException.hh"

using namespace std;


FAMBifurcations::FAMBifurcations()
  : FAMElement ()
{
  empty = true;
}

FAMBifurcations::~FAMBifurcations()
{
}

const FString& FAMBifurcations::getClassName() const
{
  static FString name("FAMBifurcations");

  return name;
}

void FAMBifurcations::getBifurcations(std::vector<FAMBifurcation>& outBif) const
{
  outBif.clear();
  outBif = bifs;
}

void FAMBifurcations::setBifurcations(const std::vector<FAMBifurcation>& inBif)
{
  bifs.clear();
  bifs = inBif;
  if (inBif.size())
    empty = false;

  change (true);
}

void FAMBifurcations::addBifurcation(const FAMBifurcation& inBif)
{
  bifs.push_back(inBif);
  empty = false;

  change (true);
}

void FAMBifurcations::addBifurcations(std::vector<FAMBifurcation>& inBif)
{
  for (positive i=0 ; i<inBif.size() ; i++)
    bifs.push_back(inBif[i]);

  if (inBif.size())
    empty = false;

  change (true);
}

//---------------------------------------------------------------------------

void FAMBifurcations::save(const FString& /*fileNameBase*/) {
  THROW_EXCEPTION( FNotImplementedException, "FAMBifurcations::save not implemented YET... ask Tom :)");
}

//---------------------------------------------------------------------------

void FAMBifurcations::load(const FString& /*fileNameBase*/) {
  THROW_EXCEPTION( FNotImplementedException, "FAMBifurcations::load not implemented YET... ask Tom :)");
}

//---------------------------------------------------------------------------
