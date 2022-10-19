//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSeparatrices.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:41 $
// Author:    $Author: garth $
// Version:   $Revision: 1.4 $
//
//--------------------------------------------------------------------------- 

#include "FAMSeparatrices.hh"
#include "FException.hh"

using namespace std;


FAMSeparatrices::FAMSeparatrices()
  : FAMElement ()
{
  empty = true;
}

//---------------------------------------------------------------------------

FAMSeparatrices::~FAMSeparatrices()
{
}

//---------------------------------------------------------------------------

const FString& FAMSeparatrices::getClassName() const
{
  static FString name("FAMSeparatrices");

  return name;
}

//---------------------------------------------------------------------------

void FAMSeparatrices::getSeparatrices(std::vector<FAMSeparatrix>& outSep) const
{
  outSep.clear();
  outSep = seps;
}

//---------------------------------------------------------------------------

void FAMSeparatrices::setSeparatrices(const std::vector<FAMSeparatrix>& inSep)
{
  seps.clear();
  seps = inSep;
  if (inSep.size())
    empty = false;

  change (true);
}

//---------------------------------------------------------------------------

void FAMSeparatrices::addSeparatrix(const FAMSeparatrix& inSep)
{
  seps.push_back(inSep);
  empty = false;
  
  change (true);
}

//---------------------------------------------------------------------------

void FAMSeparatrices::addSeparatrices(std::vector<FAMSeparatrix>& inSep)
{
  seps = inSep;
  
  if (inSep.size())
    empty = false;

  change (true);
}

//---------------------------------------------------------------------------

void FAMSeparatrices::save(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//---------------------------------------------------------------------------

void FAMSeparatrices::load(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//---------------------------------------------------------------------------
