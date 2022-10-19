//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularEdges.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:42 $
// Author:    $Author: garth $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#include "FAMSingularEdges.hh"
#include "FException.hh"

using namespace std;


FAMSingularEdges::FAMSingularEdges()
  : FAMElement ()
{
  empty = true;
}

//---------------------------------------------------------------------------

FAMSingularEdges::~FAMSingularEdges()
{
}

//---------------------------------------------------------------------------

const FString& FAMSingularEdges::getClassName() const
{
  static FString name("FAMSingularEdges");

  return name;
}

//---------------------------------------------------------------------------

void 
FAMSingularEdges::getSingularEdges(std::vector<FAMSingularEdge>& outEdges) const
{
  outEdges.clear();
  outEdges = edges;
}

//---------------------------------------------------------------------------

void 
FAMSingularEdges::setSingularEdges(const std::vector<FAMSingularEdge>& inEdges)
{
  edges.clear();
  edges = inEdges;
  if (inEdges.size())
    empty = false;
  else
    empty = true;

  change (true);
}

//---------------------------------------------------------------------------

void FAMSingularEdges::save(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//---------------------------------------------------------------------------

void FAMSingularEdges::load(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//---------------------------------------------------------------------------
