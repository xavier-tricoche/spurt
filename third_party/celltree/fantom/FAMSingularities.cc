//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularities.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:42 $
// Author:    $Author: garth $
// Version:   $Revision: 1.14 $
//
//--------------------------------------------------------------------------- 

#include "FAMSingularities.hh"
#include "FException.hh"

using namespace std;


FAMSingularities::FAMSingularities()
  : FAMElement ()
{
  empty = true;
}

//---------------------------------------------------------------------------

FAMSingularities::~FAMSingularities()
{
}

//---------------------------------------------------------------------------

const FString& FAMSingularities::getClassName() const
{
  static FString name("FAMSingularities");

  return name;
}

//---------------------------------------------------------------------------

void FAMSingularities::getSingularities(std::vector<FAMSingularPoint>& outSing) const
{
  outSing.clear();
  outSing = points;
}

//---------------------------------------------------------------------------

void FAMSingularities::setSingularities(const std::vector<FAMSingularPoint>& inSing)
{
  points.clear();
  points = inSing;
  if (inSing.size())
    empty = false;
  else
    empty = true;

  change (true);
}

//---------------------------------------------------------------------------

void FAMSingularities::addSingularity(const FAMSingularPoint& inSing)
{
  FIndex tmpIndex;
  inSing.getIncludingCellIndex(tmpIndex);
  cellToSingularitiesIndex[tmpIndex.getIndex()]=points.size();

  points.push_back(inSing);
  empty = false;
  
  change (true);
}

//---------------------------------------------------------------------------

void FAMSingularities::addSingularities(std::vector<FAMSingularPoint>& inSing)
{
  points = inSing;

  FIndex tmpIndex;

  positive i=0;
  for (i=0; i <= inSing.size(); i++)
    {
      inSing[i].getIncludingCellIndex(tmpIndex);
      cellToSingularitiesIndex[tmpIndex.getIndex()]=i;
    }
  
  if (inSing.size())
    empty = false;

  change (true);
}

//---------------------------------------------------------------------------

void FAMSingularities::save(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//---------------------------------------------------------------------------

void FAMSingularities::load(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//---------------------------------------------------------------------------

void FAMSingularities::getSingularityForCell(FIndex cellId,
                                             FAMSingularPoint& sing) const {

  map<positive, positive>:: const_iterator tmpIter;
  tmpIter = cellToSingularitiesIndex.find(cellId.getIndex());
  if (tmpIter==cellToSingularitiesIndex.end()) {
    sing = FAMSingularPoint(FAMSingularPoint::UNKNOWN);
  }
  else
    sing = points[tmpIter->second];
}

//---------------------------------------------------------------------------
