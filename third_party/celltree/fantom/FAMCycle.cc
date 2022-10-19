//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMCycle.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:30 $
// Author:    $Author: garth $
// Version:   $Revision: 1.15 $
//
//--------------------------------------------------------------------------- 

#include "FAMCycle.hh"

using namespace std;

FAMCycle::FAMCycle()
{
  aprox = 0;
  cycleCells = 0;
}

//--------------------------------------------------------------------------- 

FAMCycle::FAMCycle(const FAMCycle& c) : FObject()
{
  sing = c.sing;
  aprox = new vector<FVector>;
  vector<FVector>::iterator it, itEnd;
  it = c.aprox->begin();
  itEnd = c.aprox->end();
  
  while(it != itEnd)
    {
      aprox->push_back (*it);
      it++;
    }

  cycleCells = new list<FIndex>;
  *cycleCells = *(c.cycleCells);
}

//--------------------------------------------------------------------------- 

FAMCycle::~FAMCycle()
{
  if (aprox)
    delete aprox;
  if (cycleCells)
    delete cycleCells;
}

void FAMCycle::getSingularity (const FAMSingularPoint *s) const
{
  s = &sing;
}

//--------------------------------------------------------------------------- 

void FAMCycle::setSingularity (const FAMSingularPoint &s)
{
  sing = s;
}

//--------------------------------------------------------------------------- 

void FAMCycle::getApproximation (std::vector<FVector>* &ap) const
{
  ap = aprox;
}

//--------------------------------------------------------------------------- 

void FAMCycle::setApproximation (std::vector<FVector> *ap)
{
  if (aprox)
    delete aprox;
  aprox = ap;
}

//--------------------------------------------------------------------------- 

void FAMCycle::getCycleCells (std::list<FIndex>* &cells) const
{
  cells = cycleCells;
}

//--------------------------------------------------------------------------- 

void FAMCycle::setCycleCells (std::list<FIndex> *cells)
{
  if (cycleCells)
    delete cycleCells;
  cycleCells = cells;
}

//--------------------------------------------------------------------------- 

const FString& FAMCycle::getClassName() const 
{
  static FString className = "FAMCycle";

  return className;
}

//--------------------------------------------------------------------------- 

FAMCycle& FAMCycle::operator= (const FAMCycle& c)
{
  sing = c.sing;
  aprox->clear ();
  vector<FVector>::iterator it, itEnd;
  it = c.aprox->begin();
  itEnd = c.aprox->end();
  
  while(it != itEnd)
    {
      aprox->push_back (*it);
      it++;
    }

  return *this;
}

//--------------------------------------------------------------------------- 

