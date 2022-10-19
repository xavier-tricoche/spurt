//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellLocator.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:03 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#include "FIndex.hh"
#include "FException.hh"
#include "FCellLocator.hh"
#include "eassert.hh"

#include <iostream>

using namespace std;

FCellLocator::~FCellLocator()
{
}

// --------------------------------------------------------------------------

FCellLocator::FCellLocator( const FPositionSet* positionSet )
    : nbCells(0), pSet(positionSet)
{
}

//---------------------------------------------------------------------------
/*
void FCellLocator::save (ostream& out) const
{
  throw FException("Sorry, not implemented yet");
}
*/
//---------------------------------------------------------------------------
/*
void FCellLocator::load (istream& in)
{
  throw FException("Sorry, not implemented yet");
}
*/
//---------------------------------------------------------------------------

void FCellLocator::info(std::ostream& tmp) const
{
  tmp << "::CELLLOCATOR::" << endl
      << "============" << endl
      << endl;

  tmp << "Number of Cells: " << this->nbCells << endl
      << "============" << endl;
}

//---------------------------------------------------------------------------

positive FCellLocator::memSize() const
{
  return (positive)-1;
  FNotImplementedException e;
  e.addTraceMessage("positive FCellLocator::memSize() const");
  throw e;
}

//--------------------------------------------------------------------------- 

