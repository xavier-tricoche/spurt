//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCellDefinitions.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:01 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#include "FCellDefinitions.hh"
#include "FException.hh"

#include "eassert.hh"

#include <iostream>

//--------------------------------------------------------------------------- 

FCellDefinitions::FCellDefinitions( const std::string& newname )
  : nbCells( 0 ), nbPos( 0 ),  neighborData(0), name( newname )
{
}

//--------------------------------------------------------------------------- 

FCellDefinitions::~FCellDefinitions()
{
    // derived classes should have allocated neighborData
//    eassert( neighborData );
#ifndef NODEBUG
  if(!neighborData)
    std::cerr << "WARNING in FCellDefinitions: derived class should have allocated an FNeighborhoodData" << std::endl;
#endif

  if(neighborData)
    delete neighborData;
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions::getNbCells() const
{
    return nbCells;
}

//--------------------------------------------------------------------------- 

positive FCellDefinitions::getNbPositions() const
{
    return nbPos;
}

//--------------------------------------------------------------------------- 

const FNeighborhoodData* FCellDefinitions::getNeighborhoodData() const
{
    return neighborData;
}

//--------------------------------------------------------------------------- 

void FCellDefinitions::
getDistribution(vector<positive>& /*sizes*/,
		vector<string> & /*names*/) const
{
  THROW_EXCEPTION(FNotImplementedException," this class is not distributed! ");
}

//--------------------------------------------------------------------------- 

const std::string& FCellDefinitions::getName() const
{
  return name;
}
