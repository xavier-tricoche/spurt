//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMGridBoundary.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:41 $
// Author:    $Author: garth $
// Version:   $Revision: 1.7 $
//
//--------------------------------------------------------------------------- 

#include "FAMGridBoundary.hh"
#include "FException.hh"

using namespace std;

FAMGridBoundary::FAMGridBoundary()
	: FAMElement ()
{
  empty = true;
}

//---------------------------------------------------------------------------

FAMGridBoundary::~FAMGridBoundary()
{
}

//---------------------------------------------------------------------------

const FString& FAMGridBoundary::getClassName() const
{
  static FString name("FAMGridBoundary");

  return name;
}

//---------------------------------------------------------------------------

FIndex FAMGridBoundary::getGridBoundary( void ) const
{
  return boundaryField;
}

//---------------------------------------------------------------------------

void FAMGridBoundary::setGridBoundary(const FIndex& inBoundaryField)
{
  boundaryField = inBoundaryField;
  empty = (! boundaryField.isValid());
}

//---------------------------------------------------------------------------

void FAMGridBoundary::save(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//---------------------------------------------------------------------------

void FAMGridBoundary::load(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//---------------------------------------------------------------------------

void FAMGridBoundary::getGridBoundaryCellRelations( std::vector< std::pair< FIndex, FIndex > > &relations ) const
{
  relations=this->relations;
}

//---------------------------------------------------------------------------

void FAMGridBoundary::setGridBoundaryCellRelations( const std::vector< std::pair< FIndex, FIndex > > &relations )
{
  this->relations=relations;
}

//---------------------------------------------------------------------------

bool FAMGridBoundary::hasRelations() const
{
  return (relations.size()==0);
}

//---------------------------------------------------------------------------
