//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMInOutFlowSeparatrix.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:41 $
// Author:    $Author: garth $
// Version:   $Revision: 1.7 $
//
//--------------------------------------------------------------------------- 

#include "FAMInOutFlowSeparatrix.hh"
#include "FException.hh"

using namespace std;

FAMInOutFlowSeparatrix::FAMInOutFlowSeparatrix()
    : FAMElement () ,separatrices(), sepToCellRelations()
{
  empty = true;
}

//---------------------------------------------------------------------------

FAMInOutFlowSeparatrix::~FAMInOutFlowSeparatrix()
{
}

//---------------------------------------------------------------------------

const FString& FAMInOutFlowSeparatrix::getClassName() const
{
  static FString name("FAMInOutFlowSeparatrix");

  return name;
}

//---------------------------------------------------------------------------

void FAMInOutFlowSeparatrix::getInOutFlowSeparatrices(std::vector< std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> >&outSeparatrices ) const
{
  outSeparatrices=separatrices;
}

//---------------------------------------------------------------------------

void FAMInOutFlowSeparatrix::getInOutFlowSeparatrices(std::vector< std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> > &outSeparatrices,
                                                      std::vector< std::vector<FIndex> > &outSepToCellRelations) const
{
  outSeparatrices=separatrices;
  outSepToCellRelations=sepToCellRelations;
}

//---------------------------------------------------------------------------

void FAMInOutFlowSeparatrix::setInOutFlowSeparatrices(const std::vector< std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> > &inSeparatrices)
{
  separatrices=inSeparatrices;
  empty=separatrices.size()==0;
}
//---------------------------------------------------------------------------

void FAMInOutFlowSeparatrix::setInOutFlowSeparatrices(const std::vector< std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> > &inSeparatrices,
                                                      const std::vector< std::vector<FIndex> > &inSepToCellRelations)
{
  separatrices=inSeparatrices;
  
  sepToCellRelations=inSepToCellRelations;
  empty=separatrices.size()==0;
}

//---------------------------------------------------------------------------

void FAMInOutFlowSeparatrix::addInOutFlowSeparatrix(const std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> &inSeparatrix)
{
  separatrices.push_back(inSeparatrix);
  empty=false;
}

//---------------------------------------------------------------------------

void FAMInOutFlowSeparatrix::addInOutFlowSeparatrix(const std::vector<FAMInOutFlowSeparatrix::SeparatrixPiece> &inSeparatrix,
                                                      const std::vector< FIndex > &inSepToCellRelation)
{
  separatrices.push_back(inSeparatrix);
  sepToCellRelations.push_back(inSepToCellRelation);
  empty=false;
}

//---------------------------------------------------------------------------

bool FAMInOutFlowSeparatrix::hasRelations( void ) const
{
  return ( sepToCellRelations.size() !=0 );
}

//---------------------------------------------------------------------------

void FAMInOutFlowSeparatrix::save(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//---------------------------------------------------------------------------

void FAMInOutFlowSeparatrix::load(const FString& /*fileNameBase*/) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//---------------------------------------------------------------------------
