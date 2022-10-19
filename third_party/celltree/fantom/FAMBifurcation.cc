//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMBifurcation.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:40 $
// Author:    $Author: garth $
// Version:   $Revision: 1.11 $
//
//--------------------------------------------------------------------------- 

#include "FAMBifurcation.hh"
#include "FMath.hh"
#include "FException.hh"

#include "FVector.hh"
#include "FMatrix.hh"
#include "FTensor.hh"
#include "FCell.hh"

//---------------------------------------------------------------------------
#ifdef OUTLINE
#include "FAMBifurcation.icc"
#endif
//---------------------------------------------------------------------------
FAMBifurcation::FAMBifurcation() 
  : pos(), type(FAMBifurcation::NONE_BIF), singPathId1(), singPathId2(), cellId()
{
}

//--------------------------------------------------------------------------- 

FAMBifurcation::FAMBifurcation(const FAMBifurcation& bif) 
  : FObject(), pos(bif.pos), type(bif.type), singPathId1(bif.singPathId1), 
  singPathId2(bif.singPathId2), cellId(bif.cellId)
{
}

//--------------------------------------------------------------------------- 

FAMBifurcation::FAMBifurcation(const FPosition& _pos)
  : FObject(), pos(_pos), type(FAMBifurcation::NONE_BIF), singPathId1(), singPathId2(), cellId()
{
}

//--------------------------------------------------------------------------- 

FAMBifurcation::FAMBifurcation(const FPosition& _pos, 
			       const bifurcationType& _type,
			       const FIndex& _cellId)
  : FObject(), pos(_pos), type(_type), singPathId1(), singPathId2(), cellId(_cellId)
{
}

//--------------------------------------------------------------------------- 

FAMBifurcation::~FAMBifurcation()
{
}

//--------------------------------------------------------------------------- 

const FString& FAMBifurcation::getClassName() const 
{
  static FString className = "FAMBifurcation";

  return className;
}
      
//--------------------------------------------------------------------------- 

FAMBifurcation& FAMBifurcation::operator=(const FAMBifurcation& bif) {

  if (&bif != this) {
    this->pos = bif.pos;
    
    this->type = bif.type;

    this->singPathId1 = bif.singPathId1;

    this->singPathId2 = bif.singPathId2;

    this->cellId = bif.cellId;
  }

  return *this;
}

//---------------------------------------------------------------------------

ostream& operator<< (ostream &os, const FAMBifurcation::bifurcationType& type) {
  switch (type) {
  case FAMBifurcation::CREATION: 
    os << "CREATION";
    break;
  case FAMBifurcation::ANNIHILATION:
    os << "ANNIHILATION";
    break;
  case FAMBifurcation::HOPF:
    os << "HOPF";
    break;
  case FAMBifurcation::SWAP:
    os << "SWAP";
    break;
  case FAMBifurcation::WEDGE_SWAP:
    os << "FAMSingularPoint::WEDGE SWAP";
    break;
  case FAMBifurcation::NONE_BIF:
    os << "FAMSingularPoint::UNKNOWN/UNDEFINED TYPE";
    break;
  }
  return os;
}

//---------------------------------------------------------------------------

ostream& operator<< (ostream &os, const FAMBifurcation& bif)
{
  os << "FAMBifurcation:" << endl;
  os << "--position: " << bif.pos << endl;
  os << "--type: " << bif.type << endl;
  os << "--1st singular path involved: " << bif.singPathId1 << endl; 
  os << "--2nd singular path involved: " << bif.singPathId2 << endl; 
  os << "--cellId: " << bif.cellId << endl;

  return os;
}

