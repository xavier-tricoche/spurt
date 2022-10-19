//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSeparatrix.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:30 $
// Author:    $Author: garth $
// Version:   $Revision: 1.7 $
//
//--------------------------------------------------------------------------- 

#include "FAMSeparatrix.hh"

using namespace std;

//---------------------------------------------------------------------------

FAMSeparatrix::FAMSeparatrix() 
  : saddle(), sinksource(), last()
{
  length = 0.;
}

//--------------------------------------------------------------------------- 

FAMSeparatrix::FAMSeparatrix(const FAMSeparatrix& sep)
  : FObject(), saddle(sep.saddle), sinksource(sep.sinksource), last(sep.last)
{
  length = sep.length;
}

//--------------------------------------------------------------------------- 

FAMSeparatrix::FAMSeparatrix(const FIndex& saddleId,
			     const FIndex& sinksourceId,
			     double length) 
  : saddle(saddleId), sinksource(sinksourceId), last()
{
  this->length = length;
}

//--------------------------------------------------------------------------- 

FAMSeparatrix::~FAMSeparatrix()
{
}

//--------------------------------------------------------------------------- 

double FAMSeparatrix::getLength() const
{
  return length;
}

//--------------------------------------------------------------------------- 

void FAMSeparatrix::getSaddleIndex(FIndex& saddleId) const
{
  saddleId = saddle;
}

//--------------------------------------------------------------------------- 

void FAMSeparatrix::getSinkSourceIndex(FIndex& sinksourceId) const
{
  sinksourceId = sinksource;
}

//--------------------------------------------------------------------------- 

const FString& FAMSeparatrix::getClassName() const 
{
  static FString className = "FAMSeparatrix";

  return className;
}

//--------------------------------------------------------------------------- 

void FAMSeparatrix::setLastPoint(const FPosition& last) 
{
  this->last = last;
}

//--------------------------------------------------------------------------- 

void FAMSeparatrix::getLastPoint(FPosition& last) const
{
  last = this->last;
}

//---------------------------------------------------------------------------

ostream& operator<< (ostream &os, const FAMSeparatrix& sep)
{
  os << "FAMSeparatrix:" << endl;
  os << "\t saddle:      " << sep.saddle << endl;
  os << "\t sink/source: " << sep.sinksource << endl;
  os << "\t length:      " << sep.length << endl;
  os << "\t last:        " << sep.last << endl;

  return os;
}

//--------------------------------------------------------------------------- 

