//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCell1Din3D.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:01 $
// Author:    $Author: garth $
// Version:   $Revision: 1.4 $
//
//--------------------------------------------------------------------------- 

#include "stdAliases.hh"
#include "FCell1Din3D.hh"

#include <iostream>

//------------------------------------------------------------------------------
#if 0
class FCellEdgeInfo{

  FCellEdgeInfo(const FCell*cell);
  
  FIndex edgeToGo(double&length,
		  const FArray& start, 
		  const FArray& dir) const;

  typedef double double2[2];

  const double2*const pData;
  const int nbEdges;
  double2 edgeNormals[FCell::maxNoEdges];
  double  edgeNormalsSqlen[FCell::maxNoEdges];

  ///determines if rot.dir. of polygon is 
  ///not counter-clockwise
  bool inverted;

  friend class FCell1Din3D;

};

//------------------------------------------------------------------------------

FCellEdgeInfo::FCellEdgeInfo(const FCell*cell)
  :pData((const double2*)cell->getPositionData()),nbEdges(cell->getNumEdges())
{

  int i;

  for(i=0;i<nbEdges-1;i++)
    {
      edgeNormals[i][0] = pData[i+1][1] - pData[i][1];
      edgeNormals[i][1] = pData[i][0] - pData[i+1][0];      
    }
  edgeNormals[i][0] = pData[0][1] - pData[i][1];
  edgeNormals[i][1] = pData[i][0] - pData[0][0];


  for(i=0;i<nbEdges;i++)
    edgeNormalsSqlen[i]
      = edgeNormals[i][0] * edgeNormals[i][0]
      + edgeNormals[i][1] * edgeNormals[i][1];

  //if rot.dir. is not counter-clockwise
  double d 
    = edgeNormals[0][0] * edgeNormals[1][1]    
    - edgeNormals[1][0] * edgeNormals[0][1];

  double e
    = edgeNormals[1][0] * edgeNormals[2][1]    
    - edgeNormals[2][0] * edgeNormals[1][1];

  if(fabs(e)>fabs(d))
    inverted = e<0;
  else
    inverted = d<0;


  
}
  
//------------------------------------------------------------------------------


FIndex FCellEdgeInfo::edgeToGo(double&length,
			       const FArray& start, 
			       const FArray& dir) const
{
  
  double2 d={dir[0],dir[1]};
  double2 s={start[0],start[1]};

  length = HUGE_VAL;

  double dsqlen = d[0]*d[0]+d[1]*d[1];
  positive reti = FIndex::invalid;

  for(int i=0;i<nbEdges;i++){
    double dlen 
      = d[0]*edgeNormals[i][0]
      + d[1]*edgeNormals[i][1];

    // if d is not almost parallel to edge
    if( dlen*dlen > 1e-9 * dsqlen* edgeNormalsSqlen[i])

      // if ray leaves cell here (does not enter)
      if((dlen>0) ^ inverted ){ // " ^ " means xor

	dlen = ( (pData[i][0]-s[0])*edgeNormals[i][0] +
		 (pData[i][1]-s[1])*edgeNormals[i][1] ) / dlen;
	
	if(dlen>0){
	  if(length > dlen) 
	    { length = dlen; reti=i; }
	}
	else
	  { length = 0; reti=i; }	      
      }
  }
  
  if(reti==FIndex::invalid)
    cout<<"no aprpriate face found(in FCellEdgeInfo::edgeToGo)"<<endl;
  return reti;
}

//------------------------------------------------------------------------------

FIndex FCell1Din3D::edgeToGo(double&length,
			 const FArray& start, 
			 const FArray& dir) const
{
  if(!edgeInfo)edgeInfo = new FCellEdgeInfo(this);
  return edgeInfo->edgeToGo(length,start,dir);
}

#endif
//==============================================================================

FCell1Din3D::~FCell1Din3D()
{
//  if(edgeInfo)delete edgeInfo;
}
