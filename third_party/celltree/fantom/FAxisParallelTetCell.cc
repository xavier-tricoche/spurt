//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAxisParallelTetCell.cc,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:00 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FAxisParallelTetCell.hh"
#include "FMatrix.hh"
#include "FPosition.hh"

#include "FException.hh"
#include "FAMSingularPoint.hh"

#include <list>
#include <vector>
#include <cassert>

const FCell::geoInfo FAxisParallelTetCell::myGeoDescription=
  {
        // dimension
        3,
        // # vertices
        4,
        // # edges
        6,
        // # faces
        4,
        // edges
	{{0,1},  {1,2},  {2,0},  {0,3},  {1,3},  {2,3}},
        // face sizes
        {3, 3, 3, 3},
        // faces
        {{0,1,2},  {0,3,1},  {1,3,2},  {2,3,0}}
  };

FAxisParallelTetCell::FAxisParallelTetCell()
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
}
 
//--------------------------------------------------------------------------- 

FAxisParallelTetCell::FAxisParallelTetCell(const vector<FIndex>& vert,
					   int ls,int lm,int lb)  
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  s=ls;m=lm;b=lb;
  for(int i=0;i<4;i++) vertexIndices[i]=vert[i];
}

//--------------------------------------------------------------------------- 
FAxisParallelTetCell::FAxisParallelTetCell(const FIndex* vert,
					   int ls,int lm,int lb)  
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  s=ls;m=lm;b=lb;
  for(int i=0;i<4;i++) vertexIndices[i]=vert[i];
}

//--------------------------------------------------------------------------- 

FAxisParallelTetCell::~FAxisParallelTetCell()
{
}

//--------------------------------------------------------------------------- 

FCell* FAxisParallelTetCell::getClone() const 
{
  try 
    {
      return new FAxisParallelTetCell( vertexIndices,s,m,b);
    }
  catch (FException& e) 
    {
      e.addTraceMessage("FCell* FAxisParallelTetCell::getClone() const");
      throw e;
    }
}

//--------------------------------------------------------------------------- 

void FAxisParallelTetCell::interpolate(FTensor& result, 
				       const FPosition& position) const
{
  try{

    if(!geometryOK)
      {
	buildBoundingBox();

	ps=positions[0][s];
	pm=positions[0][m];
	pb=positions[0][b];

	invds=1/(positions[3][s]-ps);
 	invdm=1/(positions[3][m]-pm);
 	invdb=1/(positions[3][b]-pb);    	

	geometryOK=true;
      }


    double
      xs=(position[s]-ps)*invds,
      xm=(position[m]-pm)*invdm,
      xb=(position[b]-pb)*invdb;


   result
    = tensors[0] * ( 1-xb)
    + tensors[1] * ( xb-xm )
    + tensors[2] * ( xm-xs )
    + tensors[3] * xs;   
      
  }

  catch(FException&e){
    e.addTraceMessage("void FAxisParallelTetCell::"
		      "interpolate(FTensor& result, const FPosition& position) const");
    throw e;  
  }
}
 
//--------------------------------------------------------------------------- 

void FAxisParallelTetCell::derivatives(FTensor& result, 
				       const FPosition&/*p*/) const
{
  try{

    assert( !tensors[0].getOrder() || tensors[0].getDimension()==3);

    if(!geometryOK)
      {
	buildBoundingBox();

	ps=positions[0][s];
	pm=positions[0][m];
	pb=positions[0][b];

	invds=1/(positions[3][s]-ps);
 	invdm=1/(positions[3][m]-pm);
 	invdb=1/(positions[3][b]-pb);    	

	geometryOK=true;
      }

    if(!interpolationOK)
      {

	deriv.resizeTensor(3, tensors[0].getOrder()+1);
	deriv[s] = (tensors[3]-tensors[2])*invds;
	deriv[m] = (tensors[2]-tensors[1])*invdm;
	deriv[b] = (tensors[1]-tensors[0])*invdb;

	interpolationOK=true;
      }

    result=deriv;
   
 }

  catch(FException&e){
    e.addTraceMessage("void FAxisParallelTetCell::derivatives(FTensor& result, const FPosition&p) const");
    throw e;
  }
      
}

//--------------------------------------------------------------------------- 
 
bool FAxisParallelTetCell::isInside(const FPosition& position) const
{
  try{

    assert( position.getDimension() == 3 );

    if(!geometryOK)
      {
	buildBoundingBox();

	ps=positions[0][s];
	pm=positions[0][m];
	pb=positions[0][b];

	invds=1/(positions[3][s]-ps);
 	invdm=1/(positions[3][m]-pm);
 	invdb=1/(positions[3][b]-pb);    	

	geometryOK=true;
      }


    double
      xs=(position[s]-ps)*invds,
      xm=(position[m]-pm)*invdm,
      xb=(position[b]-pb)*invdb;
  
    return 
      0<=xs && xs<=xm && xm<=xb && xb<=1;

  }

  catch(FException&e){
    e.addTraceMessage("bool FAxisParallelTetCell::isInside(const FPosition& position) const");
    throw e;
  }
}  

//--------------------------------------------------------------------------- 

bool FAxisParallelTetCell::isInside(const FPosition& /*pos*/, 
				    const vector< FPosition >& /*vertices*/ ) {
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}


//--------------------------------------------------------------------------- 

FCell::CellType FAxisParallelTetCell::getCellType(void) const
{
  return FCell::AXIS_PARALLEL_TET;
}

//--------------------------------------------------------------------------- 

void FAxisParallelTetCell::getZeros( list<FAMSingularPoint>& /*result*/) const
{
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//--------------------------------------------------------------------------- 

void FAxisParallelTetCell::buildBoundingBox(void) const
{

  try{
    FPosition &p0 = positions[0], &p1= positions[3];
    bBox.setBoundingBox(p0[0],p0[1],p0[2],
			p1[0],p1[1],p1[2]);
  }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

positive FAxisParallelTetCell::
sizeOfCellType () const
{
  return 4;
}

//--------------------------------------------------------------------------- 
positive FAxisParallelTetCell::memSize() const
{
  return
    (tensorData?
     tensors[0].size()*sizeof(double)*
     (sizeOfCellType() //tensorData
      +tensors[0].getDimension() //deriv
      )
     :0)
    +
    sizeof(*this);
}
