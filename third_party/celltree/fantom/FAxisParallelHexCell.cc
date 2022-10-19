//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAxisParallelHexCell.cc,v $
// Language:  C++
// Date:      $Date: 2004/03/11 13:53:13 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FAxisParallelHexCell.hh"
#include "FMatrix.hh"
#include "FPosition.hh"

#include "FException.hh"
#include "FAMSingularPoint.hh"
#include "Fd3op.hh"
#include <iostream>
using namespace std;

#include <list>
#include <vector>

FAxisParallelHexCell::FAxisParallelHexCell()
{
}
 
//--------------------------------------------------------------------------- 

FAxisParallelHexCell::FAxisParallelHexCell(const vector<FIndex>& vert,
					   bool vtkhexcell)
  :FHexahedronCell(vert,vtkhexcell)
{
}
//--------------------------------------------------------------------------- 


FAxisParallelHexCell::FAxisParallelHexCell(const FIndex* vert)
  :FHexahedronCell(vert)
{
}

//--------------------------------------------------------------------------- 

FAxisParallelHexCell::~FAxisParallelHexCell()
{
}

//--------------------------------------------------------------------------- 

FCell* FAxisParallelHexCell::getClone() const 
{
  try 
    {
      return new FAxisParallelHexCell( vertexIndices );
    }
  catch (FException& e) 
    {
      e.addTraceMessage("FCell* FAxisParallelHexCell::getClone() const");
      throw e;
    }
}

//--------------------------------------------------------------------------- 

void FAxisParallelHexCell::interpolate(FTensor& result, 
				       const FPosition& position) const
{
  try{

    if(!tensors||!positions)
      throw FException(" positions or tensors were not set !");

    if(!interpolationOK)
      computeTensorParams();

    result.resizeTensor(tensors[0].getDimension(),tensors[0].getOrder());

    FArray localCoord
      ( (position[0]-positions[0][0])  *invLen[0],
	(position[1]-positions[0][1])  *invLen[1],
	(position[2]-positions[0][2])  *invLen[2] );
    
    tensorParams.interpolateT(localCoord,result);    

  }

  catch(FException&e){
    e.addTraceMessage("void FAxisParallelHexCell::interpolate(FTensor& result, const FPosition& position) const");
    throw e;  
  }
}
 
//--------------------------------------------------------------------------- 

void FAxisParallelHexCell::derivatives(FTensor& result, 
				       const FPosition&position) const
{
  try{

#ifndef NODEBUG
    if(tensors[0].getDimension()!=3)
      {
	throw FInvalidDimensionException();
      }
#endif

    if(!positions||!tensors)
      throw FException("positions or tensors were not set");

    //compute interpolation parameters for tensors

    if(!interpolationOK)
      computeTensorParams();

    FArray localCoord
      ( (position[0]-positions[0][0])  *invLen[0],
	(position[1]-positions[0][1])  *invLen[1],
	(position[2]-positions[0][2])  *invLen[2] );

    result.resizeTensor(3,tensors[0].getOrder()+1);

    FRefTensor a(result[0]),b(result[1]),c(result[2]);

    tensorParams.interpolatedTdx(localCoord,a);  
    a*=invLen[0];
    tensorParams.interpolatedTdy(localCoord,b);  
    b*=invLen[1];
    tensorParams.interpolatedTdz(localCoord,c);  
    c*=invLen[2];

 }

  catch(FException&e){
    e.addTraceMessage("void FAxisParallelHexCell::derivatives(FTensor& result, const FPosition&p) const");
    throw e;
  }
      
}

//--------------------------------------------------------------------------- 
 
bool FAxisParallelHexCell::isInside(const FPosition& position) const
{
  try{

    if( position.getDimension() != 3 )
      THROW_EXCEPTION( FInvalidDimensionException, "Position Dimension has to be 3" );
  
    return getBoundingBox().isInside(position);

  }

  catch(FException&e){
    e.addTraceMessage("bool FAxisParallelHexCell::isInside(const FPosition& position) const");
    throw e;
  }
}  

//--------------------------------------------------------------------------- 
 
bool FAxisParallelHexCell::isInside(const FPosition& pos,
				    const vector< FPosition >& vertices,
				    bool useVTKHexahedronEnumeration)
				    
{
  if(useVTKHexahedronEnumeration)
    return (pos[0] > vertices[0][0] - epsilon && 
	    pos[0] < vertices[6][0] + epsilon &&
	    pos[1] > vertices[0][1] - epsilon &&
	    pos[1] < vertices[6][1] + epsilon &&
	    pos[2] > vertices[0][2] - epsilon &&
	    pos[2] < vertices[6][2] + epsilon);
  else
    return (pos[0] > vertices[0][0] - epsilon && 
	    pos[0] < vertices[7][0] + epsilon &&
	    pos[1] > vertices[0][1] - epsilon &&
	    pos[1] < vertices[7][1] + epsilon &&
	    pos[2] > vertices[0][2] - epsilon &&
	    pos[2] < vertices[7][2] + epsilon);
}

//--------------------------------------------------------------------------- 

FCell::CellType FAxisParallelHexCell::getCellType(void) const
{
  return FCell::AXIS_PARALLEL_HEX;
}


//--------------------------------------------------------------------------- 

void FAxisParallelHexCell::computeTensorParams() const
{
  try{

    if(!tensors||!positions)
      throw FException("tensors or positions not set");

    invLen[0]=1.0/(positions[1][0]-positions[0][0]);
    invLen[1]=1.0/(positions[3][1]-positions[0][1]);
    invLen[2]=1.0/(positions[4][2]-positions[0][2]);
  

    vector<double> tensorP(tensors[0].size()*8);
    vector<double>::iterator tp=tensorP.begin();

    //fill tensor param buffer
    for(positive i=0;i<tensors[0].size();i++)
      for(positive j=0;j<8;j++,tp++)
	*tp=FArray(tensors[j])[i];

    tensorParams.buildParameters(tensorP);

    interpolationOK = true;

  }

  catch(FException &e){
    e.addTraceMessage
      ("void FAxisParallelHexCell::computeTensorParams() const");
    throw e;
  }  

}

//----------------------------------------------------------------------------

void FAxisParallelHexCell::buildBoundingBox(void) const
{

  try{
    FPosition &p0 = positions[0], &p1= positions[6];
    bBox.setBoundingBox(p0[0],p0[1],p0[2],
			p1[0],p1[1],p1[2]);
  }

  catch( FException e){
    e.addTraceMessage("void FAxisParallelHexCell::buildBoundingBox(void) const");
    throw e;
    }
}

positive FAxisParallelHexCell::memSize() const
{
  return
    (tensorData?
     tensors[0].size()*sizeof(double)*
     (sizeOfCellType()*2 ):0)
    +
    sizeof(*this);
}


void FAxisParallelHexCell::
computeSing(float minb[3],float maxb[3],list<FAMSingularPoint>& result) const
{
  try{

    
    if(!interpolationOK){
      computeTensorParams();      
    }


    float deltmb[3];
    //Fd3op3(deltmb,=maxb,-minb,);
    for(int i=0; i<3; i++)
		deltmb[i] = maxb[i] - minb[i];
	
	
    //    cout<<Fd3out(minb)<<' '<<d3out(maxb)<<' '<<d3out(deltmb)<<endl;


  FArray lx(3); //local coords;
  //Fd3op3(lx,=minb,+maxb,);
   for(int i=0; i<3; i++)
		lx[i] = minb[i] + maxb[i];
  lx*=0.5;

  
  FMatrix M(3,3);
  FTensor t(3,2);
  FArray x(3);

  tensorParams.interpolateT(lx,x);

  //  cout<<"local tensor"<<x<<endl;

  FVector olx(3);
  double ldelt=x.normSquare(),odelt=ldelt;
  double sqeps=epsilon*epsilon;
  int i=0;
  
  while(ldelt>sqeps){

    //if precision isn't reached with
    //the same matrix in the next step then recompute Matrix
    //( if improvement in last step is smaller than 
    //  improvement that still has to be achieved (odelt/ldelt < ldelt/sqeps)
    if(sqeps*odelt<ldelt*ldelt){
      //compute derivation of world koords after local koords
      
      FRefTensor a(t[0]),b(t[1]),c(t[2]);
      tensorParams.interpolatedTdx(lx,a);  
      tensorParams.interpolatedTdy(lx,b);  
      tensorParams.interpolatedTdz(lx,c);  
      M=t;
      //invert matrix
      M.invert();
    }
    //      cout<<"inverted"<<endl;
      //determine new point as if the derivation would be constant
    
    olx = lx;
    lx -= M * x;
    
    //      cout<<"interpolate using new trilcoords"<<endl;
    
    tensorParams.interpolateT(lx,x);

    odelt = ldelt;
    ldelt = x.normSquare();
    
    //if raw newton does not converge, use half stepsize
    while( ldelt >= odelt && i<100){
      lx=(olx+lx)*0.5;
      tensorParams.interpolateT(lx,x);
      ldelt = x.normSquare();	
      i++;
    }
    i++;
    if(i>=100){
      cout<<"reached precision:"<<sqrt(ldelt)
	  << " lx: " <<lx<<" X :" << x;
      throw FTooManyIterationsException();
    }    
  }


   if( lx[0]>minb[0]-epsilon && lx[0] < maxb[0]+epsilon
       && lx[1]>minb[1]-epsilon && lx[1] < maxb[1]+epsilon
       && lx[2]>minb[2]-epsilon && lx[2] < maxb[2]+epsilon )
    {
      //Fd3op4(x,=myPositionData,+lx,/invLen,);
	    for(int i=0; i<3; i++)
			x[i] = myPositionData[i] + lx[i] / invLen[i];
	
      FAMSingularPoint ap(x,FAMSingularPoint::NONE,1);

      FMatrix tens(3,3);
      FRefTensor 
	a(3,1,&tens(0,0)),
	b(3,1,&tens(1,0)),
	c(3,1,&tens(2,0));

      tensorParams.interpolatedTdx(lx,a);  
      tensorParams.interpolatedTdy(lx,b);  
      tensorParams.interpolatedTdz(lx,c);  


      ap.setLinearNature(tens);

      result.push_back(ap);

    }  
  }catch(FException&e)
    {
      e.addTraceMessage("FAxisParallelHexCell::computeSing");
      throw;
    }
}
