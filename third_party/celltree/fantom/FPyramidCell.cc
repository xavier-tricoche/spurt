//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPyramidCell.cc,v $
// Language:  C++
// Date:      $Date: 2004/03/11 13:53:13 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FPyramidCell.hh"
#include "FMatrix.hh"
#include "FRefArray.hh"
#include "FPosition.hh"

#include "FException.hh"

#include "FAMSingularPoint.hh"

#include <list>
#include <vector>
#include <cassert>
#include <iostream>

//macro library for 3d operations
//#include "Fd3op.hh"

// Imagine the FPyramidCell as a Pyramid with
// a base quadrilateral consisting of vertex indices 0,1,2,3
// in mathematical positive sense.
// Index 4 is the index of the top

//--------------------------------------------------------------------------- 

// the cell's geometry description
const FCell::geoInfo FPyramidCell::myGeoDescription =
  {
	// dimension
	3,
	// # vertices
	5, 
	// # edges
	8, 
	// # faces
	5, 
	// edges
	{{0,1},  {1,2},  {2,3}, {3,0},  
	 {0,4},  {1,4},  {2,4}, {3,4}},
	// face sizes
	{4,3,3,3,3},
	// faces
	{{0,1,2,3},  
	 {0,4,1}, {1,4,2}, {2,4,3}, {3,4,0}}
  };

//--------------------------------------------------------------------------- 

FPyramidCell::FPyramidCell()
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData[0],myTensors)
{
  lastPosition[0]=HUGE_VAL;
}

//--------------------------------------------------------------------------- 

FPyramidCell::FPyramidCell( const vector<FIndex>&vert )
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData[0],myTensors)
{
  for(int i=0;i<5;i++)
    vertexIndices[i] = vert[i];
  lastPosition[0]=HUGE_VAL;
}

FPyramidCell::FPyramidCell( const FIndex*vert )
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData[0],myTensors)
{
  for(int i=0;i<5;i++)
    vertexIndices[i] = vert[i];
  lastPosition[0]=HUGE_VAL;
}

//--------------------------------------------------------------------------- 

FPyramidCell::~FPyramidCell()
{
}

//--------------------------------------------------------------------------- 


positive FPyramidCell::sizeOfCellType() const
{
  return 5;
}
//--------------------------------------------------------------------------- 

FCell* FPyramidCell:: getClone() const
{
  cout<<"not implemented!"<<endl;
  assert(false);
}
//--------------------------------------------------------------------------- 

inline void 
FPyramidCell::computeLocalCoords(const double p[3]) const
{

  if((lastPosition[0]==p[0])&(lastPosition[1]==p[1])&(lastPosition[2]==p[2]))
    return;

  if(!geometryOK){
    bottom.init(myPositionData[0],myPositionData[1],
		myPositionData[2],myPositionData[3]);
    geometryOK=true;
  }

  //Fd3op2(lastPosition,=p,);
  for(int i=0; i<3; i++)
   lastPosition[i] = p[i];
   
   
  double top[3];
  
  //Fd3op2(top,=myPositionData[4],);
  for(int i=0; i<3; i++)
	top[i] = myPositionData[4][i];

  double s[3];//dir.vector

  double ret[2][3]; //return values

  //Fd3op3(s,=p,-top,);
  for(int i=0; i<3; i++)
	s[i] = p[i] - top[i];

  int n = bottom.cutWithLine(top,s,ret,false);

  if(n==0)
  {
    //Fd3op1(localCoords,=0);
	 for(int i=0; i<3; i++)
		localCoords[i] =0;
    return;
  }

  double * l = ret[0];

  if(n==2)
    if( fabs(ret[0][0]+ret[0][1]-1)+fabs(ret[0][0]-ret[0][1])
	>
	fabs(ret[1][0]+ret[1][1]-1)+fabs(ret[1][0]-ret[1][1])  
	){
	//wenn maxnorm(ret[0]-(.5,.5)) > maxnorm( ret[1]-(.5,.5)
       bottom.refineCutWithLine(top,s,ret[1]);
       l=ret[1];
       }
    else 
        bottom.refineCutWithLine(top,s,ret[0]);

  if(fabs(l[2]*epsilon)>1)
  {
    //Fd3op1(localCoords,=0);
	for(int i=0; i<3; i++)
		localCoords[i] =0;
    return;
  }


  localCoords[0]=l[0];
  localCoords[1]=l[1];
  localCoords[2]=1.0/l[2];

}

//--------------------------------------------------------------------------- 


bool 
FPyramidCell::
neighborFaceForPos(const FArray& pos,FIndex&faceId ) const
{
  computeLocalCoords(&pos[0]);
  
  //  cout<<localCoords[0]<<' '<<localCoords[1]<<' '<<localCoords[2]<<' ';

  double *l=localCoords;
  double a[3]={fabs(l[0]-0.5),fabs(l[1]-0.5),fabs(l[2]-0.5)};
 
 //if pos is inside, return false
  if((a[0]<0.5+epsilon)&(a[1]<0.5+epsilon)&(a[2]<0.5+epsilon)) 
    return false;
  
  if(fabs(l[2])<epsilon){
    faceId.setToInvalid();
    //    cout<<"not determinable"<<endl;
    return true;
  }

  a[0]*=fabs(l[2]);a[1]*=fabs(l[2]);
  
  int i=0;
  if(a[1]>a[0])i=1;
  if((a[2]>a[i])&(l[2]>0))i=2;
  
  switch(i)
    {
    case 0: 
      if( (l[0]>0.5)^(l[2]<0) ) 
	faceId=2; 
      else
	faceId=4;
      break;
    case 1:
      if( (l[1]>0.5)^(l[2]<0) )
	faceId=3;
      else
	faceId=1;
      break;
    case 2:
      faceId=0;
      break;
    }

  return true;

}

//--------------------------------------------------------------------------- 

void FPyramidCell::derivatives(FTensor& result, 
		 const FPosition& position) const
{
  try{
    positive dim=tensors[0].getDimension(),ord=tensors[0].getOrder();
    
#ifndef NODEBUG
    if(dim<3||position.size()<3)
      throw FInvalidDimensionException();
#endif
    result.resizeTensor(dim,ord+1);
    
    computeLocalCoords(&position[0]);    
    //is on top
    if(localCoords[2]==0)
    {
	cout << "local coord case!" << endl;

	localCoords[2] = 1e-6;
//       //if top of pyramid, make least squares for derivation
//       //(deriv after 4 edges here is linear)
//       double M[4][3];
//       FTensor tens[4];

//       for(int i=0;i<4;i++){
// 	Fd3op3(M[i],=myPositionData[4],-myPositionData[i],);
// 	tens[i]=(FTensor)myTensors[4]-(FTensor)myTensors[i];
//       }

//       // A * (x[i])_{0..3} = (tens[4]-tens[i])_{0..4} | A^T *()
//       //-> (x[i])_{0..3} = (A^T A)^{-1}*A^T * (tens[4]-tens[i])_{0..4}

//       FMatrix A(3,4,M[0]),AT=A.transposed();

//       A = invert(AT*A) * AT;

//       FTensor m(dim,ord);//multiplication buffer

//       for(int i=0;i<3;i++){
// 	FRefTensor t(result[i]);
// 	tens[0].mult(A(0,i,0),t);
// 	tens[1].mult(A(1,i,1),m); t+= m;
// 	tens[2].mult(A(2,i,2),m); t+= m;
// 	tens[3].mult(A(3,i,3),m); t+= m;
//       }	


    }
    
    FTensor derivTensor(dim,ord+1);
    dTensdLocal(derivTensor);
    
    double deriv[3][3];
    dPosdLocal(deriv);        
    
    FMatrix A(3,3,deriv[0]);
    
    A.invert();

    
    FTensor m(dim,ord);//multiplication buffer

    for(int i=0;i<3;i++){
      FRefTensor t(result[i]);
      derivTensor[0].mult(A(i,0),t);
      derivTensor[1].mult(A(i,1),m); t+= m;
      derivTensor[2].mult(A(i,2),m); t+= m;
    }
  }  
  catch(FException&e){
    e.addTraceMessage
      ("void FPyramidCell::derivatives"
       "(FTensor& result, const FPosition&p) const");
    throw;
  }
}

//--------------------------------------------------------------------------- 

void FPyramidCell::getZeros(list<FAMSingularPoint>& result) const
{
  cout<<"not implemented!"<<endl;
  assert(false);

  if(tensors[0].getDimension()!=3)
    THROW_EXCEPTION(FInvalidDimensionException,
		    "getZeros only works with dim 3 ord 1 tensors");

  const double
    *t0 = tensorData,
    *t1 = tensorData+3,
    *t2 = tensorData+6,
    *t3 = tensorData+9,
    *t4 = tensorData+12;

  FBilinearSurface bottom
    ( t0,t1,t2,t3 );
 

  double top[3];

  double p[3];

  
  //Fd3op2(top,=t4,);
	for(int i=0; i<3; i++)
		top[i] = t4[i];

  double s[3];//dir.vector

  double ret[2][3]; //return values

  //Fd3op2(s,=-top,);
	for(int i=0; i<3; i++)
		s[i] = - top[i];
		
  int n = bottom.cutWithLine(top,s,ret);


  for(int i=0;i<n;i++)
    {
      double *r=ret[i];
      if(r[0] > -epsilon  & r[0] < 1+epsilon && 
	 r[1] > -epsilon  & r[1] < 1+epsilon &&
	 r[2] > 1+epsilon | r[2] < -1/epsilon ){
	
	localCoords[0] = r[0];
	localCoords[1] = r[1];
	localCoords[2] = 1.0/r[2];


	double 
	  m0 = 1-localCoords[0],
	  m1 = 1-localCoords[1],
	  m1l2 = m1*localCoords[2],
	  l1l2 = localCoords[1]*localCoords[2];

	double 
	  c0 = m0*m1l2,
	  c1 = localCoords[0]* m1l2,
	  c2 = localCoords[0]* l1l2,
	  c3 = m0            * l1l2,
	  c4 = 1-localCoords[2];

	//Fd3op5(lastPosition, = c0*myPositionData[0], + c1*myPositionData[1], + c2*myPositionData[2],+ c3*myPositionData[3],);
	for(int i=0; i<3; i++)
		lastPosition[i]  = c0*myPositionData[0][i] + c1*myPositionData[1][i] + c2*myPositionData[2][i] + c3*myPositionData[3][i];
		
	//Fd3op2(lastPosition, += c4*myPositionData[4],);
	for(int i=0; i<3; i++)
		lastPosition[i] += c4*myPositionData[4][i];
	
	FArray pa(3,p);
	FAMSingularPoint ap(pa, FAMSingularPoint::NONE, 1);
	FTensor T(3,2);

	derivatives(T,pa);

	FMatrix M(T);

	cout<<M<<endl;

	ap.setLinearNature(M);
	result.push_back(ap);

	

      }
	 

    }


}



//--------------------------------------------------------------------------- 

void 
FPyramidCell::interpolate(FTensor& result, 
			       const FPosition& position) const
{
  

  computeLocalCoords(&position[0]);

  double 
    m0 = 1-localCoords[0],
    m1 = 1-localCoords[1],
    m1l2 = m1*localCoords[2],
    l1l2 = localCoords[1]*localCoords[2];

  result 
    = (m0            * m1l2) * myTensors[0] 
    + (localCoords[0]* m1l2) * myTensors[1]
    + (localCoords[0]* l1l2) * myTensors[2]
    + (m0            * l1l2) * myTensors[3] 
    + (1-localCoords[2])     * myTensors[4];

}
  
//--------------------------------------------------------------------------- 

inline void 
FPyramidCell::
dPosdLocal(double3 m[]) const
{
  //aliases
  const double3*const p = myPositionData;
  const double*const l = localCoords;

  double 
    m0 = 1-l[0],
    m1 = 1-l[1],

    m1l2 = m1*l[2],
    l1l2 = l[1]*l[2],

    m0l2 = m0*l[2],
    l0l2 = l[0]*l[2],

    l0m1 = l[0]*m1,
    m0l1 = m0*l[1],
    l0l1 = l[0]*l[1],
    m0m1 = m0*m1;
 
  
  //Fd3op5(m[0], = -m1l2 * p[0], +m1l2 * p[1], +l1l2 * p[2], -l1l2 * p[3],);
  	for(int i=0; i<3; i++)
		m[0][i] = -m1l2 * p[0][i] + m1l2 * p[1][i] + l1l2 * p[2][i] - l1l2 * p[3][i];
		
  //Fd3op5(m[1], = -m0l2 * p[0], -l0l2 * p[1], +l0l2 * p[2], +m0l2 * p[3],);
	for(int i=0; i<3; i++)
		m[1][i] = -m0l2 * p[0][i] - l0l2 * p[1][i] +l0l2 * p[2][i] +m0l2 * p[3][i];

  //Fd3op5(m[2], = +m0m1 * p[0], +l0m1 * p[1], +l0l1 * p[2], +m0l1 * p[3],);
	for(int i=0; i<3; i++)
		m[2][i] = +m0m1 * p[0][i] +l0m1 * p[1][i] +l0l1 * p[2][i] +m0l1 * p[3][i];
		
  //Fd3op2(m[2], -= p[4],);
	for(int i=0; i<3; i++)
		m[2][i] -= p[4][i];
}

//--------------------------------------------------------------------------- 

inline void 
FPyramidCell::
dTensdLocal(FTensor&result) const
{
  result.resizeTensor(myTensors[0].getDimension(),myTensors[0].getOrder()+1);

  const FTensor*const t = myTensors;
  const double*const l = localCoords;
  
  double 
    m0 = 1-l[0],
    m1 = 1-l[1],

    m1l2 = m1*l[2],
    l1l2 = l[1]*l[2],

    m0l2 = m0*l[2],
    l0l2 = l[0]*l[2],

    l0m1 = l[0]*m1,
    m0l1 = m0*l[1],
    l0l1 = l[0]*l[1],
    m0m1 = m0*m1;
  
  result[0] = -m1l2 * t[0] +m1l2 * t[1] +l1l2 * t[2] -l1l2 * t[3];
  result[1] = -m0l2 * t[0] -l0l2 * t[1] +l0l2 * t[2] +m0l2 * t[3];
  result[2] = +m0m1 * t[0] +l0m1 * t[1] +l0l1 * t[2] +m0l1 * t[3];
  result[2] -= t[4];
}
//--------------------------------------------------------------------------- 

bool FPyramidCell::isInside(const FPosition& position) const
{
    ++FCell::ntests;

  computeLocalCoords(&position[0]);

  return 
    (-epsilon < localCoords[0]) & (localCoords[0] < 1+epsilon) &
    (-epsilon < localCoords[1]) & (localCoords[1] < 1+epsilon) &
    (-epsilon < localCoords[2]) & (localCoords[2] < 1+epsilon) ;
}

//--------------------------------------------------------------------------- 

FCell::CellType FPyramidCell::getCellType(void) const
{
  return FCell::PYRAM;
}


positive FPyramidCell::memSize() const
{
  return 0;
}


