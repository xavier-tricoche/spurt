//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPrismCell.cc,v $
// Language:  C++
// Date:      $Date: 2003/09/30 14:47:49 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FPrismCell.hh"

#include "FMath.hh"
#include "FMatrix.hh"
#include "FRefArray.hh"
#include "FPosition.hh"

#include "FException.hh"
#include "FAMSingularPoint.hh"
#include <list>
#include <vector>
#include <sstream>
#include <iostream>

//macro library for 3d operations,in directory math
#include "Fd3op.hh"

// Imagine the FPrismCell as a Prism with
// a base triangle consisting of vertex indices 0,1,2
// in mathematical positive sense.
// Indices 3,4,5 are the indices of the top triangle

//--------------------------------------------------------------------------- 

// the cell's geometry description
const FCell::geoInfo FPrismCell::myGeoDescription =
  {
	// dimension
	3,
	// # vertices
	6, 
	// # edges
	9, 
	// # faces
	5, 
	// edges
	{{0,1},  {1,2},  {2,0},  
	 {5,4},  {4,3},  {3,5},  
	 {0,3},  {1,4},  {2,5}},
	// face sizes
	{3, 3, 4, 4, 4},
	// faces
	{ {0,1,2},  {5,4,3},
	  {0,3,4,1}, {1,4,5,2}, {2,5,3,0} }
  };

//--------------------------------------------------------------------------- 
FPrismCell::FPrismCell( const FIndex*vert )
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData[0],myTensors)
{
  for(int i=0;i<6;i++)
    vertexIndices[i] = vert[i];
}

//--------------------------------------------------------------------------- 

FPrismCell::FPrismCell( const vector<FIndex>&vert )
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData[0],myTensors)
{
  for(int i=0;i<6;i++)
    vertexIndices[i] = vert[i];
}


//--------------------------------------------------------------------------- 

FPrismCell::FPrismCell()
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData[0],myTensors)
{
}

//--------------------------------------------------------------------------- 
void FPrismCell::getZeros(list<FAMSingularPoint>& result) const
{
  try{

    for(int i=0;i<6;i++)
      if(tensors[i].normSquare()<1e-20)
	{cout<<'.'<<flush;return;}


    //    throw FNotImplementedException("");

    if(tensors[0].size()<3)
      throw FInvalidDimensionException("getzeros only works with tensors of dim 3 order 1");


    computeTensorParams();
    computeParams();

    const double *ze0=&te0(0),*ze1=&te1(0),*ze2=&te2(0);

    //derivations of the above after local coord 2
    const double *zde0=&tde0(0),*zde1=&tde1(0);

    //params for cubic equation
    double3 
      za,zb,zc,
      a,b,c,d,e;

    //scalar product with e2
    double ze2a,ze2b,ze2c;

    // a copy from computeParams()
    Fd3kreuz(zc,ze0,ze1);

    Fd3kreuz(a,zde0,ze1);
    Fd3kreuz(b,ze0,zde1);
    //Fd3op3(zb,=a,+b,);
	 for(int i=0; i<3; i++)
		zb[i] = a[i] + b[i];

    Fd3kreuz(za,zde0,zde1);

    ze2a = Fd3prod(ze2,za);
    ze2b = Fd3prod(ze2,zb);
    ze2c = Fd3prod(ze2,zc);

   
    complex<double> sols[3];

    FArray & tens0 = tensors[0];


    //Fd3op2(a,= -tens0,);
	for(int i=0; i<3; i++)
		a[i] = -tens0[i];
	

    // solution of the equation:
    // < p - pos0 + l2*pe2)  , ( (pe0+l2*pde0) , (pe1+l2*pde1) ) >   = 0
    // respective l2
    // with: pos0 = myPositionData[0], l2= l[2]; 
    // the rest of the vars see computeParams

    int nsols = FMath::CubicEquation
      ( -ze2a , Fd3prod(a,za)-ze2b, Fd3prod(a,zb)-ze2c, Fd3prod(a,zc) , sols);  
    
    if(nsols==-1)
      return;

    //find real solutions between 0 and 1
    int i  ;
    for(i=0;i<nsols;i++){
      complex<double> & s = sols[i];
      if( fabs(s.imag()) < 2*epsilon 
	  && -2*epsilon <= s.real() && s.real() <= 1+2*epsilon ){

	localCoords[2] = s.real();


	for(int i=0; i<3; i++)
	{
		//a = vector from first vertex of triangle to searched tensor
		//(here 0 )
		//Fd3op3(a,= -tens0,-ze2,*localCoords[2]);
		a[i] = -tens0[i] - ze2[i] * localCoords[2];
		
		//b = edge0 in height of l2
		//Fd3op3(b,=ze0,+zde0,*localCoords[2]);
		b[i] = ze0[i] + zde0[i] * localCoords[2];

		//c = edge1 in height of l2
		//Fd3op3(c,=ze1,+zde1,*localCoords[2]);
		c[i] = ze1[i] + zde1[i] * localCoords[2];
	}
	
	//d = normal vector of triangle (crossproduct of two edges)
	Fd3kreuz(d,b,c);
	double invDenom=1.0/Fd3prod(d,d);
    
	//e = normal vector to edge 1
	Fd3kreuz(e,c,d);
	localCoords[0] = Fd3prod(a,e) * invDenom;

	//e = normal vector to edge 0
	Fd3kreuz(e,d,b);
	localCoords[1] = Fd3prod(a,e) * invDenom;



	//.........................................................................
	//refine via newton iteration
  	
	static double3 
	  * mat, //jacobian
	  invMat[3]; //inverse,transposed jacobian

	
	FTensor ta(3,1);
	FTensor tMat(3,2);


	interpolateTensorsFromLocals(ta);

// 	if( localCoords[0] > -epsilon 
// 	    && localCoords[1]> -epsilon 
// 	    && 1-localCoords[0]-localCoords[1] > -epsilon )
// 	  {
// 	    cout<<"local coords"<<Fd3out(localCoords)<<endl;

// 	    cout<<"act. tens. val:"<<ta<<endl;

// 	    cout<<"te0 te1 te2 tde0 tde1 "
// 		<<Fd3out(ze0)
// 		<<Fd3out(ze1)
// 		<<Fd3out(ze2)
// 		<<Fd3out(zde0)
// 		<<Fd3out(zde1)
// 		<<endl;
// 	  }

	//Fd3op2(a,= FArray(ta),);
	for(int i=0; i<3; i++)
		a[i] = FArray(ta)[i];

	double delt = Fd3prod(a,a);
	double odelt = delt;
	double eps 
	  = (tensors[0].normSquare()+tensors[3].normSquare()) * 1e-28;

	//num.prec must be better than 1e-14
	if(delt > eps){
	  
	  dTensdLocal(tMat);
	  mat=(double3*)&tMat(0,0);
	  
	  //compute transposed inverse
	  Fd3kreuz(invMat[0],mat[1],mat[2]);
	  Fd3kreuz(invMat[1],mat[2],mat[0]);
	  Fd3kreuz(invMat[2],mat[0],mat[1]);	
    
	  double invdenom = 1.0/Fd3prod(invMat[0],mat[0]);

	  if(!isfinite(invdenom*1e5)){
	    cout<<"not a singularity _point_ "<<endl;
	    cout<<"at local coords"<<Fd3out(localCoords)<<endl;
	    cout<<"at tensor value: "<<ta<<endl;
	    FPosition pos(3);
	    interpolatePositionFromLocals(&pos[0]);
	    cout<<"and position value: "<<pos<<endl;


	    
	    continue;
	  }
	  do{
	    
	    //localCoords-= a*invMat
	    localCoords[0]-=Fd3prod(a,invMat[0])*invdenom;
	    localCoords[1]-=Fd3prod(a,invMat[1])*invdenom;
	    localCoords[2]-=Fd3prod(a,invMat[2])*invdenom;
	    
	    interpolateTensorsFromLocals(ta);
	    //Fd3op2(a,=FArray(ta),);
		for(int i=0; i<3; i++)
			a[i] = FArray(ta)[i];

	    odelt=delt;
	    delt=Fd3prod(a,a);	  
	    
	  }while(odelt>2*delt & delt>eps);

	}
		

	//.........................................................................
	


	if( localCoords[0] > -epsilon 
	    && localCoords[1]> -epsilon 
	    && 1-localCoords[0]-localCoords[1] > -epsilon )
	  {
	    FPosition pos(3);

	    interpolatePositionFromLocals(&pos[0]);


	    result.push_back(FAMSingularPoint(pos,FAMSingularPoint::NONE,1));
	  
	    //to force no recomputation of local coords
	    //in derivatives
	    bool oPC=paramsComputed;
	    paramsComputed=true;
	    bool oGOK=geometryOK;
	    geometryOK=true;
	    //Fd3op2(lastPosition,=pos,);
		for(int i=0; i<3; i++)
			lastPosition[i] = pos[i];

//  	    interpolateTensorsFromLocals(ta);
//  	    cout<<"testing getzeros:"<<ta<<endl;
	    derivatives(tMat,pos);

	    geometryOK=oGOK;
	    paramsComputed=oPC;
  	    FMatrix M(3,3,&tMat(0,0));

	    //	    cout<<M<<endl;

	    //	    result.back().setLinearNature(M);
	  	  
	  }
       
      }
    }

  }catch(FException&e)
    {
      e.addTraceMessage("FPrismCell::getZeros");
      throw;
    }

}

//--------------------------------------------------------------------------- 

FCell* FPrismCell::getClone() const
{
  throw FNotImplementedException("getClone");
}
//--------------------------------------------------------------------------- 

FPrismCell::~FPrismCell()
{
}


//--------------------------------------------------------------------------- 

bool FPrismCell::isInside(const FPosition& position) const
{
    ++FCell::ntests;

  try{
    if(position.size()<3){
      FInvalidDimensionException e;
    }
  
    if(!geometryOK)
      {
	paramsComputed=false;
	geometryOK=true;
	buildBoundingBox();
      }
    
    if(!bBox.isInside(position))
      return false;
    
    computeLocalCoords(&position[0]);
    
    return 
      localCoords[0] > -epsilon && localCoords[1] > -epsilon 
      && localCoords[0] + localCoords[1] < 1+epsilon
      && -epsilon < localCoords[2] && localCoords[2] < 1+epsilon;    
    
  }catch(FException&e)
    {
      e.addTraceMessage
	("bool FPrismCell::isInside(const FPosition& position) const");
      throw;
    }
}

//--------------------------------------------------------------------------- 

bool
FPrismCell::
neighborFaceForPos(const FArray& pos,FIndex&faceId ) const
{
  try{
    computeLocalCoords(&pos[0]);
  }catch(FInvalidPositionException){
    faceId.setToInvalid();
    return true;
  }


  double a[4]={ localCoords[0],
	        localCoords[1],
	        1-localCoords[0]-localCoords[1],
	        0.5-fabs(localCoords[2]-0.5) };

  int i=0;
  if(a[1]<a[0])i=1;
  if(a[2]<a[i])i=2;
  if(a[3]<a[i])i=3;

  if(a[i]>-epsilon)
    return false;
  
  switch(i)
    {
    case 0: faceId=4;break;
    case 1: faceId=2;break;
    case 2: faceId=3;break;
    case 3: 
      if(localCoords[2]>0.5)
	faceId=1;
      else 
	faceId=0;
      break;
    }
  
  return true;

}

//--------------------------------------------------------------------------- 


inline void FPrismCell::
interpolatePositionFromLocals(double result[3]) const
{    
    double*l=localCoords;
    //Fd3op3(a,=pe0,+l[2]*pde0,);
    //Fd3op3(b,=pe1,+l[2]*pde1,);
	//Fd3op5(result,=myPositionData[0],+ pe2,*l[2]+a,*l[0]+b,*l[1]);    
	for(int i=0; i<3; i++)
	{
		a[i] = pe0[i] + l[2]*pde0[i];
		b[i] = pe1[i] + l[2]*pde1[i];
		result[i] = myPositionData[0][i] + pe2[i] *l[2]+a[i] * l[0]+b[i] *l[1];    
	}
}
//--------------------------------------------------------------------------- 

inline void FPrismCell::
interpolateTensorsFromLocals(FTensor& result) const
{
  result 
    = tensors[0] + te2 * localCoords[2]
    + (te0 + localCoords[2]*tde0) * localCoords[0]
    + (te1 + localCoords[2]*tde1) * localCoords[1];  
}

//--------------------------------------------------------------------------- 


void FPrismCell::interpolate(FTensor& result, 
			      const FPosition& position) const
{
  try{

    if(position.size()<3)
      throw FInvalidDimensionException("wrong search pos size");    

    computeLocalCoords(&position[0]);
    
    if(!interpolationOK)
      computeTensorParams();

    interpolateTensorsFromLocals(result);

    
  }
  catch(FException &e)
    {
      e.addTraceMessage("void FPrismcell::interpolate\n"
			"(FTensor& result,const FPosition& position) const");
      throw;
    }
}


//--------------------------------------------------------------------------- 
void FPrismCell::derivatives(FTensor& result, 
  		   const FPosition& position) const
{
  try{
    positive dim=tensors[0].getDimension(),ord=tensors[0].getOrder();

#ifndef NODEBUG
    if(dim<3||position.size()<3)
      throw FInvalidDimensionException();
#endif

    computeLocalCoords(&position[0]);    
    
    FTensor derivTensor(dim,ord+1);
    dTensdLocal(derivTensor);
    
    double deriv[3][3];
    dPosdLocal(deriv);        

    FMatrix A(3,3,deriv[0]);
    
    A.invert();
    result.resizeTensor(dim,ord+1);
    
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
      ("void FPrismCell::derivatives"
       "(FTensor& result, const FPosition&p) const");
    throw;
  }


}
//--------------------------------------------------------------------------- 

void FPrismCell::computeTensorParams() const
{
  te0=tensors[1]-tensors[0];
  te1=tensors[2]-tensors[0];
  te2=tensors[3]-tensors[0];
  
  tde0=tensors[4]-tensors[3];tde0-=te0;
  tde1=tensors[5]-tensors[3];tde1-=te1;
 
  interpolationOK=true;
}



inline void FPrismCell::dTensdLocal(FTensor&result) const
{
  if(!interpolationOK)
    computeTensorParams();

  result[0] = te0 + localCoords[2] * tde0;
  result[1] = te1 + localCoords[2] * tde1;
  result[2] = te2 + localCoords[0] * tde0 + localCoords[1] * tde1;     
}

//--------------------------------------------------------------------------- 

inline void FPrismCell::dPosdLocal(double3 m[]) const
{
  //Fd3op3(m[0],=pe0,+pde0,*localCoords[2]);
  //Fd3op3(m[1],=pe1,+pde1,*localCoords[2]);
  //Fd3op4(m[2],=pde0,*localCoords[0]+pde1,*localCoords[1]+pe2,);
	for(int i=0; i<3; i++)
	{
		m[0][i] = pe0[i] + pde0[i] * localCoords[2];
		m[1][i] = pe1[i] + pde1[i] * localCoords[2];
		m[2][i] = pde0[i] * localCoords[0]+pde1[i] * localCoords[1]+pe2[i];
	}
}

//--------------------------------------------------------------------------- 
inline void FPrismCell::computeParams() const
{
  //Fd3op3(pe0,=myPositionData[1],-myPositionData[0],);
  //Fd3op3(pe1,=myPositionData[2],-myPositionData[0],);
  //Fd3op3(pe2,=myPositionData[3],-myPositionData[0],);

  //Fd3op4(pde0,=myPositionData[4],-myPositionData[3],-pe0,);
  //Fd3op4(pde1,=myPositionData[5],-myPositionData[3],-pe1,);
	for(int i=0; i<3; i++)
	{
		pe0[i] = myPositionData[1][i] - myPositionData[0][i];
		pe1[i] = myPositionData[2][i] - myPositionData[0][i];
		pe2[i] = myPositionData[3][i] - myPositionData[0][i];
		
		pde0[i] = myPositionData[4][i] - myPositionData[3][i] - pe0[i];
		pde1[i] = myPositionData[5][i] - myPositionData[3][i] - pe1[i];
	}
    
  Fd3kreuz(pc,pe0,pe1);

  Fd3kreuz(a,pde0,pe1);
  Fd3kreuz(b,pe0,pde1);
  //Fd3op3(pb,=a,+b,);
	for(int i=0; i<3; i++)
		pb[i] = a[i]+b[i];
	
  Fd3kreuz(pa,pde0,pde1);

  pe2a = Fd3prod(pe2,pa);
  pe2b = Fd3prod(pe2,pb);
  pe2c = Fd3prod(pe2,pc);

  //Fd3op1(lastPosition,=HUGE_VAL);
	for(int i=0; i<3; i++)
		lastPosition[i] = HUGE_VAL;
		
  paramsComputed=true;
}

//--------------------------------------------------------------------------- 

void FPrismCell::computeLocalCoords(const double p[3]) const
{
  if(!geometryOK)
    {
      paramsComputed=false;
      buildBoundingBox();
      geometryOK=true;
    }	
  if(!paramsComputed)
    computeParams();
  else
    if(lastPosition[0]==p[0]&&lastPosition[1]==p[1]&&lastPosition[2]==p[2])
      return;
 
  complex<double> sols[3];

  //Fd3op3(a,=p,-myPositionData[0],);
	for(int i=0; i<3; i++)
		a[i] = p[i] - myPositionData[0][i];
		
  // solution of the equation:
  // < p - pos0 + l2*pe2)  , ( (pe0+l2*pde0) , (pe1+l2*pde1) ) >   = 0
  // respective l2
  // with: pos0 = myPositionData[0], l2= localCoords[2]; 
  // the rest of the vars see computeParams

  int nsols = FMath::CubicEquation
    ( -pe2a , Fd3prod(a,pa)-pe2b, Fd3prod(a,pb)-pe2c, Fd3prod(a,pc) , sols);  
    
  if(nsols==-1)
    throw FInvalidPositionException("Cell singular");

  //find first real solution between 0 and 1
  int i  ;
  
  double mind=HUGE_VAL;
  for(i=0;i<nsols;i++){
    complex<double> & s = sols[i];
    if( fabs(s.imag()) < 2*epsilon ){
      double d=fabs(s.real()-0.5);
      if(d<mind){
	mind=d;
	localCoords[2] = s.real();
      }
    }
  }

  if( mind==HUGE_VAL ){
    ostringstream o;
    o<<"cell singular"<<endl;
    o<<"no verts:"<<geometryDescription.noVerts<<endl;
    for(unsigned int j=0;j<geometryDescription.noVerts;j++)
      o<<"vertex no "<<vertexIndices[j]<<':'<<positions[j]<<endl;
    throw FInvalidPositionException(o.str().c_str());
  }

  //.........................................................................
  //  compute other local coords

	for(int i=0; i<3; i++)
	{
		//a = vector from first vertex of triangle to searched position
		//Fd3op2(a,-= pe2,*localCoords[2]);
		a[i] -= pe2[i] * localCoords[2];
		//b = edge0 in height of l2
		//Fd3op3(b,=pe0,+pde0,*localCoords[2]);
		b[i] = pe0[i] + pde0[i] * localCoords[2];
		//c = edge1 in height of l2
		//Fd3op3(c,=pe1,+pde1,*localCoords[2]);
		c[i] = pe1[i] + pde1[i] * localCoords[2];
	}

  //d = normal vector of triangle (crossproduct of two edges)
  Fd3kreuz(d,b,c);
  double invDenom=1.0/Fd3prod(d,d);
    
  //e = normal vector to edge 1
  Fd3kreuz(e,c,d);
  localCoords[0] = Fd3prod(a,e) * invDenom;

  //e = normal vector to edge 0
  Fd3kreuz(e,d,b);
  localCoords[1] = Fd3prod(a,e) * invDenom;


  //.........................................................................
  //refine via newton iteration
  	
  static double 
    mat[3][3], //jacobian
    invMat[3][3]; //inverse,transposed jacobian

  interpolatePositionFromLocals(a);
  //Fd3op3(b,=p,-a,);
	for(int i=0; i<3; i++)
		b[i] = p[i] - a[i];
	
  double delt = Fd3prod(b,b);
  double odelt = delt;
  double eps 
    = (Fd3prod(p,p)+Fd3prod(myPositionData[0],myPositionData[0])) * 1e-28;

  //num.prec must be better than 1e-14
  if(delt > eps){
    
    dPosdLocal(mat);
    
    //compute transposed inverse
    Fd3kreuz(invMat[0],mat[1],mat[2]);
    Fd3kreuz(invMat[1],mat[2],mat[0]);
    Fd3kreuz(invMat[2],mat[0],mat[1]);	
    
    double invdenom = 1.0/Fd3prod(invMat[0],mat[0]);
    
    do{
	  
	  //localCoords+= b*invMat
	  localCoords[0]+=Fd3prod(b,invMat[0])*invdenom;
	  localCoords[1]+=Fd3prod(b,invMat[1])*invdenom;
	  localCoords[2]+=Fd3prod(b,invMat[2])*invdenom;
	  
	  //b= p-f(localCoords)
	  interpolatePositionFromLocals(a);
	  //Fd3op3(b,=p,-a,);
		for(int i=0; i<3; i++)
			b[i] = p[i] - a[i];
		
	  odelt=delt;
	  delt=Fd3prod(b,b);	  

    }while(odelt>delt && delt>eps);
  }

  //.........................................................................
  

  //Fd3op2(lastPosition,=p,);
	for(int i=0; i<3; i++)
		lastPosition[i] = p[i];
}





  
  
//--------------------------------------------------------------------------- 
  
positive FPrismCell::sizeOfCellType() const
{
  return 6;
}
 
 
//--------------------------------------------------------------------------- 
 
FCell::CellType FPrismCell::getCellType(void) const
{
  return FCell::PRISM;
}

//--------------------------------------------------------------------------- 

positive FPrismCell::memSize() const
{
  return
    (tensorData?
     tensors[0].size()*sizeof(double)*
     (sizeOfCellType()+ 5 ):0)
    +
    sizeof(*this);
}
