//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FArbitraryHexahedronCell.cc,v $
// Language:  C++
// Date:      $Date: 2004/03/11 13:53:13 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FArbitraryHexahedronCell.hh"
#include "FMatrix.hh"
#include "FPosition.hh"

#include "FException.hh"
#include "FAMSingularPoint.hh"

#include <list>
#include <vector>
#include <cassert>
#include <iostream>

//#include "Fd3op.hh"

int FArbitraryHexahedronCell::NbcomputeTril=0;
int FArbitraryHexahedronCell::Nbchecknormvect=0;

FArbitraryHexahedronCell::FArbitraryHexahedronCell()
  :trilCoords(3)
{
}
 
//--------------------------------------------------------------------------- 

FArbitraryHexahedronCell::FArbitraryHexahedronCell(const vector<FIndex>& vert,bool useVTKHexahedronEnumeration)
  :FHexahedronCell(vert,useVTKHexahedronEnumeration),trilCoords(3)
{
}

FArbitraryHexahedronCell::FArbitraryHexahedronCell(const FIndex* vert)
  :FHexahedronCell(vert),trilCoords(3)
{
}

//--------------------------------------------------------------------------- 

FArbitraryHexahedronCell::~FArbitraryHexahedronCell()
{
}

//--------------------------------------------------------------------------- 

FCell* FArbitraryHexahedronCell::getClone() const 
{
  try 
    {
      return new FArbitraryHexahedronCell( vertexIndices );
    }
  catch (FException& e) 
    {
      e.addTraceMessage("FCell* FArbitraryHexahedronCell::getClone() const");
      throw e;
    }
}

//--------------------------------------------------------------------------- 

void FArbitraryHexahedronCell::interpolate(FTensor& result, 
					   const FPosition& position) const
{
  try{

#ifndef NODEBUG
    if(!tensors||!positions)
      throw FException(" positions or tensors were not set !");
    if(position.getDimension()!=3)
      throw FInvalidDimensionException("wrong dimension of position");
#endif

    resetGeometryFlags();

    if(!pointParamsComputed)
      computePointParams();

    computeTrilinearCoords(position);

    interpolate(result);
  }

  catch(FException&e){
    e.addTraceMessage("void FArbitraryHexahedronCell::interpolate(FTensor& result, const FPosition& position) const");
    throw e;  
  }
}
 
//--------------------------------------------------------------------------- 

void FArbitraryHexahedronCell::derivatives(FTensor& result, const FPosition&p)
  const
{
  try{

#ifndef NODEBUG
    if(tensors[0].getDimension()!=3)
      throw FInvalidDimensionException();
    if(p.getDimension()!=3)
      throw FInvalidDimensionException("invalid position");
#endif

    resetGeometryFlags();

    computeTrilinearCoords(p);
  


    //compute interpolation parameters for tensors

    if(!interpolationOK)
      computeTensorParams();

    //calculate derivation of tensorfield resp. trilinear coordinates  
    FTensor derivTensor(3,tensors[0].getOrder()+1);
    computeTrilDerivTensor(derivTensor);
    
    //transform from trilinear coords to euklidian coords


    FMatrix A;
    computeTrilDeriv(A);
    A.transpose(); //because inside computeTrilderiv, matrix is also transposed
    A.invert();
    result.resizeTensor(3,tensors[0].getOrder()+1);
    FTensor m(3,tensors[0].getOrder());//multiplication buffer
    for(int i=0;i<3;i++){
      FRefTensor t(result[i]);
      derivTensor[0].mult(A(i,0),t);
      derivTensor[1].mult(A(i,1),m); t+= m;
      derivTensor[2].mult(A(i,2),m); t+= m;
    }

  }

  catch(FException&e){
    e.addTraceMessage("void FArbitraryHexahedronCell::derivatives(FTensor& result, const FPosition&p) const");
    throw e;
  }
      
}

//--------------------------------------------------------------------------- 
 
bool FArbitraryHexahedronCell::isInside(const FPosition& pos, 
					const vector< FPosition > & vertices,
					bool useVTKHexahedronEnumeration)
{
    ++FCell::ntests;

  try{

    static FArbitraryHexahedronCell cell;

    if(cell.positions[0].size()==0){	

      
      FRefArray 
	* p=cell.positions, 
	* pend = cell.positions+8;
      double*pd = cell.positionData;
      
      for(;p!=pend;p++,pd+=3)
	{
	  p->setDimension(3);
	  p->setCompPointer(pd);
	  //      cerr<<"psiz:"<<p->size()<<endl;
	}
    }



    if(!useVTKHexahedronEnumeration){
      cell.positions[0]=vertices[0];
      cell.positions[1]=vertices[1];
      //the following two are switched due to !useVTKHexahedronEnumeration numbering
      cell.positions[2]=vertices[3];
      cell.positions[3]=vertices[2];

      cell.positions[4]=vertices[4];
      cell.positions[5]=vertices[5];
      //the following two are switched due to !useVTKHexahedronEnumeration numbering
      cell.positions[6]=vertices[7];
      cell.positions[7]=vertices[6];
    }
    else{
      for (int i=0; i<8; i++)
	cell.positions[i]=vertices[i];
    }
      
    // set geometry and interpolation flag to false;
    // so that internal variables are recomputed for new positions
    cell.geometryOK = false;
    cell.interpolationOK = false;

    bool ret = cell.isInside(pos);
    //    cout << cell.trilCoords;
    return ret;
  }

  catch(FException&e){
    e.addTraceMessage("bool FArbitraryHexahedronCell::isInside(const FPosition& pos, const vector< FIndex >& verticesId, FPositionSet *posSet");
    throw e;
  }
  
}

//--------------------------------------------------------------------------- 
 
bool FArbitraryHexahedronCell::isInside(const FPosition& position) const
{
  try{

    resetGeometryFlags();

    assert( position.getDimension() == 3 );
  
    if (!bBox.isInside(position))
      return false;
  
    //calculate border planes (inner and outer ones)

    if(!normalVectorsComputed)
      computeNormalVectors();

    
    Nbchecknormvect++;

    int i;

    for(i=0;i<6;i++){

      double s=normals[i]*position;

      //if point is outside in this direction
      if(s>hih[i])
	return false;

      //if point is between outer and inner boarder plane
      if(s>loh[i])
	break;    
    }

    //in every direction point is inside
    if(i==6)
      return true;

    //point between inner and outer plane:
    //now we have to compute the trilinear coordinates
    //and check them

    //    cout<<"point inbetween:"<<endl;

    //if tril. coords should be initialized 
    //and old position is too far away from new one, do it
    trilCoords.resize(3);
    double d;
    if(!trilinearCoordsComputed
       || lastPositionEvaluated.distanceSquare(position) > ( (d=bBox.diagonal()*0.1), d*d ) ){
      
      trilCoords[0]=0.5;trilCoords[1]=0.5;trilCoords[2]=0.5;
      
      //set the coordinate normal to plane where point is inbetween
      //to 0 if the plane has tril. coord 0 or 1 if its in the plane with coord 1.
      //Ordering of planes in myGeoDescription has been optimized so
      //we can write it like this:
      trilCoords[i/2] = i%2;
    }

    //    cout<<"tril.coords set"<<endl;


    try{ 
      computeTrilinearCoords(position,i%3);       
    } 
    catch(FTooManyIterationsException){ 
      cout<<"more than 100 iterations: assuming position not inside"<<endl;
      return false;
    }


    return 
      trilCoords[0]>=-epsilon && trilCoords[0]<=1+epsilon &&
      trilCoords[1]>=-epsilon && trilCoords[1]<=1+epsilon &&
      trilCoords[2]>=-epsilon && trilCoords[2]<=1+epsilon;
  }

  catch(FException&e){
    e.addTraceMessage("bool FArbitraryHexahedronCell::isInside(const FPosition& position) const");
    throw e;
  }
}  

//--------------------------------------------------------------------------- 



void FArbitraryHexahedronCell::resetGeometryFlags() const
{
  if(!geometryOK){
    normalVectorsComputed=false;
    trilinearCoordsComputed=false;
    pointParamsComputed=false;
    
    buildBoundingBox();
    geometryOK=true;
  }
}

//--------------------------------------------------------------------------- 

FCell::CellType FArbitraryHexahedronCell::getCellType(void) const
{
  return FCell::ARBITRARY_HEX;
}



//--------------------------------------------------------------------------- 

void FArbitraryHexahedronCell::computeNormalVectors() const
{
  try{

    FVector pa(3),pb(3);

    for(int i=0;i<6;i++){

      //    if(37+i*4>=sizeof(myGeoDescription)/sizeof(int))
      //  cerr<<"wrong mygeoindex:max:"<<sizeof(myGeoDescription)/sizeof(int)<<endl;

      //cout<<"index"<<i<<":";

      //(vertexindices of surfaces begin at number  34 in myGeoDescription)
      const FPosition 
	&p0 = positions[myGeoDescription.faces[i][0]],
	&p1 = positions[myGeoDescription.faces[i][1]],
	&p2 = positions[myGeoDescription.faces[i][2]],
	&p3 = positions[myGeoDescription.faces[i][3]];
     
      //cout<<"computing ("<<p0<<'-'<<p2<<")%("<<p1<<'-'<<p3<<')'<<endl;


      //crossproduct for cross over surface
      normals[i].resize(3);
      p0.minus( p2, pa ); p1.minus( p3, pb );
      pb.crossProduct(pa,normals[i]);
      //cout<<"fertig"<<endl;

      //param. fuer ebenenglch. bestimmen
      loh[i]=p0*normals[i];
      hih[i]=p1*normals[i];

      if(loh[i]>hih[i]){
	double xch=loh[i];
	loh[i]=hih[i];
	hih[i]=xch;
      }
    
    }
  
    normalVectorsComputed=true;

  }
  catch(FException &e){
    e.addTraceMessage("void FArbitraryHexahedronCell::computeNormalVectors() const");
    throw e;
  }

}

//--------------------------------------------------------------------------- 


void FArbitraryHexahedronCell::interpolate(FTensor&t) const
{
  try{

    if(!interpolationOK)
      computeTensorParams();

    t.resizeTensor(tensors[0].getDimension(),tensors[0].getOrder());
    tensorParams.interpolateT(trilCoords,t);
  }

  catch(FException &e){
    e.addTraceMessage("void FArbitraryHexahedronCell::interpolate(FTensor&t) const");
    throw e;
  }

}

//--------------------------------------------------------------------------- 

void FArbitraryHexahedronCell::interpolate(FPosition&p) const
{
  try{

    resetGeometryFlags();

    if(!pointParamsComputed)
      computePointParams();

    p.resize(3);
    pointParams.interpolateT(trilCoords,p);

  }
  catch(FException &e){
    e.addTraceMessage("void FArbitraryHexahedronCell::interpolate(FPosition&p) const");
    throw e;
  }

}

//--------------------------------------------------------------------------- 

void FArbitraryHexahedronCell::computeTrilinearCoords(const FPosition&position,int coordToIgnore) const
{

  try{

    NbcomputeTril++;

    if(trilinearCoordsComputed && position==lastPositionEvaluated)
      return;

    //compute parameters for tril. interpolation of positions

    if(!pointParamsComputed)
      computePointParams();

    //if tril. coords should be initialized 
    //and old position is too far away from new one, do it
    double d;

    if(  
       coordToIgnore==-1 &&
       (
	! trilinearCoordsComputed ||
	lastPositionEvaluated.distanceSquare(position) > ((d=bBox.diagonal()*0.1),d*d) 
	)
       )
      {
	trilCoords[0]=0.5;trilCoords[1]=0.5;trilCoords[2]=0.5;
      }
    

    //    cout<<"start tril.coords:"<<trilCoords<<endl;
    //    cout<<"first loop"<<endl;
   
    //-------------------------------

    //     loop for raw interpolation:
    //     (it is faster in the first steps):
    //
    //     in turn for every trilinear coordinate,
    //     the derivation of positions resp. it is computed
    //     and gone in this direction to the point
    //     nearest to the point to be interpolated,
    //     which is done by solving following equation 
    //     here shown for the 1st coord:
    //
    //     < delt + dxdb1 * d , dxdb1 >   =0 
    //
    //     solving for d  gives:
    //
    //     d = - < delt ,dxdb1 > / < dxdb1, dxdb1 >
    //
    //     variables used:
    //      delt = x - p, d = b1' - b1
    //      delt' = x'-p = delt+ d* dxdb1,
    //      x,x' :old and new trilinearly interpolated position
    //      p    :position to interpolate 
    //      b1,b1':old an new 1st trilinear coordinate
                      

    FPosition x(3);
    FVector delt(3),deriv(3);
    double ldelt,ldeltend;
    int i;

    interpolate(x);

    x.minus( position ,delt );// delt=x-position
    ldelt = delt.normSquare();
    ldeltend = ldelt*1e-3; 

    double sqeps=bBox.diagonal()*epsilon;
    sqeps*=sqeps;

    if(ldeltend<sqeps*1e4)
      ldeltend=sqeps*1e4;
    
    for(i=0;ldelt>ldeltend&&i<10;i++)
      {
	 int i3=i%3;
	 if(i3!=coordToIgnore){
	   computeTrilDeriv(deriv,i3);
	   //cout<<"derivation resp. x"<<i%3<<':'<<deriv<<endl;
	   double d = -(delt*deriv)/deriv.normSquare();
	   trilCoords[i3] += d;
	   deriv *= d;
	   delt += deriv ;
	   ldelt = delt.normSquare();
	 }
       }
	 
     x=position+delt;


    //    cout<<"reached precision:"<<sqrt(ldelt)<<endl;
    //    cout<<"second(newton-raphson)- loop"<<endl;
    
    //-----------------------------

    //loop for precision via newton-raphson -method


    FMatrix M(3,3);
    FVector oTrilCoords(3);
    double odelt=ldelt;
    
    while(ldelt>sqeps){

      //if precision isn't reached with
      //the same matrix in the next step then recompute Matrix
      //( if improvement in last step is smaller than 
      //  improvement that still has to be achieved (odelt/ldelt < ldelt/sqeps)
      if(sqeps*odelt<ldelt*ldelt){
	//compute derivation of world koords after local koords
	computeTrilDeriv(M);
	//invert matrix
	M.invert();
      }
      //      cout<<"inverted"<<endl;
      //determine new point as if the derivation would be constant

      oTrilCoords=trilCoords;
      position.minus(x, delt); 
      trilCoords += M * delt;

      //      cout<<"interpolate using new trilcoords"<<endl;
  
      interpolate(x);
      odelt = ldelt;
      ldelt = x.distanceSquare(position);

      //if raw newton does not converge, use half stepsize
      while( ldelt >= odelt && i<100){
	trilCoords=(oTrilCoords+trilCoords)*0.5;
	interpolate(x);
	ldelt = x.distanceSquare(position);	
	i++;
      }
      i++;
      if(i>=100){
	cout<<"reached precision:"<<sqrt(ldelt)
	    << " trilCoords: " <<trilCoords<<" X :" << x;
	throw FTooManyIterationsException();
      }

    }


    trilinearCoordsComputed = true;
    lastPositionEvaluated   = position;

  }

  catch(FException &e){
    e.addTraceMessage("void FArbitraryHexahedronCell::computeTrilinearCoords(const FPosition&position) const");
    throw e;
  }

}


bool FArbitraryHexahedronCell::
neighborFaceForPos(const FArray& pos,FIndex&faceId ) const
{
  try{
    resetGeometryFlags();

    computeTrilinearCoords(pos);

    double b[3];
    //Fd3op2(b,=trilCoords,-0.5);
	for(int i=0; i<3; i++)
		b[i] = trilCoords[i] -0.5;
		
    double a[3]={fabs(b[0]),fabs(b[1]),fabs(b[2])};
 
    //if pos is inside, return false
    if((a[0]<0.5+epsilon)
       &(a[1]<0.5+epsilon)
       &(a[2]<0.5+epsilon)) 
      return false;
    
    positive i=0;
    if(a[1]>a[0])i=1;
    if(a[2]>a[i])i=2;
    
    //von sinn aequivalent zu:  i*2 + ( b[i]>0 ? 1 : 0 )
    faceId = (i<<1)|(b[i]>0);
    
    return true;

  }
  catch(FException &e )
    {
      cout<<e<<endl;
      faceId.setToInvalid();
      return true;
    }
  catch(exception&e)
    {
      cout<<"unknown error"<<endl;
      faceId.setToInvalid();
      return true;      
    }

}


//----------------------------------------------------------------------------


void FArbitraryHexahedronCell::computeTrilDerivTensor(FTensor&ret,int ind) const
{

  try{
    if(!interpolationOK)
      computeTensorParams();

    switch(ind)
      {
      case 0: 
	tensorParams.interpolatedTdx(trilCoords,ret);  
	break;
      case 1: 
	tensorParams.interpolatedTdy(trilCoords,ret);  
	break;
      case 2: 
	tensorParams.interpolatedTdy(trilCoords,ret);  
	break;
      default:
	throw FException("invalid koordinate Index");
      }
  }

  catch(FException &e){
    e.addTraceMessage
      ("void FArbitraryHexahedronCell::computeTrilDerivTensor(FTensor&ret,int ind) const");
    throw e;
  }


}
//----------------------------------------------------------------------------


void FArbitraryHexahedronCell::computeTrilDerivTensor(FTensor&ret) const
{

  try{

    if(!interpolationOK)
      computeTensorParams();

    
    FRefTensor a(ret[0]),b(ret[1]),c(ret[2]);
    tensorParams.interpolatedTdx(trilCoords,a);  
    tensorParams.interpolatedTdy(trilCoords,b);  
    tensorParams.interpolatedTdz(trilCoords,c);  
  }
  catch(FException &e){
    e.addTraceMessage
      ("void FArbitraryHexahedronCell::computeTrilDerivTensor(FTensor&ret) const");
    throw e;
  }  
}

//----------------------------------------------------------------------------
 

void FArbitraryHexahedronCell::computeTrilDeriv(FVector&ret,int ind) const
{
  try{

    if(!pointParamsComputed)
      computePointParams();

    switch(ind)
      {
      case 0: 
	pointParams.interpolatedTdx(trilCoords,ret);  
	break;
      case 1: 
	pointParams.interpolatedTdy(trilCoords,ret);  
	break;
      case 2: 
	pointParams.interpolatedTdz(trilCoords,ret);  
	break;
      default:
	throw FException("invalid koordinate Index");
      }
  }

  catch(FException &e){
    e.addTraceMessage
      ("void FArbitraryHexahedronCell::computeTrilDeriv(FVector&ret,int ind) const");
    throw e;
  }

}

//----------------------------------------------------------------------------
 

void FArbitraryHexahedronCell::computeTrilDeriv(FMatrix&ret) const
{
  try{

    if(!pointParamsComputed)
      computePointParams();

    FTensor t(3,2);
    FRefTensor a(t[0]),b(t[1]),c(t[2]);
    pointParams.interpolatedTdx(trilCoords,a);  
    pointParams.interpolatedTdy(trilCoords,b);  
    pointParams.interpolatedTdz(trilCoords,c);  
    ret=t;
    
  }
  catch(FException &e){
    e.addTraceMessage
      ("void FArbitraryHexahedronCell::computeTrilDeriv(FMatrix&ret) const");
    throw e;
  }  

}

//----------------------------------------------------------------------------

void FArbitraryHexahedronCell::computeTensorParams() const
{
  try{

    int tsiz=tensors[0].size(),tsiz8=tsiz*8;

    vector<double> tensorP(tsiz8);

    vector<double>::iterator tp=tensorP.begin();

    double *td = tensorData,*tdend=tensorData+tsiz;
    for(;td!=tdend;td++){
      double*td2=td,*td2end=td+tsiz8;
      for(;td2!=td2end;tp++,td2+=tsiz)
	*tp=*td2;
    }

    
    tensorParams.buildParameters(tensorP);

    interpolationOK = true;
  }
  
  catch(FException &e){
    e.addTraceMessage
      ("void FArbitraryHexahedronCell::computeTensorParams() const");
    throw e;
  }  

}

//----------------------------------------------------------------------------

void FArbitraryHexahedronCell::computePointParams() const
{

  try{
    vector<double> pointP(3*8);
    vector<double>::iterator pp=pointP.begin();
    //fill parameters
    double *pd=positionData,*pdend=pd+3;
    for(;pd!=pdend;pd++){
      double*pd2=pd,*pd2end=pd+3*8;
      for(;pd2!=pd2end;pp++,pd2+=3)
	*pp=*pd2;
    }

    //build parameters
    pointParams.buildParameters(pointP);

    pointParamsComputed=true;

  }
  catch(FException &e){
    e.addTraceMessage
      ("void FArbitraryHexahedronCell::computePointParams() const");
    throw e;
  }  
}

//----------------------------------------------------------------------------

void FArbitraryHexahedronCell::
computeSing(float minb[3],float maxb[3],list<FAMSingularPoint>& result) const
{
  try
  {
    if(!interpolationOK)
    {
      computePointParams();
      computeTensorParams();      
    }
    
    
    float deltmb[3];
    for(int i=0; i<3; i++)
      deltmb[i] = maxb[i] - minb[i];
    
  FArray lx(3); //local coords;
  for(int i=0; i<3; i++)
    lx[i] = minb[i] + maxb[i];

  lx *= 0.5;

  

  FMatrix M(3,3);
  FTensor t(3,2);
  FArray x(3);

  tensorParams.interpolateT(lx,x);
  // cout<<"local tensor"<<x<<endl;

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
      //compute derivation of world coords after local koords
      
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
      pointParams.interpolateT(lx,x);
      FAMSingularPoint ap(x,FAMSingularPoint::NONE,1);

      FMatrix tens(3,3);
      FRefTensor 
	a(3,1,&tens(0,0)),
	b(3,1,&tens(1,0)),
	c(3,1,&tens(2,0));

      tensorParams.interpolatedTdx(lx,a);  
      tensorParams.interpolatedTdy(lx,b);  
      tensorParams.interpolatedTdz(lx,c);  

      FMatrix pDeriv(3,3);

      a.setPointerOnComponents(&pDeriv(0,0));
      b.setPointerOnComponents(&pDeriv(1,0));
      c.setPointerOnComponents(&pDeriv(2,0));

      pointParams.interpolatedTdx(lx,a);  
      pointParams.interpolatedTdy(lx,b);  
      pointParams.interpolatedTdz(lx,c);  

      pDeriv.invert();

      ap.setLinearNature(tens*pDeriv);

      result.push_back(ap);

    }  
  }catch(FException&e)
    {
      e.addTraceMessage("FArbitraryHexahedronCell::computeSing");
      throw;
    }
}

//----------------------------------------------------------------------------
positive FArbitraryHexahedronCell::memSize() const
{
  return
    (
    6//normals
    +2//trilCoords,lastPositionEvaluated
    +8//pointParams
    ) *3*sizeof(double) 
    +
    (tensorData?
     tensors[0].size()*sizeof(double)*
     (sizeOfCellType()*2 //tensorData+tensorParams
      ):0)
    +
    sizeof(*this);
}
