//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FQuadrilateralCell3D.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:11 $
// Author:    $Author: garth $
// Version:   $Revision: 1.20 $
//
//--------------------------------------------------------------------------- 

#include "FQuadrilateralCell3D.hh"
#include "FException.hh"
#include "FMatrix.hh"



//---------------------------------------------------------------------------
// the cell's geometry description
const FCell::geoInfo FQuadrilateralCell3D::myGeoDescription =
  {
	// dimension
	3,
	// # vertices
	4, 
	// # edges
	4, 
	// # faces
	0, 
	// edges
	{{0,1},  {1,2},  {2,3},  {3,0}},
	// face sizes
	{},
	// faces
	{}
  };

//--------------------------------------------------------------------------- 

FQuadrilateralCell3D::FQuadrilateralCell3D(void) 
  : FCell2Din3D(myGeoDescription,myIndices,myPositions,myPositionData[0],myTensors)
{
  lastPosition[0]=HUGE_VAL;
}

//---------------------------------------------------------------------------

FQuadrilateralCell3D::FQuadrilateralCell3D( const vector< FIndex >& vert1 ) 
  : FCell2Din3D(myGeoDescription,myIndices,myPositions,myPositionData[0],myTensors)
{
  for (char i=0; i<4; i++)
    vertexIndices[i]= vert1[i];

  lastPosition[0]=HUGE_VAL;
}

FQuadrilateralCell3D::FQuadrilateralCell3D( const FIndex*vert1 ) 
  : FCell2Din3D(myGeoDescription,myIndices,myPositions,myPositionData[0],myTensors)
{
  for (char i=0; i<4; i++)
    vertexIndices[i]= vert1[i];

  lastPosition[0]=HUGE_VAL;
}

//---------------------------------------------------------------------------

FQuadrilateralCell3D::~FQuadrilateralCell3D()
{
}

//---------------------------------------------------------------------------

FCell*  FQuadrilateralCell3D::getClone() const
{
  try 
    {
      return new FQuadrilateralCell3D( vertexIndices );
    }
  catch (FException& e) 
    {
      e.addTraceMessage("FCell*  FQuadrilateralCell3D::getClone() const");
      throw;
      return (FCell *) 0;
    }
}

//---------------------------------------------------------------------------

positive FQuadrilateralCell3D::sizeOfCellType() const
{
  return 4;
}

//--------------------------------------------------------------------------- 

void FQuadrilateralCell3D::buildBoundingBox(void) const
{
  try{
#ifndef NODEBUG
    if(!positions[0].size())
      throw FException("positions not set");
#endif
    if (geometryOK)
      return;

    FPosition & p0 = positions[0];

    bBox.setBoundingBox( p0[0],p0[1],p0[2], p0[0],p0[1],p0[2] );
  
    for ( positive i = 1; i<getNumVerts(); i++ )
      bBox.resize( positions[i]);

    geometryOK = true;
  }

  catch(FException&e){
    e.addTraceMessage("void FQuadrilateralCell3D::buildBoundingBox(void) const");
    throw e;
  }
}

//---------------------------------------------------------------------------


inline void 
FQuadrilateralCell3D::computeLocalCoords(const double p[3]) const
{

  if((lastPosition[0]==p[0])&(lastPosition[1]==p[1])&(lastPosition[2]==p[2]))
    return;

  if(!geometryOK){
    surf.init(myPositionData[0],myPositionData[1],
	      myPositionData[2],myPositionData[3]);
    surf.normal(normal,0.5,0.5);
    buildBoundingBox();
  }

  //Fd3op2(lastPosition,=p,);
	for(int i=0; i<3; i++)
		lastPosition[i] = p[i];
  
  double ret[2][3]; //return values
  
  int n = surf.cutWithLine(lastPosition,normal,ret);

  if(n==0){
    //Fd3op1(localCoords,=0);
	for(int i=0; i<3; i++)
		localCoords[i] = 0;
		
    return;
  }

  double * l = ret[0];

  if(n==2)
    if( fabs(ret[0][0]-0.5)+fabs(ret[0][1]-0.5)
	>= fabs(ret[1][0]-0.5)+fabs(ret[1][1]-0.5) )
      l=ret[1];

  if(fabs(l[2]*epsilon)>1){
    //Fd3op1(localCoords,=0);
	for(int i=0; i<3; i++)
		localCoords[i] = 0;
		
    return;
  }

  //Fd3op2(localCoords,=l,);
  for(int i=0; i<3; i++)
	localCoords[i] = l[i];
  
}



//--------------------------------------------------------------------------- 

void 
FQuadrilateralCell3D::interpolate(FTensor& result, 
				  const FPosition& position) const
{
  

  computeLocalCoords(&position[0]);

  double 
    m0 = 1-localCoords[0],
    m1 = 1-localCoords[1];
  
  result 
    = (m0            * m1            ) * myTensors[0] 
    + (localCoords[0]* m1            ) * myTensors[1]
    + (localCoords[0]* localCoords[1]) * myTensors[2]
    + (m0            * localCoords[1]) * myTensors[3] ;
  
}
  
//--------------------------------------------------------------------------- 

inline void 
FQuadrilateralCell3D::
dPosdLocal(double3 m[]) const
{
  //aliases
  const double3*const p = myPositionData;
  const double*const l = localCoords;

  double 
    m0 = 1-l[0],
    m1 = 1-l[1];

  //Fd3op5(m[0], = -m1 * p[0],   +m1 * p[1], +l[1] * p[2], -l[1] * p[3],);
  //Fd3op5(m[1], = -m0 * p[0], -l[0] * p[1], +l[0] * p[2],   +m0 * p[3],);
  //Fd3op2(m[2], = normal,);
  for(int i=0; i<3; i++)
  {
	m[0][i] = -m1 * p[0][i]   +m1 * p[1][i] +l[1] * p[2][i] -l[1] * p[3][i];
	m[1][i] = -m0 * p[0][i] -l[0] * p[1][i] +l[0] * p[2][i]   +m0 * p[3][i];
	m[2][i] = normal[i];
  }
}

//--------------------------------------------------------------------------- 

inline void 
FQuadrilateralCell3D::
dTensdLocal(FTensor&result) const
{
  result.resizeTensor(myTensors[0].getDimension(),myTensors[0].getOrder()+1);

  //aliases
  const FTensor*const t = myTensors;
  const double*const l = localCoords;
  
  double 
    m0 = 1-l[0],
    m1 = 1-l[1];

  result[0] = -m1 * t[0]   +m1 * t[1] +l[1] * t[2] -l[1] * t[3];
  result[1] = -m0 * t[0] -l[0] * t[1] +l[0] * t[2]   +m0 * t[3];
  //set components to 0
  result[2].clear();
}
//--------------------------------------------------------------------------- 

void FQuadrilateralCell3D::derivatives(FTensor& result, 
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
      ("void FQuadrilateralCell3D::derivatives"
       "(FTensor& result, const FPosition&p) const");
    throw;
  }
}

//--------------------------------------------------------------------------- 


bool 
FQuadrilateralCell3D::
isInside(const FPosition& position) const
{

  computeLocalCoords(&position[0]);

  return 
    (-epsilon < localCoords[0]) & (localCoords[0] < 1+epsilon) &
    (-epsilon < localCoords[1]) & (localCoords[1] < 1+epsilon) &
    (-epsilon < localCoords[2]) & (localCoords[2] < 1+epsilon) ;
}


// -------------------------------------------------------------------------

bool FQuadrilateralCell3D::isInside( const FPosition& /*pos*/,
 				     const vector< FPosition >& /*vertices*/)
{

  THROW_DEFAULT_EXCEPTION(FNotImplementedException);
  return false;
}


// ------------------------------------------------------------------------

void FQuadrilateralCell3D::getZeros(list<FAMSingularPoint>& /*result*/) const
{ 
  THROW_DEFAULT_EXCEPTION( FNotImplementedException ); 
}

// ------------------------------------------------------------------------

FCell::CellType FQuadrilateralCell3D::getCellType(void) const
{
  return FCell::QUADRILATERAL_3D;
}

// ------------------------------------------------------------------------

positive FQuadrilateralCell3D::memSize() const
{
  return
    (tensorData?
     tensors[0].size()*sizeof(double)*
     sizeOfCellType() //tensorData
     :0)
    +
    sizeof(*this);
}
