//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCell.cc,v $
// Language:  C++
// Date:      $Date: 2004/03/11 14:32:50 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.11 $
//
//--------------------------------------------------------------------------- 

#include <cstddef>
#include <cassert>
#include <iostream>

#include "stdAliases.hh"
#include "FException.hh"
#include "FIndex.hh"
#include "FCell.hh"


// cell types: must be all included here to enable the use of their
// respective static isInside function
#include "FAxisParallelTriangleCell2D.hh"
#include "FAxisParallelQuadCell2D.hh"
#include "FAxisParallelHexCell.hh"
#include "FTriangleCell2D.hh"
#include "FQuadrilateralCell2D.hh"
#include "FTriangleCell3D.hh"
#include "FQuadrilateralCell3D.hh"
#include "FTetrahedronCell.hh"
#include "FArbitraryHexahedronCell.hh"
#include "FPrismCell.hh"
#include "FPyramidCell.hh"

//for FLocatorInfo
#include "FBilinearSurface.hh"

double FCell::zero_threshold=0.0;
double FCell::epsilon =1.0E-9;

//--------------------------------------------------------------------------- 

  
class FCellFaceInfo{
  

  FCellFaceInfo(const FCell*aCell);

  FIndex faceToGo(double&length,const FArray& start, const FArray& dir) const;


  struct face{
 
    double normal[3]; //should point to outside
    double normalSqlen;// < normal,normal >

    bool isFlat;
    FBilinearSurface bilinear;
  };

  face faces[FCell::maxNoFaces];

  const FCell::geoInfo&geo;
  
  typedef double double3[3];
  const double3*const pData;
    
  bool inverted; //if cell's faces point to the opposite direction 
  //as expected (not to outside)    

  friend class FCell;
}; 


//--------------------------------------------------------------------------- 

FCell::FCell(const geoInfo&g,FIndex*const ind,FRefArray*const pos,double *const posData,FRefTensor*const tens)
  : geometryDescription(g),
    vertexIndices(ind), 
    tensors(tens),tensorData(0),
    positions(pos), positionData(posData),
    geometryOK(false), interpolationOK(false),  
    bBox(), set_to_zero(false),faceInfo(0)
{
}

//--------------------------------------------------------------------------- 

FCell::~FCell()
{
  if(faceInfo)
    delete faceInfo;

  if (tensorData)
    delete[] tensorData;
}

//--------------------------------------------------------------------------- 

FPosition FCell::getCenterOfBoundingBox() const
{
    if (!geometryOK) 
	buildBoundingBox();

    return bBox.center();
}

//--------------------------------------------------------------------------- 

bool FCell::intersects(const FBoundingBox& inBox)
{
  if (!geometryOK) {
    buildBoundingBox();
  }

  return bBox.intersects( inBox );
}

//--------------------------------------------------------------------------- 

void FCell::setPositions(const FPositionSet* posSetPt)
{
  try{
#ifndef NODEBUG
    if(posSetPt->getDimension()!=getDimension())
      throw FInvalidDimensionException
	("positionset does not have the correct dimension");
#endif

    int size = sizeOfCellType();
    int dim = getDimension();


    if(positions[0].size()==0){

      //attention: don't put this into the constructor of FCell,
      //because the constructors of the positions
      //are invoked in the subclass after the constructor of Fcell

      FRefArray * p=positions, * pend = positions+sizeOfCellType();
      double*pd = positionData;
      for(;p!=pend;p++,pd+=dim)
	{
	  p->setDimension(dim);
	  p->setCompPointer(pd);
	  //      cerr<<"psiz:"<<p->size()<<endl;
	}
    }

    FPosition p(dim);
    
    for (int i=0; i<size; i++){
      posSetPt->getPosition(p, vertexIndices[i]);
      positions[i]=p;
    }
    
    // set geometry and interpolation flag to false;
    geometryOK = false;
    interpolationOK = false;
    if(faceInfo){delete faceInfo;faceInfo=0;}
  }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FCell::setTensors(const FTensorSet* tensSetPt)
{
  unsigned size = sizeOfCellType();
	  
  unsigned dim = tensSetPt->getDimension();
  unsigned order = tensSetPt->getOrder();

  if(order!=tensors[0].getOrder()||dim!=tensors[0].getDimension())
    {
      //resize tensors array
      if(tensorData)delete[] tensorData;
      
      int tensSiz=FTensor::pow(dim,order);
      tensorData = new double[tensSiz*size];
      
      FRefTensor *t=tensors, *tend = tensors+size;
      double * td=tensorData;

      for(;t!=tend;t++,td+=tensSiz)
	t->setDimensionOrderPointerAndSize(dim,order,td,tensSiz);            
    }
	  
  for (unsigned int i=0; i<size; i++)
    tensSetPt->getTensor(tensors[i], vertexIndices[i]);

  double norm, max=0.;
  for (unsigned int i=0; i<size; i++) {
    // check if the norm is big enough otherwise set to zero
    norm = tensors[i].norm();
    if (norm > max)
      max = norm;
  }
  if (max < getZero_Threshold() ) 
    set_to_zero = true;

  // (possible) interpolation parameters must be recomputed
  interpolationOK = false;
}

//--------------------------------------------------------------------------- 

void FCell::buildBoundingBox(void) const
{
  try{
    if(!positions)
      throw FException("positions not set");

    FPosition & p0 = positions[0];

    if( getDimension()==2)
      bBox.setBoundingBox( p0[0],p0[1], p0[0],p0[1] );
    else
      bBox.setBoundingBox( p0[0],p0[1],p0[2], p0[0],p0[1],p0[2] );
  
    for ( positive i = 1; i<getNumVerts(); i++ )
      bBox.resize( positions[i]);
  }
CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------


FIndex& FCell::operator() (positive vertex)
{
  if (vertex >= geometryDescription.noVerts) {
    THROW_EXCEPTION(FIndexOutOfBoundsException, "Tried to access non existing vertex index in a cell");
  }
  return vertexIndices[vertex];
}

//--------------------------------------------------------------------------- 

const FIndex& FCell::operator() (positive vertex) const
{
  if (vertex >= geometryDescription.noVerts) {
    THROW_EXCEPTION(FIndexOutOfBoundsException, "Tried to access non existing vertex index in a cell");
  }
  return vertexIndices[vertex];
}

//--------------------------------------------------------------------------- 

void FCell::setPrecision( const double& newepsilon )
{
  epsilon = newepsilon;
}

//--------------------------------------------------------------------------- 

const double& FCell::getPrecision()
{
  return epsilon;
}

//--------------------------------------------------------------------------- 

void FCell::setZero_Threshold( const double& value )
{
  zero_threshold = value;
}

//--------------------------------------------------------------------------- 

const double& FCell::getZero_Threshold() const
{
  return zero_threshold;
}

//--------------------------------------------------------------------------- 

void FCell::changeSetToZeroFlag() const
{
  set_to_zero = !set_to_zero;
}

//--------------------------------------------------------------------------- 

ostream& operator<<(ostream& os, FCell::CellType myCellType)
{
  static const char *tmptype[] = {    
    "TRIANGLE_2D",
    "AXIS_PARALLEL_TRI_2D",
    "QUADRILATERAL_2D",
    "AXIS_PARALLEL_QUAD_2D",
    "TRIANGLE_3D", 
    "QUADRILATERAL_3D",
    "TETRAHEDRON",
    "AXIS_PARALLEL_TET",
    "ARBITRARY_HEX", 
    "AXIS_PARALLEL_HEX",
    "PRISM", 
    "PYRAM",
    "LINE_2D",
    "LINE_3D",
    "POINT_2D",
    "POINT_3D",
    "UNDEFINED"
  };

  os << tmptype[myCellType];
  return os;
}


FCell* FCell::getCellTorso(CellType t,const std::vector<FIndex> & vertices )
{
  switch(t)
    {
    case AXIS_PARALLEL_TET:
    case TETRAHEDRON:
      return new FTetrahedronCell(vertices);

    case AXIS_PARALLEL_TRI_2D:
      return new FAxisParallelTriangleCell2D(vertices);

    case TRIANGLE_2D: 
      return new FTriangleCell2D(vertices);

    case AXIS_PARALLEL_QUAD_2D:
      return new FAxisParallelQuadCell2D(vertices);

    case QUADRILATERAL_2D: 
      return new FQuadrilateralCell2D(vertices);

    case TRIANGLE_3D:
	return new FTriangleCell3D(vertices);

    case QUADRILATERAL_3D: 
      return new FQuadrilateralCell3D(vertices);

    case AXIS_PARALLEL_HEX: 
      return new FAxisParallelHexCell(vertices);
      
    case ARBITRARY_HEX:
      return new FArbitraryHexahedronCell(vertices);

    case PRISM:
      return new FPrismCell(vertices);

    case PYRAM:
      return new FPyramidCell(vertices);

    default:
    {      
      THROW_EXCEPTION (FNotImplementedException, "unknown cell type!" );
    }
  }
  //to suppress warnings:
  return 0;
}



const FCell::geoInfo & FCell::getGeometryDescription(FCell::CellType type)
{
  switch(type)
    {
    case AXIS_PARALLEL_TET:
    case TETRAHEDRON:
      return FTetrahedronCell::myGeoDescription;

    case AXIS_PARALLEL_TRI_2D:
    case TRIANGLE_2D: 
      return FTriangleCell2D::myGeoDescription;

    case AXIS_PARALLEL_QUAD_2D:
    case QUADRILATERAL_2D: 
      return FQuadrilateralCell2D::myGeoDescription;

    case TRIANGLE_3D:
      return FTriangleCell3D::myGeoDescription;

    case QUADRILATERAL_3D: 
      return FQuadrilateralCell3D::myGeoDescription;

    case AXIS_PARALLEL_HEX: 
    case ARBITRARY_HEX:
      return FHexahedronCell::myGeoDescription;

    case PRISM:
      return FPrismCell::myGeoDescription;

    case PYRAM:
      return FPyramidCell::myGeoDescription;          

    default:
    {      
      THROW_EXCEPTION(FNotImplementedException, "unknown cell type!" );
    }
    }
}
//--------------------------------------------------------------------------- 

void FCell::getBoundingBox(FBoundingBox& result) const
{
  if (!geometryOK) {
    buildBoundingBox();
  }

  result = bBox;
}

//--------------------------------------------------------------------------- 

const FBoundingBox& FCell::getBoundingBox(void) const
{
  if (!geometryOK) {
    buildBoundingBox();
  }

  return bBox;
}

//--------------------------------------------------------------------------- 

FBoundingBox FCell::getBoundingBox(const std::vector< FIndex >& verticesId,
				   FPositionSet *posSet)
{
  try{

    positive dim = posSet->getDimension();
    positive nbPos = verticesId.size();
    vector< FPosition > pos(nbPos);
    FBoundingBox bBox;

    for (positive i=0 ; i<nbPos ; i++)
      posSet->getPosition(pos[i], i);

    if (dim == 2)
      bBox.setBoundingBox(pos[0][0], pos[0][1], 
			  pos[0][0], pos[0][1]);
    else
      bBox.setBoundingBox(pos[0][0], pos[0][1], pos[0][2], 
			  pos[0][0], pos[0][1], pos[0][2] );
  
    for (positive i=1 ; i<nbPos ; i++)
      bBox.resize(pos[i]);

    return bBox;
  }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

FBoundingBox FCell::getBoundingBox(const std::vector< FPosition >& vertices )
{
  try{
    positive dim = vertices[0].getDimension();
    positive nbPos = vertices.size();
    FBoundingBox bBox;

    if (dim == 2)
      bBox.setBoundingBox(vertices[0][0], vertices[0][1], 
			  vertices[0][0], vertices[0][1]);
    else
      bBox.setBoundingBox(vertices[0][0], vertices[0][1], vertices[0][2], 
			  vertices[0][0], vertices[0][1], vertices[0][2] );
  
    for (positive i=1 ; i<nbPos ; i++)
      bBox.resize(vertices[i]);

    return bBox;

  }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FCell::getVertexIndices(std::vector<FIndex>& result) const
{
  int siz=sizeOfCellType();
  result.resize(siz);
  for(int i=0;i!=siz;i++)
    result[i] = vertexIndices[i];
}

//--------------------------------------------------------------------------- 

void FCell::getVertices(std::vector<FPosition>& result) const
{
  if (!positions) {
    THROW_DEFAULT_EXCEPTION( FEmptyObjectException );
  }

  int siz=sizeOfCellType();
  result.resize( siz);
  
  for( int i = 0; i !=siz; i++ )
    result[i] = positions[i];
}

//--------------------------------------------------------------------------- 

void FCell::getTensor( FTensor& t, unsigned int i ) const
{
  t = tensors[ i ];
}

//--------------------------------------------------------------------------- 

void FCell::getTensors(std::vector<FTensor>& result) const
{
  try{
#ifndef NODEBUG
    if (!tensors[0].size()) {
      THROW_DEFAULT_EXCEPTION(FEmptyObjectException);
    }
#endif

    int siz=sizeOfCellType();
    result.resize( siz );

    for( int i = 0; i != siz; i++ )
      result[i] = tensors[i];
  }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 
 
void FCell::getEdges(std::vector< pair<FIndex, FIndex> >& result) const
{
  try{
    result.resize( geometryDescription.noEdges);
    
    const positive* e = &geometryDescription.edges[0][0];
    int n=geometryDescription.noEdges;
    
    for (int i=0; i<n ; i++)
      {
	result[i].first = vertexIndices[*(e++)];
	result[i].second = vertexIndices[*(e++)];
      }
  }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FCell::getEdges(std::vector< pair<FPosition, FPosition> >& result) const
{
  try{
#ifndef NODEBUG
    if (!positions[0].size()) {
      throw FEmptyObjectException();
    }
#endif
    result.resize( geometryDescription.noEdges);
    
    const positive* e = &geometryDescription.edges[0][0];
    int n=geometryDescription.noEdges;
    for (int i=0; i< n; i++)
      {
	result[i].first = positions[*(e++)];
	result[i].second = positions[*(e++)];
      }
  }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FCell::getEdges(std::vector< pair<FTensor, FTensor> >& result) const
{
  try{
#ifndef NODEBUG
    if (!tensorData) {
      throw FEmptyObjectException();
    }
#endif

    result.resize( geometryDescription.noEdges);

    const positive* e = &geometryDescription.edges[0][0];
    int n=geometryDescription.noEdges;

    for (int i=0; i<n ; i++)
      {
	result[i].first = tensors[*(e++)];
	result[i].second = tensors[*(e++)];
      }
  }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FCell::getFaces(std::vector< std::vector<FIndex> >& result) const
{

  result.resize( geometryDescription.noFaces );
  
  for (positive i=0; i< geometryDescription.noFaces; i++)
    {
      vector<FIndex> & resultI = result[i];

      positive faceSize = geometryDescription.faceSizes[i];

      resultI.resize(faceSize);

      const positive* faces = &(geometryDescription.faces[i][0]);

      for (positive j=0; j!=faceSize; j++)
	resultI[j] = vertexIndices[faces[j]];
    }
}

//--------------------------------------------------------------------------- 

void FCell::getFaces(std::vector< std::vector<FPosition> >& result) const
{
  if (!positions) {
    THROW_DEFAULT_EXCEPTION(FEmptyObjectException);
  }

  result.resize( geometryDescription.noFaces );
  
  for (positive i=0; i< geometryDescription.noFaces; i++)
    {
      vector<FArray> & resultI = result[i];

      positive faceSize = geometryDescription.faceSizes[i];

      resultI.resize(faceSize);

      const positive* faces = &geometryDescription.faces[i][0];

      for (positive j=0; j!=faceSize; j++)
	resultI[j] = positions[faces[j]];
    }
}  

//--------------------------------------------------------------------------- 

void FCell::getFaces(std::vector< std::vector<FTensor> >& result) const
{
  if (!tensors) {
    THROW_DEFAULT_EXCEPTION( FEmptyObjectException);
  }

  result.resize( geometryDescription.noFaces );
  
  for (positive i=0; i< geometryDescription.noFaces; i++)
    {
      vector<FTensor> & resultI = result[i];

      positive faceSize = geometryDescription.faceSizes[i];

      resultI.resize(faceSize);

      const positive* faces = &geometryDescription.faces[i][0];

      for (positive j=0; j!=faceSize; j++)
	resultI[j] = tensors[faces[j]];
    }

}  


//--------------------------------------------------------------------------- 

const double* FCell::getPositionData() const
{
  return positionData;
}

//--------------------------------------------------------------------------- 

const double* FCell::getTensorData() const
{
  return tensorData;
}

//--------------------------------------------------------------------------- 

positive FCell::getDimension() const
{
  return geometryDescription.dim;
}
//--------------------------------------------------------------------------- 

positive FCell::getNumVerts() const
{
  return geometryDescription.noVerts;
}

//--------------------------------------------------------------------------- 
positive FCell::getNumEdges() const
{
  return geometryDescription.noEdges;
}

//--------------------------------------------------------------------------- 
positive FCell::getNumFaces() const
{
  return geometryDescription.noFaces;
}

//--------------------------------------------------------------------------- 

positive FCell::sizeOfCellType() const
{
  return geometryDescription.noVerts;
}

//--------------------------------------------------------------------------- 

FPosition FCell::getBarycenter() const 
{
  try 
    {
      if (!positions) 
	throw FException("positions not set");

      FPosition pos = positions[0];
      for ( positive i = 1; i<getNumVerts(); i++ )
	pos += positions[i];

      return 1. / (float) getNumVerts() * pos;
    }
  CATCH_N_RETHROW( FException );
}





//--------------------------------------------------------------------------- 




FIndex FCell::faceToGo(double&length,const FArray& start, const FArray& dir) const
{

  assert(getDimension()==3);
    
  if(!faceInfo)faceInfo = new FCellFaceInfo(this);

  return faceInfo->faceToGo(length,start,dir);
 
}

#define d2inside(d,minval,maxval) \
  (d[0] > minval) & (d[0]< maxval) \
& (d[1] > minval) & (d[1]< maxval)

#define setnv(a) setNewVertex(a[0],a[1],a[2]) 

FIndex FCellFaceInfo::faceToGo(double&length,const FArray& start, const FArray& dir) const
{
    ++FCell::ntests;
  //  cout<<"faceToGo"<<endl;

  static const double eps = 1e-12;
  double bigeps=1e-6,smalleps=eps;
  
  double3 a,s,d;
  double3 ret[2];
  //Fd3op2(s,=start,);
  //Fd3op2(d,=dir,);
	for(int i=0; i<3; i++)
	{
		s[i] = start[i];
		d[i] = dir[i];
	}
	
  double dsqlen=Fd3prod(d,d);
    
  length = HUGE_VAL;    
  int reti=-1;
  
  do{

    for(unsigned int i=0;i<geo.noFaces;i++){
      const face & f =faces[i];
      
      if(!f.isFlat){ //if face is bilinear
	//	cout<<"bilin."<<endl;
	int num = f.bilinear.cutWithLine(s,d,ret);
	
	if(num==2)//sort by distance
	  if(ret[0][2]>ret[1][2])
	    {
			//Fd3op2(a,=ret[0],);
			//Fd3op2(ret[0],=ret[1],);
			//Fd3op2(ret[1],=a,);
			for(int i=0; i<3; i++)
			{
				a[i] = ret[0][i];
				ret[0][i] = ret[1][i];
				ret[1][i] = a[i];
			}
		}
	
	double dlen[2];
	
	int numToLook=0,iToLook=-1;
	for(int j=0;j<num;j++){
	  
	  double *r = ret[j];

	  if( d2inside(r,-bigeps,1+bigeps) ){
	    
	    f.bilinear.normal(a,r[0],r[1]);
	    dlen[j] = Fd3prod(a,d);
	    //if d is not almost parallel to face
	    //( <a,d>  >  smalleps*(||a||*||d||)
	    if( dlen[j]*dlen[j] > smalleps*smalleps * dsqlen * Fd3prod(a,a) ) {
	      
	      iToLook=j;
	      numToLook++;
	    
	      //if cutpoint is far inside face 
	      if( (dlen[j]>0) ^ inverted  
		  && d2inside(r,bigeps,1-bigeps) & (r[2]> bigeps) & (r[2]<1-bigeps)){
	      
		//		cout<<"direct!"<<endl;

		length = r[2];
		return i;
	      }
	    }
	    //	    else  cout<<"almost parallel"<<dlen[j]<<endl;
	  }
	  
	}  
	//only one senseful solution
	if(numToLook==1) {
	  //	  cout<<"dlen:"<<dlen[iToLook]<<" inv "
	  //<<inverted<<" newlen "<<ret[iToLook][2]<< "length "<<length<<endl;
	  // if ray exits through this face ("^"= XOR)
	  if((dlen[iToLook]>0) ^ inverted) 
	    //and if distance to this face is smaller than to the others
	    if(fabs(length)>fabs(ret[iToLook][2]))
	      {reti=i;length=ret[iToLook][2];}	  	    
	}

	//both solutions make sense
	if(numToLook==2){

	  //if the first cut is exit point
	  if( (dlen[0]>0) ^ inverted ){
	    //and if cuts go after actual point
	    if( (ret[0][2]+ret[1][2]) > 0 )
	      //and if distance to this face is smaller than to the others
	      if( fabs(length)>fabs(ret[0][2]) )
		//then set  cut face to actual face
		{ reti=i; length=ret[0][2]; }	  	    
	  }
	  //if 2nd cut is exit point
	  else
	    //if distance to this face is smaller than to the others
	    if(fabs(length)>fabs(ret[0][2]))
	      //then set  cut face to actual face
	      {reti=i;length=ret[1][2];}	  	    
	}
      
	
      }
      else{
	//	cout<< "flat "<<endl;
	double dlen=Fd3prod(f.normal,d);
      
	//if d is not almost parallel to face
	if(dlen*dlen > smalleps*smalleps * dsqlen * f.normalSqlen)
	
	  //if line goes from inside to outside  here
	  if( (dlen>0) ^ inverted ){ // "^" means xor
	    const double * f0 = pData[geo.faces[i][0]];
	  
	    //Fd3op3(a,=f0,-s,);
		for(int k=0; k<3; k++)
			a[k] = f0[k] - s[k];
	  
	    dlen = Fd3prod(f.normal,a)/dlen;
	  
	    if(fabs(length)>fabs(dlen))
	      {reti=i;length=dlen;}	  
	  }
	//	  else cout<<"flat almost parallel"<<endl;
      
      }
    }
    
    bigeps*=10;
    smalleps*=0.1;
    if(reti==-1)
      cout<<"no face found, retrying with other epsilons"<<smalleps<<endl;

  }while(reti==-1&&smalleps>1e-14);

  if(reti==-1){
    cout<<"cell with too bad proportions, exiting"<<endl;
    cout<<"start "<<start<<" dir"<<dir<<endl;
  }
    
  if(length<0)length=0;

  return reti;    
  
}

FCellFaceInfo::FCellFaceInfo(const FCell*cell)
  :
  geo(cell->geometryDescription),
  pData((const double3*)cell->getPositionData())
{
  
  double3 a,b,center;

  double lenForInv=0,d;

  //Fd3op1(center,=0);
  for(int i=0; i<3; i++)
	center[i] = 0;

  const double3*pIt=pData;
  for(unsigned int i=0;i!=geo.noVerts;i++,pIt++)
    {      
      //Fd3op2(center,+= (*pIt),);
	   for(int i=0; i<3; i++)
			center[i] += (*pIt)[i];
    }

  double invm=1.0/geo.noVerts;
  //Fd3op1(center,*=invm);
  for(int i=0; i<3; i++)
	center[i] *= invm; 


  for(unsigned int i=0;i!=geo.noFaces;i++)
    {

      face & f =faces[i];
      double *n=f.normal;
      const double 
	*f0 = pData[geo.faces[i][0]],
	*f1 = pData[geo.faces[i][1]],
	*f2 = pData[geo.faces[i][2]],
	*f3 = pData[geo.faces[i][3]];

	
      switch(geo.faceSizes[i])
	{
	  
	case 4:
	  //Fd3op3(a,=f2,-f0,);
	  //Fd3op3(b,=f3,-f1,);	  
		for(int i=0; i<3; i++)
		{
			a[i] = f2[i] - f0[i];
			b[i] = f3[i] - f1[i];	  
		}
		
	  Fd3kreuz(n,a,b);	  

	  //Fd3op3(a,=f1,-f0,);

	  //Fd3op5(b,=f0,+f1,+f2,+f3,);
	  //Fd3op3(b,=b,*0.25-center,);
		for(int i=0; i<3; i++)
		{
			a[i] = f1[i] - f0[i];
			b[i] = f0[i] + f1[i] + f2[i] + f3[i];
			b[i] = b[i] *0.25-center[i];
		}

	  f.isFlat=true;	  	  
	  
	  d=Fd3prod(b,n);
	  if(fabs(d)>fabs(lenForInv))lenForInv=d;
	  
	  if(fabs(Fd3prod(a,n)) > fabs(d) * 1e-14){
	    f.isFlat=false;
	    f.bilinear.init(f0,f1,f2,f3);	   
	  }
	  else{
	    f.normalSqlen=Fd3prod(n,n);
	  }

	  break;	    

	case 3:
	  //Fd3op3(a,=f1,-f0,);
	  //Fd3op3(b,=f2,-f1,);   
		for(int i=0; i<3; i++)
		{
			a[i] = f1[i] - f0[i];
			b[i] = f2[i] - f1[i];   
		}
	  
	  Fd3kreuz(n,a,b);	     	      

	  f.normalSqlen=Fd3prod(n,n);
	  
	  f.isFlat=true;	  

	  //Fd3op4(b,=f0,+f1,+f2,);
	  //Fd3op3(b,=b,/3-center,);
		for(int i=0; i<3; i++)
		{
			b[i] = f0[i] + f1[i] + f2[i];
			b[i] = b[i] / 3-center[i];
		}
	  
	  d = Fd3prod(b,n);      
	  if(fabs(d)>fabs(lenForInv))lenForInv=d;

	  break;
	  
	default:

	  cout<<"faces with more than 4 vertices not supported"<<endl;
	  assert(false);
	  
	}
    }
  //check if cell is inverted (i.e. the face normals point to inside instead of outside)

  inverted = (lenForInv<0);

}



  
//--------------------------------------------------------------------------- 

bool FCell::neighborFaceForPos(const FArray& /*pos*/,FIndex&/*faceId*/ ) const
{
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//--------------------------------------------------------------------------- 

unsigned int FCell::ntests = 0;
