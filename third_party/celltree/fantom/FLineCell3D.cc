//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FLineCell3D.cc,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:02 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.31 $
//
//--------------------------------------------------------------------------- 

#include "FLineCell3D.hh"

#include "eassert.hh"
#include "FException.hh"

#include "FArray.hh"
#include "FMatrix.hh"
#include "FPosition.hh"

#include <iostream>
#include <utility>
#include <vector>

//==========================================================================

#define FTRIANGLECELL3D_epsilon 1.0e-8

//==========================================================================

// first the cell's geometry description
const FCell::geoInfo FLineCell3D::myGeoDescription =
  {
    // dimension
    3,
    // # vertices
    2, 
    // # edges
    1, 
    // # faces
    0, 
    // edges
    {{0,1}},
    // face sizes
    {},
    // faces
    {}
  };

//--------------------------------------------------------------------------- 
// now comes the fun stuff ...
//--------------------------------------------------------------------------- 


//--------------------------------------------------------------------------- 

FLineCell3D::FLineCell3D()
  : FCell1Din3D(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
}

//--------------------------------------------------------------------------- 

FLineCell3D::FLineCell3D( const vector<FIndex>& vertIds )
  : FCell1Din3D(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
  eassert( vertIds.size() == 2 );
  for (char i=0; i<2; i++)
    vertexIndices[i]=vertIds[i] ;
}

//--------------------------------------------------------------------------- 

FLineCell3D::FLineCell3D( const FIndex*vertIds )
  : FCell1Din3D(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
  for (char i=0; i<2; i++)
    vertexIndices[i]=vertIds[i] ;
}

//--------------------------------------------------------------------------- 

FLineCell3D::~FLineCell3D()
{}

//--------------------------------------------------------------------------- 

FCell*  FLineCell3D::getClone() const
{
  try 
    {
      return new FLineCell3D( vertexIndices );
    }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

bool FLineCell3D::isInside( const FPosition& /*pos*/ ) const
{
  return false;
#if 0
  try
    {
      double bb[3];
      baryCoord( bb, pos );

#ifndef NODEBUG
      if ( positions[0].size()==0 )
	{
	  FException e("ERROR: missing positions !");
	  e.addTraceMessage("void FLineCell3D::isInside()");
	  throw e;
	}
#endif

      if ( !isBasisSet )
	calculateTriBasis();
      
      //---------------------------------------------------------------------
      // To check if a given position pos lies within the triangle, project
      // the vector p0pos=pos-p0 onto the local basis vectors. h in this case 
      // is the distance from the plane and alpha1and2 are the barycentric
      // coordinates of the position in lokal triangle coord. 
      //---------------------------------------------------------------------
      double h;
      FArray p0pos(3);
      
      p0pos = pos-positions[0];
      h = p0pos*basis[2];

      // well our pos is to far away from our plane ...
      if( fabs(h) >  FTRIANGLECELL3D_epsilon ) 
	{
	  cout << "normal distance was too large: " << h << endl;
	  return false;
	}
      
      double b[3];
      baryCoord( b, pos );

      // not inside the plane ...
      if( ! ( b[0] >= -epsilon && b[1] >= -epsilon && b[2] >= -epsilon ) )
	return false;
      
      return true;
    }
  catch(FException& e)
    {
      e.addTraceMessage("void FLineCell3D::isInside()");
      throw;
    }
#endif
}

//--------------------------------------------------------------------------- 

void FLineCell3D::interpolate(FTensor& /*result*/,
				  const FPosition& /*position*/) const
{
  THROW_EXCEPTION( FException, "Cannot interpolate in 3D Line Cell" );
#if 0
  try
    {
    
#ifndef NODEBUG
      if ( positions[0].size()==0 || tensors[0].size()==0 ) {
	FException e("ERROR: missing position and/or tensor information");
	throw e;
      }
#endif
      
      
      if ( !isBasisSet )
	calculateTriBasis();
      
      //-----------------------------------------------------------------------
      // The interpolation (barycentric) has the form:
      // T = b[0]*T1 + b[1]*T2 + b[2]*T3
      // T0,1 and 2 are the vectors according to the positions of the triangle
      //-----------------------------------------------------------------------
      double b[3];
      baryCoord( b, position );
      result = b[0]*tensors[0] + b[1]*tensors[1] + b[2]*tensors[2];
    }
  catch( FException& e )
    {
      e.addTraceMessage("void FLineCell3D::interpolate()");
      throw;
    }
#endif
}

//--------------------------------------------------------------------------- 
#if 0
void FLineCell3D::interpolateInPlane(FTensor& result,
					 const FPosition& position) const
{
  try
    {
#ifndef NODEBUG
      if ( positions[0].size()==0 || tensors[0].size()==0 ) {
	FException e("ERROR: missing position and/or tensor information");
	throw e;
      }
#endif
      // check to see if the tensor is an vector ..
      if( tensors[0].getOrder() != 1 )
	{
	  FException e("ERROR: Sofar we only tolerate vectors. ");
	  e.addTraceMessage("void FLineCell3D::interpolate");
	  throw e;
	}
      
      if ( !isBasisSet )
	calculateTriBasis();
      
      if ( !isProjected )
	projectTensors();
      
      //---------------------------------------------------------------------
      // The interpolation (barycentric) has the form:
      // T = alpha0*T1 + alpha1*T2 + alpha2*T3
      // T0,1 and 2 are the vectors according to the positions of the 
      // triangle
      //---------------------------------------------------------------------
      double b[3];
      baryCoord( b, position );
      
      result = ( FTensor ) ( b[0] * projTensors[0] + 
			     b[1] * projTensors[1] + 
			     b[2] * projTensors[2] );
    }
  catch( FException& e )
    {
      e.addTraceMessage("void FLineCell3D::interpolateInPlane()");
      throw;
    }
}

//--------------------------------------------------------------------------- 
#endif
void FLineCell3D::getZeros(list<FAMSingularPoint>& /*result*/) const
{
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
#if 0
  try
    {

#ifndef NODEBUG
      if ( positions[0].size()==0 || tensors[0].size()==0 ) {
	FException e("ERROR: missing position and/or tensor information");
	e.addTraceMessage("void FLineCell3D::getZero()");
	throw e;
      }
#endif

      if( !isBasisSet )
	calculateTriBasis();
      
      if ( !isProjected )
	projectTensors();

      // check to see if the tensor is an vector ..
      if( tensors[0].getOrder() != 1 )
	{
	  FException e("ERROR: Sofar we are only excepting vectors. ");
	  e.addTraceMessage("void FLineCell3D::getZero()");
	  throw e;
	}
      

      //------------------------------------------------------------------
      // Calculate the singular points which lie inside the triangle.
      // System to solve:
      // 0 = alpha1*(v1Re-v0Re)+alpha2(v2Re-v0Re)+v0Re,wherby Re^=Reduced
      // in lokal triangle coordinates.
      // Strategy: Cram'sches System !
      //------------------------------------------------------------------
      
      FMatrix tmpM(2,2);
      double det,det1,det2;
      double alpha1,alpha2;
     
      // cramer rules ...
      
      tmpM(0,0) =-1*locTensors[0](0) ; // ^= x0 of v0Re
      tmpM(1,0) =-1*locTensors[0](1) ; // ^= y0 of v0Re , etc.
      tmpM(0,1) = locTensors[2](0)-locTensors[0](0) ;
      tmpM(1,1) = locTensors[2](1)-locTensors[0](1) ;

      det1 = tmpM.detOf();
                 
      tmpM(0,1) =-1*locTensors[0](0) ; // ^= x0 of v0Re
      tmpM(1,1) =-1*locTensors[0](1) ; // ^= y0 of v0Re , etc.
      tmpM(0,0) = locTensors[1](0)-locTensors[0](0) ;
      tmpM(1,0) = locTensors[1](1)-locTensors[0](1) ;

      det2 = tmpM.detOf();

      tmpM(0,0) = locTensors[1](0)-locTensors[0](0) ;
      tmpM(1,0) = locTensors[1](1)-locTensors[0](1) ;  
      tmpM(0,1) = locTensors[2](0)-locTensors[0](0) ;
      tmpM(1,1) = locTensors[2](1)-locTensors[0](1) ;

      det = tmpM.detOf();

      if( fabs(det1) > FTRIANGLECELL3D_epsilon &&
	  fabs(det2) > FTRIANGLECELL3D_epsilon &&
	  fabs(det)  > FTRIANGLECELL3D_epsilon )
	{
		  
	  alpha1 = det1/det;
	  alpha2 = det2/det;
	  
	  // now check if position lies whithin the triangle ...
	  if( ((alpha1>=0) && 
	       (alpha2>=0) && 
	       (1-alpha1-alpha2>=0)) )
	    {
	      // calculate position in world space ...
	      FPosition singularPos(3);
	      
	      singularPos = (1-alpha1-alpha2)*positions[0] +
		alpha1*positions[1] +
		alpha2*positions[2] ;

	      FAMSingularPoint dummy( FPosition( 0., 0.) );
	      dummy.setOrder( 1 );
	      getSingularType( dummy );
	      dummy.setPosition( singularPos );

	      result.push_back( dummy );
	    }
	}
    }
  catch(FException& e){
     e.addTraceMessage("void FLineCell3D::getZero()");
    throw;
  }
#endif
}

//--------------------------------------------------------------------------- 
#if 0
void FLineCell3D::projectTensors() const
{
  try {
    
#ifndef NODEBUG

    if ( tensors[0].size()==0 ) {
      FException e("ERROR: missing position and/or tensor information");
      e.addTraceMessage("void FLineCell3D::projectTensor() const");
      throw e;
    }
//     cout << "tensors[0] has size " << tensors[0].size() << endl;
//     cout << "current cell has @" << this << endl;
//     cout << "in projectTensors: tensor values:" << endl;
//     for (positive i=0 ; i<3 ; i++)
//       cout << i << ": order = " << tensors[i].getOrder() << flush
// 	   << ", dimension = " << tensors[i].getDimension() << flush
// 	   << ", values = " << tensors[i] << endl;
//     cout << endl;
	

#endif

    if( !isBasisSet )
      calculateTriBasis();
    
    // check to see if the tensor is an vector ..
    if( tensors[0].getOrder() != 1 ){
      FException e("ERROR: Sofar we are only excepting vectors. ");
      e.addTraceMessage("void FLineCell3D::projectTensor const");
      throw e;
    }

    // project onto plane: worldcoordinates(3D)!
    FTensor e0(basis[0]),e1(basis[1]);
    projTensors.resize(3);
 
    projTensors[0].resizeTensor(3,1);
    projTensors[1].resizeTensor(3,1);
    projTensors[2].resizeTensor(3,1);

    projTensors[0] = (tensors[0]*e0)*e0 + (tensors[0]*e1)*e1;
    projTensors[1] = (tensors[1]*e0)*e0 + (tensors[1]*e1)*e1;
    projTensors[2] = (tensors[2]*e0)*e0 + (tensors[2]*e1)*e1;
    
    // projected vectors in lokal triangle coordinates(2D)!
    locTensors.resize(3);

    locTensors[0].resizeTensor(2,1);
    locTensors[1].resizeTensor(2,1);
    locTensors[2].resizeTensor(2,1);
    
    locTensors[0].setValue( 0, projTensors[0]*e0 ) ;
    locTensors[0].setValue( 1, projTensors[0]*e1 ) ;

    locTensors[1].setValue( 0, projTensors[1]*e0 ) ;
    locTensors[1].setValue( 1, projTensors[1]*e1 ) ;
    
    locTensors[2].setValue( 0, projTensors[2]*e0 ) ;
    locTensors[2].setValue( 1, projTensors[2]*e1 ) ;

    isProjected = true;
  }
  catch(FException& e){
    e.addTraceMessage("void FLineCell3D::projectTensor() const");
    throw;
  }
}

#endif
//---------------------------------------------------------------------------

FCell::CellType FLineCell3D::getCellType(void) const
{
  return FCell::LINE_3D;
}



//--------------------------------------------------------------------------- 

void FLineCell3D::derivatives(FTensor& /*result*/, const FPosition&) const
{
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//--------------------------------------------------------------------------- 
#if 0
vector<FTensor> FLineCell3D::getReducedVectors()
{
  if ( !isProjected )
    projectTensors();

  return projTensors ;
}

//--------------------------------------------------------------------------- 

void FLineCell3D::getSingularType(FAMSingularPoint &result) const
{
  try{

    if( !isBasisSet )
      calculateTriBasis();
    
    if ( !isProjected )
      projectTensors();

    FMatrix A(2,2);
    double db0dx, db0dy, db1dx, db1dy, denom;

    this->calculateLocTriPos();
    
    // computing singularity type
    denom = 1. / 
      ( (locPos[0][0]-locPos[2][0]) * (locPos[1][1]-locPos[2][1]) - 
	(locPos[0][1]-locPos[2][1]) * (locPos[1][0]-locPos[2][0]) );

    db0dx =   (locPos[1][1]-locPos[2][1]) * denom;
    db0dy = - (locPos[1][0]-locPos[2][0]) * denom;
    db1dx = - (locPos[0][1]-locPos[2][1]) * denom;
    db1dy =   (locPos[0][0]-locPos[2][0]) * denom;
    
    A(0,0) = db0dx*locTensors[0](0) + db1dx*locTensors[1](0) 
      - (db0dx+db1dx)*locTensors[2](0);
    
    A(0,1) = db0dy*locTensors[0](0) + db1dy*locTensors[1](0) 
      - (db0dy+db1dy)*locTensors[2](0);
    
    A(1,0) = db0dx*locTensors[0](1) + db1dx*locTensors[1](1) 
      - (db0dx+db1dx)*locTensors[2](1);
    
    A(1,1) = db0dy*locTensors[0](1) + db1dy*locTensors[1](1) 
      - (db0dy+db1dy)*locTensors[2](1);

    result.setLinearNature(A);

    // replace eigenvectors by their expression in global coordinates
    vector< FVector > oldEigenvectors, newEigenvectors;
    result.getEigenvectors( oldEigenvectors );
    if ( oldEigenvectors.size() )
      {
	newEigenvectors.clear();
	newEigenvectors.push_back( oldEigenvectors[0](0) * basis[0] +
				   oldEigenvectors[0](1) * basis[1] );
	newEigenvectors.push_back( oldEigenvectors[1](0) * basis[0] +
				   oldEigenvectors[1](1) * basis[1] );
	result.setEigenvectors( newEigenvectors );

//        	cout << "-----------------------" << endl;
// 	for(int i=0;i<newEigenvectors.size();i++)
// 	  cout << "newEigenvectors: " << newEigenvectors[i] << endl;
// 	cout << "-----------------------" << endl;
      }
    
    // add max norm of enclosing cell
    double tmp_norm, max_norm = locTensors[0].norm();
    tmp_norm = locTensors[1].norm();
    if (tmp_norm > max_norm)
      max_norm = tmp_norm;
    tmp_norm = locTensors[2].norm();
    if (tmp_norm > max_norm)
      max_norm = tmp_norm;
    result.setMaxNorm(max_norm);
    

  }
  catch(FException &e){
    e.addTraceMessage("void FLineCell3D::getSingularType() const");
    throw;
  }
  
}

//--------------------------------------------------------------------------- 

void FLineCell3D::calculateLocTriPos() const
{
  if( !isBasisSet )
    calculateTriBasis();

  FArray tmp(3);

  locPos.resize(3);
  locPos[0].resize(2);
  locPos[1].resize(2);
  locPos[2].resize(2);

  locPos[0](0) = 0;
  locPos[0](1) = 0;

  locPos[1](0) = basis[0] * ( positions[1] - positions[0] );
  locPos[1](1) = 0;
  
  locPos[2](0) = basis[0] * ( positions[2] - positions[0] );
  locPos[2](1) = basis[1] * ( positions[2] - positions[0] );
  
}

//--------------------------------------------------------------------------- 

void FLineCell3D::getLocTriPos(FArray &world, FArray &loc) const
{

  if( !isBasisSet )
    calculateTriBasis();
  
  loc[0] = (world-positions[0])*basis[0];
  loc[1] = (world-positions[0])*basis[1];

}

#endif
//--------------------------------------------------------------------------- 

positive FLineCell3D::memSize() const
{
  /**
   * \todo FIXME
   */
  return
    (
     3//basis
     +3//locPos
    ) *3*sizeof(double) 
    +
    (tensorData?
     tensors[0].size()*sizeof(double)*
     sizeOfCellType() //tensorData
     + (sizeof(FTensor)+tensors[0].size()*sizeof(double))
     *(3//projTensors
       +3//locTensors
       )
     :0)
    +
    sizeof(*this);
}
//----------------------------------------------------------------------------
