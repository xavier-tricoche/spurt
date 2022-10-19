//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTriangleCell3D.cc,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:02 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.31 $
//
//--------------------------------------------------------------------------- 

#include "FTriangleCell3D.hh"

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
const FCell::geoInfo FTriangleCell3D::myGeoDescription =
  {
    // dimension
    3,
    // # vertices
    3, 
    // # edges
    3, 
    // # faces
    0, 
    // edges
    {{0,1},  {1,2},  {2,0}},
    // face sizes
    {},
    // faces
    {}
  };

//--------------------------------------------------------------------------- 
// now comes the fun stuff ...
//--------------------------------------------------------------------------- 


//--------------------------------------------------------------------------- 

FTriangleCell3D::FTriangleCell3D()
  : FCell2Din3D(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
  isBasisSet = false;
  isProjected = false;

//   cout << "FTriangleCell3D::FTriangleCell3D:" << endl;
//   for ( positive i=0 ; i<3 ; i++)
//     if ( positions[i].size() )
//       cout << "ERROR: positions[" << i << "] has size "
// 	   << positions[i].size() << endl;
//   for ( positive i=0 ; i<3 ; i++)
//     if ( tensors[i].size() )
//       cout << "ERROR: tensors[" << i << "] has size "
// 	   << tensors[i].size() << endl;
}

//--------------------------------------------------------------------------- 

FTriangleCell3D::FTriangleCell3D( const vector<FIndex>& vertIds )
  : FCell2Din3D(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
  for (char i=0; i<3; i++)
    vertexIndices[i]=vertIds[i] ;
  isBasisSet = false;
  isProjected = false;

//   cout << "FTriangleCell3D::FTriangleCell3D:" << endl;
//   for ( positive i=0 ; i<3 ; i++)
//     if ( positions[i].size() )
//       cout << "ERROR: positions[" << i << "] has size "
// 	   << positions[i].size() << endl;
//   for ( positive i=0 ; i<3 ; i++)
//     if ( tensors[i].size() )
//       cout << "ERROR: tensors[" << i << "] has size "
// 	   << tensors[i].size() << endl;
}

//--------------------------------------------------------------------------- 

FTriangleCell3D::FTriangleCell3D( const FIndex*vertIds )
  : FCell2Din3D(myGeoDescription, myIndices, myPositions, myPositionData,
		myTensors)
{
  for (char i=0; i<3; i++)
    vertexIndices[i]=vertIds[i] ;
  isBasisSet = false;
  isProjected = false;

//   cout << "FTriangleCell3D::FTriangleCell3D:" << endl;
//   for ( positive i=0 ; i<3 ; i++)
//     if ( positions[i].size() )
//       cout << "ERROR: positions[" << i << "] has size "
// 	   << positions[i].size() << endl;
//   for ( positive i=0 ; i<3 ; i++)
//     if ( tensors[i].size() )
//       cout << "ERROR: tensors[" << i << "] has size "
// 	   << tensors[i].size() << endl;
}

//--------------------------------------------------------------------------- 

FTriangleCell3D::~FTriangleCell3D()
{}

//--------------------------------------------------------------------------- 

FCell*  FTriangleCell3D::getClone() const
{
  try 
    {
      return new FTriangleCell3D( vertexIndices );
    }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

bool FTriangleCell3D::isInside( const FPosition& pos ) const
{
  try
    {
      double bb[3];
      baryCoord( bb, pos );

#ifndef NODEBUG
      if ( positions[0].size()==0 )
	{
	  THROW_EXCEPTION( FException, "ERROR: missing positions !");
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
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FTriangleCell3D::baryCoord( double* bary, const FPosition& pos ) const
{
#ifndef NODEBUG
  if ( positions[0].size()==0 )
    {
      THROW_EXCEPTION( FException, "ERROR: missing positions !");
    }
#endif

  if ( !isBasisSet )
    calculateTriBasis();

  // pos in local coordinates
  FArray p0p = pos - positions[0];
  double x, y, x1, x2, y2;
  x = p0p*basis[0];
  y = p0p*basis[1];
  x1 = (positions[1] - positions[0]) * basis[0];
  x2 = (positions[2] - positions[0]) * basis[0];
  y2 = (positions[2] - positions[0]) * basis[1];
  bary[2] = y / y2;
  bary[1] = ( x - bary[2]*x2 ) / x1;
  bary[0] = 1. - bary[1] - bary[2];
}

//--------------------------------------------------------------------------- 

char FTriangleCell3D::intersectEdge( FPosition& result, const FPosition& pos,
				     const FArray& vec ) const
{
  try 
    {
      vector< FArray > pos_;
      getVertices(pos_);

      if ( !isBasisSet )
	calculateTriBasis();

#ifndef NODEBUG
      if ( !isInside( pos ) )
	{
	  THROW_DEFAULT_EXCEPTION(FInvalidPositionException );
	}
//       else if ( fabs( vec * basis[2] ) > FTRIANGLECELL3D_epsilon )
// 	{
// 	  cout << "direction vector in intersectEdge is not tangential "
// 	       << "to cell plane (" 
// 	       << fabs( vec * basis[2] ) / vec.norm() << ")" << endl;
// 	}
#endif

      // to avoid increasing roundoff error, always project onto cell
      FVector localPos = pos - positions[0];
      FPosition thePos = positions[0] + 
	localPos - localPos * basis[2] * basis[2];
      FVector theVec = vec - ( vec * basis[2] ) * basis[2];

      // determine intersected edge with cross products
      // express vectors in local coordinates.
      double v[3][2], dir[2];

      for ( unsigned int i=0 ; i<3 ; i++ )
      {	
	//project vector pointing from thepos to a vertex into triangle space and normalize it
	v[i][0] = ( positions[i] - thePos ) * basis[0];
	v[i][1] = ( positions[i] - thePos ) * basis[1];
	double length=sqrt(v[i][0]*v[i][0]+v[i][1]*v[i][1]);
	v[i][0]*=1./length;
	v[i][1]*=1./length;
      }

      {
	//project dir into triangle space and normalize it
	dir[0] = theVec * basis[0];
	dir[1] = theVec * basis[1];
	double length=sqrt(dir[0]*dir[0]+dir[1]*dir[1]);
	dir[0] *=  1./length;
	dir[1] *=  1./length;
      }

      // vector indicating if vec may intersect the corresponding edge
      bool valid[3] = { true, true, true };
      double cross;
      for ( unsigned int i=0 ; i<3 ; i++ )
	{
	  cross = v[i][0]*dir[1] - v[i][1]*dir[0];
	  if((positions[i]-thePos).norm()>epsilon)
	    {
	      if ( cross > epsilon ) 
		{
		  valid[ (i+2) % 3 ] = false;
		}
	      else 
		if ( cross < -epsilon ) 
		  {
		    valid[ i ] = false;
		  }
		else 
		  if ( v[i][0]*dir[0] + v[i][1]*dir[1] < 0. )
		    {
		      valid[ (i+2) % 3 ] = false;
		      valid[ i ] = false;
		      break;
		    }
		  else 
		    if ( cross >= 0 ) // 0 <= cross < epsilon
		      {
			valid[ (i+2) % 3 ] = true;
			valid[ i ] = false;
			valid[ (i+1) % 3 ] = false;
			break;
		      }
		    else // -epsilon <= cross < 0
		      {
			valid[ i ] = true;
			valid[ (i+1) % 3 ] = false;
			valid[ (i+2) % 3 ] = false;
			break;
		      }
	    }
	}
      
      if( !( valid[0] ||  valid[1] || valid[2] ) ) 
	return -1;
      
      int id0 = 0, id1;
      while ( !valid[ id0 ] )
	++id0;
      id1 = ( id0 + 1 ) % 3;

      // compute intersection point
      double u = 
	( crossProduct( thePos - positions[id0], theVec ) * basis[2] ) /
	( crossProduct( positions[id1] - positions[id0], theVec ) *
	  basis[2] );

      result = (1.-u)*positions[id0] + u*positions[id1];

      return id0;
    }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FTriangleCell3D::calculateTriBasis() const
{
  try
    {

#ifndef NODEBUG
      if ( positions[0].size()==0 )
	{
	  THROW_EXCEPTION( FException, "ERROR: missing positions !");
	}
#endif
      // ---------------------------------------------------------------------
      // If p0, p1, p2 are the positions of our triangle, then the local basis 
      // e0, e1 and e2 is calculated as follows:
      // e0 = p1-p0 / |p0-p1| , e1 = p0-p2 / |p0-p2| , e2 = crossProduct(e0,e1)
      // NOTE: e2 is our triangle normal and l0=|p0-p1| and  l1=|p0-p2| !
      // ---------------------------------------------------------------------
      basis.resize(3);
      FArray p0p2(3);
      
      basis[0] = positions[1]-positions[0];
      basis[0].normalize();
      p0p2 = positions[2]-positions[0];
      basis[1] = p0p2 - (p0p2 * basis[0])*basis[0];
      basis[1].normalize();

      basis[2]= crossProduct( basis[0], basis[1] );
      isBasisSet = true;
    }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

const FVector& FTriangleCell3D::normal() const
{
#ifndef NODEBUG
  if ( positions[0].size()==0 )
    {
      THROW_EXCEPTION(FException, "ERROR: missing positions !");
    }
#endif
      
  if ( !isBasisSet )
    calculateTriBasis();

  return basis[2];
}

//--------------------------------------------------------------------------- 

void FTriangleCell3D::interpolate(FTensor& result,
				  const FPosition& position) const
{
  try
    {
    
#ifndef NODEBUG
      if ( positions[0].size()==0 || tensors[0].size()==0 ) {
	THROW_EXCEPTION(FException, "ERROR: missing position and/or tensor information");
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
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FTriangleCell3D::interpolateInPlane(FTensor& result,
					 const FPosition& position) const
{
  try
    {
#ifndef NODEBUG
      if ( positions[0].size()==0 || tensors[0].size()==0 ) {
	THROW_EXCEPTION( FException, "ERROR: missing position and/or tensor information");
      }
#endif
      // check to see if the tensor is an vector ..
      if( tensors[0].getOrder() != 1 )
	{
	  THROW_EXCEPTION( FException, "ERROR: Sofar we only tolerate vectors. ");
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
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FTriangleCell3D::getZeros(list<FAMSingularPoint>& result) const
{
  try
    {

#ifndef NODEBUG
      if ( positions[0].size()==0 || tensors[0].size()==0 ) {
	THROW_EXCEPTION( FException, "ERROR: missing position and/or tensor information");
      }
#endif

      if( !isBasisSet )
	calculateTriBasis();
      
      if ( !isProjected )
	projectTensors();

      // check to see if the tensor is a vector ...
      if( tensors[0].getOrder() != 1 )
	{
	  THROW_EXCEPTION( FException, "ERROR: So far we are only excepting vectors for  FTriangleCell3D::getZeros( *** ). ");
	}
      

      //------------------------------------------------------------------
      // Calculate the singular points which lie inside the triangle.
      // System to solve:
      // 0 = alpha1*(v1Re-v0Re)+alpha2(v2Re-v0Re)+v0Re, where Re^=Reduced
      // in local triangle coordinates.
      // Strategy: Cram rule
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
//       if( fabs(det1) > FTRIANGLECELL3D_epsilon &&
// 	  fabs(det2) > FTRIANGLECELL3D_epsilon &&
// 	  fabs(det)  > FTRIANGLECELL3D_epsilon )

      if ( det != 0. )
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
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FTriangleCell3D::projectTensors() const
{
  try {
    
#ifndef NODEBUG

    if ( tensors[0].size()==0 ) {
      THROW_EXCEPTION( FException, "ERROR: missing position and/or tensor information");
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
      THROW_EXCEPTION( FException, "ERROR: Sofar we are only excepting vectors. ");
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
  CATCH_N_RETHROW( FException );
}


//---------------------------------------------------------------------------

FCell::CellType FTriangleCell3D::getCellType(void) const
{
  return FCell::TRIANGLE_3D;
}



//--------------------------------------------------------------------------- 

void FTriangleCell3D::derivatives(FTensor& , const FPosition&) const
{
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
}

//--------------------------------------------------------------------------- 

vector<FTensor> FTriangleCell3D::getReducedVectors()
{
  if ( !isProjected )
    projectTensors();

  return projTensors ;
}

vector<FTensor> FTriangleCell3D::getLocalVectors()
{
  if ( !isProjected )
    projectTensors();

  return locTensors ;
}
//--------------------------------------------------------------------------- 

void FTriangleCell3D::getSingularType(FAMSingularPoint &result) const
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
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FTriangleCell3D::calculateLocTriPos() const
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

void FTriangleCell3D::getLocTriPos(FArray &world, FArray &loc) const
{

  if( !isBasisSet )
    calculateTriBasis();
  
  loc[0] = (world-positions[0])*basis[0];
  loc[1] = (world-positions[0])*basis[1];

}


//--------------------------------------------------------------------------- 


void FTriangleCell3D::getZerosArbitrary(list<FAMSingularPoint>& result) const
{
  try
    {

#ifndef NODEBUG
      if ( positions[0].size()==0 || tensors[0].size()==0 ) {
	THROW_EXCEPTION(FException, "ERROR: missing position and/or tensor information");
      }
#endif

      if( !isBasisSet )
	calculateTriBasis();
      
      if ( !isProjected )
	projectTensors();

      // check to see if the tensor is an vector ..
      if( tensors[0].getOrder() != 1 )
	{
	  THROW_EXCEPTION(FException, "ERROR: Sofar we are only excepting vectors. ");
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
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

positive FTriangleCell3D::memSize() const
{
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
