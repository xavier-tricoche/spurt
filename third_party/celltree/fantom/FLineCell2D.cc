//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FLineCell2D.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:13 $
// Author:    $Author: garth $
// Version:   $Revision: 1.19 $
//
//--------------------------------------------------------------------------- 

#include "FLineCell2D.hh"
#include "FException.hh"

#include "FMatrix.hh"
#include "FLineCell2D.hh"
#include "FArray.hh"
#include "FPosition.hh"

#include <iostream>
#include <utility>
#include <vector>

// the cell's geometry description
const FCell::geoInfo FLineCell2D::myGeoDescription =
  {
    // dimension
    2,
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

FLineCell2D::FLineCell2D()
  : FCell1Din2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{}

//--------------------------------------------------------------------------- 

FLineCell2D::FLineCell2D( const vector<FIndex>& vertIds )
  : FCell1Din2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors),lastPos(1)
{
  for (char i=0; i<2; i++)
    vertexIndices[i]=vertIds[i] ;

//  lastBs[0] = -1.;
}

//--------------------------------------------------------------------------- 

FLineCell2D::FLineCell2D( const FIndex*vertIds )
  : FCell1Din2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors),lastPos(1)
{
  for (char i=0; i<2; i++)
    vertexIndices[i]=vertIds[i] ;

//  lastBs[0] = -1.;
}

//--------------------------------------------------------------------------- 


FLineCell2D::~FLineCell2D()
{}

//--------------------------------------------------------------------------- 

FCell*  FLineCell2D::getClone() const
{
  try 
    {
      return new FLineCell2D( vertexIndices );
    }
  CATCH_N_RETHROW( FException );
}

//--------------------------------------------------------------------------- 

void FLineCell2D::interpolate(FTensor& /*result*/, 
				  const FPosition& /*position*/) const
{
  THROW_EXCEPTION( FException, "Cannot interpolate in 2DLine Cells" );
#if 0
  try {
    if (!positions || !tensors) {
      FException e("ERROR: missing position and/or tensor information");
      throw e;
    }
    
    if (set_to_zero) {
      result = FTensor(tensors[0].getDimension(), 
		       tensors[0].getOrder());
      result = 0.;
    }
    else {
      
      double b[3];
      
      // compute barycentric coordinates
      barycentricCoordinates(b, position);
      
      // set correct dimension and order + initialize
      result = FTensor(tensors[0].getDimension(), tensors[0].getOrder(), true);
      // linear combination
      for (positive i=0 ; i<3 ; i++) 
	result += b[i]*tensors[i];
    }
  }
  catch (FException& e) {
    e.addTraceMessage("void FLineCell2D::interpolate(FTensor& result, const FPosition& position)");
    throw;
  }
#endif
} 

//--------------------------------------------------------------------------- 

void FLineCell2D::derivatives(FTensor& /*result*/, const FPosition&) const
{
  THROW_EXCEPTION( FException, "Cannot interpolate in 2DLine Cells" );
#if 0
  try {
    if (!positions || !tensors) {
      FException e("ERROR: missing position and/or tensor information");
      throw e;
    }

    if (set_to_zero) {
      result = FTensor(tensors[0].getDimension(), 
		       tensors[0].getOrder()+1);
      result = 0.;
    }
    else {
      double dbdx[3], dbdy[3];
      positive order;
	
      order = tensors[0].getOrder();
      positive dimension = tensors[0].getDimension();

      derivBaryCoord(dbdx, dbdy);
      vector<double> val1, val2, val3;

      tensors[0].getValues(val1);
      tensors[1].getValues(val2);
      tensors[2].getValues(val3);
	
      positive size = val1.size();
      vector<double> comp(2*size);
	
      for (positive i=0 ; i<size ; i++) {
	comp[i] = dbdx[0]*val1[i] + dbdx[1]*val2[i] + dbdx[2]*val3[i];
	comp[size+i] = dbdy[0]*val1[i] + dbdy[1]*val2[i] + dbdy[2]*val3[i];
      }
	
      result = FTensor(dimension, order+1, comp);
    }
  }
  catch (FException& e) {
    e.addTraceMessage("void FLineCell2D::derivatives(FTensor& result, const FPosition& position)");
    throw;
  }
#endif
}

//--------------------------------------------------------------------------- 

bool FLineCell2D::isInside(const FPosition& /*position*/) const
{
  return false;
#if 0
  try {
    if (!positions) {
      FException e("ERROR: position information missing");
      throw e;
    }
    

#ifndef NODEBUG
    if (position.getDimension() != 2)
      throw FInvalidDimensionException();
#endif

    if (position == lastPos) {

//       cout << "dynamic isInside: reuse" << endl;
      return true;
    }


    double b[2] = {0., 0.};    
    double denom = ( (positions[0][0]-positions[2][0])
		     *(positions[1][1]-positions[2][1]) -
		     (positions[0][1]-positions[2][1])
		     *(positions[1][0]-positions[2][0]) );
    double alpha = - epsilon * denom;
    
    b[0] = ( (position[0]-positions[2][0]) *
	     (positions[1][1]-positions[2][1]) -
	     (position[1]-positions[2][1]) *
	     (positions[1][0]-positions[2][0]) );
    
    if ( ( (denom > 0.) && (b[0] > alpha) ) ||
	 ( (denom < 0.) && (b[0] < alpha) ) ) {
      
      b[1] = ( (positions[0][0]-positions[2][0]) *
	       (position[1]-positions[2][1]) -
	       (positions[0][1]-positions[2][1]) *
	       (position[0]-positions[2][0]));
      
      if ( ( (denom > 0.) && (b[1] > alpha) && 
	     (denom - b[0] - b[1] > alpha) ) ||
	   ( (denom < 0.) && (b[1] < alpha) && 
	     (denom - b[0] - b[1] < alpha) ) ) {

	// save results for reuse by interpolate
	double tmp = 1. / denom;
	lastBs[0] = b[0] * tmp;
	lastBs[1] = b[1] * tmp;
	lastBs[2] = 1. - lastBs[0] - lastBs[1];
	lastPos = position;

// 	cout << "dynamic isInside said: bary. coord are " 
// 	     << lastBs[0] << ", " << lastBs[1] << ", " << lastBs[2] << endl;
	
	return true;
      }
    }


//     cout << "dynamic isInside test failed: " << b[0] << ", " << b[1] 
// 	 << ", " << b[2] << endl;

    return false;

  }
  catch (FException& e) {
    e.addTraceMessage("bool FLineCell2D::isInside(const FPosition& position) const");
    throw;
  }
#endif
}

//--------------------------------------------------------------------------- 

bool  FLineCell2D::isInside(const FPosition& /*pos*/,
				const vector< FPosition >& /*vertices*/)
{
  return false;
#if 0
  try 
    {

//       cout << "static isInside" << endl;

#ifndef NODEBUG
      // checking dimension of the given position
      if (pos.getDimension() != 2)
	throw FInvalidDimensionException();
#endif

      // compute barycentric coordinates
      
      double b[2] = {0., 0.} ;
      double denom = ( (vertices[0][0]-vertices[2][0]) *
		       (vertices[1][1]-vertices[2][1]) -
		       (vertices[0][1]-vertices[2][1]) *
		       (vertices[1][0]-vertices[2][0]) );
      double alpha = - epsilon * denom;
      
      b[0] = ( (pos[0]-vertices[2][0]) *
	       (vertices[1][1]-vertices[2][1]) -
	       (pos[1]-vertices[2][1]) *
	       (vertices[1][0]-vertices[2][0]));

      if ( (alpha < 0 && b[0] > alpha) ||
	   (alpha > 0 && b[0] < alpha) ) {
      
	b[1] = ( (vertices[0][0]-vertices[2][0]) *
		 (pos[1]-vertices[2][1]) -
		 (vertices[0][1]-vertices[2][1]) *
		 (pos[0]-vertices[2][0]));

	if ( (alpha < 0 && b[1] > alpha && denom - b[0] - b[1] > alpha) ||
	     (alpha > 0 && b[1] < alpha && denom - b[0] - b[1] < alpha) ) {

// 	  cout << "static isInside: pos = " << pos 
// 	       << ", vertices: " << endl
// 	       << " #1 = " << vertices[0] << ", " << endl
// 	       << " #2 = " << vertices[1] << ", " << endl
// 	       << " #3 = " << vertices[2] << ", " << endl
// 	       << ", bary. coord. = " << b[0]/denom << ", "
// 	       << b[1]/denom << ", " << 1. - b[0]/denom - b[1]/denom
// 	       << endl << "-> TRUE"
// 	       << endl;

	  return true;
	}
      }

//       cout << "static isInside: pos = " << pos
// 	   << ", vertices: " << endl
// 	   << " #1 = " << vertices[0] << ", " << endl
// 	   << " #2 = " << vertices[1] << ", " << endl
// 	   << " #3 = " << vertices[2] << ", " << endl
// 	   << ", b[0] = " << b[0]/denom << ", b[1] = " 
// 	   << b[1]/denom << ", b[2] = "
// 	   << 1. - b[1]/denom - b[0]/denom
// 	   << endl << "FALSE" << endl;

      return false;
    }
  catch (FException& e) 
    {
      e.addTraceMessage("bool FLineCell2D::isInside(const FPosition& pos, const vector< FPosition >& vertices)");

      cout << e << endl;

      throw;
    }
#endif
}



//--------------------------------------------------------------------------- 

void FLineCell2D::buildBoundingBox(void) const
{
  try {
    
    if (geometryOK)
      return;

    if (!positions) {
      FException e("ERROR: position information missing");
      throw e;
    }
    
    double minX, minY, maxX, maxY;
    minX = positions[0][0];
    maxX = minX;
    minY = positions[0][1];
    maxY = minY;

      if (positions[1][0] < minX) 
	minX = positions[1][0];
      else if (positions[1][0] > maxX) 
	maxX = positions[1][0];

      if (positions[1][1] < minY) 
	minY = positions[1][1];
      else if (positions[1][1] > maxY) 
	maxY = positions[1][1];


    bBox = FBoundingBox(minX, minY, maxX, maxY);
  }
  catch (FException& e) {
    e.addTraceMessage("void FLineCell2D::buildBoundingBox(void)");
    throw;
  }
}

//--------------------------------------------------------------------------- 
#if 0 
void FLineCell2D::barycentricCoordinates(double b[3], 
					     const FPosition& position) const
{
  try {

    if (!positions) {
      FException e("ERROR: position information missing");
      throw e;
    }


    // check for possible reuse of precomputed values
    if (lastBs[0] >= -.5 && lastPos == position) {
      b[0] = lastBs[0];
      b[1] = lastBs[1];
      b[2] = lastBs[2];

//       cout << "position " << position << " has barycentric coord. " 
// 	   << b[0] << ", " << b[1] << ", " << b[2] << " (reuse)" << endl;

      return;
    }
      

    // barycentric coordinates
    double denom = 1. / 
      ( (positions[0][0]-positions[2][0])
	*(positions[1][1]-positions[2][1]) -
	(positions[0][1]-positions[2][1])
	*(positions[1][0]-positions[2][0]) );
    
    b[0] = ((position[0]-positions[2][0])
	    *(positions[1][1]-positions[2][1]) -
	    (position[1]-positions[2][1])
	    *(positions[1][0]-positions[2][0])) * denom;
    
    b[1] = ((positions[0][0]-positions[2][0])
	    *(position[1]-positions[2][1]) -
	    (positions[0][1]-positions[2][1])
	    *(position[0]-positions[2][0])) * denom;
  
    b[2] = 1.0 - b[0] - b[1];


//     cout << "position " << position << " has barycentric coord. " 
// 	 << b[0] << ", " << b[1] << ", " << b[2] << endl;

//     if (! ( b[0] > -epsilon && b[0] < 1. + epsilon &&
// 	    b[1] > -epsilon && b[1] < 1. + epsilon &&
// 	    b[2] > -epsilon && b[2] < 1. + epsilon ) )

//       cout << "WARNING WARNING WARNING" << endl;

  }
  catch (FException& e) {
    e.addTraceMessage("void FLineCell2D::barycentricCoordinates(double b[3], const FPosition& position) const");
    throw;
  }
}

//--------------------------------------------------------------------------- 

void FLineCell2D::derivBaryCoord(double dbdx[3], double dbdy[3]) const
{
  try {

    if (!positions) {
      FException e("ERROR: position information missing");
      throw e;
    }
  
    // barycentric coordinates derivatives
    double denom = 1. /
      ( (positions[0][0]-positions[2][0])*(positions[1][1]-positions[2][1]) -
	(positions[0][1]-positions[2][1])*(positions[1][0]-positions[2][0]) );
    
    dbdx[0] = (positions[1][1]-positions[2][1]) * denom;
    dbdx[1] = (positions[2][1]-positions[0][1]) * denom;
    dbdx[2] = (positions[0][1]-positions[1][1]) * denom;
    
    dbdy[0] = (positions[2][0]-positions[1][0]) * denom;
    dbdy[1] = (positions[0][0]-positions[2][0]) * denom;
    dbdy[2] = (positions[1][0]-positions[0][0]) * denom;

  }
  catch (FException& e) {
    e.addTraceMessage("void FLineCell2D::derivBaryCoord(double dbdx[3], double dbdy[3]) const");
    throw;
  }
}
#endif
//---------------------------------------------------------------------------

FCell::CellType FLineCell2D::getCellType(void) const
{
  return FCell::LINE_2D;
}

//--------------------------------------------------------------------------- 
 
void FLineCell2D::getZeros(list<FAMSingularPoint>& /*result*/) const
{
  THROW_DEFAULT_EXCEPTION( FNotImplementedException );
#if 0
  try {

    if (!positions || !tensors) {
      FException e("ERROR: missing position and/or tensor information");
      throw e;
    }
  
    if (set_to_zero) 
      return;
  
    FTensor testTensor;
    
    // check tensor type: 2D vectors or tensors expected
    if ((tensors[0].getOrder() != 1 && tensors[0].getOrder() != 2) || 
	tensors[0].getDimension() != 2) {
      FInvalidDimensionException e;
      throw e;
    }
  
    double b1, b2, denom;
    FTensor *T = this->tensors; // simplify notations

    // solve (b0*T0+b1*T1+(1-b0-b1)*T2) = 0 in b0, b1.
    if (T[0].getOrder() == 1) { // vector case
      b1    = (-T[2](0) * (T[1](1)-T[2](1))) - (-T[2](1) * (T[1](0)-T[2](0)));
      b2    = (-T[2](1) * (T[0](0)-T[2](0))) - (-T[2](0) * (T[0](1)-T[2](1)));
      denom = ((T[0](0)-T[2](0)) * (T[1](1)-T[2](1))) - 
	((T[0](1)-T[2](1)) * (T[1](0)-T[2](0)));
    } 
    else { // tensor case
      b1    = (T[2](1,1)-T[2](0,0))*(T[1](0,1)-T[2](0,1))
	- (-T[2](0,1))*(T[1](0,0)-T[1](1,1)-T[2](0,0)+T[2](1,1));
      b2    = (T[0](0,0)-T[0](1,1)-T[2](0,0)+T[2](1,1))*(-T[2](0,1))
	- (T[0](0,1)-T[2](0,1))*(T[2](1,1)-T[2](0,0));
      denom = (T[0](0,0)-T[0](1,1)-T[2](0,0)+T[2](1,1))*(T[1](0,1)-T[2](0,1)) -
	(T[1](0,0)-T[1](1,1)-T[2](0,0)+T[2](1,1))*(T[0](0,1)-T[2](0,1));
    }
  
    if (denom != 0.) {
      denom = 1. / denom;
      b1 *= denom;
      b2 *= denom;
    }
    else {
      cerr << "singular system encountered in void FLineCell2D::getZeros()"
	   << endl;
      return;
    }
  
    // simplify notations
    FPosition *P = this->positions;

    if ((b1>=-epsilon) && (b2>=-epsilon) && (b1+b2<=1+epsilon)) {

      FPosition pos = P[0] * b1 + P[1] * b2 + P[2] * (1.-b1-b2);

      if (T[0].getOrder() == 1) { // vector case
	// add the found zero to the list
	FAMSingularPoint tmp = FAMSingularPoint(pos, FAMSingularPoint::NONE, 1);
	// computing singularity type
	FMatrix A(2,2);
	double db1dx, db1dy, db2dx, db2dy, denom;
	denom = 1. / 
	  ( (P[0][0]-P[2][0]) * (P[1][1]-P[2][1]) - 
	    (P[0][1]-P[2][1]) * (P[1][0]-P[2][0]) );
	db1dx =   (P[1][1]-P[2][1]) * denom;
	db1dy = - (P[1][0]-P[2][0]) * denom;
	db2dx = - (P[0][1]-P[2][1]) * denom;
	db2dy =   (P[0][0]-P[2][0]) * denom;
	
	A(0,0) = db1dx*T[0](0) + db2dx*T[1](0) - (db1dx+db2dx)*T[2](0);
	A(0,1) = db1dy*T[0](0) + db2dy*T[1](0) - (db1dy+db2dy)*T[2](0);
	A(1,0) = db1dx*T[0](1) + db2dx*T[1](1) - (db1dx+db2dx)*T[2](1);
	A(1,1) = db1dy*T[0](1) + db2dy*T[1](1) - (db1dy+db2dy)*T[2](1);
	
	tmp.setLinearNature(A);

	// add max norm of enclosing cell
	double tmp_norm, max_norm = T[0].norm();
	tmp_norm = T[1].norm();
	if (tmp_norm > max_norm)
	  max_norm = tmp_norm;
	tmp_norm = T[2].norm();
	if (tmp_norm > max_norm)
	  max_norm = tmp_norm;
	tmp.setMaxNorm(max_norm);
	
	result.push_back(tmp);
      } 
      else {          // tensor case
	FAMSingularPoint tmp(pos, FAMSingularPoint::DEGENERATE, FAMSingularPoint::NONE, 1);
      
	double a, b, c, d, dbdx[3], dbdy[3];
	derivBaryCoord(dbdx,dbdy);
	a = 0.5*dbdx[0]*(T[0](0,0)-T[0](1,1))
	  + 0.5*dbdx[1]*(T[1](0,0)-T[1](1,1))
	  + 0.5*dbdx[2]*(T[2](0,0)-T[2](1,1));
	b = 0.5*dbdy[0]*(T[0](0,0)-T[0](1,1))
	  + 0.5*dbdy[1]*(T[1](0,0)-T[1](1,1))
	  + 0.5*dbdy[2]*(T[2](0,0)-T[2](1,1));
	c = 0.5*dbdx[0]*T[0](0,1)
	  + 0.5*dbdx[1]*T[1](0,1)
	  + 0.5*dbdx[2]*T[2](0,1);
	d = 0.5*dbdy[0]*T[0](0,1)
	  + 0.5*dbdy[1]*T[1](0,1)
	  + 0.5*dbdy[2]*T[2](0,1);
	
	tmp.setDegenerateType(a,b,c,d, this);

	result.push_back(tmp);
      }
    }
  }
  catch (FException& e) {
    e.addTraceMessage("void FLineCell2D::getZeros(list<FAMSingularPoint>& result) const");
    throw;
  }
#endif
}

//---------------------------------------------------------------------------

positive FLineCell2D::memSize() const
{
  return
    (
     +2//lastPos
    ) *2*sizeof(double) 
    +
    (tensorData?
     tensors[0].size()*sizeof(double)*
     sizeOfCellType() //tensorData
     :0)
    +
    sizeof(*this);
}
//----------------------------------------------------------------------------
