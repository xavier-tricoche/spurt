//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FQuadrilateralCell2D.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:29 $
// Author:    $Author: garth $
// Version:   $Revision: 1.21 $
//
//--------------------------------------------------------------------------- 

#include "FQuadrilateralCell2D.hh"
#include "FTensorSet.hh"
#include "FPositionSet.hh"
#include "FTriangleCell2D.hh"
#include "FException.hh"
#include "FMath.hh"
#include "FMatrix.hh"
#include <complex>
#include <utility>

#ifdef DEBUG_RESOLUTION
#include <iomanip>
#endif

#ifndef NDEBUG
#include <iostream>
#include "eassert.hh"
#endif

using namespace std;

//---------------------------------------------------------------------------
// the cell's geometry description
const FCell::geoInfo FQuadrilateralCell2D::myGeoDescription =
  {
	// dimension
	2,
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

FQuadrilateralCell2D::FQuadrilateralCell2D(void) 
  : FCell2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
}

//---------------------------------------------------------------------------

FQuadrilateralCell2D::FQuadrilateralCell2D( const vector< FIndex >& vert1 ) 
  : FCell2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  for (char i=0; i<4; i++)
    vertexIndices[i]= vert1[i];

}

FQuadrilateralCell2D::FQuadrilateralCell2D( const FIndex*vert1 ) 
  : FCell2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  for (char i=0; i<4; i++)
    vertexIndices[i]= vert1[i];

}

//---------------------------------------------------------------------------

FQuadrilateralCell2D::~FQuadrilateralCell2D()
{
}

//---------------------------------------------------------------------------

FCell*  FQuadrilateralCell2D::getClone() const
{
  try 
    {
      return new FQuadrilateralCell2D( vertexIndices );
    }
  catch (FException& e) 
    {
      e.addTraceMessage("FCell*  FQuadrilateralCell2D::getClone() const");
      throw;
      return (FCell *) 0;
    }
}

//---------------------------------------------------------------------------

positive FQuadrilateralCell2D::sizeOfCellType() const
{
  return 4;
}

//--------------------------------------------------------------------------- 

void FQuadrilateralCell2D::buildBoundingBox(void) const
{
  try{
#ifndef NODEBUG
    if(!positions[0].size())
      throw FException("positions not set");
#endif
    if (geometryOK)
      return;

    FPosition & p0 = positions[0];

    if(getDimension()==2)
      bBox.setBoundingBox( p0[0],p0[1], p0[0],p0[1] );
    else
      bBox.setBoundingBox( p0[0],p0[1],p0[2], p0[0],p0[1],p0[2] );
  
    for ( positive i = 1; i<getNumVerts(); i++ )
      bBox.resize( positions[i]);

    geometryOK = true;
  }

  catch(FException&e){
    e.addTraceMessage("void FQuadrilateralCell2D::buildBoundingBox(void) const");
    throw e;
  }
}

//---------------------------------------------------------------------------

void FQuadrilateralCell2D::interpolate(FTensor& result, 
                                       const FPosition& position) const
{
  // No precomputation is required for the interpolation.
  // Therefore, positionOK and interpolationOK are of no use here.

  try {
#ifndef NODEBUG
    if (!positions[0].size() || !tensors[0].size()) {
      FException e("ERROR: missing position and/or tensor information");
      throw e;
    }
#endif    
    if (set_to_zero) {
      result = FTensor(tensors[0].getDimension(), 
		       tensors[0].getOrder());
      result = 0.;
    }

    // local coordinates' computation to proceed interpolation
    double local[2]; // local coordinates
    double pos[2] = {position[0], position[1]};
    computeLocalCoord(local, pos);

    // interpolate in local coordinates
    result = 
      (1.-local[1]) * ((1.-local[0]) * tensors[0] + local[0] * tensors[1]) +
      local[1] * ((1.-local[0]) * tensors[3] + local[0] * tensors[2]);
  }
  catch (FException& e) {
    e.addTraceMessage("void FQuadrilateralCell2D::interpolate(FTensor& result, const FPosition& position) const");
    throw;
  }
}

//---------------------------------------------------------------------------

void FQuadrilateralCell2D::derivatives(FTensor& result, 
                                       const FPosition& position) const
{ 
  // No precomputation is required for the derivation.
  // Therefore, positionOK and interpolationOK are of no use here.

  try {
#ifndef NODEBUG
    if (!positions[0].size() || !tensors[0].size()) {
      FException e("ERROR: missing position and/or tensor information");
      throw e;
    }
#endif
    double local[2];
    double pos[2] = {position[0], position[1]};
    double deriv[4];

    vector<double> comps[4];
    for (int i=0; i<4; i++)
      comps[i].resize(sz);

    vector<double> comp(2*sz);
    
    // computation of local coordinates
    computeLocalCoord(local, pos);
    computeLocalCoordDeriv(deriv, local);
    for (positive i=0 ; i<4 ; i++)
      tensors[i].getValues(comps[i]);

    for (positive i=0; i<sz ; i++) {
      comp[i] = deriv[0]*(comps[1][i]-comps[0][i]) + 
	deriv[1]*(comps[3][i]-comps[0][i]) +
	(deriv[0]*local[1]+deriv[1]*local[0])*
	(comps[0][i]-comps[1][i]+comps[2][i]-comps[3][i]);
      comp[i+sz] = deriv[2]*(comps[1][i]-comps[0][i]) + 
	deriv[3]*(comps[3][i]-comps[0][i]) +
	(deriv[2]*local[1]+deriv[3]*local[0])*
	(comps[0][i]-comps[1][i]+comps[2][i]-comps[3][i]);
    }
      
    result = FTensor(2, order+1, comp);
  }
  catch (FException& e) {
    e.addTraceMessage("void FQuadrilateralCell2D::derivatives(FTensor& result, const FPosition& position) const");
    throw;
  }
}

// ------------------------------------------------------------------------

bool FQuadrilateralCell2D::isInside(const FPosition& position) const
{
  try {

    buildBoundingBox();
    
    if (!bBox.isInside(position))
      return false;

    // quadrilateral decomposition into 2 triangles
    return ( isInsideSubTriangle(positions[0], positions[1], 
				 positions[3], position) ||
	     isInsideSubTriangle(positions[1], positions[2],
				 positions[3], position) );
  }
  catch (FException& e) {
    e.addTraceMessage("bool FQuadrilateralCell2D::isInside(const FPosition& position) const");
    throw;
  }
  
  // for the compiler...
  return false;
}

// -------------------------------------------------------------------------

bool FQuadrilateralCell2D::isInside( const FPosition& pos,
 				     const vector< FPosition >& vertices)
{
  double denom, alpha, b0, b1;

  // 1st subtriangle
  denom = (vertices[0][0]-vertices[3][0])*(vertices[1][1]-vertices[3][1]) -
    (vertices[0][1]-vertices[3][1])*(vertices[1][0]-vertices[3][0]);

  alpha = - denom * epsilon;
  
  b0 = (pos[0]-vertices[3][0])*(vertices[1][1]-vertices[3][1]) -
    (pos[1]-vertices[3][1])*(vertices[1][0]-vertices[3][0]);

  if ( (alpha < 0 && b0 > alpha) || 
       (alpha > 0 && b0 < alpha) ) {
      
      b1 = (vertices[0][0]-vertices[3][0])*(pos[1]-vertices[3][1]) -
	(vertices[0][1]-vertices[3][1])*(pos[0]-vertices[3][0]);
      
      if ( (alpha < 0 && b1 > alpha && denom - b0 - b1 > alpha) ||
	   (alpha > 0 && b1 < alpha && denom - b0 - b1 < alpha) )
	return true;
  }  
  

  // 2nd subtriangle
  denom = (vertices[1][0]-vertices[3][0])*(vertices[2][1]-vertices[3][1]) -
    (vertices[1][1]-vertices[3][1])*(vertices[2][0]-vertices[3][0]);
  
  alpha = - denom * epsilon;
  
  b0 = (pos[0]-vertices[3][0])*(vertices[2][1]-vertices[3][1]) -
    (pos[1]-vertices[3][1])*(vertices[2][0]-vertices[3][0]);
  
  if ( (alpha < 0 && b0 > alpha) ||
       (alpha > 0 && b0 < alpha) ) {
    
    b1 = (vertices[1][0]-vertices[3][0])*(pos[1]-vertices[3][1]) -
      (vertices[1][1]-vertices[3][1])*(pos[0]-vertices[3][0]);
    
    if ( (alpha < 0 && b1 > alpha && denom - b0 - b1 > alpha) ||
	 (alpha > 0 && b1 < alpha && denom - b0 - b1 < alpha) )
      return true;
  }
  
  return false;
}

// ------------------------------------------------------------------------

bool FQuadrilateralCell2D::isInsideSubTriangle(const FPosition& pos1, 
                                               const FPosition& pos2, 
                                               const FPosition& pos3,
                                               const FPosition& posX) const 
{
  // barycentric coordinates
  double denom = (pos1[0]-pos3[0])*(pos2[1]-pos3[1]) -
    (pos1[1]-pos3[1])*(pos2[0]-pos3[0]);

  double alpha = - denom * epsilon;
  
  double b0 = (posX[0]-pos3[0])*(pos2[1]-pos3[1]) -
    (posX[1]-pos3[1])*(pos2[0]-pos3[0]); 

  if ( (alpha < 0 && b0 > alpha) || 
       (alpha > 0 && b0 < alpha) ) {
    
    double b1 = (pos1[0]-pos3[0])*(posX[1]-pos3[1]) -
      (pos1[1]-pos3[1])*(posX[0]-pos3[0]);
    
    if ( (alpha < 0 && b1 > alpha && denom - b0 - b1 > alpha) ||
	 (alpha > 0 && b1 < alpha && denom - b0 - b1 < alpha) ) {

      return true;
    }
  }

  return false;
}

// ------------------------------------------------------------------------

void FQuadrilateralCell2D::getZeros(list<FAMSingularPoint>& result) const
{ 
  // no precomputation is required by this function

  try {
#ifndef NODEBUG
    if (!positions[0].size() || !tensors[0].size()) {
      FException e("ERROR: missing position and/or tensor information");
      throw e;
    }
#endif
    // if cell has been set to zero, we're done
    if (set_to_zero) {
      return;
    }

    if (tensors[0].getDimension() != 2 || 
	(tensors[0].getOrder() != 1 && tensors[0].getOrder() != 2)) {
      FInvalidDimensionException e;
      throw e;
    }
  
    // first compute zero position in local coordinates
    double a, b, c; // coefficients of the quadratic equation to solve
    double local[2];
    double d1, d2; 
    FAMSingularPoint sing;
    FTensor sum = tensors[0]-tensors[1]+tensors[2]-tensors[3];
    a = sum(0)*(tensors[3]-tensors[0])(1) -
      sum(1)*(tensors[3]-tensors[0])(0);
    b = (tensors[1]-tensors[0])(0)*(tensors[3]-tensors[0])(1) -
      (tensors[1]-tensors[0])(1)*(tensors[3]-tensors[0])(0) +
      sum(0)*tensors[0](1) -
      sum(1)*tensors[0](0);;
    c = (tensors[1]-tensors[0])(0)*tensors[0](1) -
      (tensors[1]-tensors[0])(1)*tensors[0](0);
  
    complex<double> roots[2];
    int nbRoots = FMath::QuadraticEquation(a, b, c, roots);
    
    for (int i=0 ; i<nbRoots ; i++) {
      
      local[1] = -1.; // dummy
      
      if (roots[i].imag() == 0. &&  
	  roots[i].real() > -epsilon && 
	  roots[i].real() < 1.+epsilon) {
	
	local[1] = roots[i].real() ;
	d1 = (tensors[1](0) - tensors[0](0)) + local[1]*sum(0);
	d2 = (tensors[1](1) - tensors[0](1)) + local[1]*sum(1);
	
	// check for numerical stability
	if (fabs(d1) > fabs(d2)) {
	  
#ifndef NODEBUG
	  if (fabs(d1) <= epsilon)
	    cout << "in FQuadrilateralCell2D::getZeros: fabs(d1) = "
		 << fabs(d1) << " (< epsilon = " << epsilon << ")" 
		 << endl;
#endif      
	  
	  local[0] = - ( tensors[0](0) + 
			 local[1] * (tensors[3](0)-tensors[0](0)) ) / d1;
	}
	else {
	  
#ifndef NODEBUG
	  if (fabs(d2) <= epsilon)
	    cout << "in FQuadrilateralCell2D::getZeros: fabs(d2) = "
		 << fabs(d2) << " (< epsilon = " << epsilon << ")" 
		 << endl;
#endif   
	  
	  local[0] = - ( tensors[0](1) + 
			 local[1] * (tensors[3](1)-tensors[0](1)) ) / d2;
	}
	
	if (local[0] > -epsilon && local[0] < 1.+epsilon) {

	  // convert position to global coordinates
	  FPosition pos = positions[0] +
	    local[0] * (positions[1]-positions[0]) +
	    local[1] * (positions[3]-positions[0]) +
	    local[0] * local[1] * 
	    (positions[0]-positions[1]+positions[2]-positions[3]);

// 	  cout << "checking zero location: ";
// 	  FTensor dum;
// 	  interpolate(dum, pos);
// 	  cout << dum << endl; 
	  
	  sing = FAMSingularPoint(pos, FAMSingularPoint::NONE, 1);
	  
	  // singularity type computation
	  FMatrix A(2,2);
	  double deriv[4];
	  // compute local coordinates' derivatives
	  computeLocalCoordDeriv(deriv, local);
	  A(0,0) = deriv[0]*(tensors[1](0) - tensors[0](0)) +
	    deriv[1]*(tensors[3](0) - tensors[0](0)) +
	    (deriv[0]*local[1] + deriv[1]*local[0])*
	    (tensors[0](0)-tensors[1](0)+tensors[2](0)-tensors[3](0));
	  A(0,1) = deriv[0]*(tensors[1](1) - tensors[0](1)) +
	    deriv[1]*(tensors[3](1) - tensors[0](1)) +
	    (deriv[0]*local[1] + deriv[1]*local[0])*
	    (tensors[0](1)-tensors[1](1)+tensors[2](1)-tensors[3](1));
	  A(1,0) = deriv[2]*(tensors[1](0) - tensors[0](0)) +
	    deriv[3]*(tensors[3](0) - tensors[0](0)) +
	    (deriv[2]*local[1] + deriv[3]*local[0])*
	    (tensors[0](0)-tensors[1](0)+tensors[2](0)-tensors[3](0));
	  A(1,1) = deriv[2]*(tensors[1](1) - tensors[0](1)) +
	    deriv[3]*(tensors[3](1) - tensors[0](1)) +
	    (deriv[2]*local[1] + deriv[3]*local[0])*
	    (tensors[0](1)-tensors[1](1)+tensors[2](1)-tensors[3](1));

	    sing.setLinearNature(A);
	    result.push_back(sing);
	}
      }
    }
  }
  catch (FException& e) {
    e.addTraceMessage("void FQuadrilateralCell2D::getZeros(list<FAMSingularPoint>& result) const");
    throw;
  }
}

// ------------------------------------------------------------------------

FCell::CellType FQuadrilateralCell2D::getCellType(void) const
{
  return FCell::QUADRILATERAL_2D;
}

// ------------------------------------------------------------------------

void FQuadrilateralCell2D::localCoordToInitialize(void) const
{

  if (localCoordInitialized && geometryOK)
    return;

  p1[0] = positions[1][0] - positions[0][0];
  p2[0] = positions[2][0] - positions[0][0];
  p3[0] = positions[3][0] - positions[0][0];
  p4[0] = -p1[0] + p2[0] - p3[0];

  p1[1] = positions[1][1] - positions[0][1];
  p2[1] = positions[2][1] - positions[0][1];
  p3[1] = positions[3][1] - positions[0][1];
  p4[1] = -p1[1] + p2[1] - p3[1];

  // precomputation of useful determinants
  deter1 = p3[0]*p4[1] - p3[1]*p4[0]; // |P0P3 Delta|
  deter2 = p3[0]*p1[1] - p3[1]*p1[0]; // |P0P3 P0P1|
  deter3 = p1[0]*p4[1] - p1[1]*p4[0]; // |P0P1 Delta|

  localCoordInitialized = true;

  if (!geometryOK)
    buildBoundingBox(); // afterwards, geometryOK = true.
}

// ------------------------------------------------------------------------

void FQuadrilateralCell2D::computeLocalCoord(double *local, 
					     double *pos) const
{
  localCoordToInitialize();
  
  double a, b, c; // coefficients of the quadratic equation to solve
  double p[2]; // relative position of current position with respect to P0
  // the quadratic equation in s is:
  // |P0P3 Delta| s*s + (|P0P3 P0P1| + |Delta P0P1|) s + |P0P1 P0P| = 0.
  p[0] = pos[0] - positions[0][0];
  p[1] = pos[1] - positions[0][1];
  a = deter1;
  b = deter2 + p4[0]*p[1] - p4[1]*p[0];
  c = p1[0]*p[1] - p1[1]*p[0];
  
  complex<double> roots[2];
  positive nbRoots = FMath::QuadraticEquation(a, b, c, roots);
  local[1] = -1.; // dummy
  

  for (positive i=0 ; i<nbRoots ; i++) {

#ifndef NODEBUG
    eassert (roots[i].imag() == 0.); // the cell is convex
#endif
    
    // increase approximation tolerance to avoid inconsistency 
    // with triangle based isInside test
    if (roots[i].real() > -3.*epsilon && 
	roots[i].real() < 1.+ 3.*epsilon) { 
      local[1] = roots[i].real();
      break;
    }
  }

#ifndef NODEBUG
  if (local[1] == -1.) {
//     cout << "position " << "(" << pos[0] << ", " << pos[1] << ")"
// 	 << " lies outside the cell " 
// 	 << "in FQuad2D::computeLocalCoord"
// 	 << endl;
//     cout << "current cell has positions: "
// 	 << "#1 : " << positions[0] 
// 	 << ", #2 : " << positions[1] 
// 	 << ", #3 : " << positions[2]
// 	 << ", #4 : " << positions[3]
// 	 << endl;
//     cout << "... and this is quite surprising since the isInside function "
// 	 << "answered : " << isInside(FPosition(pos[0], pos[1]))
// 	 << endl;

//     cout << "roots of quadratic equation were: " << endl;
//     for (positive i=0 ; i<nbRoots ; i++)
//       cout << roots[i] << endl;
//     for (positive i=0 ; i<nbRoots ; i++)
//       cout << "1.-roots[i].real() = " << 1. - roots[i].real() << endl;

//     cout << endl;
    
#ifdef DEBUG_RESOLUTION 
    std::cerr << setprecision(12);
    std::cerr << "a=" << a << " b=" << b << " c=" << c << " roots=" << roots[0] << "; " << roots[1] << std::endl;
    std::cerr << "pos=" << pos[0] << ", " << pos[1] << std::endl;
#endif
    throw FException("illegal local coord computation in current cell");
  }
  // eassert (local[1] != -1.); // the position lies inside the quadrilateral
#endif

  double d1, d2; 
  d1 = p1[0] + local[1]*p4[0];
  d2 = p1[1] + local[1]*p4[1];

  // check for numerical stability
  if (fabs(d1) > fabs(d2)) {

#ifndef NODEBUG
    eassert (fabs(d1) >= epsilon); // degenerate case
#endif      

    local[0] = (p[0] - local[1]*p3[0]) / d1;
  }
  else {

#ifndef NODEBUG
    eassert (fabs(d2) >= epsilon); // degenerate case
#endif   

    local[0] = (p[1] - local[1]*p3[1]) / d2;
  }


//   cout << "checking local coordinates' computation:" << endl
//        << "pos = (" << pos[0] << ", " << pos[1] << ")" << ", phi(pos) = "
//        << positions[0] + local[0] * (positions[1]-positions[0]) +
//     local[1] * (positions[3]-positions[0]) +
//     local[0]*local[1] * 
//     (positions[0] - positions[1] + positions[2] - positions[3])
//        << endl;
}

// ------------------------------------------------------------------------

void FQuadrilateralCell2D::computeLocalCoordDeriv(double *deriv, 
						  double *local) const
{
  localCoordToInitialize();
  
  double tmp1 = 2.*deter3*local[0] + p4[0]*local[1] - p4[1]*local[0] - deter2;
  double tmp2 = 2.*deter1*local[1] + deter2 + p4[0]*local[1] - p4[1]*local[0];

  deriv[0] = (p4[1]*local[0] + p3[1]) / tmp1;  // dr/dx
  deriv[2] = -(p4[0]*local[0] + p3[0]) / tmp1; // dr/dy
  deriv[1] = (p4[1]*local[1] + p1[1]) / tmp2;  // ds/dx
  deriv[3] = -(p4[0]*local[1] + p1[0]) / tmp2; // ds/dy
}

//===========================================================================
positive FQuadrilateralCell2D::memSize() const
{
  return
    (tensorData?
     tensors[0].size()*sizeof(double)*
     sizeOfCellType() //tensorData
     :0)
    +
    sizeof(*this);
}
