//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAxisParallelQuadCell2D.cc,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:00 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#include "FAxisParallelQuadCell2D.hh"
#include "FTensorSet.hh"
#include "FPositionSet.hh"
#include "FException.hh"
#include "FMath.hh"
#include "FMatrix.hh"
#include <complex>
#include <utility>
#include <iostream>

//---------------------------------------------------------------------------

// the cell's geometry description
const FCell::geoInfo FAxisParallelQuadCell2D::myGeoDescription =
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

FAxisParallelQuadCell2D::FAxisParallelQuadCell2D(void) 
  : FCell2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{}

//--------------------------------------------------------------------------- 

FAxisParallelQuadCell2D::
FAxisParallelQuadCell2D(const vector<FIndex>& vertIds) 
  : FCell2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  for (positive i=0 ; i<4 ; i++)
    myIndices[i] = vertIds[i];
}

FAxisParallelQuadCell2D::
FAxisParallelQuadCell2D(const FIndex* vertIds) 
  : FCell2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  for (positive i=0 ; i<4 ; i++)
    myIndices[i] = vertIds[i];
}

//--------------------------------------------------------------------------- 

FAxisParallelQuadCell2D::~FAxisParallelQuadCell2D()
{}

//--------------------------------------------------------------------------- 

FCell*  FAxisParallelQuadCell2D::getClone() const
{
  try 
    {
      return new FAxisParallelQuadCell2D(vertexIndices);
    }
  catch (FException& e) 
    {
      e.addTraceMessage("FCell*  FAxisParallelQuadCell2D::getClone() const");
      throw;
      return (FCell *) 0;
    }
}

//--------------------------------------------------------------------------- 

positive FAxisParallelQuadCell2D::sizeOfCellType() const
{
  return 4;
}

//--------------------------------------------------------------------------- 

void FAxisParallelQuadCell2D::interpolate(FTensor& result, 
					  const FPosition& position) const
{
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
      double u, v;
      u = (position[0]-positions[0][0]) / (positions[1][0]-positions[0][0]);
      v = (position[1]-positions[0][1]) / (positions[3][1]-positions[0][1]);
      
      result = (1.-v) * ((1.-u)*tensors[0] + u*tensors[1]) +
	v * ((1.-u)*tensors[3] + u*tensors[2]);
    }
  }
  catch (FException& e) {
    e.addTraceMessage("void FAxisParallelQuadCell2D::interpolate(FTensor& result, const FPosition& position)");
    throw e;
  }
}

//--------------------------------------------------------------------------- 

void FAxisParallelQuadCell2D::derivatives(FTensor& result, 
					  const FPosition& position) const
{ 
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
      // interpolant derivatives computation
      double u, v;
      u = (position[0]-positions[0][0]) / (positions[1][0]-positions[0][0]);
      v = (position[1]-positions[0][1]) / (positions[3][1]-positions[0][1]);
      double dudx, dvdy;
      dudx = 1. / (positions[1][0] - positions[0][0]);
      dvdy = 1. / (positions[3][1] - positions[0][1]);
      
      vector< double > comp; // result's components
      vector< double > comps[4]; // components of the cell's tensors
      for (positive i=0 ; i<4 ; i++)
	tensors[i].getValues(comps[i]);
      
      positive sz, derivSz;
      sz = comps[0].size();
      derivSz = 2*sz;
      comp.resize(derivSz);

      for (positive i=0 ; i<sz ; i++) {
	comp[i] = 
	  dudx * (comps[1][i]-comps[0][i] +
		  (comps[0][i]-comps[1][i]+comps[2][i]-comps[3][i])*v);
	comp[sz+i] = 
	  dvdy * (comps[3][i]-comps[0][i] +
		  (comps[0][i]-comps[1][i]+comps[2][i]-comps[3][i])*u);
      }

      result = FTensor( tensors[0].getOrder() 
			? tensors[0].getDimension() 
			: comp.size(), 
			tensors[0].getOrder()+1,
			comp);
    }
  }
  catch (FException& e) {
    e.addTraceMessage("void FAxisParallelQuadCell2D::derivatives(FTensor& result, const FPosition& position)");
    throw;
  }  
}

//--------------------------------------------------------------------------- 

bool FAxisParallelQuadCell2D::isInside(const FPosition& position) const
{
//   cout << "getting in isInside()" << endl;
  try {
    if (!positions) {
      FException e("ERROR: position information missing");
      throw e;
    }

    if (position.getDimension() != 2)
      throw FInvalidDimensionException();
    
    return (position[0] > positions[0][0]-epsilon &&
	    position[0] < positions[1][0]+epsilon &&
	    position[1] > positions[0][1]-epsilon &&
	    position[1] < positions[3][1]+epsilon);

  }
  catch (FException& e) {
    e.addTraceMessage("bool FAxisParallelQuadCell2D::isInside(const FPosition& position) const");
    throw;
  }
  
  // for the compiler...
  return false;
}

//--------------------------------------------------------------------------- 

void FAxisParallelQuadCell2D::getZeros(list<FAMSingularPoint>& result) const
{ 
  try {

    if (!positions || !tensors) {
      FException e("ERROR: missing position and/or tensor information");
      throw e;
    }

    // if cell has been set to zero, we're done
    if (set_to_zero)
      {
	return;
      }

    // check applicability of the function
    if (tensors[0].getDimension() != 2 || 
	(tensors[0].getOrder() != 1 && tensors[0].getOrder() != 2) ) {
      FInvalidDimensionException e;
      throw e;
    }

    double alpha, beta, gamma, u, v;

    if (tensors[0].getOrder() == 1) { // VECTORS
    
      // interpolation parameters (in local coordinates)
      double a[2], b[2], c[2], d[2];
      a[0] = tensors[0](0);
      a[1] = tensors[0](1); 
      b[0] = tensors[1](0)-tensors[0](0);
      b[1] = tensors[1](1)-tensors[0](1); 
      c[0] = tensors[3](0)-tensors[0](0);
      c[1] = tensors[3](1)-tensors[0](1); 
      d[0] = tensors[0](0)-tensors[1](0)+tensors[2](0)-tensors[3](0);
      d[1] = tensors[0](1)-tensors[1](1)+tensors[2](1)-tensors[3](1);   

      gamma = a[1]*c[0] - c[1]*a[0];
      beta  = b[1]*c[0] + a[1]*d[0] - c[1]*b[0] - a[0]*d[1];
      alpha = b[1]*d[0] - d[1]*b[0];
    
      complex<double> roots[2];
      FAMSingularPoint tmp;
      int nbRoots = FMath::QuadraticEquation(alpha, beta, gamma, roots);

      for (int i=0 ; i<nbRoots ; i++) {

	if (!roots[i].imag()) {
	  // real root found

	  if (c[0]+d[0]*roots[i].real() != 0.) {
	    u = roots[i].real();
	    v = -(a[0]+b[0]*u) / (c[0]+d[0]*u);

	    if (u>-epsilon && u<1.+epsilon && v>-epsilon && v<1.+epsilon) {

	      // compute physical coordinates
	      double x = (1.-u) * positions[0][0] + u * positions[1][0];
	      double y = (1.-v) * positions[0][1] + v * positions[3][1];
	      FPosition pos(x,y);

	      tmp = FAMSingularPoint(pos, FAMSingularPoint::NONE, 1);
	    
	      // singularity type computation
	      FMatrix A(2,2);
	      FTensor derivs;
	      derivatives(derivs, pos);
	      for (positive i=0 ; i<2 ; i++)
		for (positive j=0 ; j<2 ; j++)
		  A(i,j) = derivs(i,j);
	      tmp.setLinearNature(A);
	      result.push_back(tmp);
	    }
	  }
	
	}
      }
    } 
    else { // TENSORS

      // interpolation parameters in local coordinates
      double a[4], b[4], c[4], d[4];

      a[0] = tensors[0](0,0);
      a[1] = tensors[0](0,1);
      a[2] = tensors[0](1,0);
      a[3] = tensors[0](1,1);

      b[0] = tensors[1](0,0)-tensors[0](0,0);
      b[1] = tensors[1](0,1)-tensors[0](0,1);
      b[2] = tensors[1](1,0)-tensors[0](1,0);
      b[3] = tensors[1](1,1)-tensors[0](1,1);

      c[0] = tensors[3](0,0)-tensors[0](0,0);
      c[1] = tensors[3](0,1)-tensors[0](0,1);
      c[2] = tensors[3](1,0)-tensors[0](1,0);
      c[3] = tensors[3](1,1)-tensors[0](1,1);

      d[0] = tensors[0](0,0)-tensors[1](0,0)+tensors[2](0,0)-tensors[3](0,0);
      d[1] = tensors[0](0,1)-tensors[1](0,1)+tensors[2](0,1)-tensors[3](0,1);
      d[2] = tensors[0](1,0)-tensors[1](1,0)+tensors[2](1,0)-tensors[3](1,0);
      d[3] = tensors[0](1,1)-tensors[1](1,1)+tensors[2](1,1)-tensors[3](1,1);
    
    
      double adiff, bdiff, cdiff, ddiff;
      adiff=a[0]-a[3];
      bdiff=b[0]-b[3];
      cdiff=c[0]-c[3];
      ddiff=d[0]-d[3];
    
      gamma = adiff*c[1] - cdiff*a[1];
      beta  = adiff*d[1] + bdiff*c[1] - cdiff*b[1] -ddiff*a[1];
      alpha = bdiff*d[1] - ddiff*b[1];
    
      complex<double> roots[2];
      double u, v;
      int nbRoots = FMath::QuadraticEquation(alpha, beta, gamma, roots);

      for (int i=0; i<nbRoots; i++) {

	if (!roots[i].imag()) {  //real root found

	  if (c[1] +d[1]*roots[i].real() ) {
	    u = roots[i].real();
	    v = - (a[1]+b[1]*u) / (c[1]+d[1]*u);
	  
	    if (u>-epsilon && u<1.+epsilon && v>-epsilon && v<1.+epsilon) {
	      double x = (1.-u) * positions[0][0] + u * positions[1][0];
	      double y = (1.-v) * positions[0][1] + v * positions[3][1];
	      FPosition pos=FPosition(x,y);
	      FAMSingularPoint tmp(pos, FAMSingularPoint::DEGENERATE, FAMSingularPoint::NONE, 1);

#ifndef NODEBUG
	      FTensor tens = tensors[0] + 
		u * (tensors[1]-tensors[0]) +
		v * (tensors[3]-tensors[0]) + 
		u*v * (tensors[0]-tensors[1]+tensors[2]-tensors[3]);

	      if (fabs( tens(0,0) - tens(1,1) ) > epsilon ||
		  fabs( tens(0,1) ) > epsilon ||
		  fabs( tens(1,0) ) > epsilon )
		cout << "TEST FAILED at: "<< pos << endl
		     << "   with tensor: "<< tens << endl;
#endif

	      FTensor derivs;
	      derivatives(derivs, pos);
	      double aa,bb,cc,dd;
	      aa=0.5*(derivs(0,0,0)-derivs(1,1,0));
	      bb=0.5*(derivs(0,0,1)-derivs(1,1,1));
	      cc=derivs(0,1,0);
	      dd=derivs(0,1,1);
	      tmp.setDegenerateType(aa, bb, cc, dd, this);
	      result.push_back(tmp);
	    }  
	  }
	}
      }
    }
  
  }
  catch (FException& e) {
    e.addTraceMessage("void FAxisParallelQuadCell2D::getZeros(list<*FSingularPoint>& result)");
    throw;
  }
}

//--------------------------------------------------------------------------- 

void FAxisParallelQuadCell2D::buildBoundingBox(void) const {
  
  if (geometryOK)
    return;

  if (!positions) {
    FException e("ERROR: position information missing");
    throw e;
  }

  try {
    bBox = FBoundingBox(positions[0][0], positions[0][1], 
			positions[1][0], positions[3][1]);
  }
  catch (FException& e) {
    e.addTraceMessage("void FAxisParallelQuadCell2D::buildBoundingBox(void)");
    throw;
  }
  
  // no additional cached information depends on geometry
  geometryOK = true;
}

//--------------------------------------------------------------------------- 

FCell::CellType FAxisParallelQuadCell2D::getCellType(void) const
{
  return FCell::AXIS_PARALLEL_QUAD_2D;
}

//===========================================================================
positive FAxisParallelQuadCell2D::memSize() const
{
  return
    (tensorData?
     tensors[0].size()*sizeof(double)*
     sizeOfCellType() //tensorData
     :0)
    +
    sizeof(*this);
}
