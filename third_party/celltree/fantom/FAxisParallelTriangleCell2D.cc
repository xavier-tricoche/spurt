//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAxisParallelTriangleCell2D.cc,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:00 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#include "FAxisParallelTriangleCell2D.hh"
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
const FCell::geoInfo FAxisParallelTriangleCell2D::myGeoDescription =
  {
    // dimension
    2,
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

FAxisParallelTriangleCell2D::FAxisParallelTriangleCell2D(void) 
  : FCell2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{}

//--------------------------------------------------------------------------- 

FAxisParallelTriangleCell2D::
FAxisParallelTriangleCell2D(const vector<FIndex>& vertIds) 
  : FCell2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  for (positive i=0 ; i<3 ; i++)
    vertexIndices[i]=vertIds[i];

  // is the triangle above the rectangle's splitting diagonal?
  if ( (positive)(vertIds[1]-vertIds[0]) == 1)
    up = false;
  else
    up = true;
}

FAxisParallelTriangleCell2D::
FAxisParallelTriangleCell2D(const FIndex*vertIds) 
  : FCell2D(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  for (positive i=0 ; i<3 ; i++)
    vertexIndices[i]=vertIds[i];

  // is the triangle above the rectangle's splitting diagonal?
  if ( (positive)(vertIds[1]-vertIds[0]) == 1)
    up = false;
  else
    up = true;
}

//--------------------------------------------------------------------------- 

FAxisParallelTriangleCell2D::~FAxisParallelTriangleCell2D()
{}

//--------------------------------------------------------------------------- 

FCell*  FAxisParallelTriangleCell2D::getClone() const
{
  try 
    {
      return new FAxisParallelTriangleCell2D(vertexIndices);
    }
  catch (FException& e) 
    {
      e.addTraceMessage("FCell* FAxisParallelTriangleCell2D::getClone() const");
      throw;
      return (FCell *) 0;
    }
}

//--------------------------------------------------------------------------- 

positive FAxisParallelTriangleCell2D::sizeOfCellType() const
{
  return 3;
}

//--------------------------------------------------------------------------- 

void FAxisParallelTriangleCell2D::interpolate(FTensor& result, 
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
      // local coordinates
      double u, v;
      // barycentric coordinates
      double b0, b1, b2;
      if (up) {
	u = (position[0]-positions[2][0]) / (positions[1][0]-positions[2][0]);
	v = (position[1]-positions[0][1]) / (positions[2][1]-positions[0][1]);
	b0 = 1.-v;
	b1 = u+v-1.;
	b2 = 1.-u;
      }
      else {
	u = (position[0]-positions[0][0]) / (positions[1][0]-positions[0][0]);
	v = (position[1]-positions[0][1]) / (positions[2][1]-positions[0][1]);
	b0 = 1.-u-v;
	b1 = u;
	b2 = v;
      }

      // TEST
      if (b0 < -epsilon || b0 > 1.+epsilon || 
	  b1 < -epsilon || b1 > 1.+epsilon ||
	  b2 < -epsilon || b2 > 1.+epsilon) {
	cout << "WARNING: asked interpolate for a position that lies outside"
	     << endl;
	cout << "local coordinates are: (" << u << ", " << v << ")" << endl;
      }

      result = b0*tensors[0] + b1*tensors[1] + b2*tensors[2];
    }
  }
  catch (FException& e) {
    e.addTraceMessage("void FAxisParallelTriangleCell2D::interpolate(FTensor& result, const FPosition& position)");
    throw e;
  }
}

//--------------------------------------------------------------------------- 

void FAxisParallelTriangleCell2D::derivatives(FTensor& result, 
					      const FPosition& /*position*/) const
{ 
  try {
    if (!positions || !tensors) {
      THROW_EXCEPTION(FException, "ERROR: missing position and/or tensor information");
    }

    if (set_to_zero) {
      result = FTensor(tensors[0].getDimension(), 
		       tensors[0].getOrder()+1);
      result = 0.;
    }
    else {
      // local coordinates
      double dudx, dvdy;
      // barycentric coordinates
      double db0[2], db1[2], db2[2];
      if (up) {
	dudx = 1. / (positions[1][0]-positions[2][0]);
	dvdy = 1. / (positions[2][1]-positions[0][1]);
	db0[0] = 0.; db0[1] = -dvdy;
	db1[0] = dudx; db1[1] = dvdy;
	db2[0] = -dudx; db2[1] = 0.;
      }
      else {
	dudx = 1. / (positions[1][0]-positions[0][0]);
	dvdy = 1. / (positions[2][1]-positions[0][1]);
	db0[0] = -dudx; db0[1] = -dvdy;
	db1[0] = dudx; db1[1] = 0.;
	db2[0] = 0.; db2[1] = dvdy;
      }

      vector< double > comp; // result's components
      vector< double > cell_comp[3]; // components of the cell's tensors
      for (positive i=0 ; i<3 ; i++)
	tensors[i].getValues(cell_comp[i]);
      
      positive sz, derivSz;
      sz = cell_comp[0].size();
      derivSz = 2*sz;
      comp.resize(derivSz);

      for (positive i=0 ; i<sz ; i++) {
	comp[i] = db0[0]*cell_comp[0][i] + db1[0] * cell_comp[1][i] + 
	  db2[0] * cell_comp[2][i];
	comp[sz+i] = db0[1]*cell_comp[0][i] + db1[1] * cell_comp[1][i] + 
	  db2[1] * cell_comp[2][i];
      }

      result = FTensor( 2, tensors[0].getOrder()+1,
			comp);
    }
  }
  catch (FException& e) {
    e.addTraceMessage("void FAxisParallelTriangleCell2D::derivatives(FTensor& result, const FPosition& position)");
    throw;
  }  
}

//--------------------------------------------------------------------------- 

bool FAxisParallelTriangleCell2D::isInside(const FPosition& position) const
{
  try {
    if (!positions) {
      FException e("ERROR: position information missing");
      throw e;
    }

    if (position.getDimension() != 2)
      throw FInvalidDimensionException();
    
    // local coordinates
    double u, v;
    u = (position[0]-positions[2][0]) / (positions[1][0]-positions[2][0]);
    v = (position[1]-positions[0][1]) / (positions[2][1]-positions[0][1]);
    
    if (u > -epsilon && u < 1.+epsilon && 
	v > -epsilon && v < 1.+epsilon &&
	((up && (u+v > 1.-epsilon)) ||
	 (!up && (u+v < 1.+epsilon))))
	return true;
  }
  catch (FException& e) {
    e.addTraceMessage("bool FAxisParallelTriangleCell2D::isInside(const FPosition& position) const");
    throw;
  }
  
  return false;
}

//--------------------------------------------------------------------------- 

void FAxisParallelTriangleCell2D::getZeros(list<FAMSingularPoint>& result) 
const
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

    // compute singularity position in barycentric coordinates 
    double b1, b2, denom;
    FTensor *T = this->tensors; // simplify notations
    
    // solve (b0*T0+b1*T1+(1-b0-b1)*T2) = 0 in b0, b1.
    if (T[0].getOrder() == 1) {    // vector case
      b1    = (-T[2](0) * (T[1](1)-T[2](1))) - (-T[2](1) * (T[1](0)-T[2](0)));
      b2    = (-T[2](1) * (T[0](0)-T[2](0))) - (-T[2](0) * (T[0](1)-T[2](1)));
      denom = ((T[0](0)-T[2](0)) * (T[1](1)-T[2](1))) - 
	((T[0](1)-T[2](1)) * (T[1](0)-T[2](0)));
    } 
    else {                       // tensor case
      b1    = (T[2](1,1)-T[2](0,0))*(T[1](0,1)-T[2](0,1))
	- (-T[2](0,1))*(T[1](0,0)-T[1](1,1)-T[2](0,0)+T[2](1,1));
      b2    = (T[0](0,0)-T[0](1,1)-T[2](0,0)+T[2](1,1))*(-T[2](0,1))
	- (T[0](0,1)-T[2](0,1))*(T[2](1,1)-T[2](0,0));
      denom = (T[0](0,0)-T[0](1,1)-T[2](0,0)+T[2](1,1))*(T[1](0,1)-T[2](0,1))
	- (T[1](0,0)-T[1](1,1)-T[2](0,0)+T[2](1,1))*(T[0](0,1)-T[2](0,1));
    }
  
    if (denom != 0.) {
      denom = 1. / denom;
      b1 *= denom;
      b2 *= denom;
    }
    else
      b1 = -1000.;

    if (b1 > -epsilon && b2 > -epsilon && 1.-b1-b2 > -epsilon) {
      // singularity found lies inside the cell
      FPosition pos = b1 * positions[0] + b2 * positions[1] +
	(1.-b1-b2) * positions[2];
      FAMSingularPoint tmp;

      FTensor dummy;
      interpolate(dummy, pos);
      cout << "interpolate at singular point" << dummy.norm() << endl;

      if (T[0].getOrder() == 1) { // vector case
	tmp = FAMSingularPoint(pos, FAMSingularPoint::NONE, 1);
	FTensor derivs;
	derivatives(derivs, pos);
	FMatrix A(2,2);
	for (positive i=0 ; i<2 ; i++)
	  for (positive j=0 ; j<2 ; j++)
	    A(i,j) = derivs(i,j);
	tmp.setLinearNature(A);
	result.push_back(tmp);
      }
      else { // tensor case
	FAMSingularPoint tmp(pos, FAMSingularPoint::DEGENERATE, FAMSingularPoint::NONE, 1);
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
  catch (FException& e) {
    e.addTraceMessage("void FAxisParallelTriangleCell2D::getZeros(list<*FSingularPoint>& result)");
    throw;
  }
}

//--------------------------------------------------------------------------- 

void FAxisParallelTriangleCell2D::buildBoundingBox(void) const {
  
  if (geometryOK) 
    return;

  if (!positions) {
    FException e("ERROR: position information missing");
    throw e;
  }

  try {
    bBox = FBoundingBox(positions[2][0], positions[0][1], 
			positions[1][0], positions[2][1]);
  }
  catch (FException& e) {
    e.addTraceMessage("void FAxisParallelTriangleCell2D::buildBoundingBox(void)");
    throw;
  }
  
  // no additional cached information depends on geometry
  geometryOK = true;
}

//--------------------------------------------------------------------------- 

FCell::CellType FAxisParallelTriangleCell2D::getCellType(void) const
{
  return FCell::AXIS_PARALLEL_TRI_2D;
}

//===========================================================================
positive FAxisParallelTriangleCell2D::memSize() const
{
  return
    (tensorData?
     tensors[0].size()*sizeof(double)*
     sizeOfCellType() //tensorData
     :0)
    +
    sizeof(*this);
}
