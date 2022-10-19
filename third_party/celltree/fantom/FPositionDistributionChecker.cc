//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPositionDistributionChecker.cc,v $
// Language:  C++
// Date:      $Date: 2003/09/10 13:13:47 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//---------------------------------------------------------------------------

#include "FPositionDistributionChecker.hh"
#include "FBoundingBox.hh"
#include "FException.hh"

#include <cmath>
#include <cstdio>
#include <cassert>
#include <cstdlib>
#include <iostream>

using namespace std;

//---------------------------------------------------------------------------

FPositionDistributionChecker::
FPositionDistributionChecker()
  : gCount(0), nx(0), ny(0), nz(0), theBBox(3)
{
}

//---------------------------------------------------------------------------

FPositionDistributionChecker::
~FPositionDistributionChecker()
{
}

//---------------------------------------------------------------------------

void
FPositionDistributionChecker::
initializeDimensions ( positive nx, positive ny, positive nz )
{
  this->nx = nx;
  this->ny = ny;
  this->nz = nz;
  size[0] = nx; size[1] = ny; size[2] = nz;

  // initialize structure information 

  // coordinates' container
  for (positive i=0 ; i<3 ; i++) {
    // at the beginning, we suppose a one-dimensional dependency
    // of the coordinates wrt the corresponding dimension. 
    // Since a dependency in an additional dimension is possible
    // too, we set a priori the index of this addition dimension
    // to a dummy value.
    coord[i].first = 3; // dummy value not in {0, 1, 2}
    coord[i].second.resize(size[i]);
  }

  // dependencies
  for (positive i=0 ; i<3 ; i++)
    for (positive j=0 ; j<3 ; j++)
      if (i==j)
	depend[i][j] = true; // linearity in each dim at the beginning
      else
	depend[i][j] = false;

  // uniformity
  uniform[0] = uniform[1] = uniform[2] = true;
  distance[0] = distance[1] = distance[2] = 0.;
  
  // local coordinates
  discreteCoord[0] = discreteCoord[1] = discreteCoord[2] = 0; 

  gCount = 0;

  have_to_check = true;

}

//---------------------------------------------------------------------------

void
FPositionDistributionChecker::
checkStep( double comp[] )
{

  if (gCount == 0) {

    // initialize coordinates' arrays with first position's coordinates
    for (positive i=0 ; i<3 ; i++) 
      coord[i].second[0] = comp[i];

    gCount++;

    return;
  }

  // only check for bounding box extension if considered position lies
  // on the grid boundary
  positive i,j,k,dim;
  bool check_bb = false;
  if (nz == 1) { // 2D
    
    dim = 2;

    j = (gCount+1) / nx;
    if (j==0 || j==ny-1)
      check_bb = true;
    else {
      i = (gCount+1) - j*nx;
      if (i==0 || i==nx-1)
	check_bb = true;
    }
  }
  else { // 3D

    dim = 3;

    k = (gCount+1) / (nx*ny);
    if (k==0 || k==nz-1)
      check_bb = true;
    else {
      j = ((gCount+1) - k*nx*ny) / nx;
      if (j==0 || j==ny-1)
	check_bb = true;
      else {
	i = (gCount+1) - k*nx*ny - j*nx;
	if (i==0 || i==nx-1)
	  check_bb = true;
      }
    }
  }

  double compD[3];
  compD[0] = comp[0];
  compD[1] = comp[1];
  if ( dim == 3 )
    compD[2] = comp[2];
  
  if (check_bb) 
    theBBox.resize(FPosition(dim, compD));

  // process only if no unstructuredness has been encountered so far
  if (!have_to_check)
    return;


  // update discrete coordinates
  discreteCoord[0]++;
  if (discreteCoord[0] == nx) {
    discreteCoord[0] = 0;
    //    cout << endl << "\tmoving to next line..." << endl;
    discreteCoord[1]++;
    if (discreteCoord[1] == ny) {
      discreteCoord[1] = 0;
      //      cout << endl << "\t moving to next plane..." << endl;
      discreteCoord[2]++;
    }
  }

  // check every single coordinate
  for (positive i=0 ; i<3 ; i++) {
    checkCoord(comp[i], i);
  }

  gCount++;
}

//---------------------------------------------------------------------------

void
FPositionDistributionChecker::
checkCoord(double a_coord, positive dim) {
  
  // determine which value to check by consulting dependencies
  positive checkId, i[4]; // size 4 is a ugly trick... (xa4)
  positive dep2;

  for (positive l=0 ; l<3 ; l++) {
    if (depend[dim][l]) {
      //      cout << "dimension " << dim 
      //	   << " depends on discrete dimension #" << l << endl;
      i[l] = discreteCoord[l];
    }
    else {
      // a 0 is used to express the independency of the coordinate
      // with respect to the considered dimension.
      i[l] = 0; 
      //      cout << "dimension " << dim 
      //	   << " does not depend on discrete dimension #" << l << endl;
    }
  }

  // corresponding index
  checkId = i[2]*nx*ny + i[1]*nx + i[0];
  //  cout << "checkId = " << checkId << " and gCount = " << gCount << endl;

  dep2 = coord[dim].first; // possible second dependency

  // does the reference value exist?
  if ( checkId < gCount ) { // we met this index already

    double reference = coord[dim].second[index(dim, i[dim], dep2, i[dep2])];

    //cout << "reference value = " << reference << endl;
    //cout << "current value = " << a_coord << endl;
    // check the corresponding value
    if ( reference != a_coord ) {

      //      cout << "new value found" << endl;
      
      // update coordinates: we infered a rectilinear structure
      // for the dimension j (such that i[j] has been set to 0).
      // Now, a new value has been encountered and we must 
      // complete a posteriori the contents of the coordinates'
      // array by setting all empty coordinates to the constant
      // value observed in the previous steps.
      for (positive l=0 ; l<3 ; l++) {
	if (l == dim) {
	  
// 	  cout << "dependency wrt its own dimension (" << i[l] << ")" << endl;

	}
	else if (discreteCoord[l] != 0 && i[l] == 0) {

// 	  cout << "new dependency wrt a new dimension (" << l << ")" << endl;

	  // we do not inspect this discrete coordinate for the first
	  // time but get a new value: An additional dependency has
	  // been exhibited.
	  depend[dim][l] = true;
	  if (coord[dim].first < 3) 
	    // a valid additional dimension has been set already
	    // no need to proceed: no optimization will be possible
	    have_to_check = false;
	  else {
	    dep2 = coord[dim].first = l;
	    coord[dim].second.resize(size[dim]*size[l]);

	    // complete the array by setting all missing values
	    // to the old coordinate (that was constant so far)
	    for (positive k=0 ; k<size[l] ; k++)
	      coord[dim].second[index(dim, i[dim], l, k)] = reference; 
	  }
	}
	else {

	  //  cout << "pre-existing dependency (" << l << ")" << endl;

	}
      }

      // insert new value

      //      cout << "insert a new reference value" << endl;

      coord[dim].second[index(dim, i[dim], dep2, discreteCoord[dep2])] = 
	a_coord;
    }
    // otherwise we are done: structure check succeeds
  }
  else {
    
    //    cout << "first time we encounter this position" << endl;

    // update coord
    coord[dim].second[index(dim, i[dim], dep2, i[dep2])] = a_coord;
  }


  return;
}

//---------------------------------------------------------------------------

void
FPositionDistributionChecker::
getDependencies(bool depend[3][3]) const
{
  for (positive i=0 ; i<3 ; i++)
    for (positive j=0 ; j<3 ; j++)
      depend[i][j] = this->depend[i][j]; 
}

//---------------------------------------------------------------------------

void
FPositionDistributionChecker::
getBoundingBox ( FBoundingBox& bbox ) const
{
  bbox = this->theBBox;
}

//---------------------------------------------------------------------------

void
FPositionDistributionChecker::
getCoordinates1D(vector< double >& coord, positive dim) const
{

  if (dim > 2 || this->coord[dim].first < 3 ) {
    
    FException e("ERROR: wrong parameters");
    e.addTraceMessage("void FPositionDistributionChecker::getCoordinates1D(vector< double >& coord, positive dim) const");
    throw e;
  }
  
  coord.resize(size[dim]);
  for (positive i=0 ; i<size[dim] ; i++)
    coord[i] = this->coord[dim].second[i];
}

//---------------------------------------------------------------------------

void
FPositionDistributionChecker::
getCoordinates2D(vector<double> coord[2], positive dim1, positive dim2) const
{
  if (dim1 > 2 || dim2 > 2 ||
      this->coord[dim1].first != dim2 ||
      this->coord[dim2].first != dim1 ||
      dim1 > dim2) {
    
    FException e("ERROR: wrong dependency parameters");
    e.addTraceMessage("void FPositionDistributionChecker::getCoordinates2D(vector<double> coord[2], positive dim1, positive dim2) const");
    throw e;
  }  
  
  coord[0].resize(size[dim1]*size[dim2]);
  coord[1].resize(size[dim1]*size[dim2]);
  for (positive i=0 ; i<size[dim1] ; i++)
    for (positive j=0 ; j<size[dim2] ; j++) {
      coord[0][j*dim1+i] = this->coord[dim1].second[j*size[dim1]+i];
      coord[1][j*dim1+i] = this->coord[dim2].second[i*size[dim2]+j];
    }
}

//---------------------------------------------------------------------------

positive
FPositionDistributionChecker::
index(positive dim1, positive id1, positive dim2, positive id2)
{
  if (dim2 == 3) // dummy value
    return id1;
  else 
    return id2*size[dim1]+id1;
}

//---------------------------------------------------------------------------

bool
FPositionDistributionChecker::
isStructured() const
{
  for (positive i=0 ; i<3 ; i++)
    if (!depend[i][(i+1)%3] && !depend[i][(i+2)%3] &&
	!depend[(i+1)%3][i] && !depend[(i+2)%3][i])
      // at least one dimension is 1D structured
      return true;

  for (positive i=0 ; i<3 ; i++)
    if ((depend[i][(i+1)%3] && !depend[i][(i+2)%3] &&
	 depend[(i+1)%3][i] && !depend[(i+1)%3][(i+2)%3] &&
	 !depend[(i+2)%3][i] && !depend[(i+2)%3][(i+1)%3]) ||
	(depend[i][(i+2)%3] && !depend[i][(i+1)%3] &&
	 depend[(i+2)%3][i] && !depend[(i+2)%3][(i+1)%3] &&
	 !depend[(i+1)%3][i] && !depend[(i+1)%3][(i+2)%3]))
      // a 2D structured subspace exists
      return true;
  

  return false;
}
