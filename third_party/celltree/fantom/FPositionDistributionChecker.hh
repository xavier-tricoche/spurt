//---------------------------------------------------------------------------
//
// Project:   TOP
// Module:    $RCSfile: FPositionDistributionChecker.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:08 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//---------------------------------------------------------------------------  

#ifndef FPositionDistributionChecker_hh
#define FPositionDistributionChecker_hh

#include "FObject.hh"
#include "stdAliases.hh"
#include "FBoundingBox.hh"
#include <vector>

using std::vector;

/**
 * This Class is used to characterize the structure of a structured
 * positionset whose dimensions are known (the VTK format i.e.).
 * Note that the positions are expected to be sorted according to
 * the usual X-Y-Z ordering of 3D data.
 * Usage: create an object of this type, initialize with the dimensions
 * of the grid and feed it with all positions successively.
 * At the end you can request the information on structure of the 
 * positionset.
 */
class FPositionDistributionChecker: public FObject
{
public:
  // constructor
  FPositionDistributionChecker();

  // destructor
  ~FPositionDistributionChecker();

  // intitalize the checker
  void initializeDimensions ( positive nx, positive ny, positive nz );

  // check the consecutive steps, considering that they are enumerated
  // in x, y, and z direction ...
  void checkStep( double comp[] );
  
  // gets the bounding box of all the checked positions
  void getBoundingBox ( FBoundingBox& theBBox ) const;

  // returns the matrix of boolean indicating for each coordinate
  // the dependency wrt every dimension
  void getDependencies(bool depend[3][3]) const;

  // returns the coordinates of a particular 'independent' dimension
  void getCoordinates1D(vector< double >& coord, positive dim) const;

  // returns the coordinates depending on two interrelated dimensions 
  void getCoordinates2D(vector< double > coord[2], positive dim1, 
			positive dim2) const;

  // returns true if a (any) structure was detected
  bool isStructured() const;

private:  

  // global counter to track the current position
  positive gCount;

  // dimension of the source data
  positive nx, ny, nz, size[3];

  // counter that tell the adress of the currently checked Position in
  // the structured grid
  positive discreteCoord[3];

  // for each dimension, characterize the identified structure:
  // dimensions' dependencies 
  bool depend[3][3];
  // uniform distribution of the coordinates for each dimension
  bool uniform[3];

  // the distances that coordinates in those directions covered
  double distance[3]; 

  FBoundingBox theBBox;

  void checkCoord( double , positive);

  // structure to save positions' coordinates
  std::pair< positive, std::vector< double > > coord[3];

  // compute the index corresponding to 2D coordinates according to the
  // structure above
  positive index(positive dim1, positive id1, positive dim2, positive id2);

  std::vector< double > coordinates;

  bool have_to_check;
};

#endif // FPositionDistributionChecker_hh

