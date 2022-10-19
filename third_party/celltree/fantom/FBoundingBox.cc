//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FBoundingBox.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.13 $
//
//--------------------------------------------------------------------------- 

#include <cmath>
#include <cstdio>
#include <algorithm>

#include <util.hh>

#include "FBoundingBox.hh"

#ifdef OUTLINE
#include "FBoundingBox.icc"
#endif

#include <iostream>

using namespace std;

//---------------------------------------------------------------------------

FBoundingBox::~FBoundingBox()
{
  dim = 0; 
  if (range != NULL)
    delete[] range;
}

//---------------------------------------------------------------------------

bool FBoundingBox::isInside(const FPosition& position) const 
{
#ifndef NODEBUG
  if (position.getDimension() != dim || dim == 0) {
    FInvalidDimensionException e;
    e.addTraceMessage("bool FBoundingBox::isInside(const FPosition& position) const");
    throw e;
  }
#endif
  

  if (dim == 2)
    return (position(0)>=range[0]-epsilon && 
	    position(0)<=range[2]+epsilon &&
            position(1)>=range[1]-epsilon && 
	    position(1)<=range[3]+epsilon);
  else
    return (position(0)>=range[0]-epsilon && 
	    position(0)<=range[3]+epsilon &&
            position(1)>=range[1]-epsilon && 
	    position(1)<=range[4]+epsilon &&
            position(2)>=range[2]-epsilon && 
	    position(2)<=range[5]+epsilon);
}

//---------------------------------------------------------------------------

bool FBoundingBox::isInside(const FBoundingBox& inBBox) const
{
  //test if inBBox is inside of 'this'

#ifndef NODEBUG
  if (inBBox.getDimension() != dim){
    FInvalidDimensionException e("ERROR: Dimensions of both BBoxes have to be equal");
    e.addTraceMessage("bool FBoundingBox::isInside(const FBoundingBox&)");
    throw e;
  }
#endif

  // got it from CGAL

  if (getDimension() == 2) {
    if (range[2] < inBBox.range[0] || inBBox.range[2] < range[0])
      return false;
    if (range[3] < inBBox.range[1] || inBBox.range[3] < range[1])
      return false;
  } 
  else if (getDimension() == 3) {
    if (range[3] < inBBox.range[0] || inBBox.range[3] < range[0])
      return false;
    if (range[4] < inBBox.range[1] || inBBox.range[4] < range[1])
      return false;
    if (range[5] < inBBox.range[2] || inBBox.range[5] < range[2])
      return false;
  } 
  else {
    FInvalidDimensionException e("ERROR: Unkown dimension");
    e.addTraceMessage("bool FBoundingBox::isInside(const FBoundingBox&)");
    throw e;  
  }
  
  return true;
}

//---------------------------------------------------------------------------

bool FBoundingBox::includes( const FBoundingBox& inBBox ) const
{
  //test if inBBox is included in 'this'

#ifndef NODEBUG
  if (inBBox.getDimension() != dim){
    FInvalidDimensionException e("ERROR: Dimensions of both BBoxes have to be equal");
    e.addTraceMessage("bool FBoundingBox::includes(const FBoundingBox&)");
    throw e;
  }
#endif
  
  for ( positive i=0;i<dim; i++)
    if ( (range[i]>=inBBox.range[i]) || 
	 (range[dim+i]<=inBBox.range[dim+i]) ) 
      return (false);

  return (true);
}

//---------------------------------------------------------------------------

double FBoundingBox::diagonal() const {

  double res=0;

  for (positive i=0 ; i<dim ; i++)
    res += (range[i+dim]-range[i])*(range[i+dim]-range[i]);
  
  return sqrt(res);
}

//---------------------------------------------------------------------------

double FBoundingBox::computeBoundingSphereRadius() const
{
  if (dim != 3) 
    {
      FInvalidDimensionException e;
      e.addTraceMessage("FBoundingBox::computeBoundingSphereRadius: dimension != 3");
      cout << "FBoundingBox::computeBoundingSphereRadius: dimension != 3" << endl;
      throw e;
    }  

  // find longest vertex
  double vertex1 = fabs(range[1] - range[0]);
  double vertex2 = fabs(range[3] - range[2]);
  double vertex3 = fabs(range[5] - range[4]);
  
  return 0.5*std::max( vertex1, std::max(vertex2, vertex3) );
}

//===========================================================================

void FBoundingBox::resize(const FPosition& inPos)
{
  positive i, dim = inPos.getDimension();

  for (i=0; i<dim ; i++){
    if (range[i] > inPos[i])
      range[i] = inPos[i];
    else if (range[i + dim]  < inPos[i])
      range[i + dim] = inPos[i];
  }

}

//===========================================================================

ostream& operator << (ostream& os, const FBoundingBox& inBB) 
{
  
  os << "(";

  if (inBB.dim==2){
    os << inBB.range[0] << ", " << inBB.range[1] << ") - (";
    os << inBB.range[2] << ", " << inBB.range[3] << ")";
  } else if (inBB.dim==3) {
    os << inBB.range[0] << ", " << inBB.range[1] << ", " <<
      inBB.range[2] << ") - (";
    os << inBB.range[3] << ", " << inBB.range[4] << ", " <<
      inBB.range[5] << ")";
  } else {
    os << "Bounding Box has bad dimension " << inBB.dim << ".";
  }

  return os;
  
}

//===========================================================================

istream& operator>> (istream& is, FBoundingBox& inBB) 
{
    string tmp;

  getline( is, tmp );

  if( inBB.range )
    delete inBB.range;

  if( count( tmp.begin(), tmp.end(), ',' ) == 4 )
  {
      inBB.dim = 3;
      inBB.range = new double[6];
      
      const char *t = tmp.c_str();

      for(;(*t!='\0') && (*t!='('); t++);
      t++;
      sscanf (t, "%lf",&(inBB.range[0]));
      for(;(*t!='\0') && (*t!=','); t++);
      t++;
      sscanf (t, "%lf",&(inBB.range[1]));
      for(;(*t!='\0') && (*t!=','); t++);
      t++;
      sscanf (t, "%lf",&(inBB.range[2]));
      for(;(*t!='\0') && (*t!='('); t++);
      t++;
      sscanf (t, "%lf",&(inBB.range[3]));
      for(;(*t!='\0') && (*t!=','); t++);
      t++;
      sscanf (t, "%lf",&(inBB.range[4]));
      for(;(*t!='\0') && (*t!=','); t++);
      t++;
      sscanf (t, "%lf",&(inBB.range[5]));
    }
  else
    {
      inBB.dim = 2;
      inBB.range = new double[4];

      const char *t = tmp.c_str();

      for(;(*t!='\0') && (*t!='('); t++);
      t++;
      sscanf (t, "%lf",&(inBB.range[0]));
      for(;(*t!='\0') && (*t!=','); t++);
      t++;
      sscanf (t, "%lf",&(inBB.range[1]));
      for(;(*t!='\0') && (*t!='('); t++);
      t++;
      sscanf (t, "%lf",&(inBB.range[2]));
      for(;(*t!='\0') && (*t!=','); t++);
      t++;
      sscanf (t, "%lf",&(inBB.range[3]));
   }

  return is;
}

//===========================================================================

void FBoundingBox::getRange( double& xmin, double& xmax, 
			     double& ymin, double& ymax, 
			     double& zmin, double& zmax ) const 
{
    if( dim == 2 )
    {
	xmin = range[0]; ymin = range[1]; zmin = 0.0;
	xmax = range[2]; ymax = range[3]; zmax = 0.0;
    }
    else
    {
	xmin = range[0]; ymin = range[1]; zmin = range[2];
	xmax = range[3]; ymax = range[4]; zmax = range[5];
    }
}        

//--------------------------------------------------------------------------- 

FArray FBoundingBox::random() const
{
    if( dim == 2 )
    {
	double x = drand48(), y = drand48();

	return FArray( (1-x)*range[0] + x*range[2],
		       (1-y)*range[1] + y*range[3], 0 );
    }
    else
    {
	double x = drand48(), y = drand48(), z = drand48();

	return FArray( (1-x)*range[0] + x*range[3],
		       (1-y)*range[1] + y*range[4], 
		       (1-z)*range[2] + z*range[5]);
    }
}

//--------------------------------------------------------------------------- 

FArray FBoundingBox::min() const
{
    if( dim == 2 )
	return FArray( range[0], range[1], 0 );
    else
	return FArray( range[0], range[1], range[2] );
}

//--------------------------------------------------------------------------- 

FArray FBoundingBox::max() const
{
    if( dim == 2 )
	return FArray( range[2], range[3], 0 );
    else
	return FArray( range[3], range[4], range[5] );
}

//--------------------------------------------------------------------------- 
  
FBoundingBox FBoundingBox::intersection( const FBoundingBox& rhs ) const
{
  double range[ 6 ];
  for ( unsigned int i = 0; i< dim; ++i )
    range[ i ] = std::max( this->range[ i ], rhs.range[ i ] );

  for ( unsigned int i = 0; i< dim; ++i )
    range[ dim + i ] = std::min( this->range[ dim + i ], rhs.range[ dim + i ] );
  if ( dim == 2 )
    return FBoundingBox( range[ 0 ], range[ 1 ], range[ 2 ], range[ 3 ] );
  else if ( dim == 3 )
    return FBoundingBox( range[ 0 ], range[ 1 ], range[ 2 ], range[ 3 ], range[ 4 ], range[ 5 ] );
  return FBoundingBox();
}
