//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularPath.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:42 $
// Author:    $Author: garth $
// Version:   $Revision: 1.18 $
//
//--------------------------------------------------------------------------- 

#include "FAMSingularPath.hh"
#include "FException.hh"

#include "FVector.hh"

#include <utility>
#include <iostream>

using namespace std;


//---------------------------------------------------------------------------
// Standard Constructor
FAMSingularPath::FAMSingularPath() 
  : pos(), cellIds(), type(FAMSingularPath::NONE_PATH), sepStart(), start(), stop(), dimension(0)
{
}

//--------------------------------------------------------------------------- 
// Copy Constructor
FAMSingularPath::FAMSingularPath(const FAMSingularPath& path) 
  : FObject(), pos(path.pos.size()), cellIds(path.cellIds.size()), type(path.type), sepStart(), start(path.start), stop(path.stop), dimension(path.dimension)
{
   for ( positive i=0; i<path.pos.size(); i++ ) 
     pos[i] = path.pos[i];
   for ( positive i=0; i<path.cellIds.size(); i++ ) 
     cellIds[i] = path.cellIds[i];
   for ( positive i=0; i<path.sepStart.size(); i++ ) 
     sepStart.push_back(path.sepStart[i]);
   
}

//--------------------------------------------------------------------------- 
// Constructor initialized by a vector of singular points forming the path
FAMSingularPath::FAMSingularPath(const FAMSingularPath::pathType& _type, 
				 const std::vector<FVector>& _pos)
  : pos(_pos.size()), cellIds(), type(_type), sepStart(), start(), stop()
{

   for (positive i=0; i<_pos.size(); i++) pos[i] = _pos[i];
  
   if (pos.size() > 0) dimension = pos[0].size();
   else dimension = 0;
}

//--------------------------------------------------------------------------- 
// Constructor initialized by a vector of singular points forming the path
FAMSingularPath::FAMSingularPath(const FAMSingularPath::pathType& _type, 
				 const std::vector<FVector>& _pos,
				 const std::vector< positive >& _cellIds)
  : pos(_pos.size()), cellIds( _cellIds.size() ), type(_type), sepStart(), start(), stop()
{

   for (positive i=0; i<_pos.size(); i++) pos[i] = _pos[i];
   for (positive i=0; i<_cellIds.size(); i++) cellIds[i] = _cellIds[i];
  
   if (pos.size() > 0) dimension = pos[0].size();
   else dimension = 0;
}

//--------------------------------------------------------------------------- 
//
positive FAMSingularPath::getDimension() const
{   
  return this->dimension;
}


//--------------------------------------------------------------------------- 

void FAMSingularPath::getSeparatricesStart(std::vector< FVector >& _sepStart) const
{
  if (type != FAMSingularPath::SADDLE_PATH && 
      type != FAMSingularPath::ATTACH_PATH && 
      type != FAMSingularPath::SEPAR_PATH && 
      type != FAMSingularPath::TRISECTOR_PATH && 
      type != FAMSingularPath::WEDGE_PATH && 
      type != FAMSingularPath::WEDGE_SINGLE_PATH) 
  {
    THROW_EXCEPTION( FException, "ERROR: not of correct type");
  }

  if (!sepStart.size()) THROW_DEFAULT_EXCEPTION( FEmptyObjectException );
  
  _sepStart.clear();
  for (positive i=0; i<sepStart.size(); i++) _sepStart.push_back(sepStart[i]);
  
}

//--------------------------------------------------------------------------- 

void FAMSingularPath::setSeparatricesStart(const std::vector<FVector>& _sepStart)
{
  if (type != FAMSingularPath::SADDLE_PATH && 
      type != FAMSingularPath::ATTACH_PATH && 
      type != FAMSingularPath::SEPAR_PATH && 
      type != FAMSingularPath::TRISECTOR_PATH && 
      type != FAMSingularPath::WEDGE_PATH && 
      type != FAMSingularPath::WEDGE_SINGLE_PATH) 
  {
    THROW_EXCEPTION( FException, "ERROR: not of correct type");
  }

  if (!pos.size()) THROW_DEFAULT_EXCEPTION( FEmptyObjectException );
  
  // apparently the following code is used to check if there is the correct number of
  // starting points given for a separatrix according to the type of the singular point
  if ( ( type == FAMSingularPath::SADDLE_PATH &&
	 _sepStart.size() != 4*pos.size() ) ||
       ( type == FAMSingularPath::ATTACH_PATH &&
	 _sepStart.size() != pos.size() ) ||
       ( type == FAMSingularPath::SEPAR_PATH &&
	 _sepStart.size() != pos.size() ) ||
       ( type == FAMSingularPath::TRISECTOR_PATH &&
	 _sepStart.size() != 3*pos.size() ) ||
       ( type == FAMSingularPath::WEDGE_PATH &&
	 _sepStart.size() != 2*pos.size() ) ||
       ( type == FAMSingularPath::WEDGE_SINGLE_PATH && 
	 _sepStart.size() != pos.size() ) ) 
  {
    THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
  }

  sepStart.clear();
  for (positive i=0 ; i<_sepStart.size() ; i++) sepStart.push_back(_sepStart[i]);
  
}

//---------------------------------------------------------------------------  
// Returns a vector containing the positions of the singular points forming the path
void FAMSingularPath::getPositions(std::vector<FVector>& result) const 
{

  if (!pos.size()) THROW_DEFAULT_EXCEPTION( FEmptyObjectException );
  
  result.clear();
  for (positive i=0 ; i<pos.size() ; i++)  
    result.push_back(pos[i]);
}


//---------------------------------------------------------------------------  
// Returns a vector containing the cells containing the positions of the singular points forming the path
void FAMSingularPath::getCellIds( std::vector< positive >& result ) const 
{
  if (!cellIds.size()) THROW_DEFAULT_EXCEPTION( FEmptyObjectException );
  
  result.clear();

  for (positive i=0 ; i<cellIds.size() ; i++)  
    result.push_back( cellIds[i] );
}


//--------------------------------------------------------------------------- 

FAMSingularPath::~FAMSingularPath()
{
}


//--------------------------------------------------------------------------- 

const FString& FAMSingularPath::getClassName() const 
{
  static FString className = "FAMSingularPath";

  return className;
}

      
//--------------------------------------------------------------------------- 
// assignment operator
FAMSingularPath& FAMSingularPath::operator=(const FAMSingularPath& path) 
{
  
  if (&path != this) 
  {
    pos.clear();
    sepStart.clear();
    cellIds.clear();
    for ( positive i=0; i<path.pos.size(); i++ )
      pos.push_back( path.pos[i] );
    for ( positive i=0; i<path.sepStart.size(); i++ )
      sepStart.push_back( path.sepStart[i] );
    for ( positive i=0; i<path.cellIds.size(); i++ )
      cellIds.push_back( path.cellIds[i] );

    type = path.type;
  }
  
  return *this;
}

//--------------------------------------------------------------------------- 

void FAMSingularPath::setStartBifId(const FIndex& start) 
{
   this->start = start;
}

//--------------------------------------------------------------------------- 

void FAMSingularPath::getStartBifId(FIndex& start) const {
  
  start = this->start;
}

//--------------------------------------------------------------------------- 

void FAMSingularPath::setStopBifId(const FIndex& stop) {
  
  this->stop = stop;
}

//--------------------------------------------------------------------------- 

void FAMSingularPath::getStopBifId(FIndex& stop) const {
  
  stop = this->stop;
}

//--------------------------------------------------------------------------- 
// if there is no path stored the result is undefined
void FAMSingularPath::getStartPosition(FPosition& pos) {

  if (this->pos.size()) 
    pos = FPosition(this->pos.front());
  else
    pos = FPosition();
}

//--------------------------------------------------------------------------- 
// if there is no path stored the result is undefined
void FAMSingularPath::getStopPosition(FPosition& pos) {

  if (this->pos.size()) 
    pos = FPosition(this->pos.back());
  else
    pos = FPosition();
}

positive FAMSingularPath::getNbPos() 
{
  return pos.size();
}

double FAMSingularPath::getArcLength() 
{
  double result = 0;
  for( positive i = 0; i < pos.size()-1; i++)
    result += ( pos[i] - pos[i+1]).norm();
  return result;
}


//--------------------------------------------------------------------------- 
// checks whether a given point lies on the path or not
bool FAMSingularPath::isOnPath(FVector& aPos, double tol) const{

  if(aPos.size() != dimension) THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
  
  double t = aPos[dimension-1];

//   cout << "in isOnPath: t = " << t << ", start = " << pos.front()[2] 
//        << ", stop = " << pos.back()[2] << endl;
  
  if (!pos.size() || (float) t < (float) pos.front()[dimension-1] || (float) t > (float) pos.back()[dimension-1])
    return false;

//   cout << "came through" << endl;

  positive old, first, last, current;
  old = 1234567890; // dummy
  first = 0;
  last = pos.size()-1;
  current = (first + last) / 2;
  
  while( true ) 
  {
    if ((float) pos[current][dimension-1] < (float) t)
      first = current;
    else if ((float) pos[current][dimension-1] > (float) t)
      last = current;
    else
      break;
    old = current;
    current = (first + last) / 2;
    if ((float) old == (float) current)
      break;
  }

  FVector best;
  if (pos[current][dimension-1] != t) 
  {   // linear interpolate best position on path
      double u = (t-pos[first][dimension-1]) / (pos[last][dimension-1]-pos[first][dimension-1]);
      best = (1.-u)*pos[first] + u*pos[last];
  }
  else
    best = pos[current];

//   cout << "aPos = " << aPos << ", best = " << best 
//        << "found distance = " << sqrt((aPos[0]-best[0])*(aPos[0]-best[0]) +
//    				      (aPos[1]-best[1])*(aPos[1]-best[1]) +
//    				      (aPos[2]-best[2])*(aPos[2]-best[2])) 
//        << endl;

  // we have found the best input position approximation
  // along the time axis: Compute the actual time/space distance
  double glitch = 0.0;
  for(positive i=0; i<dimension-1; i++) glitch += (aPos[i]-best[i])*(aPos[i]-best[i]);
  
  if(glitch < tol*tol) 
  {
    //replace input position with (satisfying) best approximation
    aPos = best;
    return true;
  }

  return false;
}

//---------------------------------------------------------------------------

ostream& operator<< (ostream &os, const FAMSingularPath::pathType& type) {
  switch (type) {
  case FAMSingularPath::SADDLE_PATH: 
    os << "SADDLE PATH";
    break;
  case FAMSingularPath::SOURCE_PATH:
    os << "SOURCE PATH";
    break;
  case FAMSingularPath::SINK_PATH:
    os << "SINK PATH";
    break;
  case FAMSingularPath::ATTACH_PATH: 
    os << "ATTACH PATH";
    break;
  case FAMSingularPath::SEPAR_PATH: 
    os << "SEPAR PATH";
    break;
  case FAMSingularPath::TRISECTOR_PATH:
    os << "TRISECTOR PATH";
    break;
  case FAMSingularPath::WEDGE_PATH:
    os << "WEDGE PATH";
    break;
  case FAMSingularPath::WEDGE_SINGLE_PATH:
    os << "WEDGE SINGLE PATH";
    break;
  case FAMSingularPath::NONE_PATH:
    os << "UNKNOWN/UNDEFINED PATH";
    break;
  }
  return os;
}

//---------------------------------------------------------------------------

ostream& operator<< (ostream &os, const FAMSingularPath& path)
{
  os << "FAMSingularPath:" << endl;
  os << "--positions: ";
  for ( positive i=0 ; i<path.pos.size(); i++ ) 
    os << path.pos[i] << " ";
  os << "--cellsIds: ";
  for ( positive i=0 ; i<path.cellIds.size(); i++ ) 
    os << path.cellIds[i] << " ";
  os << endl;
  os << "--type: " << path.type << endl;
  if ( path.type == FAMSingularPath::SADDLE_PATH || 
       path.type == FAMSingularPath::ATTACH_PATH || 
       path.type == FAMSingularPath::SEPAR_PATH || 
       path.type == FAMSingularPath::TRISECTOR_PATH || 
       path.type == FAMSingularPath::WEDGE_PATH || 
       path.type == FAMSingularPath::WEDGE_SINGLE_PATH ) {
    os << "--separatrices start: ";
    
    for (positive i=0 ; i<path.sepStart.size(); i++) 
      os << path.sepStart[i] << endl;
  }
  os << "-- start bifurcation has index " << path.start << endl;
  os << "-- stop bifurcation has index " << path.stop << endl;
   
  return os;
}


//---------------------------------------------------------------------------
// Checks if two paths are connected or not
bool connected(const FAMSingularPath &path1, const FAMSingularPath &path2) {

  positive dim1 = path1.getDimension();
  positive dim2 = path2.getDimension();

  if(dim1 != dim2) THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );  
    
  if (path1.pos.front()[dim1-1] > path1.pos.back()[dim1-1] ||
      path2.pos.front()[dim1-1] > path2.pos.back()[dim1-1])
    cout << "WARNING: reverse order detected" << endl;
      
  bool conn = false;
  for(positive i=0; i<dim1-1; i++) conn = conn && ((float) path1.pos.front()[i] == (float) path2.pos.back()[i]);
  if(!conn) 
    for(positive i=0; i<dim1-1; i++) conn = conn && ((float) path1.pos.back()[i] == (float) path2.pos.front()[i]);
  
  
  return (path1.type == path2.type && conn);
/*	  
	  ((
	    (float) path1.pos.front()[0] == (float) path2.pos.back()[0] &&
	    (float) path1.pos.front()[1] == (float) path2.pos.back()[1] &&
	    (float) path1.pos.front()[2] == (float) path2.pos.back()[2]
	    ) 
	   ||
	   (
	    (float) path1.pos.back()[0] == (float) path2.pos.front()[0] &&
	    (float) path1.pos.back()[1] == (float) path2.pos.front()[1] &&
	    (float) path1.pos.back()[2] == (float) path2.pos.front()[2]
	    )
	   )
	  );
*/
}

//---------------------------------------------------------------------------
// Concatenates two given paths
void concatenatePaths( FAMSingularPath& path, const FAMSingularPath& toAppend) 
{
  positive dim1 = path.getDimension();
  positive dim2 = toAppend.getDimension();

  if(dim1 != dim2) THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );  


  vector< FVector > pos;
  vector< positive > cellIds;
  vector< FVector > sepStart;

  FAMSingularPath::pathType type1, type2;
  path.getType(type1);
  toAppend.getType(type2);

//   cout << "going to append " << type1 << "("
//        << path.pos.size() << ", " << path.sepStart.size()
//        << ") and " << type2 << "(" 
//        << toAppend.pos.size() << ", " 
//        << toAppend.sepStart.size() << ")" << endl;

  if ((float) path.pos.front()[dim1-1] == (float) toAppend.pos.back()[dim1-1]) 
  {
//     cout << "path1 has first position " << path.pos.front() << endl
// 	 << "path2 has last position " << toAppend.pos.back() << endl;

    for (positive i=0 ; i<toAppend.pos.size() ; i++)
      pos.push_back(toAppend.pos[i]);
    for (positive i=1 ; i<path.pos.size() ; i++) 
      pos.push_back(path.pos[i]);
    for (positive i=0 ; i<toAppend.cellIds.size() ; i++)
      cellIds.push_back(toAppend.cellIds[i]);
    for (positive i=1 ; i<path.cellIds.size() ; i++) 
      cellIds.push_back(path.cellIds[i]);

    if (path.type == FAMSingularPath::SADDLE_PATH) {

      for (positive j=0 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);

      if ( sepStart.size() && path.sepStart.size() && !sepStart.back().size()) {
// 	cout << sepStart[sepStart.size()-4] 
// 	     << " " << sepStart[sepStart.size()-3] 
// 	     << " " << sepStart[sepStart.size()-2]
// 	     << " " << sepStart[sepStart.size()-1]
// 	     << " " << path.sepStart[0] 
// 	     << " " << path.sepStart[1] 
// 	     << " " << path.sepStart[2] 
// 	     << " " << path.sepStart[3] << endl;

	sepStart[sepStart.size()-3] = path.sepStart[0];
	sepStart.back() = path.sepStart[2];
      }

      for (positive j=4 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);
    }
    else if (path.type == FAMSingularPath::ATTACH_PATH) {

      for (positive j=0 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);

      for (positive j=1 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);
    }
    else if (path.type == FAMSingularPath::SEPAR_PATH) {

      for (positive j=0 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);

      for (positive j=1 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);
    }
    else if (path.type == FAMSingularPath::TRISECTOR_PATH) {

      for (positive j=0 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);

      for (positive j=3 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);
    }
    else if (path.type == FAMSingularPath::WEDGE_PATH) {

      for (positive j=0 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);

      for (positive j=2 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);
    }
    else if (path.type == FAMSingularPath::WEDGE_SINGLE_PATH) {

      for (positive j=0 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);

      for (positive j=1 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);
    }
  }
  else {

//     cout << "path1 has first position " << toAppend.pos.front() << endl
// 	 << "path2 has last position " << path.pos.back() << endl;

    for (positive i=0 ; i<path.pos.size() ; i++) 
      pos.push_back(path.pos[i]);
    for (positive i=1 ; i<toAppend.pos.size() ; i++) 
      pos.push_back(toAppend.pos[i]);
    for (positive i=0 ; i<path.cellIds.size() ; i++) 
      cellIds.push_back(path.cellIds[i]);
    for (positive i=1 ; i<toAppend.cellIds.size() ; i++) 
      cellIds.push_back(toAppend.cellIds[i]);

    if (path.type == FAMSingularPath::SADDLE_PATH) {
      for (positive j=0 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);

      if (!sepStart.back().size()) {
// 	cout << sepStart[sepStart.size()-4] 
// 	     << " " << sepStart[sepStart.size()-3] 
// 	     << " " << sepStart[sepStart.size()-2]
// 	     << " " << sepStart[sepStart.size()-1]
// 	     << " " << toAppend.sepStart[0] 
// 	     << " " << toAppend.sepStart[1] 
// 	     << " " << toAppend.sepStart[2] 
// 	     << " " << toAppend.sepStart[3] << endl;
	sepStart[sepStart.size()-3] = toAppend.sepStart[0];
	sepStart.back() = toAppend.sepStart[2];
      }

      for (positive j=4 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);
    }
    else if (path.type == FAMSingularPath::ATTACH_PATH) {

      for (positive j=0 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);

      for (positive j=1 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);
    }       
    else if (path.type == FAMSingularPath::SEPAR_PATH) {

      for (positive j=0 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);

      for (positive j=1 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);
    }       
    else if (path.type == FAMSingularPath::TRISECTOR_PATH) {

      for (positive j=0 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);

      for (positive j=3 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);
    }
    else if (path.type == FAMSingularPath::WEDGE_PATH) {

      for (positive j=0 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);

      for (positive j=2 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);
    }   
    else if (path.type == FAMSingularPath::WEDGE_SINGLE_PATH) {

      for (positive j=0 ; j<path.sepStart.size() ; j++)
	sepStart.push_back(path.sepStart[j]);

      for (positive j=1 ; j<toAppend.sepStart.size() ; j++)
	sepStart.push_back(toAppend.sepStart[j]);
    }       
  }
    
  path.pos = pos;
  path.cellIds = cellIds;
  path.sepStart = sepStart;
}

//--------------------------------------------------------------------------- 
