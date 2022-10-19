//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularPath.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:42 $
// Author:    $Author: garth $
// Version:   $Revision: 1.21 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMSingularPath_hh
#define __FAMSingularPath_hh

#include "FObject.hh"
#include "FPosition.hh"
#include "FIndex.hh"

#include "FException.hh"
#include "FAMSingularPoint.hh"

#include "FVector.hh"

class FMatrix;
class FCell;

//===========================================================================

/** 
 *Structure to handle type and successive locations of a singular
 *point in a 2D unsteady vector field.
 */
class FAMSingularPath : public FObject
{
public:

  //===========================================================================
  
  /**
   * Generic path types
   */
  enum pathType { SADDLE_PATH, SOURCE_PATH, SINK_PATH, ATTACH_PATH,
		  SEPAR_PATH, // vector
		  TRISECTOR_PATH, WEDGE_PATH, WEDGE_SINGLE_PATH, // tensor
		  NONE_PATH };
  
  /**
   * Conversion from a singularity type to its corresponding path 
   * in a time dependent case.
   */
  static pathType main2path( FAMSingularPoint::mainType type);

  //===========================================================================
  
  /**
   * Undocumented ... tricoche
   */
  friend std::ostream& operator<< (std::ostream &os, const pathType& type);

  /** 
   *{\bf Description:}\\
   *Constructor: returns an empty FAMSingularPath.
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   */
  FAMSingularPath();

  /** 
   *{\bf Description:}\\
   *copy constructor
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *\param
   *{\bf sing}: singularity to copy
   */
  FAMSingularPath(const FAMSingularPath& path);
  
  /** 
   *{\bf Description:}\\
   *Constructor: returns a FAMSingularPath with given type and positions 
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *\param
   *{\bf type}: singularity type
   *\param
   *{\bf pos}: successive positions of the singularity.
   */
  FAMSingularPath(const pathType& _type, const std::vector<FVector>& _pos);

  /** 
   *{\bf Description:}\\
   *Constructor: returns a FAMSingularPath with given type and positions and cellIds
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *\param
   *{\bf type}: singularity type
   *\param
   *{\bf pos}: successive positions of the singularity.
   */
  FAMSingularPath(const pathType& _type, 
		  const std::vector< FVector >& _pos,
		  const std::vector< positive >& _cellIds);

  /** 
   *{\bf Description:}\\
   *Sets the separatrices start positions corresponding to the successive 
   *positions
   *\\{\bf Precondition:}\\
   *Positions have been set.
   *Singular path is of type SADDLE_PATH, TRISECTOR_PATH or WEDGE_PATH
   *The number of start positions is 4x the number of positions 
   *(FAMSingularPoint::SADDLE), 1x (ATTACH_PATH or SEPAR_PATH)
   *3x (FAMSingularPoint::TRISECTOR), 2x (FAMSingularPoint::WEDGE)
   *\\{\bf Postcondition:}\\
   *start positions have been set
   *\exception
   *FEmptyObjectException
   *@exception
   *FException: not of correct type
   *@exception
   *FInvalidDimensionException
   *\param
   *{\bf _sepStart}: separatrices start positions
   */
  void setSeparatricesStart(const std::vector< FVector >& _sepStart);
  /// Undocumented.
  void addSeparatricesStart(const std::vector< FVector >& _sepStart);

  /** 
   *{\bf Description:}\\
   *Gets the separatrices start position corresponding to the successive 
   *positions
   *\\{\bf Precondition:}\\
   *Singular path is of type SADDLE_PATH, ATTACH_PATH, SEPAR_PATH,
   *TRISECTOR_PATH or WEDGE_PATH
   *Start positions have been set.
   *of positions.
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *FEmptyObjectException
   *@exception
   *FException: not of correct type
   *@param
   *{\bf _sepStart}: returned start positions
   */
  void getSeparatricesStart(std::vector< FVector >& _sepStart) const;

  /** 
   *{\bf Description:}\\
   *returns the positions of the FAMSingularPath.
   *\\{\bf Precondition:}\\
   *Positions have been set.
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *FEmptyObjectException
   *\param
   *{\bf result:} returned positions.
   */
  void getPositions(std::vector<FVector>& result) const;
  
//  /// Undocumented.
//  void addPositions(std::vector<FVector>& steps);

  /** 
   *{\bf Description:}\\
   *returns the ids of the cells containing positions of the FAMSingularPath.
   *\\{\bf Precondition:}\\
   *Cell ids have been set.
   *\exception
   *FEmptyObjectException
   *\param
   *{\bf result:} returned cell ids.
   */
  void getCellIds(std::vector< positive >& result) const;

  positive getDimension() const;

  /** 
   *{\bf Description:}\\
   *returns the type of the FAMSingularPath.
   *\\{\bf Precondition:}\\
   *a type has been set.
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *FException: no type has been set.
   *\param
   *{\bf result:} returned type.
   */
  void getType(pathType& result) const;

  /** 
   *{\bf Description:}\\
   *Destructor.
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   */
  ~FAMSingularPath();

  /**
   *{\bf Description:}\\
   *Returns the class name.
   */
  virtual const FString& getClassName() const;

  /// Assignment operator.
  FAMSingularPath& operator=(const FAMSingularPath& path);


  /** 
   *{\bf Description:}\\
   *prints the contents of the singular path
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *none
   *\exception
   *none
   *\param
   *{\bf path}: path to print
   */
  friend std::ostream& operator<< (std::ostream &os, const FAMSingularPath& path);

  /// Undocumented.
  friend bool connected(const FAMSingularPath& path1, const FAMSingularPath& path2);

  /// Concatenate two paths.
  friend void concatenatePaths(FAMSingularPath& path, const FAMSingularPath& toAppend);

  /// Undocumented.
  void setStartBifId(const FIndex& start);
  /// Undocumented.
  void getStartBifId(FIndex& start) const;

  /// Undocumented.
  void setStopBifId(const FIndex& stop);
  /// Undocumented.
  void getStopBifId(FIndex& stop) const;

  /// Undocumented.
  void getStartPosition(FPosition& pos);

  /// Undocumented.
  void getStopPosition(FPosition& pos);

  /// returns the number of nodes of the path
  positive getNbPos();
  
  /// returns the arc length of the path
  double getArcLength() ;


  /** 
   *{\bf Description:}\\
   *Checks if input position lies on the path (with a given tolerance)
   *\\{\bf Precondition:}\\
   *none
   *\\{\bf Postcondition:}\\
   *returned TRUE iff input position "almost" lies on the path.
   *\exception
   *none
   *\param
   *{\bf pos}: input position replaced by its approximation
   *on the curve if TRUE has been returned, unchanged otherwise.
   *\param
   *{\bf tol}: numerical tolerance
   */
  bool isOnPath(FVector& aPos, double tol = 1.0E-3) const;
  
private:
  
  // bifurcation location
  std::vector< FVector > pos;
  std::vector< positive > cellIds;

  // singular path type 
  pathType type;

  // in the saddle case, separatrices start positions corresponding to 
  // the successive positions
  std::vector< FVector > sepStart; 

  // start/stop bifurcation
  FIndex start, stop;

  positive dimension;  // the overall dimension (e.g. 4 = 3 space + 1 time)
};

//--------------------------------------------------------------------------- 

inline void FAMSingularPath::getType(FAMSingularPath::pathType& result) const 
{
  if (type == FAMSingularPath::NONE_PATH) {
    THROW_EXCEPTION( FException, "ERROR: no type has been set");
  }
  result = type;
}

//--------------------------------------------------------------------------- 

inline FAMSingularPath::pathType FAMSingularPath::main2path(FAMSingularPoint::mainType type) 
{  
  switch (type) {
  case FAMSingularPoint::SADDLE: {
    return FAMSingularPath::SADDLE_PATH;
  }
  case FAMSingularPoint::SINK: {
    return FAMSingularPath::SINK_PATH;
  }
  case FAMSingularPoint::SOURCE: {
    return FAMSingularPath::SOURCE_PATH;
  }
  case FAMSingularPoint::TRISECTOR: {
    return FAMSingularPath::TRISECTOR_PATH;
  }
  case FAMSingularPoint::WEDGE: {
    return FAMSingularPath::WEDGE_PATH;
  }
  case FAMSingularPoint::WEDGE_SINGLE: {
    return FAMSingularPath::WEDGE_SINGLE_PATH;
  }
  case FAMSingularPoint::CENTER: 
  case FAMSingularPoint::OTHER:
    return FAMSingularPath::NONE_PATH;
  }

  return FAMSingularPath::NONE_PATH;
}


#endif // __FAMSingularPath_hh
