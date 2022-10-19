//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularPoint.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:33:33 $
// Author:    $Author: garth $
// Version:   $Revision: 1.30 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMSingularPoint_hh
#define __FAMSingularPoint_hh

#include "FObject.hh"
#include "FPosition.hh"
#include "FIndex.hh"

#include "FException.hh"

#include <list>
#include <complex>

#include "FVector.hh"
class FMatrix;
class FCell;

/** 
 * Class to handle the singular points that occur in std::vector-
 * and tensorfields. Not only are position and the containing cell stored for
 * the singularity, but also its nature, type and main type, as well as
 * references (FIndex) to separatrices starting or ending in it.
 */
class FAMSingularPoint : public FObject
{
public:

//===========================================================================

  /// Singularity nature enumeration
  enum singularityNature { CRITICAL, DEGENERATE, UNKNOWN };
  
  /// Singularity type enumeration
  enum singularityType {
    // std::vectors:
    SADDLE_2D, ATTRACT_NODE, REPELL_NODE, ATTRACT_FOCUS,
    REPELL_FOCUS, CENTER_2D, NONE, HIGHER, BOUNDARY,
    ATTRACT_NODE_3D, REPELL_NODE_3D, CENTER_3D,
    ATTRACT_FOCUS_3D, REPELL_FOCUS_3D,
    SADDLE_3D, SPIRAL_SADDLE_3D,
    //Tensor
    TRISECTOR_POINT, WEDGE_POINT, WEDGE_SINGLE_POINT, HIGHER_DEGENERATE};  
  
  /// Main type enumeration
  enum mainType { SADDLE, SINK, SOURCE, CENTER, TRISECTOR, WEDGE,
		  WEDGE_SINGLE, OTHER };
  
  static mainType whichType(singularityType type);

//===========================================================================

  
  /** 
   *Description:
   *Constructor: returns an empty FAMSingularPoint.
   */
  FAMSingularPoint();

  
  /** 
   *Description:
   *Constructor: returns an empty FAMSingularPoint with given nature.
   *\param
   * nature: singularity nature (critical or degenerate)
   */
  FAMSingularPoint(const singularityNature& nature);

  /** 
   *Description:
   *copy constructor
   *\param
   * sing: singularity to copy
   */
  FAMSingularPoint(const FAMSingularPoint& sing);
  
  /** 
   *Description:
   *Constructor: returns a FAMSingularPoint located at the given position.
   *\param
   * pos: position where the zero is.
   */
  FAMSingularPoint(const FPosition& pos);

  /** 
   *Description:
   *Constructor: returns a FAMSingularPoint located at the given position with the
   *given nature.
   *\param
   * pos: position where the zero is.
   *\param
   * nature: singularity nature (critical or degenerate)
   */
  FAMSingularPoint(const FPosition& pos, const singularityNature& nature);  

  /** 
   *Description:
   *Constructor: returns a FAMSingularPoint located at the given position 
   *and having the given type with given order.
   *\exception
   *none
   *\param
   * pos: position where the zero is.
   *\param
   * type: singularity type. 
   *\param
   * order: singularity order
   */
  FAMSingularPoint(const FPosition& pos, const singularityType& type, 
		   positive order);

  /** 
   *Description:
   *Constructor: returns a FAMSingularPoint located at the given position 
   *and having the given type with given order and nature.
   *\param
   * pos: position where the zero is.
   *\param
   * nature: singularity nature (critical or degenerate)
   *\param
   * type: singularity type. 
   *\param
   * order: singularity order   
   */
  FAMSingularPoint(const FPosition& pos, const singularityNature& nature,
		   const singularityType& type, positive order); 

  /** 
   *Description:
   *returns the index of the cell containing the singularity.
   * \pre 
   *a cell index has been set
   * \post
   *none
   *\exception
   *FException: no cell index has been set.
   *\param
   * cellId: returned cell index
   */
  void getIncludingCellIndex(FIndex& cellId) const;
  /**
   * Undocumented 
   */
  void setIncludingCellIndex(const FIndex& cellId);

  /** 
   *Description:
   *returns the position of the FAMSingularPoint.
   * \pre
   *a position has been set.
   * \post
   *none
   *\exception
   *FEmptyObjectException
   *\param
   * result: returned position of the singularity.
   */
  void getPosition(FPosition& result) const;
  /**
   * Undocumented 
   */
  void setPosition(const FPosition& result);

  /** 
   *Description:
   *returns the type of the FAMSingularPoint.
   * \pre
   *a type has been set.
   * \post
   *none
   *\exception
   *FException: no type has been set.
   *\param
   * result: returned type of the singularity.
   */
  void getType(singularityType& result) const;
  void setType(singularityType val);

  /** 
   *Description:
   *returns the nature of the FAMSingularPoint.
   * \pre
   *a type has been set.
   * \post
   *none
   *\exception
   *FException: no type has been set.
   *\param
   * result: returned nature of the singularity.
   */
  void getNature(singularityNature& result) const;
  void setNature(singularityNature val);

  /** 
   *Description:
   *returns the order of the FAMSingularPoint.
   * \pre
   *an order has been set.
   * \post
   *none
   *\exception
   *FException: no order has been set.
   */
  positive getOrder() const;
  void setOrder( positive order);


  /** 
   *Description:
   *Destructor.
   * \pre
   *none
   * \post
   *none
   *\exception
   *none
   */
  virtual ~FAMSingularPoint();

  /**
   *Description:
   *Returns the class name.
   */
  virtual const FString& getClassName() const;

  /** 
   *Description:
   *Finds out the nature of the given singular point with the linear 
   *approximation given as matrix in input.
   * \pre
   *The field has order 1.
   * \post
   *The nature of the singular point has been set. 
   *If the singular point turns out to be a saddle point, the
   *eigenstd::vectors are set.
   *\exception
   *FInvalidDimensionException
   *\exception
   *FException: degenerate or higher order singularity found! 
   *\param
   * approx: linear approximation of the field around the
   *singular point.
   *\param
   * zero: considered singular point.
   */
  void setLinearNature(const FMatrix& approx);

  /**
   *{\Description:
   *Finds out the type of the degenerate point and computes the
   *angles of the separatrices
   * \pre
   *The singular point is a degenerate point and has dimension 2.
   * \post
   *The type of the degenerate point has been set.
   *\exception
   *FInvalidDimensionException
   *\param
   * a,b,c,d: values to compute the rotation invariant
   */
  void setDegenerateType(double a, double b, double c, double d,
			 const FCell *cell);

  /**
   *Description:
   * gets the Eigenvector of the std::vectorfield in linear case
   * \pre
   * eigenstd::vectors have been set
   * \post
   * eigenVec contains copy of eigenstd::vectors
   *\exception
   *none
   *\param
   * eigenVec the eigenvalues to be copied
   */
  void getEigenvectors(std::vector<FVector>& eigenVec) const;

  /**
   *Description:
   * sets the Eigenvector of the std::vectorfield in linear case
   * \pre
   * none
   * \post
   * eigenstd::vectors have been set
   *\exception
   *none
   *\param
   * eigenVec the eigenvalues to be copied
   */
  void setEigenvectors(const std::vector<FVector>& eigenVec);

  /**
   *Description:
   * gets the Eigenvalues of the std::vectorfield in linear case
   * \pre
   * eigenvalues have been set
   * \post
   * eigenVal contains copy of eigenvalues
   *\exception
   *none
   *\param
   * eigenVal the eigenvalues to be copied
   */
  void getEigenvalues(std::vector<std::complex<double> >& eigenVal) const;

  /**
   *Description:
   * Sets the Eigenvalues of the std::vectorfield in linear case
   * \post
   *eigenvalues set
   *\param
   * eigenVal the eigenvalues to be copied
   */
  void setEigenvalues(const std::vector<std::complex<double> >& eigenVal);

  /**
   *Description:
   *Adds a separatrix index to the std::list of indices of separatrices
   *starting at the singular point.
   *\param
   * start: index of a separatrix starting at the singular point
   */
  void addStartingSeparatrix(const FIndex& start);

  /**
   * Undocumented 
   */
  void getStartingSeparatrices(std::list<FIndex>& starts) const;

  /**
   *Description:
   *Adds a separatrix index to the std::list of indices of separatrices
   *ending at the singular point.
   *\param
   * stop: index of a separatrix ending at the singular point
   */
  void addStoppingSeparatrix(const FIndex& stop);
  
  /**
   * Undocumented 
   */
  void getStoppingSeparatrices(std::list<FIndex>& stops) const;

  /**
   * Undocumented 
   */
  void getSeparatrixAngles(std::vector<double>& sAngles);

  /**
   * Undocumented 
   */
  void setSeparatrixAngles(const std::vector<double>& sAngles);

  /**
   * Undocumented 
   */
  void getSeparatrixAngle(unsigned int sepIndex, double& angle);

  /**
   * Undocumented 
   */
  double getSeparatrixAngle(unsigned int sepIndex);

  /**
   * Undocumented 
   */
  FAMSingularPoint& operator=(const FAMSingularPoint& sing);

  /**
   * Set connections of a singular point with another one: 
   * indices of existing FAMSeparatrix objects must be provided
   */
  void setConnections(const std::vector< FIndex >& connections);
  void addConnection(const FIndex& connection);

  /**
   * Get connections of a singular point with another one: 
   * indices of existing FAMSeparatrix objects are returned
   */
  void getConnections(std::vector< FIndex >& connections) const;


  /** 
   *Description:
   *prints the contents of the singular point
   *\param
   * singularity: singular point to print.
   */
  friend std::ostream& operator<< (std::ostream &os, const FAMSingularPoint& singularity);

  /** 
   *Description:
   *return (in the std::vector case) the enstrophy, that is the absolute value
   *of the vorticity
   * \pre
   *none (saddle point has negative enstrophy)
   */
  double getEnstrophy() const;


  void setMaxNorm(double max_norm);
  double getMaxNorm() const;

private:
  
  // singular point position
  FPosition pos;

  // singular point nature (critical,degenerate)
  singularityNature nature;

  // singular point type
  singularityType type;

  // singularity order
  positive order;

  // eigenvector::vectors
  std::vector<FVector> eigenvector;
  // eigenvalues
  std::vector< std::complex<double> > eigenvalues;

  // separation angles
  std::vector<double> sepAngles;

  //-----------------------------------------------------------------

  // indices of the separatrices that are starting at the singular point
  std::list<FIndex> sepStart;

  // indices of the separatrices that are ending at the singular point
  std::list<FIndex> sepStop;

  //-----------------------------------------------------------------

  // index of the containing cell
  FIndex cellId;

  // indices of separatrices
  std::vector< FIndex > connections;

  double enstrophy;

  double max_norm;
};


/**
 * Undocumented ... 
 */
std::ostream& operator<< (std::ostream &os, const FAMSingularPoint::singularityNature& type);

/**
 * Undocumented ... 
 */
std::ostream& operator<< (std::ostream &os, const FAMSingularPoint::singularityType& type);


/**
 * Undocumented ... 
 */
std::ostream& operator<< (std::ostream &os, const FAMSingularPoint::mainType& type);




//===========================================================================
#ifndef OUTLINE
#include "FAMSingularPoint.icc"
#endif
//=========================================================================== 

#endif // __FAMSingularPoint_hh
