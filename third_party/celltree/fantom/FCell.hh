//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FCell.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:36:59 $
// Author:    $Author: garth $
// Version:   $Revision: 1.24 $
//
//--------------------------------------------------------------------------- 

#ifndef __FCell_hh
#define __FCell_hh

#include "FObject.hh"

#include "stdAliases.hh"
#include "FIndex.hh"
#include "FBoundingBox.hh"

#include "FRefArray.hh"
#include "FPositionSet.hh"
#include "FRefTensor.hh"
#include "FTensorSet.hh"

#include <list>
#include <vector>
#include <utility>

class FAMSingularPoint;

//this class is only present in FCell.cc
//contains data for faceToGo
class FCellFaceInfo;

//===========================================================================

/**
 * The FCell class provides an abstract geometric cell.
 *
 * ATTENTION! every class derived from FCell MUST set the pointer
 * GeometryDescription in its constructor to a static array in a format
 * described below (section above geometryDescription)
 *
 * When introducing a newcell type, you should also correct the constant
 *      maxNoEdges, maxNoVertices, maxNoVertsPerFace, maxNoFaces
 * here, if their value is too small.
 *
 * And you should add your cell in the getGeometryDestription Function in
 * the FCell.cc file.
 *
 * \ingroup Cells
 */
class FCell : public FObject
{
public:

  /**
   * Enumeration listing all the types of cells that appear in FAnToM
   */
  // ATTENTION !!!! WHEN INSERTING OR CHANGING THINGS HERE, CHANGE THEM
  // ACCORDINGLY IN THE .CC FILE for the ostream operator !!!!!!!!!!!!!
  // and the Marching Cubes algorithm in visKernel!
  // and in FTensorFieldReader* (?!)
  // and maybe much more ==> pay attention!!!!
  enum CellType { 
    TRIANGLE_2D,             
    AXIS_PARALLEL_TRI_2D,    
    QUADRILATERAL_2D,        
    AXIS_PARALLEL_QUAD_2D,   
    TRIANGLE_3D,             
    QUADRILATERAL_3D,        
    TETRAHEDRON,             
    AXIS_PARALLEL_TET,       
    ARBITRARY_HEX,           
    AXIS_PARALLEL_HEX,       
    PRISM,                   
    PYRAM,                   
    LINE_2D,
    LINE_3D,
    POINT_2D,
    POINT_3D,
    UNDEFINED                
  };

  
  /** 
   *\struct geoInfo
   * Structure holding the geometry information of the cell, this will be a
   * static variable in the respective derived class.
   */
  struct geoInfo;

  /** 
   *\par Description:
   *Constructor; multiple pointers/references are set to array sin the derived class
   *\param g:
   *the static,const geoInfo structure of the derived cell class
   *\param aVertexIndices:
   *should point on a FIndex array with sizeOfCellType() entries
   *\param aPositions:
   *should point on a FRefArray array with sizeOfCellType() entries
   *\param posData:
   *should point on a double array with sizeofCellType()*getDimension() entries
   *\param  aTensors
   *should point on a FRefTensor array with sizeOfCellType() entries
   */
  FCell( const geoInfo&g,
	 FIndex * const indices,
	 FRefArray * const positions,double*const posData,
	 FRefTensor * const tensors );

  /** 
   *\par Description:
   *Destructor
   */
  virtual ~FCell();

  /** 
   *\par Description:
   * returns a pointer on a copy of the cell.
   */
  virtual FCell* getClone() const = 0;



  /** 
   *\par Description:
   * returns a pointer on a new cell torso of type t
   * with vertex numbers vertices.
   *\param t: type of the new cell
   *\param vertices: vertex indices of the new cell
   *\exception: FNotImplementedException
   * if no entry for the given cell type is there
   */
  static FCell* getCellTorso(CellType t,const std::vector<FIndex> & vertices );


  /** 
   *\par Description:
   * returns the size of the FIndex array that is needed by the class
   * Default value: geometryDescription.noVerts
   *\return
   * number of needed FIndexes.
   */
  virtual positive sizeOfCellType() const;

  /**
   *\par Description:
   * returns the cell dimension.
   */
  virtual positive getDimension() const;

  /**
   *\par Description:
   * getter for the indices of the cell vertices
   *\param
   * result: returned vertices indices
   */
  void getVertexIndices(std::vector<FIndex>& result) const;
  void getVertices(std::vector<FPosition>& result) const;
  void getTensors(std::vector<FTensor>& result) const;
  void getTensor(FTensor& t, unsigned int i ) const;

  
  /** 
   *\par Description:
   *Returns the cell edges.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned edges formated as pairs of vertices indices.
   */
  void getEdges(std::vector< pair<FIndex, FIndex> >& result) const;
  void getEdges(std::vector< pair<FPosition, FPosition> >& result) const;
  void getEdges(std::vector< pair<FTensor, FTensor> >& result) const;

 /** 
   *\par Description:
   *Returns the cell faces.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned faces formated as lists of vertices indices.
   */
  void getFaces(std::vector< std::vector<FIndex> >& result) const;
  void getFaces(std::vector< std::vector<FArray> >& result) const;
  void getFaces(std::vector< std::vector<FTensor> >& result) const; 
  
  /** 
   *\par Description:
   *Returns the interpolated tensor value at \b position, for the
   *tensorSet and positionSet that were defined by positionsToUpdate()
   *and tensorsToUpdate()
   *\pre
   *the position cache and tensor cache must be valid, and the position
   *should (must ?) be inside the cell.
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned interpolated tensor.
   *\param
   *position: position at which to compute the interpolant
   */
  virtual void interpolate(FTensor& result,
			   const FPosition& position ) const = 0;

  /** 
   *\par Description:
   *Returns the interpolated derivative value at \b position, for the
   *tensorSet and positionSet that were defined by positionsToUpdate()
   *and tensorsToUpdate()
   *\pre
   *the position cache and tensor cache must be valid, and the position
   *should (must ?) be inside the cell.
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned interpolated tensor.
   *\param
   *position: position at which to compute the interpolant
   */
  virtual void derivatives(FTensor& result,
			   const FPosition& position ) const = 0;

  /** 
   *\par Description:
   *Returns zeros lying in the cell by adding them to the given list.
   *\pre
   *the position cache and tensor cache must be valid
   *\post
   *if a zero is inside the cell, it is added to the list.
   *\exception
   *none
   *\param
   *result: returned zeros.
   */
  virtual void getZeros(std::list<FAMSingularPoint>& result) const = 0;

  /** 
   *\par Description:
   *Indicates if the given \b position lies in the cell.
   *\pre
   *the position cache and tensor cache must be valid
   *\post
   *none
   *\exception
   *none
   *\param
   *position: position to test.
   */
  virtual bool isInside(const FPosition& position) const = 0;

  /** 
   *\par Description:
   *Returns the cell bounding box.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned cell bounding box.
   */
  void getBoundingBox(FBoundingBox& result) const;
  const FBoundingBox& getBoundingBox(void) const;		       


  /** 
   *\par Description:
   *Constructs the bounding box of the given positions
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   *\param
   *result: returned bounding box.
   */
  static FBoundingBox getBoundingBox(const std::vector< FIndex >& verticesId,
				     FPositionSet *posSet);

  static FBoundingBox getBoundingBox(const std::vector< FPosition >& vertices );

  
  /** 
   *\par Description:
   *Returns barycenter
   *\pre
   *positions have been set
   *\post
   *none
   *\exception
   *FException: positions not set
   *\param
   *result: returned barycenter
   */
  virtual FPosition getBarycenter() const;

  
  /** 
   *\par Description:
   *returns the type of the cell. (similar to getClassName, but more efficient
   */
  virtual CellType getCellType(void) const = 0;

//===========================================================================
// VIRTUAL METHODS
//===========================================================================
  
  /** 
   *\par Description:
   *Returns the center of the bounding box.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  virtual FPosition getCenterOfBoundingBox() const;

  /** 
   *\par Description:
   *Indicates if the given bounding box intersects the bounding box of
   *the cell.
   *\pre
   *the position cache cache must be valid
   *\post
   *none
   *\exception
   *none
   *\param
   *position: the bounding box which we test against
   */
  virtual bool intersects(const FBoundingBox & inBox);

  /** 
   *\par Description:
   * set new positions from the given position set
   *\pre
   * the indices of the vertices have been set
   *\post
   * new position values have been set
   *\exception
   *none
   */
  void setPositions( const FPositionSet* posSetPt );

  /** 
   *\par Description:
   * set new tensors from the given tensor set
   *\pre
   * the indices of the vertices have been set
   *\post
   * new tensor values have been set
   *\exception
   *none
   */
  void setTensors( const FTensorSet* tensSetPt );



  /**
   *\par Description:
   * gives a pointer to the internal array holding the
   * coordinates of the cell vertices.
   *( stored like this: {x0,y0,z0,x1,y1,z1,x2,y2,z3,x4,y4 ...} )
   * Pay attention that the memory area the returned pointer
   * points to gets invalid when the cell is deleted
   * \return 
   * pointer to an array of doubles
   */
  const double* getPositionData() const;

  /**
   *\par Description:
   * gives a pointer to the internal array holding the 
   * tensor components.
   * Pay attention that the memory area the returned pointer
   * points to gets invalid when 
   * either setTensors is invoked again or the cell deleted.
   *\pre
   * setTensors has been invoked
   *\return
   * pointer to an array of doubles
   */
  const double* getTensorData() const;

			    
  /**
   * help function for cell locator:
   * get face where line through start
   * in direction dir exits the cell
   *\pre
   * start is inside cell
   * \param start
   * start position in actual cell
   * \param dir
   * direction in which line goes
   *\retval length
   * parameter for line : start + length * dir = cutpoint with face
   *\return
   * Id of next face where line leaves the cell
   */
  virtual FIndex faceToGo(double&length,const FArray& start, const FArray& dir) const;

  /**
   * help function for cell locator:
   * get face id where neighbor cell lies
   * where pos could lie in
   * by computing the local coordinates
   *\pre
   * none
   * \param pos
   * position to search for
   * \retval faceId: face where better cell could lie, or invalid
   * if local coords couldn't be computed
   *\return
   * false if pos is inside,
   * true otherwise
   */
  virtual bool neighborFaceForPos(const FArray& pos,FIndex&faceId ) const;


protected:

  /** 
   *\par Description:
   *Builds the cell bounding box.
   *\pre
   *the position cache must be valid, !AND! this function assumes, that
   *the bounding box of the position cache represents the boundingbox for the
   *whole cell, if this is not the case, the concrete cells have to overload
   *this function!
   *\post
   *none
   *\exception
   *none
   */
  virtual void buildBoundingBox(void) const;

//===========================================================================
// ARBITRARY ACCESS METHODS
//===========================================================================

public:
  
  /** 
   *\par Description:
   *Grants access to the defining vertex indices.
   *Example: Tetrahedron cell = 4 vertices, (*aCell)(0) = vertex entry 0
   *         After changing the values, positionsToUpdate() and
   *         tensorsToUpdate() must be called to set the cell to a valid
   *         state for integration.
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  FIndex& operator() (positive vertexId);

  /** 
   *\par Description:
   *Grants read access to the defining vertex indices.
   *Example: Tetrahedron cell = 4 vertices, (*aCell)(0) = vertex entry 0
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  const FIndex& operator() (positive vertexId) const;
  
  /** 
   *\par Description:
   * set the numerical precision for the inside test.
   */
  void setPrecision (const double& eps);

  /** 
   *\par Description:
   * get the numerical precision for the inside test.
   */
  const double& getPrecision ();


  /** 
   *\par Description:
   *none
   *\pre
   *none
   *\post
   *none
   *\exception
   *none
   */
  void setZero_Threshold( const double& value );

  /** 
   *\par Description:
   *none
   *\pre
   *none
    *\post
   *none
   *\exception
   *none
   */
  const double& getZero_Threshold(void) const;

  void changeSetToZeroFlag() const;


  /**
   * returns approximate size of cell in memory in bytes
   */
  virtual positive memSize() const=0;

  /**
   *\par Description
   * functions for access to geometryDescription
   */  

  positive  getNumVerts( ) const;
  positive  getNumEdges( ) const;
  positive  getNumFaces( ) const;

  //to use static mem allocation(in some classes) these constants are used:

  /// maximum number of vertices per face of all cell types
  static const int maxNoVertsPerFace=4;

  ///maximum number of vertices of all cell types
  static const int maxNoVerts=8;

  ///maximum number of faces    of all cell types
  static const int maxNoFaces=6;

  ///maximum number of edges of all cell types
  static const int maxNoEdges=12;

  
  struct geoInfo{ 
    positive dim,                 //dimension of the cell vertices
      noVerts,                    //number of vertices in cell
      noEdges,                    //number of edges in cell
      noFaces,                    //number of faces in cell
    
      edges[maxNoEdges][2],       //vertex numbers of edges
      faceSizes[maxNoFaces],      //number of vertices of faces
      faces[maxNoFaces][maxNoVertsPerFace];//vertex numbers of faces
  };

  static const geoInfo & getGeometryDescription(CellType type);

  const geoInfo & geometryDescription;

protected:

    
  // the vector that contains the actual cell definition
  /**
   *\par Description
   * the vector containing the indices of the cell vertices
   * NB: This array MUST be provided in the (not empty) constructor
   * of any cell class.
   */
  FIndex * const vertexIndices;

  /** 
   *\par Description
   * pointer to cached Tensors.
   */
  FRefTensor *tensors;
  ///point on array where the data for this tensor is stored
  double*tensorData;

  /** 
   *\par Description
   * pointer to cached Positions.
   */
  FRefArray * const positions;
  ///pointer on array where the data for the FRefArrays is stored
  double*const positionData;


  /** 
   *\par Description
   * set to false when setPositions() or setTensors() is called. 
   * Used to indicate (if needed) that some intern parameters of
   * the cell (related to its geometry and/or interpolation) may
   * require a new computation. Therefore, it is NECESSARY to handle
   * them properly in the case of cell classes that need to 
   * precompute part of their information.
   */
  mutable bool geometryOK;
  mutable bool interpolationOK;

  /** 
   *\par Description
   * precision for the inside test
   */
  static double epsilon;

  /** 
   *\par Description
   * defines the lower bound for treating a tensor value as zero
   */
  static double zero_threshold;

  /**
   *\par Description
   * bounding box of the cell
   */
  mutable FBoundingBox bBox;

  /**
   *\par Description
   * if this flag is set (usually when updating the tensors),
   * no singularity is calculated, and interpolate returns a ZERO tensor
   */
  mutable bool set_to_zero;


  //------------------------------------------------------------------------

  ///data for cutWithLine/faceToGo:
  mutable FCellFaceInfo * faceInfo;

public:
  static unsigned int ntests;
  
};

//===========================================================================

ostream& operator<<(ostream& os, FCell::CellType);



#endif // __FCell_hh
/*! 
  @ingroup DataSet
  \defgroup Cells Cells (The different cell types)

  \brief
  This submodule contains the FCell class and all it subclasses that consitute the different cell types.

  \par WARNING
  Be aware that this module consists of many more classes that are not listed here because noone has added them yet.
*/

