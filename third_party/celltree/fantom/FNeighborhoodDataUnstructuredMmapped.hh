
//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FNeighborhoodDataUnstructured.hh,v $
// Language:  C++
// Date:      $Date: 2003/10/22 15:27:40 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef __FNeighborhoodDataUnstructuredMmapped_hh
#define __FNeighborhoodDataUnstructuredMmapped_hh


#include "FNeighborhoodData.hh"
#include "FDataSet.hh"
#include "FCell.hh"
#include<vector>
#include "FanyArray.hh"

class FCell;
class FCellDefinitions;

class FNeighborhoodDataUnstructuredMmapped
  :public FNeighborhoodData
{

public:

  //optimized one
  FNeighborhoodDataUnstructuredMmapped
  (const FCellDefinitions*d, 
   shared_ptr< FanyArray< positive > > indices,
   shared_ptr< FanyArray< positive > > offsets);
  

  /**
   * Get all the cells connected to a given face of a given cell.
   * \param srcCellId
   * Id of the cell which owns the face
   * \param faceId
   * Id of the face in the given cell
   *\param dstCellId
   * Id of the cell lying on the other side of the face or invalid
   * If face is on the border
   */
  void getCellFaceNeighbor(const FIndex& srcCellId,
                                   const FIndex& faceId,
                                   FIndex& dstCellId) const;
  /// Undocumented.
  void getCellEdgeNeighbors(const FIndex& srcCellId,
                                   const FIndex& edgeId,
                                   vector<FIndex>& dstCellId) const;
  /// Undocumented.
  void getCellVertexNeighbors(const FIndex& srcCellId,
                                     const FIndex& vertexId,
                                     vector<FIndex>& dstCellId) const;
  /**
   * Get all the cells connected to a given position.
   * \param srcPosId
   * Id of the position to check for neighbors.
   * \param neighborCells
   * Storage for neighboring cells found.
   */
  void 
  getPositionCellNeighbors(const FIndex& srcPosId,
			   vector<FIndex> &neighborCells) const; 
  /**
   * Get all the cells connected to the given edge.
   * \param firstVertex
   * Id of the first vertex of the edge to check for neighbors.
   * \param secondVertex
   * Id of the second vertex of the edge to check for neighbors.
   * \param neighborCells
   * Storage for neighboring cells found.
   * \pre
   * You have to ensure (or check afterwards) that the two vertices given
   * really build up an edge of the cells found!
   */
  void 
  getEdgeCellNeighbors(const FIndex& firstVertex, 
		       const FIndex& secondVertex, 
		       vector<FIndex> &neighborCells) const;

  /**
   * Get all the points connected to a given position
   * by cell edge
   * \param srcPosId
   * Id of the position to check for neighbors.
   * \param neighborCells
   * Storage for neighboring positions found.
   */
  void
  getPositionPointNeighbors(const FIndex& srcPosId,
			   vector<FIndex> &neighborPos) const; 





  /**
   *look where the different part of the points are stored
   *hint: FIntervalArray can be used to 
   *determine the interval for a point index
   *
   *\param sizes 
   *number of points in one file
   *
   *\param filenames
   *filename for one block
   */
  void getDistribution(vector<positive>& sizes,vector<string> & names) const;


  ~FNeighborhoodDataUnstructuredMmapped();

  positive memSize() const;

private:

  /**
   * Get all the cells connected to a given position.
   * \param srcPosId
   * Id of the position to check for neighbors.
   * \param neighborCells
   * Storage for neighboring cells found.
   */
  inline void 
  getPositionCellNeighbors(positive srcPosId,
			   vector<positive> &neighborCells) const; 

  const FCellDefinitions*cellDef;

  //vector containing cell indices referenced by cellsAtP
  shared_ptr< FanyArray< positive > > data;
  
  //array of offsets in data on beginning of list of cells 
  //in data for each point
  shared_ptr< FanyArray< positive > > cellsAtP; 

  //buffer variables
  mutable vector<FIndex> vertices;
  mutable vector<positive> setbuf1,setbuf2,setbuf3;


  //vars for positionPointNeighbors

  typedef vector<positive> vIndex;

  struct ppCacheEntry {
    ppCacheEntry();
    template<class citer_positive>
    void init( citer_positive cbegin,
	       citer_positive cend);
    ~ppCacheEntry();
    positive*cells,*cellsEnd;
    vIndex points;
    positive pointId;    
  };

  typedef list<ppCacheEntry>::iterator ppCacheIter;
  static list<ppCacheEntry> EmptyList;

  //maximum entries in ppCache
  mutable positive ppMaxCacheEntries;  

  mutable positive ppNbCacheEntries;
  mutable list<ppCacheEntry> ppCache;
  mutable vector<ppCacheIter>  ppArray;




};

#endif
