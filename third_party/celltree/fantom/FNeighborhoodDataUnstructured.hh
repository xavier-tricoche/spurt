
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

#ifndef __FNeighborhoodDataUnstructured_hh
#define __FNeighborhoodDataUnstructured_hh

#include "FNeighborhoodData.hh"
#include "FDataSet.hh"
#include "FCell.hh"
#include<vector>

class FCell;
class FCellDefinitions;

class FNeighborhoodDataUnstructured
  :public FNeighborhoodData
{

public:

  ///constructor
  FNeighborhoodDataUnstructured(const FCellDefinitions*d);
  virtual ~FNeighborhoodDataUnstructured();

  //optimized one
  FNeighborhoodDataUnstructured 
  (const FCellDefinitions*d, 
   const vector<pair<FCell::CellType,unsigned int> >& cellTypes,
   const vector<FIndex> & cellIndices);


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
  /**
   * Get all the cells connected to a given edge of a given cell.
   * \param srcCellId
   * Id of the cell which owns the edge
   * \param faceId
   * Id of the edge in the given cell
   *\param dstCellId
   * Ids of the cells lying on the other side of the edge or invalid
   * If there are no edge neighbors
   */
  void getCellEdgeNeighbors(const FIndex& srcCellId,
                                   const FIndex& edgeId,
                                   vector<FIndex>& dstCellId) const;
  /**
   * Get all the cells connected to a given vertex of a given cell.
   * \param srcCellId
   * Id of the cell which owns the vertex
   * \param faceId
   * Id of the vertex in the given cell
   *\param dstCellId
   * Ids of the cells lying on the other side of the vertex or invalid
   * If there are no vertex neighbors
   */
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






  positive memSize() const;

private:

  const FCellDefinitions*cellDef;

  //vector containing cell indices referenced by cellsAtP
  vector< positive > data;
  
  //array of pointers on beginning of list of cells 
  //in data for each point
  vector< positive * > cellsAtP; 

  //buffer variables
  mutable vector<FIndex> vertices;
  mutable vector<positive> setbuf1,setbuf2,setbuf3;


  //vars for positionPointNeighbors

  typedef vector<positive> vIndex;

  struct ppCacheEntry {
    ppCacheEntry();
    void init(const positive*cbeg,const positive*cend);
    ~ppCacheEntry();
    positive*cells,*cellsEnd;
    vIndex points;
    positive pointId;    
  };

  typedef list<ppCacheEntry>::iterator ppCacheIter;
  //typedef ppCacheEntry *ppCacheIter;

  //maximum entries in ppCache
  mutable positive ppMaxCacheEntries;  

  mutable positive ppNbCacheEntries;
  mutable list<ppCacheEntry> ppCache;
  mutable vector<ppCacheIter>  ppArray;

  static std::list<ppCacheEntry> EmptyList;


};

#endif
