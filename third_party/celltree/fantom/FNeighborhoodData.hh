//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FNeighborhoodData.hh,v $
// Language:  C++
// Date:      $Date: 2003/09/19 19:34:15 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FNeighborhoodData_hh
#define __FNeighborhoodData_hh
 
#include "FIndex.hh"

#include <vector>

#include "FObject.hh"

/** 
 * Abstract class that handles the neighborhood information of any grid type
 */
class FNeighborhoodData : public FObject
{
  
public:
  /// Undocumented.
  FNeighborhoodData();
  /// Undocumented.
  virtual ~FNeighborhoodData();
  
  /** 
   * get the cell which is connected over the face No. faceId
   * with the actual cell,if possible
   *\param srcCellId: 
   * id of actual cell
   *\param faceId: 
   * id of face of actual Cell
   *\retval dstCellId: 
   * id of the cell connected with actual cell with face faceid
   * or invalid if no other cell is at this face
   */

  virtual void getCellFaceNeighbor(const FIndex& srcCellId,
                                   const FIndex& faceId,
                                   FIndex& dstCellId) const = 0;

  /** 
   * get the cells which are connected over the edge No. edgeId
   * with the actual cell and which have no common face with actual cell
   *\param srcCellId: 
   * id of actual cell
   *\param edgeId: 
   * id of edge of actual Cell
   *\retval dstCellId: 
   * ids of the cells connected with actual cell with edge edgeid
   */
  virtual void getCellEdgeNeighbors(const FIndex& srcCellId,
                                   const FIndex& edgeId,
                                   std::vector<FIndex>& dstCellId) const = 0;
  /** 
   * get the cells which are connected over the vertex No.vertexId
   * with the actual cell and which have no common face or edge with actual cell
   *\param srcCellId: 
   * id of actual cell
   *\param vertexId: 
   * id of vertex of actual Cell
   *\retval dstCellId: 
   * ids of the cells connected with actual cell with vertex vertexid
   */
  virtual void getCellVertexNeighbors(const FIndex& srcCellId,
                                     const FIndex& vertexId,
                                     std::vector<FIndex>& dstCellId) const = 0;
  /**
   * Get all the cells connected to a given position.
   * \param srcPosId
   * Id of the position to check for neighbors.
   * \retval neighborCells
   * Storage for neighboring cells found.
   */
  virtual void 
  getPositionCellNeighbors(const FIndex& srcPosId,
      std::vector<FIndex> &neighborCells) const = 0; 
  /**
   * Get all the cells connected to the given edge.
   * \param firstVertex
   * Id of the first vertex of the edge to check for neighbors.
   * \param secondVertex
   * Id of the second vertex of the edge to check for neighbors.
   * \retval neighborCells
   * Storage for neighboring cells found.
   * \pre
   * You have to ensure (or check afterwards) that the two vertices given
   * really build up an edge of the cells found!
   */
  virtual void 
  getEdgeCellNeighbors(const FIndex& firstVertex, 
		       const FIndex& secondVertex, 
           std::vector<FIndex> &neighborCells) const=0;


  /**
   * Get all the points connected to a given position
   * by cell edge
   * \param srcPosId
   * Id of the position to check for neighbors.
   * \param neighborCells
   * Storage for neighboring positions found.
   */
  virtual void
  getPositionPointNeighbors(const FIndex& srcPosId,
      std::vector<FIndex> &neighborPos) const=0;


  virtual positive memSize() const = 0;
};

//===========================================================================

#endif // __FNeighborhoodData_hh
