//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FNeighborhoodDataStructured3D.hh,v $
// Language:  C++
// Date:      $Date: 2003/09/19 19:34:15 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FNeighborhoodDataStructured3D_hh
#define __FNeighborhoodDataStructured3D_hh
 
#include "FString.hh"
#include "FIndex.hh"
#include "FNeighborhoodData.hh"
#include <vector>

/** 
 * Specialization of the abstract class FNeighborhoodData in the case of 
 * structured 3D grids (i.e. rectilinear or curvilinear). The triangulated 
 * case is handled too. The structuredness of the grid enables an on-the-fly
 * computation of the connectivity information and thus permits to optimize
 * memory costs.
 */
class FNeighborhoodDataStructured3D : public FNeighborhoodData
{
  
public:
  /// Undocumented.
  FNeighborhoodDataStructured3D(positive iDim, positive jDim, positive kDim,
				bool tetrahedrized);

  /// Undocumented.
  ~FNeighborhoodDataStructured3D();
  
  /// Undocumented.
  virtual void getCellFaceNeighbor(const FIndex& srcCellId,
                                   const FIndex& faceId,
                                   FIndex& dstCellId) const;
  /// Undocumented.
  virtual void getCellEdgeNeighbors(const FIndex& srcCellId,
				    const FIndex& edgeId,
            std::vector<FIndex>& dstCellId) const;
  /// Undocumented.
  virtual void getCellVertexNeighbors(const FIndex& srcCellId,
				      const FIndex& vertexId,
              std::vector<FIndex>& dstCellId) const;
  /**
   * Get all the cells connected to a given position.
   * \param srcPosId
   * Id of the position to check for neighbors.
   * \param neighborCells
   * Storage for neighboring cells found.
   */
  virtual void 
  getPositionCellNeighbors(const FIndex& srcPosId,
      std::vector<FIndex> &neighborCells) const; 
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
  virtual void 
  getEdgeCellNeighbors(const FIndex& firstVertex, 
		       const FIndex& secondVertex, 
           std::vector<FIndex> &neighborCells) const;

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
      std::vector<FIndex> &neighborPos) const;

  positive memSize() const;

protected:

  // flag to indicate if the structured grid is handled as a tetrahedrization
  bool tetrahedrized;

  // grid dimensions (number of positions in each dimension)
  positive iDim, jDim, kDim;
};

//===========================================================================

#endif // __FNeighborhoodDataStructured3D_hh
