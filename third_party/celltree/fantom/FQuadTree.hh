//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FQuadTree.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:11 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//---------------------------------------------------------------------------

#include "FBoundingBox.hh"
#include "FNode.hh"
#include "FIndex.hh"
#include "stdAliases.hh"
#include "FCell.hh"
#include "FPosition.hh"

#ifndef FQuadTree_hh
#define FQuadTree_hh

/*!
  This is a frontend for FNode. This class is used in 2D unstructured grid 
  types like FGrid2DArbitrary or FGrid2DCurvilinear.
*/

class FQuadTree
{
public:
  /*!
   \par Description:
    default constructor
   */
  FQuadTree();

  /*! 
    \par Description:
    value constructor
   */
  FQuadTree( const FBoundingBox& mainBBox,
	     const positive& newmaxDepth,
	     const positive& newmaxCells );

  /*!
    \par Description:
    value constructor
   */
  FQuadTree( const FBoundingBox& mainBBox,
	     const positive& newmaxDepth,
	     const positive& newmaxCells,
	     const vector< pair< FIndex, FPosition > >& newcells );

  /*!
    \par Description:
    value constructor
   */
  FQuadTree( const FBoundingBox& mainBBox,
	     const positive& newmaxDepth,
	     const positive& newmaxCells,
	     const vector<FPosition>& newcells );

  /*! 
   \par Description:
    destructor
  */
  ~FQuadTree();

  /*! 
    \par Description:
    adds a cell to the quadtree
   */
  void addCell(const pair< FIndex, FPosition >& cell );

  /*! 
   \par Description:
   returns the cell indices at position pos
   */
  vector<FIndex> search( const FPosition& pos );

  /*! 
   \par Description:
   returns the cell indices and the matching bounding boxes at position pos
   */
  vector< pair< FIndex, FPosition > > searchboxes( const FPosition& pos );

  /*! 
   \par Description:
   prints debugging information
   */
  void print();

  /*!
    \par Description:
    returns maximum bucket size
  */
  unsigned int maxBucketSize();
    
  /*!
    \par Description:
    Saves the tree into the given filename.
  */
  void save( char* filename );

  /*!
    \par Description:
    Loads the tree from the given filename.
  */
  void load( char* filename );


private:
  //! pointer to the root node
  FNode* root;

  //! maximum depth of the quadtree
  positive maxDepth;

  //! maximum number of cells per bucket
  positive maxCells;

  //! dimension
  positive dimension;
};

#endif
