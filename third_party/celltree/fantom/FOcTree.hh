//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FOcTree.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:08 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//---------------------------------------------------------------------------

#include "FPosition.hh"
#include "FNode.hh"
#include "FIndex.hh"
#include "stdAliases.hh"
#include "FCell.hh"
#include "FPosition.hh"

#ifndef FOcTree_hh
#define FOcTree_hh

/*!
  This is a frontend for FNode for use as a three-dimensional cell locator.
  This class is used in 3D unstructured grid types.
*/

class FOcTree
{
public:
  /*!
   \par Description:
    default constructor
   */
  FOcTree();

  /*! 
    \par Description:
    value constructor
   */
  FOcTree( const FBoundingBox& mainBBox,
	   const positive& newmaxDepth,
	   const positive& newmaxCells );

  /*! 
    \par Description:
    value constructor
   */
  FOcTree( const FBoundingBox& mainBBox,
	   const positive& newmaxDepth,
	   const positive& newmaxCells,
	   const vector< pair< FIndex, FPosition > >& newcells );

  /*! 
    \par Description:
    value constructor
   */
  FOcTree( const FBoundingBox& mainBBox,
	   const positive& newmaxDepth,
	   const positive& newmaxCells,
	   const vector<FPosition>& newcells );


  /*!
    \par Description:
    destructor
   */
  ~FOcTree();

  /*!
    \par Description:
    adds a cell to the octree
   */
  void addCell(const pair< FIndex, FPosition >& cell );

  /*!
    \par Description:
    returns the cell index at position pos
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

  //! maximum depth of the octree
  positive maxDepth;

  //! maximum number of cells per bucket
  positive maxCells;

  //! dimension
  positive dimension;
};

#endif
