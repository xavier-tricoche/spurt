//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FNode.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:08 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//---------------------------------------------------------------------------

#include "FObject.hh"
#include <vector>
#include "FIndex.hh"
#include "stdAliases.hh"
#include <utility>
#include "FCell.hh"
#include "FBoundingBox.hh"
#include "FPosition.hh"

#ifndef FNode_hh
#define FNode_hh

/*!
  This is a node for FQuadTree or FOcTree, saving the indices and
  bounding boxes of the cells contained. Furthermore, there is a
  set of pointers to the childs of the node saved. You should not use
  this class, use the frontends FQuadTree and FOctree. \n
  Definitions: \n
  A node is called a leaf if it contains no pointers to other nodes. \n
  A node is called a root node if its parent is the null pointer.
*/

class FNode : public FObject
{
public:


  /*!
    \par Description:
    Standard constructor, never used.
   */
  FNode();

  /*!
    \par Description
    Constructor for a node. \n
    newdimension = 2 : A node for a FQuadTree is generated \n 
    newdimension = 3 : A node for a FOctree is generated \n
    This method is used by addCell and the standard constructors of
    FQuadTree and FOcTree.

    \pre
    newdimension must be two or three.
    \pre
    Dimension of nodeBorders must be equal to newdimension.

    \exception
    FInvalidDimensionException : newdimension non-equal to two or three.
    \exception
    FInvalidDimensionException : nodeBorders dimension non-equal to
    newdimension.

    \param
    newdimension : Dimension of this node, must be two for a FQuadtree or
    three for a FOctree.
    \param
    maxNumberCells : Maximum Number of cells per bucket.
    \param
    newdepth : Actual depth of this node.
    \param
    newmaxdepth : Maximum depth of the FQuadTree or FOctree. If the depth of
    this node is equal to or greater than this value, the tree isn't splitted 
    anymore and all cells are stored inside this bucket.
    \param
    nodeBorders : Bounding box of this bucket.
  */
  FNode( const unsigned int& newdimension, 
	 const unsigned int& maxNumberCells,
	 const unsigned int& newdepth,
	 const unsigned int& newmaxdepth,
	 const FBoundingBox& nodeBorders );	 

  /*! 
    \par Description:
    Constructor for a whole grid. The cells are added in the order
    they appear in the newcells vector. This method is used in
    FQuadTree and FOcTree to build the cell locator.

    \pre
    newdimension must be two or three.
    \pre
    Dimension of nodeBorders must be equal to newdimension.
    \pre
    Dimension of each FBoundingBox contained in newcells must be equal
    to newdimension.

    \exception
    FInvalidDimensionException : newdimension non-equal to two or three.
    \exception
    FInvalidDimensionException : nodeBorders dimension non-equal to
    newdimension.
    \exception
    FInvalidDimensionException : One of the cells contained in the 
    newcells vector doesn't match the dimension precondition.

    \param
    newdimension : Dimension of the nodes, must be two for a FQuadtree or
    three for a FOctree.
    \param
    maxNumberCells : Maximum Number of cells per bucket.
    \param
    newmaxdepth : Maximum depth of the FQuadTree or FOctree. If the depth of
    one of the build nodes is equal to this value, the tree isn't splitted 
    anymore and all cells are stored inside this bucket.
    \param
    newcells : Cell indices and their suitable bounding boxes. 
    \param
    nodeBorders : Main bounding box of the structure.
   */
  FNode( const unsigned int& newdimension, 
	 const unsigned int& maxNumberCells,
	 const unsigned int& newmaxdepth,
	 const vector< pair< FIndex, FPosition > >& newcells,
	 const FBoundingBox& nodeBorders);

  /*!
    \par Description:
    constructor for load
  */
  FNode( char* filename );

  /*! 
    \par Description:
    Destructor for the whole underlying subtree.
   */
  ~FNode();

  /*! 
    \par Description:
    Adds a cell to the tree. An overflow is handled by splitting the cell,
    see split.

    \pre
    Dimension of the FBoundingBox must be equal to the dimension of the 
    tree.
    \pre 
    The bounding box of the newcell has to intersect the main bounding
    of the tree.

    \exception
    FInvalidDimensionException : Dimension of newcells bounding box is
    non-equal to the tree dimension.
    \exception
    FException : addCell was called from outside, cell lies outside
    the main bounding box of the tree

    \param
    newcell : cell to be added to the tree
   */
  void addCell( const pair< FIndex, FPosition >& newcell );
  
  /*! 
    \par Description:
    The tree is searched for a cell.

    \param
    pos : The FPosition where the corresponding bucket is wanted.

    \return
    A vector which contains the indices and bounding boxes of all cells
    inside the bucket containing pos.
   */
  vector< pair< FIndex, FPosition > > search( const FPosition& pos ) const;
  
  /*! 
    \par Description:
    Prints recursively informations about the tree. Each nodes (buckets)
    bounding box, depth, isLeaf status (boolean), number of childs and number 
    of cells is printed. This method is made for debugging reasons
    only.
    */
  void print() const;


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

  /*!
    \par Description:
    Returns the maximum number of cells of the tree.
  */
  void maxBucketSize( unsigned int& size );

private:
  /*! 
    \par Description:
    Checks if pos is inside the nodes bbox.
    \return
    true : pos is inside bounding box of the node
    \return
    false : pos is outside bounding box of the node
   */
  bool isInsideBBox( const FPosition& pos ) const;

  /*! 
    \par Description:
    Checks if cell is inside nodes bbox.
    \return
    true : cell intersects bounding box of the node
    \return
    false : cell doesn't intersect bounding box of the node
   */
  bool isInsideBBox( const pair< FIndex, FPosition >& cell ) const;

  /*! 
    \par Description:
    returns all contained cells in this node
    \return
    A vector containing cell indices and their bounding boxes.
   */
  vector< pair< FIndex, FPosition > > returnCells() const;

  /*! 
    \par Description:
    Checks if node is a leaf.
    \return
    true : node is a leaf.
    \return
    false : node is not a leaf (inner node).
   */
  bool isLeaf() const;

  /*! 
    \par Description:
    Checks if the maximum number of cells is reached inside the bucket.
    \return
    true : Number of cells inside the bucket is equal or greater than the 
    value given to the constructor. Adding cells will cause a split or
    overflow if the maximum depth of the tree is reached.
    \return
    false : Number of cells inside the bucket is lesser than the value given
    to the constructor. Cells may be added to this node.
   */
  bool maximumReached() const;

  /*! 
    \par Description:
    Access to the parent of the bucket. If a null pointer is returned,
    the node is a root node.
    \return
    A pointer to the parent of the node.
   */
  FNode* getParent() const;

  /*! 
    \par Description:
    Changes the pointer to the parent of the actual node. This method
    is used within split.
    
    \pre
    newparent is non-equal to a null pointer. 

    \exception
    FNullPointerAssignmentException : newparent is a null pointer.

    \param
    newparent : pointer the new parent node.
   */
  void setParent( FNode* newparent );

  // adds a child 
  void addChild( FNode* newchild );

  /*! 
    \par Description:
    Splits the bucket if a cell is added and the maximum number of cells
    inside a node is reached. Depending on the dimension of the tree
    a node is split into 4 (two-dimensional case) or 8 (three-dimensional
    case) subnodes. split is not called if the maximum tree depth is
    reached, in this case the cells are added to the lowest level nodes.
   */
  void split();

  /*!
    \par Description:
    dumps tree information to a file
  */
  void dumpTree( ofstream& out );


  // ATTRIBUTES

  //! dimension of the node
  unsigned int dimension;

  //! maximum Number of Cells
  unsigned int maxCells;

  //! actual depth
  unsigned int depth;

  //! maximal depth
  unsigned int maxDepth;

  //! bounding box of the node
  FBoundingBox nodeBox;

  //! bounding boxes and indices of cells in the bucket
  vector< pair< FIndex, FPosition > > cells;
  
  // Quadtree : 4 childs 
  // 
  //        ----------- 
  //      /  0  |  3  / 
  //      -----------   
  //    /  1  |  2  /   
  //    -----------     
  // 
  // Octree   : 8 childs 
  // 
  //        ----------- 
  //      /  4  |  7  / 
  //      -----------   
  //    /  5  |  6  /   
  //    -----------     
  //        ----------- 
  //      /  0  |  3  / 
  //      -----------   
  //    /  1  |  2  /   
  //    -----------     
  //

  //! pointer to the nodes parent  
  FNode* parent;

  //! pointers to childs
  vector<FNode*> childs;
};

#endif
