#ifndef __FkdTree_hh
#define __FkdTree_hh

#include "stdAliases.hh"
#include<boost/shared_ptr.hpp>
#include "FArray.hh"
#include "stdAliases.hh"
#include "FanyArray.hh"
#include <math.h>
#include "FBoundingBox.hh"

using namespace boost;

class FPositionSet;
class FNeighborhoodData;
class FIndex;
//class FBoundingBox;
//#define HUGE_VAL 0xffffff
/*
 * this class provides a kdtree 
 * (at tree where the positions are recursively split at the median
 * resp. the different coordinates in the positionset)
 * where the nodes are NOT elements of the positionset,
 * but only the leaves.
 * if there is more than one element (position index) in a leaf,
 * the first element is taken as leaf entry.
 * also, if multiple elements have one coordinate equal,
 * there is no guarantee, that one element is excatly re-found.
 * so, it is ONLY guaranteed, that an element near the searched element is found
 */

class FkdTree
{ 
public:
  //forward declaration
  struct nodeBlock;

  /**
   * parameters for use with distriubuted mmapped array:
   *
   *\param Atree
   * array which will  contain the inner nodes of the tree,
   * if not zero, it should have at least
   * ((2^levelsPerBlock)^ (ceil(numlevels/levelsPerBlock)) - 1 )/(2^levelsPerBlock -1)
   *( = \sum_{i=0}^{\lceil numlevels/levelsPerBlock \rceil -1 } (2^levelsPerBlock)^i )
   * entries.
   *
   *\param Aleaves
   * array which will contain 
   * the indices of the points contained in the kdtree leaves, if set,
   * it should have 2^numLevels entries 
   *
   *parameters for build kdtree only in part of the above arrays:
   *\param nLevels
   * levels the should be buiild below the root block
   *
   *\param rBlock
   * the root block will not be written into
   * Aleaves[0], but into rBlock.
   *
   *\param startBlockId :
   * blockId where the first block will be written to 
   *(excluding the root block)
   *
   *\param Abmin \param Abmax :
   * kdtree box from upper processes
   *
   */
  FkdTree( const FPositionSet* pos,
		    shared_ptr< FAssignableAnyArray<nodeBlock> > Atree,	   
		    shared_ptr< FAssignableAnyArray<positive> > Aleaves,
		    shared_ptr< FAssignableAnyArray<double> >   AposList,
		    shared_ptr< FAssignableAnyArray<positive> > AposListIndices,
		    positive startBlockId,
		    nodeBlock & rBlock,
		    const double * Abmin,
		    const double * Abmax,	  
		    positive nLevels
		    );

  /**
   * constructor that builds 
   * kdtree with floor(log2(number of positions))-lessLevels
   * levels
   *\pre: 
   * positionset is valid
   *
   * \param pos:
   * positionSet from which kdTree is built
   *
   * \param lesslevels:
   * nb of leaves = 2^( ceil(log2(nbPos)) - lesslevels )
   *
   *\param roundBottom:
   * determines if number of leaves should be
   * max{x: x=2^levels, x<=nbPos} (if true)
   * or 
   * min{x: x=2^levels, x>=nbPos} (if false)
   *
   *\param AposList
   * array for optimization
   * which describes the coordinates in the positionset
   * in the order x0,y0,z0,x1,y1,z1,... in 3d
   * and x0,y0,x1,y1,... in 2d 
   */

  FkdTree( const FPositionSet* pos, 
	   const FNeighborhoodData* dat = 0,
	   positive lessLevels = 0,
	   bool roundBottom = false,
	   shared_ptr< FAssignableAnyArray<double> > 
	   AposList  = shared_ptr< FAssignableAnyArray<double> >(),	   
	   shared_ptr< FAssignableAnyArray<positive> > 
	   AposListIndices = shared_ptr< FAssignableAnyArray<positive> >()
	   );

  /**
   *constructor for taking an already existing kdtree,
   * not building it
   *\param pos:
   * positionset for which kdtree was built
   * and on which the indices in leaves point
   *
   *\param tree
   * array with blocks containing the inner nodes of the tree
   * ( first: the highest level, then:
   * all blocks of the next level, and so on
   * so it has 1+ 2^levelsPerBlock + 2^(2*levelsPerBlock)+ ... entries
   * = ( (2^levelsPerBlock)^( 1 + ceil(numlevels/levelsPerBlock) ) - 1 )/(2^levelsPerBlock -1)
   * entries.
   *
   *\param leaves
   * indices of the points contained in the kdtree leaves
   * ( 2^numLevels  entries )
   *
   *\param numlevels
   * number of recursion levels of the tree,
   * should be <= ceil(log2(number of positions))
   *
  */
  FkdTree( const FPositionSet* pos, 
	   shared_ptr< FAssignableAnyArray<nodeBlock> > tree,	   
	   shared_ptr< FAssignableAnyArray<positive> > leaves,
	   positive nLevels
	  );
  //variables in search structure
  /** 
   * destructor
   */
  ~FkdTree();

  /**
   * gets entry of kdtree leaf which contains the position p
   * if inside bounding box, invalid index otherwise
   *\param p
   * position to search for
   *\return
   * entry of kdtree leaf which contains p
   */
  FIndex getEntry( const FPosition&p ) const;


  /**
   * gets entry of kdtree leaf which contains the position p
   * if inside bounding box, invalid index otherwise
   * and the box of the found kdtree leaf
   *\param p
   * position to search for
   *\param minbox:
   * coords of lower left back corner of box
   *\param maxbox:
   * coords of upper top front  corner of box
   *\return
   * entry of kdtree leaf which contains p
   */
  FIndex getEntryAndBox( const FPosition&p ,
  			 double*minbox,double*maxbox) const;

  /**
   * gets number of leaf which contains the position p
   *\param p
   * position to search for
   *\return
   * entry of kdtree leaf which contains p
   */

  positive getLeafNumber( const double*p ) const;


  /**
   *gets the box of kdtree leaf with number i
   *\retval minbox:
   *minimum box coords
   *\retval maxbox:
   *maximum box coords
   *\return
   * box of kdtree leaf with number i
   */
  void getLeafBox( positive i,double minbox[3],double maxbox[3],positive level=0) const;



  /**
   * gets the entries of those kdtree leaves
   * which intersect the box b
   *\param b
   * bounding box to find leaves intersecting it 
   *\retval indices
   * entries of the found leaves
   */
  void getEntries( const FBoundingBox &b, vector<FIndex>&indices ) const;


  /**
   * gets the entries of those kdtree leaves
   * which intersect the 
   * implicit surface formed by f(x)=0
   *\param f 
   * function object which has the signature
   * double f(double*), e.g. an instance of plane3d
   * described below
   * 
   *\param maxdist
   * maximum size of boxes after split
   *( should be smaller than the minimum  
   * size of structures of the surface )
   *\param minlevels
   *minimum levels gone down before lloking at the splits
   *\retval indices
   * entries of the kdtree leaves found
   * which intersect the surface
   */
  template< class functionScalarKD >
  void getEntriesIntersecting( const functionScalarKD & f,vector<FIndex>&indices ,
			       double epsilon=1e-9, int maxlevels=18 , 
			       double maxdist=HUGE_VAL ) const;

  //plane for use with getEntriesIntersecting
  struct plane3d{
    //distance to zero
    double d;
    //normal of plane
    double n[3];
    inline double operator() (const double*x) const
    { return x[0]*n[0]+x[1]*n[1]+x[2]*n[2]-d; }    
  };


  /**
   * replaces the position indices in the leaves
   * by the first neighbor cell index returned by the neighborhood
   * \pre this is the first invocation
   * \post 
   * the indices returned by the getEntry or getEntries functions
   * are now cell indices and not position indices
   * \param dat
   * neighborhood used to return the neighboring indices to position
   */
  void replacePosIndByCellInd( const FNeighborhoodData*dat );  


  /**
   * for ShowkdTree visalgo
   */
  shared_ptr< const FAssignableAnyArray<nodeBlock> >  getTree(int&levels) const;

  /**
   * for cellLocatorNeighbors
   */
  shared_ptr< const FAssignableAnyArray<positive> > getLeaves() const;

  positive memSize() const;

  /**
   *gets number of blocks in the array tree
   *when the kdtree has levels levels
   *
   *\param levels
   *number of levels in a kdtree
   *
   *\return
   *number of blocks needed in array tree
   */
  static positive numberOfBlocks(unsigned levels);

  static const positive levelsPerBlock = 6 ;

  struct nodeBlock{
    //for memory alignment, one entry more than needed is stored
    char coordIds[(1<<levelsPerBlock)];
    double coordVals[1<<levelsPerBlock];    
  };

  /// class for going through the tree; defined below
  class kdTreeIter;

private:



  //temporary variables used for computeTreeData::computeTree

  const FPositionSet* posSet;

  mutable double* aTreePos;
  mutable char* aTreeCoord;
  mutable FAssignableAnyArray<positive>::iterator aLeafPos;
  mutable int aLevels;
  mutable FAssignableAnyArray<nodeBlock>::iterator aBlockId;
  
  //helper classes for constructor
  struct computeTreeData;

  //help function for computetree:
  //gets bbox of 1000 random positions in [left,right)
  //writes minbbox,maxbbox
  void getRandomBBox(positive left,positive right);


  //splits list into two parts (at middle) 
  //with (in coordinate coord):
  // [left,middle)<=pivot<=[middle,right)
  //return: value of pivot 
  double split(positive left,positive middle,positive right,positive coord);






  //-----------------------

  //buffer variables used for getBBoxIndices

  mutable vector<FIndex>*indexList;
  mutable double minbbox[3],maxbbox[3];//also for getRandomBBox

  /* getBBoxIndices is used in getEntries
   * \pre
   * - minbbox,maxbbox is set to the bbox to search
   * - indexList is set to the vector to fill with indices
   * - indexList points to an empty Vector
   */
  void getBBoxIndices(kdTreeIter & it) const;

  template< class functionScalarKD>
  friend class dataForGetEntries;


  //variables used in constructor

  FAssignableAnyArray<double> * posList;
  FAssignableAnyArray<positive> * posListIndices;



  //--------------------------------------

  FBoundingBox *box;

  int levels; //2^levels entries in kdTree
  int dim;

  //variables in search structure
  shared_ptr< FAssignableAnyArray<nodeBlock> > tree;

  shared_ptr< FAssignableAnyArray<positive> > leaves;

  //vars for progress indication
  mutable double numcomp,gescomp;
  mutable positive oprogress;

};


class FkdTree::kdTreeIter
{

  static const positive lpb
  = FkdTree::levelsPerBlock;

  const FkdTree * tree;

  positive aHeight,
    aLeafId,
    aBlockId,
    blocksPerSubtree,
    //IB means "In Block" ,
    //so these are block-local variables
    subtreeSizeIB,       
    nodeIdIB, leafIdIB;

  inline static unsigned 
  countbits1(register unsigned long x);
    
public:

   kdTreeIter(const FkdTree*tree);    

  /** 
   *go to the left sub -node or -leaf
   *\pre top() != true
   */
   void left();

  /**
   *\pre top() != true
   *go to the right sub -node or -leaf
   */
   void right();

  /**
   *\pre bottom() != true
   * to the the father node
   */
   void up();

  /**
   *get coord value of node
   */
   double coordVal() const;

  /**
   *get coord number at which the split goes
   */
   positive coordId() const;

  /**
   * get actual leaf id;
   * if iter is not at bottom,
   * the smallest leaf id in the subtree
   * is returned
   */
   positive leafId() const;


  /**get actual height in the tree
   *(height=0 -> is at a leaf,
   *(height=tree->levels -> is at root )
   */
   positive height() const;

  ///is iterator at the root node ?
   bool top() const;

  /// is iterator at a leaf ?
   bool bottom() const;

   bool operator != (const kdTreeIter &x) const;

        
};

template< class functionScalarKD>
class dataForGetEntries{


  const FkdTree&kd;

  FkdTree::kdTreeIter iter;

  const functionScalarKD & f;
    
  vector<FIndex> & indexList;
    
  double eps;
  int minAL;
  double maxD;
    
  void getLeavesIntersecting(double*minbb,double*maxbb,
			     positive plus,positive minus);
public:
    
  dataForGetEntries(const FkdTree&k,const functionScalarKD & func,vector<FIndex> & indices,
		    double epsilon,int maxlev,double maxdist);
    
};
  
template< class functionScalarKD >
void FkdTree::getEntriesIntersecting( const functionScalarKD & f,
				      vector<FIndex>&indices ,double epsilon,
				      int maxlevels , double maxdist ) const
{  
  dataForGetEntries < functionScalarKD > 
    (*this,f,indices,epsilon,maxlevels,maxdist);  
}
  
  
// helper class for getEntries 
template< class functionScalarKD >
dataForGetEntries<functionScalarKD>::
dataForGetEntries(const FkdTree&k,
		  const functionScalarKD & func,vector<FIndex>&indices,
		  double epsilon,int maxlev,double maxdist)
  : kd(k),iter(&k),f(func),indexList(indices),eps(epsilon) ,maxD(maxdist)
{ 
  double minbbox[3],maxbbox[3];
  kd.box->getRange(minbbox[0],maxbbox[0],
		   minbbox[1],maxbbox[1],
		   minbbox[2],maxbbox[2]);
    
    
  minAL= iter.height()-maxlev;
    
  indices.clear();    
    
  positive plus=0,minus=0;
  positive pdim = 1<<kd.dim;
    
  for(positive i=0 ; i < pdim ; i++){
    double coords[3]={ (i&1)?maxbbox[0]:minbbox[0],
		       (i&2)?maxbbox[1]:minbbox[1],
		       (i&4)?maxbbox[2]:minbbox[2] };
    double wert=f(coords);
    plus|=positive(wert > eps)<<i;
    minus|=positive(wert < -eps)<<i;
      
  }

  getLeavesIntersecting(minbbox,maxbbox,minus,plus);
}
  

#define infnorm(a,b)\
 ( fabs(b[0]-a[0])+fabs(b[1]-a[1])+fabs(b[2]-a[2]) )
  
template< class functionScalarKD>
void dataForGetEntries<functionScalarKD>::getLeavesIntersecting
( double * minbb, double * maxbb ,
  positive minus,positive plus
  )
{
    
  positive coord = iter.coordId();
  double coordVal= iter.coordVal();
    
    
  positive 
    powc = 1<<coord, // 0tes, 1stes oder 2tes bit 1
    mpowd = (1<<(1<<kd.dim))-1, // 2^(2^dim) -1, also 2^dim bits auf 1
    pdim = 1<<(kd.dim);
    
  positive nminus=0,nplus=0;

  //calculate signs of values at new cutpoints (of splitDim with bbox)
  for(positive i=0 ; i < pdim ; i++)
      
    if( (i&powc) == 0 ){
	
      double coords[3]={ i&1?maxbb[0]:minbb[0],
			 i&2?maxbb[1]:minbb[1],
			 i&4?maxbb[2]:minbb[2] };
	
      coords[coord]=coordVal;
	
      double wert = f(coords);
	
      nplus|=positive(wert > eps)<< i;
	
      nminus|=positive(wert < -eps)<< i;
	
    }
    
    
  //if 1st box part contains surface
    
  static const positive
    upbitsc[3]={0xAA,0xCC,0xF0};
    
  positive 
    upbits = upbitsc[coord],
    downbits = ~ upbits;
    
  bool dontCheckBits = iter.height()>1 && int(iter.height()) > minAL && infnorm(minbb,maxbb) > maxD ;

  positive mA,pA,mB,pB;

  //  cout<<"first"<<powc<<' '<<upbits<<' '<<(downbits&0xff)<<endl;
  //cout<<minus<<' '<<nminus<<' '<<plus<<' '<<nplus<<endl;

  mA = (minus & downbits) | (nminus << powc);
  pA = (plus  & downbits) | (nplus  << powc);    

  //cout<<mA<<' '<<pA<<endl;
    
  //if 1st box part contains surface
  bool bitCheckA = dontCheckBits || 
    ( (mA == 0) & (pA == 0) ) || 
    ( (mA != 0) & (mA != mpowd) ) || 
    ( (pA != 0) & (pA != mpowd) ) ;

  mB = (minus & upbits) | nminus ;
  pB = (plus  & upbits) | nplus ;    
  //if 2nd box part contains surface
  bool bitCheckB = dontCheckBits || 
    ( (mB == 0) & (pB == 0) ) || 
    ( (mB != 0) & (mB != mpowd) ) || 
    ( (pB != 0) & (pB != mpowd) ) ;



  if( iter.height() > 1){

    double newbb[3];

    //cout<<" 1st"<<endl;
    if(bitCheckA)
      {
	memcpy(newbb,maxbb,sizeof(newbb));
	newbb[coord]=coordVal;       
	
	iter.left();
	getLeavesIntersecting(minbb,newbb,mA,pA);
	iter.up();
      }
    //cout<<" 2nd"<<endl;
    if(bitCheckB)
      {
	memcpy(newbb,minbb,sizeof(newbb));
	newbb[coord]=coordVal;       

	iter.right();
	getLeavesIntersecting(newbb,maxbb,mB,pB);
	iter.up();
      }
  }
    
  else{    

    positive a=(*kd.leaves)[iter.leafId()];
    positive b=(*kd.leaves)[iter.leafId()+1];

    if(bitCheckA)
      indexList.push_back(FIndex(a));

    if( b!=a && bitCheckB)
      indexList.push_back(FIndex(b));

  }

  
}

#ifndef OUTLINE
#include "FkdTree.icc"
#endif


#endif //__FkdTree_hh
