#include "FkdTree.hh"

#include "FPositionSet.hh"
#include "FNeighborhoodData.hh"
#include "FStatusReporter.hh"
 
#include <algorithm>

#include <iostream>
#include <cassert>
#include <cmath>

#ifdef OUTLINE
#include "FkdTree.icc"
#endif

#include "FIndexPositionIter.hh"

#ifdef WIN32
// drand48()
#include "util.hh"
#endif

// #include <timer.hpp>
// #include <memstat.hpp>
#include <cstdio>

struct FkdTree::computeTreeData{
  
  nodeBlock actBlock;
  
  FkdTree * tree;
  
  positive restLevels;    
  
  //iterators for data in actBlock
  double* aTreePos;
  char* aTreeCoord;
  
  //preconditions of global  vars:
  // - posList is valid and indices 
  //   are set increasingly from 0 to (number of positions)-1
  // - actualTreePos is set to tree
  // - actualLeafPos is set to leaves
  // - aLevels is set to 0
  
  // makes kdtree of interval [left,right) in posList
  // and enters it into the lists tree and leaves
  // at the correct position
  void computeTree(positive left,positive right,const double *bbmin,const double*bbmax);
  
};



positive 
FkdTree::numberOfBlocks(unsigned l)
{


  positive 
    n=levelsPerBlock,
    nn=1<<n,
    r = l%n;

  // 1 rootBlock 
  // + 2^(l \mod n) bloecke mit \sum_{i=0}^{l/n-1} nn^i 
  // einzelbloecken
  return 
    unsigned(r!=0)   
    + ((( (1<<(l-r)) -1 ) / (nn-1) ) << r);			     
  

}

//-----------------------------------------------------------------------------

FkdTree::FkdTree( const FPositionSet* pos, 
		  shared_ptr< FAssignableAnyArray<nodeBlock> > Atree,	   
		  shared_ptr< FAssignableAnyArray<positive> > Aleaves,
		  positive nLevels
		  )
  : posSet(pos),levels(nLevels),tree(Atree),leaves(Aleaves)
{
  dim = pos->getDimension();
  box = new FBoundingBox(pos->getBoundingBox());
    

  if( Atree->size() != numberOfBlocks(levels) ){
    cerr<<" levels"<<levels<<" given size:"<<Atree->size()
	<<" expected size:"<<numberOfBlocks(levels)<<endl;
    THROW_EXCEPTION(FException,"invalid array size for Atree");
  }
  
  if( Aleaves->size() != 1LU<< levels )
    THROW_EXCEPTION(FException,"invalid array size for ALeaves");
  


}

//-----------------------------------------------------------------------------

FkdTree::FkdTree( const FPositionSet* pos, const FNeighborhoodData* dat,
		  positive lessLevels,bool roundBottom,
		  shared_ptr< FAssignableAnyArray<double> >   AposList,
		  shared_ptr< FAssignableAnyArray<positive> >   AposListIndices
		  )
  : posSet(pos)
{
    theStatusRep->progressBegin( "creating kdtree" );

  dim = pos->getDimension();
  box = new FBoundingBox(pos->getBoundingBox());
  

  positive noPos = pos->getNbPositions();

  if(!AposListIndices)
    AposListIndices.reset(new FanyVector<positive>(noPos));

  if(!dat)
    for(positive i=0;i<noPos;i++)
      (*AposListIndices)[i]=i;

    

  if(!dat && !AposList){
    AposList.reset(new FanyVector<double>(noPos*dim));
    FPosition p(dim);
    positive k=0;
    for(positive i=0;i<noPos;i++){
      posSet->getPosition(p,i);
      for(int j=0;j<dim;++j,++k)
	(*AposList)[k]=p[j];
      
    }
  }
  

  if(dat){
    //if neighborhoodData given,
    //include only positions which are vertices of cells

    positive posConnected=0;

    vector<FIndex> cells;
    

    for(positive i=0;i!=noPos;i++)
      {
	dat->getPositionCellNeighbors(i,cells);

	if(cells.size()!=0){
	  (*AposListIndices)[posConnected]=i;
	  posConnected++;
	}
	  
      }

    noPos = posConnected;

    if(!AposList)
      AposList.reset(new FanyVector<double>(noPos*dim));
	
    FPosition p;
    positive k=0;

    for(positive i=0;i!=noPos;i++){
      posSet->getPosition(p,(*AposListIndices)[i]);

      for(int j=0;j<dim;j++,k++)
	(*AposList)[k]=p[j];
    }
      
  }

  posList = AposList.get(); posListIndices = AposListIndices.get();
  
  double bmin[3],bmax[3];

  //determine levels so that 2^levels >= noPos
  levels=0;
  while((1LU<<levels)<noPos)levels++;
  
  if(roundBottom){
    //determine levels so that 2^levels <= noPos
    if((1LU<<levels)>noPos)
      levels--;
  }
  
  levels -= lessLevels;
  
  tree.reset( new FanyVector<nodeBlock> ( numberOfBlocks(levels) ));

  vector<positive> leavesV(1LU<<levels,~0);
  leaves.reset( new FanyVector<positive> (leavesV) );  
  
  box->getRange(bmin[0],bmax[0],bmin[1],bmax[1],bmin[2],bmax[2]);
  

  //for progress indication
  gescomp = noPos*double(levels); // gescomp = gesamtaufwand
  numcomp = 0;
  oprogress=0;

  theStatusRep->progress( 0, 1 );
  
  computeTreeData rootBlock;
  
  rootBlock.tree = this;
  rootBlock.aTreeCoord = rootBlock.actBlock.coordIds ;
  rootBlock.aTreePos = rootBlock.actBlock.coordVals ;
  rootBlock.restLevels = (levels + levelsPerBlock - 1)% levelsPerBlock +1;

  aLeafPos = leaves->begin();

  aBlockId = tree->begin()+1;

  aLevels = 0;

  rootBlock.computeTree(0,noPos,bmin,bmax);

  (*tree)[0] = rootBlock.actBlock;

  theStatusRep->progressEnd();
}






FkdTree::FkdTree( const FPositionSet* pos,
		  shared_ptr< FAssignableAnyArray<nodeBlock> > Atree,	   
		  shared_ptr< FAssignableAnyArray<positive> > Aleaves,
		  shared_ptr< FAssignableAnyArray<double> >   AposList,
		  shared_ptr< FAssignableAnyArray<positive> > AposListIndices,
		  positive startBlockId,
		  nodeBlock & rBlock,
		  const double * Abmin,
		  const double * Abmax,	  
		  positive nLevels
		  )
  : posSet(pos),levels(nLevels),tree(Atree),leaves(Aleaves)
{
  theStatusRep->progressBegin( "creating kdtree" );

  dim = pos->getDimension();
  box = new FBoundingBox(pos->getBoundingBox());
  

  positive noPos;

  noPos = pos->getNbPositions();  
  
  posList = AposList.get(); posListIndices = AposListIndices.get();
  

  //for progress indication
  gescomp = noPos*double(levels); // gescomp = gesamtaufwand
  numcomp = 0;
  oprogress=0;

  theStatusRep->progress( 0, 1 );

  computeTreeData rootBlock;

  rootBlock.tree = this;
  rootBlock.aTreeCoord = rootBlock.actBlock.coordIds ;
  rootBlock.aTreePos = rootBlock.actBlock.coordVals ;
  rootBlock.restLevels = levels % levelsPerBlock ;

  if(rootBlock.restLevels == 0){ 
    //if levels fit exactly into local block structure,
    // restlevels would be 0
    rootBlock.restLevels = levelsPerBlock;
  }

  aLeafPos = leaves->begin();

  aBlockId = tree->begin() + startBlockId;
  
  aLevels = 0;
  
  rootBlock.computeTree(0,noPos,Abmin,Abmax);    
  rBlock = rootBlock.actBlock;

  theStatusRep->progressEnd();
}


FkdTree::~FkdTree()
{
  delete box;
}





//-----------------------------------------------------------------------------




void FkdTree::replacePosIndByCellInd(const FNeighborhoodData*dat)
{
  vector<FIndex> neighbors;
  for(positive i=0;i < positive(1<<levels);i++){

    dat->getPositionCellNeighbors((*leaves)[i],neighbors);
    (*leaves)[i]=neighbors[0];
  }
}


//-----------------------------------------------------------------------------

void FkdTree::getRandomBBox(positive left,positive right)
{

  for(int i=0;i<dim;i++){
    minbbox[i]=maxbbox[i]=(*posList)[left*dim+i];
  }

  //if <= 1 element
  if(left+1>=right)
    {
      for(int i=0;i<dim;i++)
	maxbbox[i]=minbbox[i]+1;
      return;
    }
  
  if(right-left<1000){

    if(dim==3)
      for(positive i=left;i!=right;i++)
	{
	  positive k=i*dim;
	  double a[]={(*posList)[k],(*posList)[k+1],(*posList)[k+2]};
	  if(a[0]<minbbox[0])minbbox[0]=a[0];
	  if(a[1]<minbbox[1])minbbox[1]=a[1];
	  if(a[2]<minbbox[2])minbbox[2]=a[2];
	  if(a[0]>maxbbox[0])maxbbox[0]=a[0];
	  if(a[1]>maxbbox[1])maxbbox[1]=a[1];
	  if(a[2]>maxbbox[2])maxbbox[2]=a[2];
	}
    else
      for(positive i=left;i!=right;i++)
	{
	  positive k=i*dim;
	  double a[]={(*posList)[k],(*posList)[k+1]};
	  if(a[0]<minbbox[0])minbbox[0]=a[0];
	  if(a[1]<minbbox[1])minbbox[1]=a[1];
	  if(a[0]>maxbbox[0])maxbbox[0]=a[0];
	  if(a[1]>maxbbox[1])maxbbox[1]=a[1];
	}    
  }
  else{
    int rl=right-left;

    if(dim==3)
      for(positive i=0;i<1000;i++)
	{
	  positive k=dim*(left+positive(rl*drand48()));
	  double a[]={(*posList)[k],(*posList)[k+1],(*posList)[k+2]};
	  if(a[0]<minbbox[0])minbbox[0]=a[0];
	  if(a[1]<minbbox[1])minbbox[1]=a[1];
	  if(a[2]<minbbox[2])minbbox[2]=a[2];
	  if(a[0]>maxbbox[0])maxbbox[0]=a[0];
	  if(a[1]>maxbbox[1])maxbbox[1]=a[1];
	  if(a[2]>maxbbox[2])maxbbox[2]=a[2];
	}
    else
      for(positive i=0;i<1000;i++)
	{
	  positive k=dim*(left+positive(rl*drand48()));
	  double a[]={(*posList)[k],(*posList)[k+1]};
	  if(a[0]<minbbox[0])minbbox[0]=a[0];
	  if(a[1]<minbbox[1])minbbox[1]=a[1];
	  if(a[0]>maxbbox[0])maxbbox[0]=a[0];
	  if(a[1]>maxbbox[1])maxbbox[1]=a[1];
	}    
  }
}

//-----------------------------------------------------------------------------

double FkdTree::split(positive left,positive middle,positive right,positive coord)
{
  
  FIndexPositionRef aref
    = { posList, posListIndices,
	dim,coord,0 };
  FIndexPositionIter  posListIter = aref;

  
  if(left+1==right)//only one element
    {
      return posListIter[left];      
    }

  //then,try to split it at middle
  

  positive ll=left,rr=right-1;
  double pivot;
  positive itNo=0;

  do{
    //    cerr<<"1 ll:"<<ll<<" rr:"<<rr<<endl;
    
    positive l=ll,r=rr;
    //determine pivot element
    //every 4th time by the means between two random elements
    //else by the means between left and right element
    // itNo & 3  is equivalent to:  itNo % 4
    itNo++;
    if( itNo & 3 == 0 && l+5 < r ) 
      pivot = 
	( posListIter[l + positive(drand48()*( middle-l ))]
	  + posListIter[middle + positive(drand48()*( r+1-middle ))] 
	  )/2;
    else 
      pivot = (posListIter[l]+posListIter[r])/2;
    
// Can't compile with msvc_2005 because of error c2784
	/*
D:\Visual Studio\VC\INCLUDE\algorithm(1551) : error C2784: 'bool std::operator !
=(const std::list<_Ty,_Ax> &,const std::list<_Ty,_Ax> &)' : could not deduce tem
plate argument for 'const std::list<_Ty,_Ax> &' from 'FIndexPositionIter'
        D:\Visual Studio\VC\INCLUDE\list(1266) : see declaration of 'std::operat
or !='
        D:\Visual Studio\VC\INCLUDE\algorithm(1570) : see reference to function
template instantiation '_BidIt std::_Partition<FIndexPositionIter,_Pr>(_BidIt,_B
idIt,_Pr)' being compiled
        with
        [
            _BidIt=FIndexPositionIter,
            _Pr=std::binder2nd<std::less_equal<double>>
        ]
        .\FkdTree.cc(426) : see reference to function template instantiation '_B
idIt std::partition<FIndexPositionIter,std::binder2nd<_Fn2>>(_BidIt,_BidIt,_Pr)'
 being compiled
        with
        [
            _BidIt=FIndexPositionIter,
            _Fn2=std::less_equal<double>,
            _Pr=std::binder2nd<std::less_equal<double>>
        ]

*/

    r =  partition(posListIter+ll,posListIter+rr+1,
		   bind2nd(less_equal<double>(),pivot)) 
      - posListIter;


    if(itNo>right-left) throw FException("sth wrong in ktree:too many it");

    //if everything <=pivot,
    //split it into the parts <pivot,=pivot
    if(r==rr+1){
      
      r = partition(posListIter+ll,posListIter+rr+1,
		    bind2nd(less<double>(),pivot)) 
		    - posListIter;      

      //if middle is in equal array ,break
      if((r<=middle) & (middle<=rr))
	break;

    }


    if (r<=middle) ll = r;
    else rr = r-1;

    //cerr<<" ll:"<<ll<<" rr:"<<rr<<endl;
    
  }while(ll!=middle);

  if(! ( posListIter[left] <= pivot && posListIter[middle] >= pivot ) ){

    cout<<left<<' '<<middle<<' '<<posListIter[left]<<' '<<pivot<<' '<<posListIter[middle]<<endl;
    throw FException();
  }
  
  return pivot;
}

//-----------------------------------------------------------------------------


void FkdTree::computeTreeData::
computeTree(positive left,positive right,const double*bmin,const double*bmax)
{

  //compute bbox of contained positions
  
    
  tree->getRandomBBox(left,right);
  
  //determine coordinate in which to split
  
  static double maxd[3];
  for(int i=0;i<tree->dim;i++)
    {
      maxd[i]=(bmax[i]-bmin[i])*(tree->maxbbox[i]-tree->minbbox[i]);
    }  
  
  //  maxd[coords[i]] is sorted list, greatest first
  
  int mbb=(maxd[0]<maxd[1]);
  if(tree->dim==3)
    mbb+= 2*(maxd[1]<maxd[2]) + 4*(maxd[0]<maxd[2]);
  
  static const int coordssort[][3]={{0,1,2},{1,0,2}, //* 0 0
				    {0,2,1},   {-1}, //* 1 0 
				    {-1},   {1,2,0}, //* 0 1
				    {2,0,1},{2,1,0}};//* 1 1    


  const int *coords= coordssort[mbb];
  
  assert(coords[0]!=-1);

  assert(maxd[coords[0]]>=maxd[coords[1]] 
	 &&  (tree->dim==2 || maxd[coords[1]]>=maxd[coords[2]] ) );

    
    
  positive middle=(left+right)/2;
  
  double pivot;
  int coord=coords[0];

  //more than 1 element
  pivot = tree->split(left,middle,right, coord );    
  
  //set tree elements

  *(aTreePos++)=pivot;
  *(aTreeCoord++)=coord;

  //loop until non-null split is reached
//   for(int i=0;i<dim;i++){

//     coord = coords[i];
        

 //    double v = (pivot-bmin[coord])/(bmax[coord]-pivot);
    
//     if(finite(v)&&finite(1/v))
//    break;    
//  }
  
  //indicate progress
  
  tree->numcomp += right-left;
  if(tree->levels - tree->aLevels>10) 
      theStatusRep->progress( (unsigned)tree->numcomp, (unsigned)tree->gescomp );

  //make recursive invocations

  //increase distance from top of tree
  ++ tree->aLevels;

  //decrease distance from bottom of block
  -- restLevels;


  if(tree->aLevels==tree->levels)
    {
      //      cerr<<"inserted"<<endl;
      *(tree->aLeafPos++)=(*tree->posListIndices)[left];
      *(tree->aLeafPos++)=(*tree->posListIndices)[middle];
    }
  else{
    double a[3];

    //if there are still levels in actual block:
    if(restLevels){
      
      
      memcpy(a,bmax,sizeof(a));
      a[coord]=pivot;
      computeTree(left,middle,bmin,a);

      memcpy(a,bmin,sizeof(a));
      a[coord]=pivot;
      computeTree(middle,right,a,bmax);    

    }
    else{

      computeTreeData dat;

      dat.tree = tree;

      //...................compute left subtree
      memcpy(a,bmax,sizeof(a));
      a[coord]=pivot;
	
      dat.aTreePos = dat.actBlock.coordVals ;
      dat.aTreeCoord = dat.actBlock.coordIds ;
      dat.restLevels = FkdTree::levelsPerBlock;


      FAssignableAnyArray<nodeBlock>::iterator 
	oBlockId = tree->aBlockId++;

      dat.computeTree(left,middle,bmin,a);

      *oBlockId = dat.actBlock;

      //....................compute right subtree
      memcpy(a,bmin,sizeof(a));
      a[coord]=pivot;

      dat.aTreePos = dat.actBlock.coordVals ;
      dat.aTreeCoord = dat.actBlock.coordIds ;

      oBlockId = tree->aBlockId++;

      dat.computeTree(middle,right,a,bmax);

      *oBlockId = dat.actBlock;

    }
  }

  ++ restLevels;
  -- tree->aLevels;

}

//-----------------------------------------------------------------------------

void FkdTree::getEntries(const FBoundingBox&b,vector<FIndex>&indices) const
{

  try{

    if(!b.intersects(*box))
      return;

    indexList = &indices;
    b.getRange(minbbox[0],maxbbox[0],
	       minbbox[1],maxbbox[1],
	       minbbox[2],maxbbox[2]);

    indices.clear();

    kdTreeIter it(this);
    getBBoxIndices(it);
    
  }
  CATCH_N_RETHROW(FException);
}


//-----------------------------------------------------------------------------

void FkdTree::getBBoxIndices(kdTreeIter&it) const
{
  double s  = it.coordVal();
  positive c = it.coordId();

  bool goleft  = minbbox[ c ] <= s;
  bool goright = maxbbox[ c ] >  s;
 
  if(it.height()<=1){
    positive i = it.leafId();

    if(goleft)
      indexList->push_back((*leaves)[i]);

    if(goright)
      indexList->push_back((*leaves)[i+1]);
  }
  else{
    if(goleft){
      it.left();
      getBBoxIndices(it);
      it.up();
    }
    if(goright){
      it.right();
      getBBoxIndices(it);
      it.up();    
    }

  }


}

//-----------------------------------------------------------------------------





void FkdTree::getLeafBox( positive index,double minbox[3],double maxbox[3],positive lesslevels) const
{
  box->getRange(minbox[0],maxbox[0],
		minbox[1],maxbox[1],
		minbox[2],maxbox[2]);

  kdTreeIter it(this);

  unsigned int powlevel=levels;
     
  while(powlevel>lesslevels){
    powlevel--;

    if( (index>>powlevel) & 1 ){
      minbox[it.coordId()]=it.coordVal();
      it.right();
    }
    else{
      maxbox[it.coordId()]=it.coordVal();
      it.left();
    }
  }
             
}



//-----------------------------------------------------------------------------

FIndex FkdTree::getEntry(const FPosition&p) const
{

  double*pos=new double[dim];
  for(int i=0;i!=dim;i++)pos[i]=p[i];

  if(!box->isInside(p))
    return FIndex::invalid;

  return (*leaves)[getLeafNumber(pos)];
}
//-----------------------------------------------------------------------------

positive FkdTree::getLeafNumber( const double*pos ) const

{
 
#ifdef OUTLINE
#define register
#endif

  register positive 
    powlevel, numtogo, 
    bigger, bigp, 
    aBlockId, aLeafPos, aLeafId, subBlockSize;
  register const char * aTreeCoord;
  register const double * aTreePos;
  register const nodeBlock * b;

  aLeafId = 0;
  aBlockId = 0;

  powlevel =  (levels + levelsPerBlock - 1) % levelsPerBlock + 1;

  numtogo  =  levels - powlevel + levelsPerBlock;

  //number of blocks in sub-tree part
  subBlockSize 
    = ( (1<<numtogo) - 1 ) 
    / ( (1<<levelsPerBlock) - 1 ); 
             
  do{
    numtogo -= levelsPerBlock;
    subBlockSize >>= levelsPerBlock;

    b = &((*tree)(aBlockId));
    aTreePos   = b -> coordVals ;
    aTreeCoord = b -> coordIds ;

    aLeafPos = 0;

    do{
      --powlevel;

      bigger = (*aTreePos < pos[*aTreeCoord]) << powlevel;
      bigp = !bigger + bigger;
      aTreePos += bigp; 
      aTreeCoord += bigp;
      aLeafPos += bigger;

    }while(powlevel);
    //aLeafPos is now the id of the leaf of the block,
    //and now: the same multiplied by 2^numtogo
    aLeafId += aLeafPos<<numtogo;
    aBlockId += 1 + aLeafPos*subBlockSize;

    powlevel = levelsPerBlock;

  }while(numtogo);            
  
  return aLeafId;

#ifdef OUTLINE
#undef register
#endif
}


//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

























//-----------------------------------------------------------------------------

FIndex FkdTree::getEntryAndBox(const FPosition&pos,
			       double*minbox,double*maxbox) const

{
#ifdef OUTLINE
#define register
#endif

  register positive 
    powlevel, numtogo, 
    bigger, bigp, 
    aBlockId, aLeafPos, aLeafId, subBlockSize;
  register const char * aTreeCoord;
  register const double * aTreePos;
  register const nodeBlock * b;

  if(!box->isInside(pos))
    return FIndex::invalid;

  box->getRange(minbox[0],maxbox[0],
		minbox[1],maxbox[1],
		minbox[2],maxbox[2]);

  aLeafId = 0;
  aBlockId = 0;

  powlevel =  (levels + levelsPerBlock - 1) % levelsPerBlock + 1;

  numtogo  =  levels - powlevel + levelsPerBlock;

  //number of blocks in sub-tree part
  subBlockSize 
    = ( (1<<numtogo) - 1 ) 
    / ( (1<<levelsPerBlock) - 1 ); 
             
  do{
    numtogo -= levelsPerBlock;
    subBlockSize >>= levelsPerBlock;

    b = &((*tree)(aBlockId));
    aTreePos   = b -> coordVals ;
    aTreeCoord = b -> coordIds ;

    aLeafPos = 0;

    do{
      --powlevel;

      bigger = (*aTreePos < pos[*aTreeCoord]) << powlevel;
      if(bigger) 
	minbox[*aTreeCoord]=*aTreePos;
      else
	maxbox[*aTreeCoord]=*aTreePos;

      bigp = !bigger + bigger;
      aTreePos += bigp; 
      aTreeCoord += bigp;
      aLeafPos += bigger;

    }while(powlevel);
    //aLeafPos is now the id of the leaf of the block,
    //and now: the same multiplied by 2^numtogo
    aLeafId += aLeafPos<<numtogo;
    aBlockId += 1 + aLeafPos*subBlockSize;

    powlevel = levelsPerBlock;

  }while(numtogo);            
  
  return (*leaves)[aLeafId];

#ifdef OUTLINE
#undef register
#endif
 
}

//-----------------------------------------------------------------------------


shared_ptr< const FAssignableAnyArray<FkdTree::nodeBlock> > FkdTree::getTree(int&levels) const
{
  levels=this->levels;
  return tree;
}
//-----------------------------------------------------------------------------
shared_ptr< const FAssignableAnyArray<positive> > 
FkdTree::getLeaves() const
{
  return leaves;
}
//-----------------------------------------------------------------------------

positive FkdTree::memSize() const
{
//   cout << "Size of FkdTree (" 
//        << ((1<<levels)-1) << " levels) : " << endl;
//   cout << "tree: " << ((1<<levels)-1) * sizeof( double ) << endl
//        << "treeCoord: " << ((1<<levels)-1) * sizeof( char ) << endl
//        << "leaves: " << ((1<<levels)-1) * sizeof( unsigned int )
//        << endl << endl;

  return 
    ((1<<levels)-1) * sizeof( double ) // tree
    +
    ((1<<levels)-1) * sizeof( char )  // treeCoord
    +
    ((1<<levels)) * sizeof( unsigned int ) // leaves
    + 
    sizeof(*this);
}


