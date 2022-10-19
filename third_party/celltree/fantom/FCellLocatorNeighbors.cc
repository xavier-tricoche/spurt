 
#include "FCellLocatorNeighbors.hh"
#include "FTetrahedronCell.hh"
#include "FGrid.hh"
#include "FNeighborhoodData.hh"
#include "FkdTree.hh"
#include "Fd3op.hh"

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <map> 

//#define startAtCellCenter
//#define testIsInsideNeighbors
//#define dontUseLastCell


FCellLocatorNeighbors::FCellLocatorNeighbors( const FGrid* agrid )
    : FCellLocator( agrid->getPositionSet().get() ), grid(agrid)
{
    neighborData = grid->getCellDefinitions()->getNeighborhoodData();

    assert( neighborData );
    assert( pSet );

    tree = pSet->getKdTree(neighborData);

    //lastGonePositions=0;lastGoneCells=0;
}

FCellLocatorNeighbors::~FCellLocatorNeighbors()
{
}


bool FCellLocatorNeighbors::searchCell( FIndex& aIndex, const FPosition& aPosition ) const
{
  return searchCell( aIndex, aPosition, &mydata );
#if 0
  if(!pSet->getBoundingBox().isInside(aPosition))
    return false;


  // check distance searched position - old cell center
  // if it is smaller than max. distance of cell center 
  // to its vertices, squared,
  // use last cell as start cell

  FArray actPos;
#ifndef dontUseLastCell
  if(lastCell.isValid() &&
     aPosition.distanceSquare( cellCenter )
     < radiusBSphereSq )
    {      
      aIndex=lastCell;
      actPos=cellCenter;

      if(grid->getNextCell(cell,aIndex,actPos,aPosition,
			   lastGoneCells,lastGonePositions))
	{
	
	  if(aIndex!=lastCell){
	    lastCell=aIndex;
	    computeCellBall();
	  }
	  return true;
	}
      
    }
#endif
  vector<FIndex> indices;

  double minbox[3],maxbox[3];

  FIndex startPos
    = tree->getEntryAndBox(aPosition,minbox,maxbox);
    
  neighborData->getPositionCellNeighbors(startPos,indices);
    
  aIndex = indices[0];
  cell = grid->getCell(aIndex);
    
  pSet->getPosition(actPos,startPos);
#ifdef startAtCellCenter
  actPos=cell->getBarycenter();
#endif
  
  if( grid->getNextCell(cell,aIndex,actPos,aPosition,
			lastGoneCells,lastGonePositions) ){
    lastCell=aIndex;

#ifdef testIsInsideNeighbors
    aIndex=isInsideNeighbors(aIndex,aPosition,2,true);	  

    if(aIndex.isValid()&&lastCell!=aIndex){
      cell=grid->getCell(aIndex);
      lastCell=aIndex;
    }
#endif

    computeCellBall();
    return true;
  }
    
  
  //ifcell has not been found by beginning at 
  //current start position,try kdtree leaf neighbors
    
    
  double delt[3];
  double start[3];
    
    
  double d,eps;
    
  int dim=pSet->getDimension();
  static const double extension=0.001;
    
  if(dim==2){
    start[0]= actPos[0];
    start[1]= actPos[1];
    delt[0]=aPosition[0]-start[0];
    delt[1]=aPosition[1]-start[1];      

    eps=( maxbox[0]-minbox[0] +
	  maxbox[1]-minbox[1] )*extension;
  }
  else{
    //Fd3op2(start,=actPos,);
    //Fd3op3(delt,=aPosition,-start,);
	for(int i=0; i-3;i++){
		start[i]=actPos[i];
		delt[i] =aPosition[i] -start[i];
	}
		
    eps=( maxbox[0]-minbox[0] +
	  maxbox[1]-minbox[1] +
	  maxbox[2]-minbox[2] )*extension;
  }
    
    
  //determine coordinate where ray hits kdtree leaf border
  double dmin=HUGE_VAL;
  int imin=-1;
  for(int i=0;i<dim;i++){

    if(delt[i]>0)
      d=fabs((maxbox[i]-start[i])/delt[i]);
    else
      d=fabs((minbox[i]-start[i])/delt[i]);

    if(d<dmin)
      {dmin=d;imin=i;}	

    minbox[i]-=eps;maxbox[i]+=eps;
  }

  // set box to boarder
  if(delt[imin]>0)
    minbox[imin] = maxbox[imin]-2*eps;
  else
    maxbox[imin] = minbox[imin]+2*eps;


  vector<FIndex> leafIndices;
    
  if(dim==2)
    tree->getEntries
      (FBoundingBox(minbox[0],minbox[1],maxbox[0],maxbox[1]),
       leafIndices);
  else
    tree->getEntries
      (FBoundingBox(minbox[0],minbox[1],minbox[2],
		    maxbox[0],maxbox[1],maxbox[2]),
       leafIndices);
    
  //now try all kdtree leaves in the border
  //if a ray from there does find the right cell

  vector<FIndex>::iterator indIt;
    
  //  cout<<"searching kdtree neighbors("<<indices.size()<<')'<<endl;

  for(indIt=leafIndices.begin();
      indIt!=leafIndices.end();
      ++indIt)
    {

      FIndex startPos=*indIt;
	
      neighborData->getPositionCellNeighbors(startPos,indices);
	
      aIndex = indices[0];
      cell = grid->getCell(aIndex);
	
      pSet->getPosition(actPos,startPos);
      
      //      cout<<"trying"<<actPos<<endl;

#ifdef startAtCellCenter
      actPos=cell->getBarycenter();
#endif
	
      if(lastGonePositions)
	lastGonePositions->push_back(FPosition(HUGE_VAL,HUGE_VAL,HUGE_VAL));

      if(lastGoneCells)
	lastGoneCells->push_back(FIndex::invalid);

      if( grid->getNextCell(cell,aIndex,actPos,aPosition
			    ,lastGoneCells,lastGonePositions) )
	{
	  lastCell=aIndex;
	  
#ifdef testIsInsideNeighbors
	  aIndex=isInsideNeighbors(aIndex,aPosition,2,true);	  

	  if(aIndex.isValid()&&lastCell!=aIndex){
	    cell=grid->getCell(aIndex);
	    lastCell=aIndex;
	  }
#endif	  
	  computeCellBall();
	  return true;
	}            
    }
  
  computeCellBall();
  lastCell=aIndex;

  aIndex.setToInvalid();

//  cout.precision(16);
//   cout<<"not found:"<<aPosition[0]<<' ';
//   cout.precision(16);
//   cout<<aPosition[1]<<endl;
  return false;
#endif
}



void FCellLocatorNeighbors::computeCellBall() const
{
  mydata.computeCellBall();
#if 0
  cellCenter=cell->getBarycenter();
  vector<FPosition> cellPositions;
  cell->getVertices(cellPositions);
  vector<FPosition>::iterator posIt;
  
  radiusBSphereSq=0;
  for(posIt=cellPositions.begin();
      posIt!=cellPositions.end();
      posIt++){
    
    double d=posIt->distanceSquare(cellCenter);
    
    if(d > radiusBSphereSq) radiusBSphereSq = d;
    
  }	  	
  radiusBSphereSq*=(1+2*cell->getPrecision());
#endif
}




#if 0
const vector<FArray>& FCellLocatorNeighbors::getLastGonePositions() const
{
  return *lastGonePositions;
}
const vector<FIndex>& FCellLocatorNeighbors::getLastGoneCells() const
{
  return *lastGoneCells;
}

void FCellLocatorNeighbors::setFillLastGone(bool b) const
{
  if(b){
    if(!lastGoneCells)
      lastGoneCells=new vector<FIndex>(); 
    else
      lastGoneCells->clear();
    if(!lastGonePositions)
      lastGonePositions=new vector<FArray>();
    else
      lastGonePositions->clear();
  }
  else{
    if(lastGoneCells)
      {delete lastGoneCells;lastGoneCells=0;}
    if(lastGonePositions)
      {delete lastGonePositions;lastGonePositions=0;}
  }
}
#endif

positive FCellLocatorNeighbors::memSize() const
{
  return 
    tree->memSize()
    +
    sizeof(*this);
}


FIndex 
FCellLocatorNeighbors::
isInsideNeighbors(FIndex start,const FPosition&searchPos,int maxlayers,bool startIsCell) const
{

  typedef set<positive,less<positive> > aset;

  aset 
    nc,oc,//new cells,old cells
    np,op;//new positions,old positions

  aset::iterator i;
  vector<FIndex> buf;


  if(startIsCell){
    oc.insert(start);
  }
  else{
    op.insert(start);    
    neighborData->getPositionCellNeighbors(start,buf);
    oc.insert(buf.begin(),buf.end());    
  }
    
  for(int its=0;its<maxlayers;its++,np.swap(op),nc.swap(oc))
    {
      
      for(i=oc.begin();i!=oc.end();++i)
      {	
	  shared_ptr<FCell> cell = grid->getCell(*i);

	  if(cell->isInside(searchPos))
	    return *i;

	  cell->getVertexIndices(buf);
	  np.insert(buf.begin(),buf.end());			       	
      }

      for(i=op.begin();i!=op.end();++i)np.erase(*i);

      for(i=np.begin();i!=np.end();++i){

	neighborData->getPositionCellNeighbors(*i,buf);
	nc.insert(buf.begin(),buf.end());              	
      }

      for(i=oc.begin();i!=oc.end();i++)nc.erase(*i);
            
    }

  cout<<"seems not to be reachable(more than "<<maxlayers<<" layers): startcell:"<<start<<" searchpos "<<searchPos<<endl;

  return FIndex::invalid;
}

shared_ptr<const FkdTree> FCellLocatorNeighbors::getKdTree() const
{
  return tree;
}

//---------------------------------------------------------------------------

void FCellLocatorNeighbors::info(ostream &out) const
{
  int levels;

  tree->getTree(levels);

  out
    <<" cell locator using kdtree for approximation of start cell"<<endl
    <<" then goes throgh neighbors to find correct cell "<<endl
    <<" kdtree has:"<<levels<<" levels"<<endl;
}

//---------------------------------------------------------------------------
// Test for multithreaded version
//---------------------------------------------------------------------------

void FCellLocatorNeighbors::AccelData::computeCellBall()
{
  cellCenter=cell->getBarycenter();
  vector<FPosition> cellPositions;
  cell->getVertices(cellPositions);
  vector<FPosition>::iterator posIt;
  
  radiusBSphereSq=0;
  for(posIt=cellPositions.begin();
      posIt!=cellPositions.end();
      posIt++){
    
    double d=posIt->distanceSquare(cellCenter);
    
    if(d > radiusBSphereSq) radiusBSphereSq = d;
    
  }	  	
  radiusBSphereSq*=(1+2*cell->getPrecision());
}


bool FCellLocatorNeighbors::searchCell( FIndex& aIndex, const FPosition& aPosition, AccelData*data ) const
{
  AccelData &dd = *data;

  if(!pSet->getBoundingBox().isInside(aPosition))
    return false;


  // check distance searched position - old cell center
  // if it is smaller than max. distance of cell center 
  // to its vertices, squared,
  // use last cell as start cell

  FArray actPos;
#ifndef dontUseLastCell
  if(dd.lastCell.isValid() &&
     aPosition.distanceSquare( dd.cellCenter )
     < dd.radiusBSphereSq )
    {      
      aIndex=dd.lastCell;
      actPos=dd.cellCenter;

      if(grid->getNextCell(dd.cell,aIndex,actPos,aPosition,
			   dd.lastGoneCells,dd.lastGonePositions))
	{
	
	  if(aIndex!=dd.lastCell){
	    dd.lastCell=aIndex;
	    dd.computeCellBall();
	  }
	  return true;
	}
      
    }
#endif
  vector<FIndex> indices;

  double minbox[3],maxbox[3];

  FIndex startPos
    = tree->getEntryAndBox(aPosition,minbox,maxbox);
    
  neighborData->getPositionCellNeighbors(startPos,indices);
    
  aIndex = indices[0];
  dd.cell = grid->getCell(aIndex);
    
  pSet->getPosition(actPos,startPos);
#ifdef startAtCellCenter
  actPos=dd.cell->getBarycenter();
#endif
  
  if( grid->getNextCell(dd.cell,aIndex,actPos,aPosition,
			dd.lastGoneCells,dd.lastGonePositions) ){
    dd.lastCell=aIndex;

#ifdef testIsInsideNeighbors
    aIndex=isInsideNeighbors(aIndex,aPosition,2,true);	  

    if(aIndex.isValid()&&dd.lastCell!=aIndex){
      dd.cell=grid->getCell(aIndex);
      dd.lastCell=aIndex;
    }
#endif

    computeCellBall();
    return true;
  }
    
  
  //ifcell has not been found by beginning at 
  //current start position,try kdtree leaf neighbors
    
    
  double delt[3];
  double start[3];
    
    
  double d,eps;
    
  int dim=pSet->getDimension();
  static const double extension=0.001;
    
  if(dim==2){
    start[0]= actPos[0];
    start[1]= actPos[1];
    delt[0]=aPosition[0]-start[0];
    delt[1]=aPosition[1]-start[1];      

    eps=( maxbox[0]-minbox[0] +
	  maxbox[1]-minbox[1] )*extension;
  }
  else{
    //Fd3op2(start,=actPos,);
    //Fd3op3(delt,=aPosition,-start,);
	for(int i=0; i-3;i++){
		start[i]=actPos[i];
		delt[i] =aPosition[i] -start[i];
	}
		
    eps=( maxbox[0]-minbox[0] +
	  maxbox[1]-minbox[1] +
	  maxbox[2]-minbox[2] )*extension;
  }
    
    
  //determine coordinate where ray hits kdtree leaf border
  double dmin=HUGE_VAL;
  int imin=-1;
  for(int i=0;i<dim;i++){

    if(delt[i]>0)
      d=fabs((maxbox[i]-start[i])/delt[i]);
    else
      d=fabs((minbox[i]-start[i])/delt[i]);

    if(d<dmin)
      {dmin=d;imin=i;}	

    minbox[i]-=eps;maxbox[i]+=eps;
  }

  // set box to boarder
  if(delt[imin]>0)
    minbox[imin] = maxbox[imin]-2*eps;
  else
    maxbox[imin] = minbox[imin]+2*eps;


  vector<FIndex> leafIndices;
    
  if(dim==2)
    tree->getEntries
      (FBoundingBox(minbox[0],minbox[1],maxbox[0],maxbox[1]),
       leafIndices);
  else
    tree->getEntries
      (FBoundingBox(minbox[0],minbox[1],minbox[2],
		    maxbox[0],maxbox[1],maxbox[2]),
       leafIndices);
    
  //now try all kdtree leaves in the border
  //if a ray from there does find the right cell

  vector<FIndex>::iterator indIt;
    
  //  cout<<"searching kdtree neighbors("<<indices.size()<<')'<<endl;

  for(indIt=leafIndices.begin();
      indIt!=leafIndices.end();
      ++indIt)
    {

      FIndex startPos=*indIt;
	
      neighborData->getPositionCellNeighbors(startPos,indices);
	
      aIndex = indices[0];
      dd.cell = grid->getCell(aIndex);
	
      pSet->getPosition(actPos,startPos);
      
      //      cout<<"trying"<<actPos<<endl;

#ifdef startAtCellCenter
      actPos=cell->getBarycenter();
#endif
	
      if(dd.lastGonePositions)
	dd.lastGonePositions->push_back(FPosition(HUGE_VAL,HUGE_VAL,HUGE_VAL));

      if(dd.lastGoneCells)
	dd.lastGoneCells->push_back(FIndex::invalid);

      if( grid->getNextCell(dd.cell,aIndex,actPos,aPosition
			    ,dd.lastGoneCells,dd.lastGonePositions) )
	{
	  dd.lastCell=aIndex;
	  
#ifdef testIsInsideNeighbors
	  aIndex=isInsideNeighbors(aIndex,aPosition,2,true);	  

	  if(aIndex.isValid()&&lastCell!=aIndex){
	    cell=grid->getCell(aIndex);
	    lastCell=aIndex;
	  }
#endif	  
	  computeCellBall();
	  return true;
	}            
    }
  
  computeCellBall();
  dd.lastCell=aIndex;

  aIndex.setToInvalid();

//  cout.precision(16);
//   cout<<"not found:"<<aPosition[0]<<' ';
//   cout.precision(16);
//   cout<<aPosition[1]<<endl;
  return false;
}

