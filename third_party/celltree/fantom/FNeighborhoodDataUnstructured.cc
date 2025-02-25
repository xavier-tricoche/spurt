#include "FNeighborhoodDataUnstructured.hh"
#include "FGrid.hh"
#include "FCellDefinitions.hh"
#include "FStatusReporter.hh"

#include <algorithm>
#include <iomanip>
#include <cassert>
#include <iostream>

#include <wall_timer.hpp>
// #include <memstat.hpp>

//--------------------------------------------------------------------------- 

std::list<FNeighborhoodDataUnstructured::ppCacheEntry> FNeighborhoodDataUnstructured::EmptyList;

/// Undocumented.
void  FNeighborhoodDataUnstructured::
getCellFaceNeighbor(const FIndex& srcCellId,
		    const FIndex& faceId,
		    FIndex& dstCellId) const
{

  FCell::CellType t;
  cellDef->getCellType(srcCellId,t);
  const FCell::geoInfo &geo = FCell::getGeometryDescription(t);
  
  {
    shared_ptr<FCell> cell;
    cell=cellDef->getCellTorso(srcCellId);
    cell->getVertexIndices(vertices);
  }

  assert(faceId < geo.noFaces);
  
  positive aface[FCell::maxNoVertsPerFace];
  
  positive vertsPerFace = geo.faceSizes[faceId];
  const positive * vertsOfFace = geo.faces[faceId];

  for(positive i=0;i!=vertsPerFace;i++)
    aface[i]=vertices[vertsOfFace[i]];

  setbuf1.resize(cellsAtP[aface[0]+1]-cellsAtP[aface[0]]);
  setbuf2.resize(setbuf1.size());
  vector<positive>::iterator end;

  positive a=aface[0],b=aface[1];
  
  end=set_intersection(cellsAtP[a],cellsAtP[a+1],
		       cellsAtP[b],cellsAtP[b+1],
		       setbuf1.begin());

  for(positive i=2;i<vertsPerFace;i++)
    {
      positive c=aface[i];
      end=set_intersection(cellsAtP[c],cellsAtP[c+1],
			   setbuf1.begin(),end,
			   setbuf2.begin());
      setbuf1.swap(setbuf2);
    }

  if(setbuf1[0]!=srcCellId)
    dstCellId=setbuf1[0];
  else if( end-setbuf1.begin() >= 2)
    dstCellId=setbuf1[1];
  else dstCellId.setToInvalid();

}

//--------------------------------------------------------------------------- 

/// Undocumented.
void FNeighborhoodDataUnstructured::
getCellEdgeNeighbors(const FIndex& srcCellId,
		     const FIndex& edgeId,
		     vector<FIndex>& dstCellId) const
{
  FCell::CellType t;
  cellDef->getCellType(srcCellId,t);
  const FCell::geoInfo &geo = FCell::getGeometryDescription(t);
  {
    shared_ptr<FCell> cell;
    cell=cellDef->getCellTorso(srcCellId);
    cell->getVertexIndices(vertices);
  }
  
  positive 
    pa=geo.edges[edgeId][0],
    pb=geo.edges[edgeId][1],
    a=vertices[pa],
    b=vertices[pb];
  
  setbuf1.resize(cellsAtP[a+1]-cellsAtP[a]);
  setbuf2.resize(setbuf1.size());
  setbuf3.resize(setbuf1.size());
  
  vector<positive>::iterator end,end1;

  //get all cells connected to vertices a AND b 
  end1 = set_intersection( cellsAtP[a],cellsAtP[a+1],
			   cellsAtP[b],cellsAtP[b+1],
			   setbuf1.begin() );
   
  setbuf2=setbuf1;
  end=setbuf2.end();
    
  for(positive i=0;i<vertices.size();i++)
  {
    positive c=vertices[i];
    end=set_intersection(setbuf2.begin(),end,
			 cellsAtP[c],cellsAtP[c+1],			   
			 setbuf3.begin());
    setbuf3.swap(setbuf2);
  }
  
  
  end=set_difference(setbuf1.begin(),end1,
		     setbuf2.begin(),end,			   
		     setbuf3.begin());
  dstCellId.resize(end-setbuf3.begin());
  copy( setbuf3.begin(),end, dstCellId.begin());
}

//--------------------------------------------------------------------------- 

void FNeighborhoodDataUnstructured::getCellVertexNeighbors( const FIndex& srcCellId,
							    const FIndex& vertexId,
							    vector<FIndex>& dstCellId ) const
{
  cellDef->getCellVerticesIndices(srcCellId,vertices);

    positive a = vertices[vertexId];

    setbuf1.clear();
    setbuf1.insert(setbuf1.begin(),cellsAtP[a],cellsAtP[a+1]);
    setbuf2.resize(setbuf1.size());

    vector<positive>::iterator end=setbuf1.end();

    for(positive i=0;i<vertices.size();i++)
    {
	if(i!=vertexId)
	{
	    positive c=vertices[i];
	    end=set_difference(setbuf1.begin(),end,
			       cellsAtP[c],cellsAtP[c+1],
			       setbuf2.begin());
	    setbuf1.swap(setbuf2);
	}
    }

    dstCellId.clear();
    dstCellId.insert( dstCellId.begin(),setbuf1.begin(),end );
}

//--------------------------------------------------------------------------- 

/**
   * Get all the cells connected to a given position.
   * \param srcPosId
   * Id of the position to check for neighbors.
   * \param neighborCells
   * Storage for neighboring cells found.
   */

void 
FNeighborhoodDataUnstructured::getPositionCellNeighbors(const FIndex& srcPosId,
							vector<FIndex> &neighborCells) const
{
    neighborCells.resize( cellsAtP[positive(srcPosId)+1] - 
			  cellsAtP[srcPosId]);

    copy( cellsAtP[srcPosId],
	  cellsAtP[(positive)srcPosId+1],
	  neighborCells.begin() );
} 

//--------------------------------------------------------------------------- 

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
FNeighborhoodDataUnstructured::getEdgeCellNeighbors(const FIndex& firstVertex, 
						    const FIndex& secondVertex, 
						    vector<FIndex> &neighborCells) const
{

  neighborCells.resize(cellsAtP[(positive)firstVertex+1]-cellsAtP[firstVertex]);

  neighborCells.erase
    (
     set_intersection(cellsAtP[firstVertex],cellsAtP[(positive)firstVertex+1],
		      cellsAtP[secondVertex],cellsAtP[(positive)secondVertex+1],
		      neighborCells.begin()),
     neighborCells.end()
     );		 
  
}

//--------------------------------------------------------------------------- 

FNeighborhoodDataUnstructured::~FNeighborhoodDataUnstructured()
{
}


//---------------------------------------------------------------------------

FNeighborhoodDataUnstructured::
FNeighborhoodDataUnstructured (const FCellDefinitions*d)
  :cellDef(d),ppMaxCacheEntries(10000)
{
  try{

    positive nPos=d->getNbPositions();
    //initializing arrays

  positive nCells=d->getNbCells();
  //  cout<<" npos"<<nPos<<endl;

  vector<positive> nCellsPerP(nPos,0);//number of cells at point i

  positive nCellVect=0;

  //init nCellsPerP

  for(positive i=0;i<nCells;i++)
    {

      cellDef->getCellVerticesIndices(i,vertices);

      nCellVect+=vertices.size();

      vector<FIndex>::iterator vertIt=vertices.begin();

      for(;vertIt!=vertices.end();vertIt++)
	nCellsPerP[ *vertIt ]++;	    	  
    }

  data.resize(nCellVect);
  cellsAtP.resize(nPos+1);

  //initializes cellsAtP vector

  cellsAtP[0] = &data[0] + nCellsPerP[0];

  for(positive i=1;i<nPos;i++)
    {
      //      cout<<nCellsPerP[i]<<' '<<flush;
      cellsAtP[i] = cellsAtP[i-1] + nCellsPerP[i]; 
    }

  cellsAtP[nPos]=cellsAtP[nPos-1];

 
  //fill data vector

  for(positive i=0;i<nCells;i++)
    {

      cellDef->getCellVerticesIndices(i,vertices);

      vector<FIndex>::iterator vertIt=vertices.begin();

      for(;vertIt!=vertices.end();vertIt++)
	*( --cellsAtP[*vertIt] )  = i ;	
    }
      


  //sort cell Vector with stl-sort

  for(positive i=0;i<nPos;i++){
    //cout<<"sorting range"<<cellsAtP[i]-data.begin()<<" - "<<cellsAtP[i+1]-data.begin()<<endl;
    sort( cellsAtP[i],cellsAtP[i+1] );
  }
  } 
  catch(FException&e)
    {
      e.addTraceMessage("FNeighborhoodDataUnstructured::"
			"FNeighborhoodDataUnstructured (const FCellDefinitions*d)");
      throw;
    }

}

//---------------------------------------------------------------------------

FNeighborhoodDataUnstructured::
FNeighborhoodDataUnstructured 
(const FCellDefinitions*d, 
 const vector<pair<FCell::CellType,unsigned int> >& cellTypes,
 const vector<FIndex> & cellIndices)
  :cellDef(d),ppMaxCacheEntries(10000)
{
  try{

      // memstat_reset();

    positive nPos=d->getNbPositions();
    //initializing arrays

    positive nCells=d->getNbCells();
    //  cout<<" npos"<<nPos<<endl;
    
    vector<positive> nCellsPerP(nPos,0);//number of cells at point i

    positive nCellVect=0;

    //init nCellsPerP

    theStatusRep->progressBegin( "initializing unstructured neighborhood info" );

    wall_timer timer;

    positive oprogress=0;
    for(positive i=0;i<nCells;i++)
      {
	if((i&0xfff)==0){
	  positive progress=40*i/nCells;
	  theStatusRep->progress( progress );
	}
	
	nCellVect+= cellTypes[i+1].second - cellTypes[i].second;

	vector<FIndex>::const_iterator 
	  vertIt = cellIndices.begin()+cellTypes[i].second,
	vertEnd = cellIndices.begin()+cellTypes[i+1].second;
	
      for(;vertIt!=vertEnd;vertIt++)
	nCellsPerP[ *vertIt ]++;	    	  
    }
    
    data.resize(nCellVect);
    cellsAtP.resize(nPos+1);
    
    //initializes cellsAtP vector

    cellsAtP[0] = &data[0] + nCellsPerP[0];
    
    for(positive i=1;i<nPos;i++)
      {
	//      cout<<nCellsPerP[i]<<' '<<flush;
	cellsAtP[i] = cellsAtP[i-1] + nCellsPerP[i]; 
      }
    
    cellsAtP[nPos]=cellsAtP[nPos-1];
    
    
    //fill data vector
    oprogress=0;
    
    for(positive i=0;i<nCells;i++)
      {
	
	if((i&0xfff)==0){
	  positive progress=50*i/nCells+40;
	  theStatusRep->progress( progress );
	}
	
	vector<FIndex>::const_iterator 
	  vertIt = cellIndices.begin()+cellTypes[i].second,
	  vertEnd = cellIndices.begin()+cellTypes[i+1].second;
	
	for(;vertIt!=vertEnd;vertIt++)
	*( --cellsAtP[*vertIt] )  = i ;	
      }
    
    

    //sort cell Vector with stl-sort
    
    for(positive i=0;i<nPos;i++){
      if((i&0xfff)==0){
	positive progress=10*i/nPos+90;
	theStatusRep->progress( progress );
	}
      
      //cout<<"sorting range"<<cellsAtP[i]-data.begin()<<" - "<<cellsAtP[i+1]-data.begin()<<endl;
      sort( cellsAtP[i],cellsAtP[i+1] );
    }
    theStatusRep->progressEnd();

    // printf( "construction: %.2fs, %.2f %.2f\n", timer.elapsed(), memstat_max()/1048576.0, memstat_current()/1048576.0 ); 
  } 
  catch(FException&e)
    {
      e.addTraceMessage("FNeighborhoodDataUnstructured::"
			"FNeighborhoodDataUnstructured (const FCellDefinitions*d)");
      throw;
    }

}

//---------------------------------------------------------------------------

positive FNeighborhoodDataUnstructured::memSize() const 
{

//   cout << "Size of FNeighborhoodDataUnstructured: " << endl;
//   cout << "data (" << data.capacity() << ") : "
//        << data.capacity() * sizeof( unsigned int ) << endl
//        << "cellsAtP (" << cellsAtP.capacity() << ") : "
//        << cellsAtP.capacity() * sizeof( void * ) << endl;
//   cout << "total: " 
//        << ( data.capacity() * sizeof( unsigned int ) +
// 	    cellsAtP.capacity() * sizeof( void * ) )
//        << endl;
  return ( data.capacity() * sizeof( unsigned int ) +
	   cellsAtP.capacity() * sizeof( void * ) );
    
}

//---------------------------------------------------------------------------


/**
 *help fcn: inserts value val into sorted array vec.
 *\return
 *true if val was not existing
 */

inline bool insertValSortedVec
( vector<positive> & vec,
  positive val )
{

  vector<positive>::iterator 
    insIt = lower_bound(vec.begin(),vec.end(),val);

  if(insIt!=vec.end() && *insIt==val)return false;

  vec.insert(insIt,val);
  return true;
    
}


//---------------------------------------------------------------------------

FNeighborhoodDataUnstructured::  
ppCacheEntry::ppCacheEntry()
{
  cells=0;
}

void
FNeighborhoodDataUnstructured::  
ppCacheEntry::
init(const positive*cbeg,const positive*cend)
{

  positive nCells = cend-cbeg;

  //if 3d grid and with cells that have 3 edges per point ,
  //the number of edge-connected points
  //is like that:
  points.reserve(nCells/2+2);
  
  cells = new positive[nCells];
  cellsEnd = cells+nCells;
  //stl copy
  copy(cbeg,cend,cells);
}

FNeighborhoodDataUnstructured::  
ppCacheEntry::
~ppCacheEntry()
{
  if(cells)delete[]cells;
}

void
FNeighborhoodDataUnstructured::  
getPositionPointNeighbors(const FIndex& srcPosId,
			  vector<FIndex> &neighborPos) const
{

  positive scid=srcPosId;

  if(ppArray.empty()){
    ppArray.resize(cellsAtP.size()-1,EmptyList.begin());
    ppNbCacheEntries=0;
  }
 
  ppCacheIter cacheIt = ppArray[scid];

  if(cacheIt==EmptyList.begin()){

    ppCache.push_front( ppCacheEntry() );
    ppNbCacheEntries++;

    cacheIt = ppCache.begin();

    cacheIt->init( cellsAtP[scid],cellsAtP[scid+1]);
    cacheIt->pointId = scid;


    ppArray[scid] = cacheIt;
  }
  else{
    //move actual pos to front of actual list
    ppCache.splice(ppCache.begin(),ppCache,cacheIt);
  }

  //if point info is not ready for this point
  if( cacheIt->cells!=0 ){    

    positive 
      * cellIt    = cacheIt->cells, 
      * cellItEnd = cacheIt->cellsEnd;

    for(;cellIt!=cellItEnd;++cellIt){

      FCell::CellType t;

      cellDef->getCellType( *cellIt, t);
      
      {
	shared_ptr<FCell> cell;
	cell=cellDef->getCellTorso(*cellIt);
	cell->getVertexIndices(vertices);
      }

      //initialize cache Entries for cell vectors if not yet done

      for(positive i=0;i<vertices.size();i++){
	positive avid=vertices[i];
	ppCacheIter acIt=ppArray[avid];		

	if(acIt==EmptyList.begin()){
	  ppCache.push_front( ppCacheEntry() );
	  ppNbCacheEntries++;

	  acIt = ppCache.begin();

	  acIt->init( cellsAtP[avid],cellsAtP[avid+1]);
	  acIt->pointId = avid;

	  ppArray[avid] = acIt;
	  
	}

	if(acIt->cells && acIt!=cacheIt){
	  //stl remove
	  acIt->cellsEnd = remove(acIt->cells,acIt->cellsEnd,*cellIt);
	  
	  if(acIt->cellsEnd==acIt->cells){
	    delete [] acIt->cells;
	    acIt->cells = 0;
	  }
	}


      }
      	
      const FCell::geoInfo & g = FCell::getGeometryDescription(t);
		
      for(positive i=0;i<g.noEdges;i++){
	
	positive 
	  aId = vertices[g.edges[i][0]],
	  bId = vertices[g.edges[i][1]];

	ppCacheIter 
	  aIt = ppArray[aId],
	  bIt = ppArray[bId];

	insertValSortedVec( aIt->points,bId );
	insertValSortedVec( bIt->points,aId );	  
		
      }
                     
    }
    
    delete [] cacheIt->cells;
    cacheIt->cells = 0;

  }


  neighborPos.resize( cacheIt->points.size() ); 
  
  copy(  cacheIt->points.begin(),
	 cacheIt->points.end(),
	 neighborPos.begin()    );


  // remove overstanding list entries at end
  
  while( ppNbCacheEntries > ppMaxCacheEntries ){

    ppNbCacheEntries -- ;

    ppArray[ppCache.back().pointId]=EmptyList.begin();

    ppCache.pop_back();
    
  }

}

