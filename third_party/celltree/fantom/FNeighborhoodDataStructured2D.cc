//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FNeighborhoodDataStructured2D.cc,v $
// Language:  C++
// Date:      $Date: 2003/09/19 19:34:15 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FNeighborhoodDataStructured2D.hh"
#include "FException.hh"
#include <vector>
#include <set>
#include <algorithm>

using namespace std;

FNeighborhoodDataStructured2D::
FNeighborhoodDataStructured2D(positive iDim, positive jDim, 
			      bool triangulated)
{
  if (!iDim*jDim) {
    THROW_DEFAULT_EXCEPTION(FInvalidDimensionException);
  }

  this->iDim = iDim;
  this->jDim = jDim;
  this->triangulated = triangulated;
}

//--------------------------------------------------------------------------- 

FNeighborhoodDataStructured2D::
~FNeighborhoodDataStructured2D()
{}

//---------------------------------------------------------------------------
 
void
FNeighborhoodDataStructured2D::
getCellFaceNeighbor(const FIndex& /*srcCellId*/,
		    const FIndex& /*faceId*/,
		    FIndex& /*dstCellId*/) const
{
THROW_EXCEPTION( FException, "ERROR: Asking for face neighborhood information in 2D");
}

//---------------------------------------------------------------------------
 
void
FNeighborhoodDataStructured2D::
getCellEdgeNeighbors(const FIndex& srcCellId,
		     const FIndex& edgeId,
		     std::vector< FIndex >& dstCellId) const
{
  try {

    if (srcCellId >= (iDim-1)*(jDim-1)*(triangulated?2:1)) {
      FIndexOutOfBoundsException e("invalid cell id");
      throw e;
    }
  
    dstCellId.clear();
    
    
    if (!triangulated) {
      
      // quad cells

      // convert cellId
      
      ldiv_t d=ldiv(srcCellId,iDim-1);

      switch (edgeId.getIndex()) {
      case 0: {
	if (d.quot > 0)
	  dstCellId.push_back((positive)srcCellId-(iDim-1));
	break;
      }
      case 1: {
	if ((positive)d.rem < iDim-2)
	  dstCellId.push_back((positive)srcCellId+1);
	break;
      }
      case 2: {
	if ((positive)d.quot < jDim-2)
	  dstCellId.push_back((positive)srcCellId+iDim-1);
	break;
      }
      case 3: {
	if (d.rem > 0)
	  dstCellId.push_back((positive)srcCellId-1);
	break;
      }
      default: {
	FIndexOutOfBoundsException e("invalid edge id");
	throw e;
      }
      }
     
    }
    else {
      
      // triangles

      // convert cellId
      ldiv_t d=ldiv(srcCellId>>1,iDim-1);

      //upper triangle
      if (srcCellId&1)
	switch (edgeId) 
	  {
	  case 0: 
	    if ((positive)d.rem < iDim-2) 
	      dstCellId.push_back((positive)srcCellId+1);
	    break;
	
	  case 1: 
	    if ((positive)d.quot < jDim-2)
	      dstCellId.push_back((positive)srcCellId + (iDim-1)*2 -1);
	    break;
	    
	  case 2: 
	    dstCellId.push_back((positive)srcCellId-1);
	    break;
	    
	  default: 
	  FIndexOutOfBoundsException e("invalid edge id");
	  throw e;	  
	  }

      //lower triangle
      else  	
	switch (edgeId) 
	  {
	  case 0: 
	    if (d.quot > 0)
	      dstCellId.push_back((positive)srcCellId- (iDim-1)*2 +1);
	    break; 	
    
	  case 1: 
	    dstCellId.push_back((positive)srcCellId+1);
	    break;

	  case 2: 
	    if(d.rem>0)
	      dstCellId.push_back((positive)srcCellId-1);
	    break;

	  default: 
	    FIndexOutOfBoundsException e("invalid edge id");
	    throw e;
	  }
    }
  }
  catch (FException& e) {
    
    e.addTraceMessage(" void FNeighborhoodDataStructured2D::getCellEdgeNeighbors(const FIndex& srcCellId, const FIndex& edgeId, vector< FIndex >& dstCellId) const");
    throw;
  }
}

//---------------------------------------------------------------------------
 
void FNeighborhoodDataStructured2D::
getCellVertexNeighbors(const FIndex& srcCellId,
		       const FIndex& vertexId,
		       std::vector<FIndex>& dstCellId) const {

  try {

    if (srcCellId >= (iDim-1)*(jDim-1)*(triangulated?2:1)) {
      FIndexOutOfBoundsException e("invalid cell id");
      throw e;
    }
  
    dstCellId.clear();

    if(!triangulated){      

      ldiv_t d =  ldiv(srcCellId,iDim-1);

      switch(positive(vertexId))
	{
	case 0: //lower left
	  if(d.rem>0&&d.quot>0)
	    dstCellId.push_back((positive)srcCellId -(iDim-1) - 1);
	  break;
	case 1: //lower right
	  if((positive)d.rem<iDim-2&&d.quot>0)
	    dstCellId.push_back((positive)srcCellId -(iDim-1) + 1);
	  break;
	case 2: //upper right
	  if((positive)d.rem<iDim-2&(positive)d.quot<jDim-2)
	    dstCellId.push_back((positive)srcCellId + iDim-1 + 1);
	  break;
	case 3://upper left
	  if(d.rem>0&&(positive)d.quot<jDim-2)
	    dstCellId.push_back((positive)srcCellId + (iDim-1) - 1);
	  break;	  

	default:
	  throw FException("invalid vertex id");
	}    
    }
    else{

      positive scid=srcCellId;
      positive stepY=(iDim-1)*2;
      bool upper = scid&1;

      ldiv_t d =  ldiv(scid>>1,iDim-1);

      //upper triangle
      if(upper)
	switch(positive(vertexId))
	  {
	  case 0: //lower right
	    if(d.quot>0){
	      positive ind0=scid-stepY;
	      dstCellId.push_back(ind0);

	      if((positive)d.rem<iDim-2){
		dstCellId.push_back(ind0+1);
		dstCellId.push_back(ind0+2);			   
	      }
	    }
	    
	    break;

	  case 1: //upper right

	    if((positive)d.rem<iDim-2)
	      dstCellId.push_back(scid+2);

	    if((positive)d.quot<jDim-2){
	      positive ind0=scid+stepY;
	      dstCellId.push_back(ind0);

	      if((positive)d.rem<iDim-2){
		dstCellId.push_back(ind0+1);
	      }
	    }

	    
	    break;

	  case 2: //upper left

	    if(d.rem>0){

	      dstCellId.push_back(scid-2);

	      if((positive)d.quot<jDim-2){
		positive ind0=scid+stepY;
		dstCellId.push_back(ind0-3);
		dstCellId.push_back(ind0-2);
	      }


	    }
	    
	    break;

	  default:
	    throw FException("invalid vertex id");
	    
	    
	  }    
      
      //lower triangle
      else
	switch(positive(vertexId))
	  {
	  case 0: //lower left
	    if(d.quot>0){
	      positive ind0=scid-stepY;

	      if(d.rem>0)
		dstCellId.push_back(ind0-1);

	      dstCellId.push_back(ind0);
	    }

	    if(d.rem>0)
	      dstCellId.push_back(scid-2);			   
	    
	    break;

	  case 1: //lower right

	    if((positive)d.rem<iDim-2){

	      if(d.quot>0){
		positive ind0=scid-stepY;
		dstCellId.push_back(ind0+2);
		dstCellId.push_back(ind0+3);
	      }

	      dstCellId.push_back(scid+2);
	    }
	    
	    break;

	  case 2: //upper left

	    if((positive)d.quot<jDim-2){
	      positive ind0=scid+stepY;


	      if(d.rem>0){
		dstCellId.push_back(ind0-2);
		dstCellId.push_back(ind0-1);
	      }

	      dstCellId.push_back(ind0);
	    }
	    
	    break;

	  default:
	    throw FException("invalid vertex id");
	    	    
	  }    
    }
  
  }
  catch (FException& e) {
    e.addTraceMessage("void FNeighborhoodDataStructured2D::");
    e.addTraceMessage("getCellVertexNeighbors(const FIndex& srcCellId,");
    e.addTraceMessage("const FIndex& edgeId,");
    e.addTraceMessage("vector<FIndex>& dstCellId) const");
    throw;
  }

}

//---------------------------------------------------------------------------
 
void FNeighborhoodDataStructured2D::
getEdgeCellNeighbors(const FIndex& firstVertex, 
		     const FIndex& secondVertex, 
		     std::vector<FIndex> &neighborCells) const {

  neighborCells.clear();
  vector< FIndex > tmpCellIds1, tmpCellIds2;
  
  // get cell neighbors of both given vertices
  getPositionCellNeighbors(firstVertex, tmpCellIds1);
  getPositionCellNeighbors(secondVertex, tmpCellIds2);

  // return intersection of both vectors... (STL)
  set< positive > tmp;
  insert_iterator< set< positive > > insert_iter(tmp, 
						 tmp.begin());
  set_intersection(tmpCellIds1.begin(), tmpCellIds1.end(),
		   tmpCellIds2.begin(), tmpCellIds2.end(), insert_iter);

  for (set< positive >::const_iterator iter = tmp.begin() ;
       iter != tmp.end() ; iter++)
    neighborCells.push_back(FIndex(*iter));
}

//---------------------------------------------------------------------------

void FNeighborhoodDataStructured2D::
getPositionCellNeighbors(const FIndex& srcPosId,
			 std::vector<FIndex> &neighborCells) const
{
  
  neighborCells.clear();
  positive sz = iDim - 1; // number of quad cells per row 

  // position coordinates
  ldiv_t d = ldiv(srcPosId,iDim);
  positive i=(positive)d.rem, j=(positive)d.quot;

  //cell index of right upper quad cell
  positive cellInd = sz*j+i;
  
  // quad cells    
  if(!triangulated){
  
    //lower
    if(j>0){
      
      int a= cellInd - sz;
      
      //left
      if(i>0)
	neighborCells.push_back(a-1);
      
      //right
      if(i<sz)
	neighborCells.push_back(a);
    }
  
    //upper
    if(j<jDim-1){
      
      //left
      if(i>0)
	neighborCells.push_back(cellInd-1);
      
      //right
      if(i<sz)
	neighborCells.push_back(cellInd);     
    }
  }


  //triangle cells
  else{

    if(j>0){
      
      int a= cellInd - sz;

      //left lower cell
      if(i>0)
	neighborCells.push_back((a-1)*2+1);
      
      //right lower cells
      if(i<sz){
	neighborCells.push_back(a*2);
	neighborCells.push_back(a*2+1);
      }

    }
  
    if(j<jDim-1){
      
      //left upper cells
      if(i>0){
	neighborCells.push_back((cellInd-1)*2);
	neighborCells.push_back((cellInd-1)*2+1);
      }

      //right upper cell
      if(i<sz)
	neighborCells.push_back(cellInd*2);     
    }    
    
  }  

}
//---------------------------------------------------------------------------
 /**
   * Get all the points connected to a given position
   * by cell edge
   * \param srcPosId
   * Id of the position to check for neighbors.
   * \param neighborCells
   * Storage for neighboring positions found.
   */
void
FNeighborhoodDataStructured2D::
getPositionPointNeighbors(const FIndex& srcPosId,
			  std::vector<FIndex> &neighborPos) const
{

  positive scid = srcPosId;
  // position coordinates
  ldiv_t d = ldiv(scid,iDim);
  positive i=(positive)d.rem, j=(positive)d.quot;

  
  neighborPos.clear();

  if(j>0){

    neighborPos.push_back( scid - iDim );        

    if(triangulated & i < iDim-1)
      neighborPos.push_back( scid - iDim + 1 );
  }

  if(i>0)
    neighborPos.push_back( scid - 1 );        

  if(i<iDim-1){
    neighborPos.push_back( scid + 1 );
  }

  if(j<jDim-1){

    if(triangulated & i>0 )
      neighborPos.push_back( scid + iDim - 1 );

    neighborPos.push_back( scid + iDim );

  }
    
}

//---------------------------------------------------------------------------

positive FNeighborhoodDataStructured2D::memSize() const
{
    return 0;
}
