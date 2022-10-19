//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FNeighborhoodDataStructured3D.cc,v $
// Language:  C++
// Date:      $Date: 2003/09/19 19:34:15 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FNeighborhoodDataStructured3D.hh"
#include "FException.hh"
#include "FTetrahedronCell.hh"
#include<stdlib.h>
#include<algorithm>

#include "FNeighborhoodDataStructuredTetTables.hh"




FNeighborhoodDataStructured3D::
FNeighborhoodDataStructured3D(positive iDim, positive jDim,  positive kDim,
			      bool tetrahedrized)
{
  if (!iDim*jDim*kDim) {

    FInvalidDimensionException e;
    e.addTraceMessage("FNeighborhoodDataStructured3D::FNeighborhoodDataStructured3D(positive iDim, positive jDim, positive kDim, bool tetrahedrized)");
    throw e;
  }

  this->iDim = iDim;
  this->jDim = jDim;
  this->kDim = kDim;
  this->tetrahedrized = tetrahedrized;
}

//--------------------------------------------------------------------------- 

FNeighborhoodDataStructured3D::
~FNeighborhoodDataStructured3D()
{}

//---------------------------------------------------------------------------

 
void
FNeighborhoodDataStructured3D::
getCellFaceNeighbor(const FIndex& srcCellId,
		    const FIndex& faceId,
		    FIndex& dstCellId) const
{
  try {
    if (srcCellId >= (iDim-1)*(jDim-1)*(kDim-1)*(tetrahedrized?6:1)) {
      FIndexOutOfBoundsException e("invalid cell id");
      throw e;
    }
    
    positive scid=srcCellId;
    dstCellId.setToInvalid();

    if(tetrahedrized){
      
      ldiv_t d =ldiv(scid,6);
      positive tetInd=d.rem;
      d = ldiv(d.quot,iDim-1);//d.rem= srcPosid % iDim, d.quot= srcPosId / iDim
      positive i=d.rem; 
      d = ldiv(d.quot,jDim-1);
      positive j=d.rem,k=d.quot;

      int *tab = 
	FNeighborhoodDataStructuredTetTables::
	tetCellFaceNeighbor[tetInd][faceId];

      i+=tab[0];j+=tab[1];k+=tab[2];

      if(i<iDim-1&&j<jDim-1&&k<kDim-1)
	dstCellId
	  =( ( k*(jDim-1) + j )*(iDim-1) + i ) * 6  +  tab[3];
	
    }      

    else {

      ldiv_t d;
      d = ldiv(scid,iDim-1);//d.rem= srcPosid % iDim, d.quot= srcPosId / iDim
      positive i=d.rem; 
      d = ldiv(d.quot,jDim-1);
      positive j=d.rem,k=d.quot;

      switch(positive(faceId))
	{
	case 0: if ( i>0      ) dstCellId=scid-1; break;
	case 1: if ( i<iDim-2 ) dstCellId=scid+1; break;
	case 2: if ( j>0      ) dstCellId=scid-(iDim-1); break;
	case 3: if ( j<jDim-2 ) dstCellId=scid+(iDim-1); break;
	case 4: if ( k>0      ) dstCellId=scid-(iDim-1)*(jDim-1); break;
	case 5: if ( k<kDim-2 ) dstCellId=scid+(iDim-1)*(jDim-1); break;
	default:
	throw FIndexOutOfBoundsException("wrong face id");
	}
    }

  }
  catch (FException& e) {
    e.addTraceMessage("void FNeighborhoodDataStructured3D::getCellFaceNeighbor(const FIndex& srcCellId, const FIndex& faceId, FIndex& dstCellId) const");
    throw e;
  }
}

//---------------------------------------------------------------------------
 
void
FNeighborhoodDataStructured3D::
getCellEdgeNeighbors(const FIndex& srcCellId,
		    const FIndex& edgeId,
		    std::vector< FIndex >& dstCellId) const
{
  try{

    if (srcCellId >= (iDim-1)*(jDim-1)*(kDim-1) *(tetrahedrized ? 6 : 1)) {
      FIndexOutOfBoundsException e("invalid cell id");
      throw e;
    }


    positive scid=positive(srcCellId);

    if(tetrahedrized){

      ldiv_t d=  ldiv(scid,6);
      positive tetInd=d.rem;

      d = ldiv(d.quot,iDim-1);
      //d.rem= srcPosid % iDim, d.quot= srcPosId / iDim
      positive i=d.rem; 
      d = ldiv(d.quot,jDim-1);
      positive j=d.rem,k=d.quot;

      typedef int int7[7];
      int7 * tabvector 
	= FNeighborhoodDataStructuredTetTables::
	tetCellEdgeNeighbors[tetInd][edgeId];

      dstCellId.clear();

      for(;(*tabvector)[3]!=-1;tabvector++)
	{
	  int*tab = *tabvector;
	  positive ii=i+tab[0],jj=j+tab[1],kk=k+tab[2];
	  
	  if(ii<iDim-1&&jj<jDim-1&&kk<kDim-1){

	    //if hexcell is inside grid, insert all tetrahedron indices inside
	    //the hexcell which are vertex neighbors
	    
	    positive baseIndex =( ( kk*(jDim-1) + jj )*(iDim-1) + ii ) * 6 ;
	    
	    for(int l=3;tab[l]!=-1;l++)
	      dstCellId.push_back( baseIndex + tab[l]  );
	    
	  }	  
      }
    }

    else{

      ldiv_t d = ldiv(scid,iDim-1);
      //d.rem= srcPosid % iDim, d.quot= srcPosId / iDim
      positive i, i1, i2;
      i = i1 = i2 = d.rem; 
      d = ldiv(d.quot,jDim-1);
      positive j, k;
      j = d.rem;
      k = d.quot;
      
      switch(positive(edgeId))
	{
	case 0:  j--;k--; break;
	case 1:  j++;k--; break;
	case 2:  j--;k++; break;
	case 3:  j++;k++; break;

	case 4:  k--;i--; break;
	case 5:  k++;i--; break;
	case 6:  k--;i++; break;
	case 7:  k++;i++; break;

	case 8:  i--;j--; break;
	case 9:  i++;j--; break;
	case 10: i--;j++; break;
	case 11: i++;j++; break;
	
	default:
	  throw FIndexOutOfBoundsException("wrong edge Index");
	}
      
      if( i<iDim-1 && j<jDim-1 && k<kDim-1 ) //if cell inside grid
	{
	  dstCellId.resize(1);
	  dstCellId[0] = ( k *(jDim-1) + j  ) * (iDim-1) + i;
	}
      else
	dstCellId.resize(0);

    }
  }
  catch (FException& e) {
    
    e.addTraceMessage(" void FNeighborhoodDataStructured3D::getCellEdgeNeighbor(const FIndex& srcCellId, const FIndex& edgeId, vector< FIndex >& dstCellId) const");
    throw e;
  }
}

//---------------------------------------------------------------------------
 
void FNeighborhoodDataStructured3D::
getCellVertexNeighbors(const FIndex& srcCellId,
		       const FIndex& vertexId,
		       std::vector<FIndex>& dstCellId) const {

  try {


    ldiv_t d;

    positive scid=positive(srcCellId);

    if(tetrahedrized){

      ldiv_t d=  ldiv(scid,6);
      positive tetInd=d.rem;

      d = ldiv(d.quot,iDim-1);
      //d.rem= srcPosid % iDim, d.quot= srcPosId / iDim
      positive i=d.rem; 
      d = ldiv(d.quot,jDim-1);
      positive j=d.rem,k=d.quot;

      typedef int int10[10];
      int10 * tabvector 
	= FNeighborhoodDataStructuredTetTables::
	tetCellVertexNeighbors[tetInd][vertexId];

      dstCellId.clear();

      for(;(*tabvector)[3]!=-1;tabvector++)
	{
	  int*tab = *tabvector;
	  positive ii=i+tab[0],jj=j+tab[1],kk=k+tab[2];
	  
	  if(ii<iDim-1&&jj<jDim-1&&kk<kDim-1){

	    //if hexcell is inside grid, insert all tetrahedron indices inside
	    //the hexcell which are vertex neighbors
	    
	    positive baseIndex =( ( kk*(jDim-1) + jj )*(iDim-1) + ii ) * 6 ;
	    
	    for(int l=3;tab[l]!=-1;l++)
	      dstCellId.push_back( baseIndex + tab[l]  );
	    
	  }	  
	}
      
    }

    else{

      positive vid=positive(vertexId);

      //vid is now  id in voxel-enumerated hex cell
        
      d = ldiv(scid,iDim-1);//d.rem= srcPosid % iDim, d.quot= srcPosId / iDim
      positive i=d.rem; 
      d = ldiv(d.quot,jDim-1);
      positive j=d.rem,k=d.quot;
          
      if(vid&1) i++; else i--;
      if(vid&2) j++; else j--;
      if(vid&4) k++; else k--;

      if(i<iDim-1&&j<jDim-1&&k<kDim-1) //if neighbor cell inside grid
	{
	  dstCellId.resize(1);
	  dstCellId[0] = ( k*(jDim-1) + j ) * (iDim-1) + i ;
	}
      else 
	dstCellId.resize(0);
    }
	
	
  }
  catch (FException& e) {
    e.addTraceMessage("void FNeighborhoodDataStructured3D::");
    e.addTraceMessage("getCellVertexNeighbors(const FIndex& srcCellId,");
    e.addTraceMessage("const FIndex& edgeId,");
    e.addTraceMessage("vector<FIndex>& dstCellId) const");
    throw e;
  }

}

//---------------------------------------------------------------------------
 
void FNeighborhoodDataStructured3D::
getEdgeCellNeighbors(const FIndex& firstVertex, 
		     const FIndex& secondVertex, 
		     std::vector<FIndex> &neighborCells) const {

  vector< FIndex > tmpCellIds1, tmpCellIds2;
  
  // get cell neighbors of both given vertices;
  // vertices should be sorted
  getPositionCellNeighbors(firstVertex, tmpCellIds1);
  getPositionCellNeighbors(secondVertex, tmpCellIds2);

  // return intersection of both vectors... (STL)
  neighborCells.resize(24);
  vector<FIndex>::iterator resEnd;
  resEnd = set_intersection(
			    tmpCellIds1.begin(), tmpCellIds1.end(),
			    tmpCellIds2.begin(), tmpCellIds2.end(), 
			    neighborCells.begin()
			    );

  neighborCells.erase(resEnd,neighborCells.end());
}

//---------------------------------------------------------------------------



void FNeighborhoodDataStructured3D::
getPositionCellNeighbors(const FIndex& srcPosId,
			 std::vector<FIndex> &neighborCells) const
{

  try{

    
    ldiv_t d = ldiv(srcPosId,iDim);//d.rem= srcPosid % iDim, d.quot= srcPosId / iDim
    positive i=d.rem; 
    d=ldiv(d.quot,jDim);
    positive j=d.rem,k=d.quot;    
    int stepJ=(iDim-1),stepK=stepJ*(jDim-1);
    positive cellInd = k*stepK + j*stepJ + i; //cell Index of right upper front cell

    FIndex v[8];//vector to store possible hexahedroncell indices in it


    if(k>0){

      int a= cellInd-stepK;

      if(j>0){

	int b= a-stepJ;

	if(i>0)
	  v[0]=(b-1);

	if(i < iDim-1)
	  v[1]=b;

      }

      if(j<jDim-1){

	if(i>0)
	  v[2]=a-1;

	if(i<iDim-1)
	  v[3]=a;

      }
    }

    if(k<kDim-1){

      if(j>0){

	int b= cellInd-stepJ;

	if(i>0)
	  v[4]=b-1;

	if(i<iDim-1)
	  v[5]=b;

      }

      if(j<jDim-1){

	if(i>0)
	  v[6]=cellInd-1;

	if(i<iDim-1)
	  v[7]=cellInd;

      }
    }
    

    neighborCells.clear();

    if(tetrahedrized){      

      neighborCells.reserve(24);

      for(int l=0;l<8;l++)
	if(v[l].isValid()){
	  int*tab 
	    = FNeighborhoodDataStructuredTetTables::
	    tetPositionCellNeighbors[l];
	  int baseIndex=(positive)v[l]*6;
	  for(;*tab!=-1;tab++)
	    neighborCells.push_back(baseIndex+*tab);
	}
	    
    }	

    else{

      neighborCells.reserve(8);

      for(int l=0;l<8;l++)
	if(v[l].isValid())neighborCells.push_back(v[l]);

    }
    
  }
  
  catch (FException& e) {
    e.addTraceMessage("void FNeighborhoodDataStructured3D::");
    e.addTraceMessage("getPositionCellNeighbors(const FIndex& srcPosId,");
    e.addTraceMessage("vector<FIndex>& neighborCells) const");
    throw e;
  }
 
}

void FNeighborhoodDataStructured3D::
getPositionPointNeighbors(const FIndex& srcPosId,
			 std::vector<FIndex> &neighborPos) const
{

  try{

    positive scid = srcPosId;
    ldiv_t d = ldiv(scid,iDim);//d.rem= srcPosid % iDim, d.quot= srcPosId / iDim
    positive i=d.rem; 
    d=ldiv(d.quot,jDim);
    positive j=d.rem,k=d.quot;    
    positive stepJ=iDim,stepK=stepJ*jDim;

    neighborPos.clear();

    positive a,b;

    if(k>0){
      a = scid - stepK;
      if( tetrahedrized ){
	if(j>0 ){
	  b = a-stepJ;
	  if( i>0 )
	    neighborPos.push_back( b - 1 );
	  neighborPos.push_back( b );        
	}
	if(i>0)
	  neighborPos.push_back( a - 1);
      }
      neighborPos.push_back( a );        
    }

    if(j>0){
      a = scid - stepJ;
      if( tetrahedrized & i>0 )
	neighborPos.push_back( a - 1 );
      neighborPos.push_back( a );              
    }

    if(i>0)
      neighborPos.push_back( scid - 1 );        

    if(i<iDim-1){
      neighborPos.push_back( scid + 1 );
    }

    if(j<jDim-1){
      a = scid + stepJ;
      neighborPos.push_back( a );
      if( tetrahedrized & i < iDim-1 )
	neighborPos.push_back( a + 1 );
    }

    if(k<kDim-1){
      a = scid + stepK;
      neighborPos.push_back( a );
      
      if( tetrahedrized){

	if(i < iDim-1)
	  neighborPos.push_back( a + 1 );
	  
	if(j < jDim-1 ){
	  b = a + stepJ;
	  neighborPos.push_back( b );
	  if( i < iDim - 1 )
	    neighborPos.push_back(  b + 1 );
	}
      }
      
    }
    
  }
  CATCH_N_RETHROW(FException);
}

positive FNeighborhoodDataStructured3D::memSize() const
{
    return sizeof(*this);
}
