//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FHexahedronCell.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:06 $
// Author:    $Author: garth $
// Version:   $Revision: 1.15 $
//
//--------------------------------------------------------------------------- 

#include "FHexahedronCell.hh"
#include "FException.hh"
#include <vector>
#include "Fd3op.hh"


//--------------------------------------------------------------------------- 

// the cell's geometry description:

// !!! don't change This:               !!!
// !!! the algorithms inside this cell  !!!
// !!! depend on it                     !!!
// logic: 
// edges: for every koord.axis 
//        4 edges with voxel enumeration follow each other
// faces: for every koord.axis
//        2 faces normal to it and with rotation
//        direction right (as shown in hexahedron.fig)
const FCell::geoInfo FHexahedronCell::myGeoDescription =
  {
    //dimension
    3,
    // # vertices
    8, 
    // # edges
    12, 
    // # faces
    6, 
    // edges
    {{0,1},  {3,2},  {4,5},  {7,6},
     {0,3},  {4,7},  {1,2},  {5,6},
     {0,4},  {1,5},  {3,7},  {2,6}},  
    // face sizes
    {4, 4, 4, 4, 4, 4},
    // faces
    {{3,7,4,0},  {1,5,6,2},  
     {0,4,5,1},  {2,6,7,3},
     {0,1,2,3},  {7,6,5,4}}
  };

//--------------------------------------------------------------------------- 

positive FHexahedronCell::sizeOfCellType() const
{ return 8; }


FHexahedronCell::FHexahedronCell()
  : FCell(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
}
 
//--------------------------------------------------------------------------- 

FHexahedronCell::FHexahedronCell(const vector<FIndex>& vert,bool useVTKHexahedronEnumeration)
  : FCell(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  try{

    for(int i=0;i<8;i++)
      vertexIndices[i]=vert[i];

    if(!useVTKHexahedronEnumeration)
      fromToVoxelEnum(vertexIndices);

  }
  CATCH_N_RETHROW( FException );
}

FHexahedronCell::FHexahedronCell(const FIndex*vert)
  : FCell(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  try{

    for(int i=0;i<8;i++)
      vertexIndices[i]=vert[i];
    // to be consitent with old behavior of this function
    /*
     *\todo is this correct?
     */
    fromToVoxelEnum(vertexIndices);

  }
  CATCH_N_RETHROW( FException );
}


void FHexahedronCell::
splitZ(const double * const t[8],//3*8 floats 
       //(tensor values of edges)
       float minb[3],//3 ints
       float maxb[3],
       int level,
       byte signsx, //signs of 1st,2nd, 3rd tensor comp.
       byte signsy,
       byte signsz, 
       list<FAMSingularPoint>& result
       ) const
{
  try{

    level--;

    if(level==0){
      //sing.found
      computeSing(minb,maxb,result);
      return;
    }

    static const byte
      bitshi=0xF0,bitslo=0x0F;  //bits for left,right signs
  
  
    //interpolate tensor values in the middle

    double s[4][3];

    for(int i=0; i<3; i++)
    {	
      s[0][i] = t[0][i] + t[4][i];
      s[1][i] = t[1][i] + t[5][i];
      s[2][i] = t[3][i] + t[7][i];
      s[3][i] = t[2][i] + t[6][i];
    }
    for(int i=0; i<3; i++)
    {
      s[0][i] /= 2;
      s[1][i] /= 2;
      s[2][i] /= 2;
      s[3][i] /= 2;
    }
    

    //compute signs of tensor values in the middle into correct positions of bit vector

    byte sxb = ((s[3][0]>0)<<3) | ((s[2][0]>0)<<2) | ((s[1][0]>0)<<1) | (s[0][0]>0);
    byte syb = ((s[3][1]>0)<<3) | ((s[2][1]>0)<<2) | ((s[1][1]>0)<<1) | (s[0][1]>0);
    byte szb = ((s[3][2]>0)<<3) | ((s[2][2]>0)<<2) | ((s[1][2]>0)<<1) | (s[0][2]>0);

  
    byte nsxb,nsyb,nszb;
  
    //compute new bitvector of back part

    nsxb= (signsx&bitslo) | (sxb<<4);
    nsyb= (signsy&bitslo) | (syb<<4);
    nszb= (signsz&bitslo) | (szb<<4);
  
    if( (nsxb!=0) & (nsxb!=0xff) & (nsyb!=0) & (nsyb!=0xff) & (nszb!=0) & (nszb !=0xff) )
      {
	const double * const newT[8] = {t[0],t[1],t[2],t[3],s[0],s[1],s[3],s[2]};       
	float newMaxb[3] = {maxb[0],maxb[1],(minb[2]+maxb[2])/2};	  
	splitX(newT,minb,newMaxb,level,nsxb,nsyb,nszb,result);	
      }
  
    //compute new bitvector of front part

    nsxb= (signsx&bitshi)|sxb;
    nsyb= (signsy&bitshi)|syb;
    nszb= (signsz&bitshi)|szb;
  
    if( (nsxb!=0) & (nsxb!=0xff) & (nsyb!=0) & (nsyb!=0xff) & (nszb!=0) & (nszb !=0xff) )
      {
	const double * const newT[8] = {s[0],s[1],s[3],s[2],t[4],t[5],t[6],t[7]};       
	float newMinb[3] = {minb[0],minb[1],(minb[2]+maxb[2])/2};	  
	splitX(newT,newMinb,maxb,level,nsxb,nsyb,nszb,result);	
      }      
  }
  CATCH_N_RETHROW( FException );
}

void FHexahedronCell::
splitY(const double * const t[8],//3*8 floats 
       //(tensor values of edges)
       float minb[3],//3 ints
       float maxb[3],
       int level,
       byte signsx, //signs of 1st,2nd, 3rd tensor comp.
       byte signsy,
       byte signsz, 
       list<FAMSingularPoint>& result
       ) const
{
  try{

    level--;

    if(level==0){
      //sing.found
      computeSing(minb,maxb,result);
      return;
    }

    static const byte
      bitshi = 0xCC,bitslo = 0x33;  //bits for up,down signs
  
  
    //interpolate tensor values in the middle

    double s[4][3];
	
    for(int i = 0; i<3; i++)
    {
      s[0][i] = t[0][i] + t[3][i];
      s[1][i] = t[1][i] + t[2][i];
      s[2][i] = t[4][i] + t[7][i];
      s[3][i] = t[5][i] + t[6][i];
    }
    
    for(int i = 0; i<3; i++)
    {
      s[0][i] /= 2;
      s[1][i] /= 2;
      s[2][i] /= 2;
      s[3][i] /= 2;
    }

    //compute signs of tensor values in the middle into correct positions of bit vector

    byte sxb = ((s[3][0]>0)<<5) | ((s[2][0]>0)<<4) | ((s[1][0]>0)<<1) | (s[0][0]>0);
    byte syb = ((s[3][1]>0)<<5) | ((s[2][1]>0)<<4) | ((s[1][1]>0)<<1) | (s[0][1]>0);
    byte szb = ((s[3][2]>0)<<5) | ((s[2][2]>0)<<4) | ((s[1][2]>0)<<1) | (s[0][2]>0);

  
    byte nsxb,nsyb,nszb;
  
    //compute new bitvector of down part

    nsxb= (signsx&bitslo) | (sxb<<2);
    nsyb= (signsy&bitslo) | (syb<<2);
    nszb= (signsz&bitslo) | (szb<<2);
  
    if( (nsxb!=0) & (nsxb!=0xff) & (nsyb!=0) & (nsyb!=0xff) & (nszb!=0) & (nszb !=0xff) )
      {
	const double * const newT[8] = {t[0],t[1],s[1],s[0],t[4],t[5],s[3],s[2]};        
	float newMaxb[3] = {maxb[0],(minb[1]+maxb[1])/2,maxb[2]};	  
	splitZ(newT,minb,newMaxb,level,nsxb,nsyb,nszb,result);	
      }
  
    //compute new bitvector of up part

    nsxb= (signsx&bitshi)|sxb;
    nsyb= (signsy&bitshi)|syb;
    nszb= (signsz&bitshi)|szb;
  
    if( (nsxb!=0) & (nsxb!=0xff) & (nsyb!=0) & (nsyb!=0xff) & (nszb!=0) & (nszb !=0xff) )
      {
	const double * const newT[8] = {s[0],s[1],t[2],t[3],s[2],s[3],t[6],t[7]};       
	float newMinb[3] = {minb[0],(minb[1]+maxb[1])/2,minb[2]};	  
	splitZ(newT,newMinb,maxb,level,nsxb,nsyb,nszb,result);	
      }      

  }
  CATCH_N_RETHROW( FException );
}

void FHexahedronCell::
splitX(const double * const t[8],//3*8 floats 
       //(tensor values of edges)
       float minb[3],//3 ints
       float maxb[3],
       int level,
       byte signsx, //signs of 1st,2nd, 3rd tensor comp.
       byte signsy,
       byte signsz, 
       list<FAMSingularPoint>& result
       ) const
{
  try{

    level--;

    if(level==0){
      //sing.found
      computeSing(minb,maxb,result);
      return;
    }

    static const byte
      bitshi = 0xAA,bitslo = 0x55;  //bits for left,right signs 0x55=01010101
  
  
    //interpolate tensor values in the middle

    double s[4][3];

    for(int i = 0; i<3; i++)
    {
      s[0][i] = t[0][i] + t[1][i];
      s[1][i] = t[3][i] + t[2][i];
      s[2][i] = t[4][i] + t[5][i];
      s[3][i] = t[7][i] + t[6][i];
    }
    for(int i = 0; i<3; i++)
    {
      s[0][i] /= 2;
      s[1][i] /= 2;
      s[2][i] /= 2;
      s[3][i] /= 2;
    }

    //compute signs of tensor values in the middle into correct positions of bit vector

    byte sxb = ((s[3][0]>0)<<6) | ((s[2][0]>0)<<4) | ((s[1][0]>0)<<2) | (s[0][0]>0);
    byte syb = ((s[3][1]>0)<<6) | ((s[2][1]>0)<<4) | ((s[1][1]>0)<<2) | (s[0][1]>0);
    byte szb = ((s[3][2]>0)<<6) | ((s[2][2]>0)<<4) | ((s[1][2]>0)<<2) | (s[0][2]>0);

  
    byte nsxb,nsyb,nszb;
  
    //compute new bitvector of left part

    nsxb= (signsx&bitslo) | (sxb<<1);
    nsyb= (signsy&bitslo) | (syb<<1);
    nszb= (signsz&bitslo) | (szb<<1);
  
    if( (nsxb!=0) & (nsxb!=0xff) & (nsyb!=0) & (nsyb!=0xff) & (nszb!=0) & (nszb !=0xff) )
      {
	const double * const newT[8]={t[0],s[0],s[1],t[3],t[4],s[2],s[3],t[7]};      
	float newMaxb[3]={(minb[0]+maxb[0])/2,maxb[1],maxb[2]};	  
	splitY(newT,minb,newMaxb,level,nsxb,nsyb,nszb,result);	
      }
  
    //compute new bitvector of right part

    nsxb= (signsx&bitshi)|sxb;
    nsyb= (signsy&bitshi)|syb;
    nszb= (signsz&bitshi)|szb;
  
    if( (nsxb!=0) & (nsxb!=0xff) & (nsyb!=0) & (nsyb!=0xff) & (nszb!=0) & (nszb !=0xff) )
      {
	const double * const newT[8]={s[0],t[1],t[2],s[1],s[2],t[5],t[6],s[3]};       
	float newMinb[3]={(minb[0]+maxb[0])/2,minb[1],minb[2]};	  
	splitY(newT,newMinb,maxb,level,nsxb,nsyb,nszb,result);	
      }      
  }
  CATCH_N_RETHROW( FException );
}
				     


void FHexahedronCell::getZeros( list<FAMSingularPoint>& result) const
{
  try{
    
#ifndef NODEBUG
    if(!tensorData||!positions[0].size())
      THROW_EXCEPTION( FException, "tensors or positions not set" );
    
    if(tensors[0].getOrder()!=1)
      THROW_EXCEPTION( FException, "getzeros only works with order 1 tensors" );
#endif

  
    //check if zero point is possible
    //also check if cell is zero at all
  
    set_to_zero = true; 
    byte s[3]={0,0,0};
  
    for(int i=0;i<3;i++){
    
      bool plus=false;//is there a value >0 ?
      bool minus=false;//is there a value <0 ?

      for(byte j=0;j<8;j++){
	double d= tensorData[j*3+i];
	s[i] |= (d>0) <<j;
	plus  |= (d > zero_threshold);
	minus |= (d < -zero_threshold);
      }
    
      //if at least one tensor value is !=0, set set_to_zero is false
      set_to_zero &=  ! ( plus | minus ) ;
    
      //if no sign change is possible inside of cell, there is no zero inside
      //(but could be on border)
      if( plus ^ minus)
	return;    
    }


    if(set_to_zero)return;
    
    float minb[3]={0,0,0},maxb[3]={1,1,1};
    double * tens[8] = {tensorData,tensorData+3,tensorData+6,tensorData+9,
			tensorData+12,tensorData+15,tensorData+18,tensorData+21};

    splitX(tens,minb,maxb,34,s[0],s[1],s[2],result);

  }
  CATCH_N_RETHROW( FException );

}
