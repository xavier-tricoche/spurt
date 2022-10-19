//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTetrahedronCell.cc,v $
// Language:  C++
// Date:      $Date: 2003/11/05 20:57:09 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.22 $
//
//--------------------------------------------------------------------------- 

#include "FTetrahedronCell.hh"
#include "FMatrix.hh"
#include "FRefArray.hh"
#include "FPosition.hh"

#include "FException.hh"

#include "FAMSingularPoint.hh"
#include <list>
#include <vector>
//macro library for 3d operations
#include "Fd3op.hh"

// Imagine the FTetrahedronCell as a Tetrahedron with
// a base triangle consisting of vertex indices 0,1,2
// in mathematical positive sense.
// Index 3 is the index of the top.

//--------------------------------------------------------------------------- 

// the cell's geometry description
const FCell::geoInfo FTetrahedronCell::myGeoDescription =
  {
	// dimension
	3,
	// # vertices
	4, 
	// # edges
	6, 
	// # faces
	4, 
	// edges
	{{0,1},  {1,2},  {2,0},  {0,3},  {1,3},  {2,3}},
	// face sizes
	{3, 3, 3, 3},
	// faces
	{{0,1,2},  {0,3,1},  {1,3,2},  {2,3,0}}
  };

//--------------------------------------------------------------------------- 

FTetrahedronCell::FTetrahedronCell()
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  reuseInit = false;
  normalVectorsComputed = false;
}

//--------------------------------------------------------------------------- 

FTetrahedronCell::FTetrahedronCell( const vector<FIndex>&vert )
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  for(int i=0;i<4;i++)
    vertexIndices[i] = vert[i];

  reuseInit = false;
  normalVectorsComputed = false;
}

FTetrahedronCell::FTetrahedronCell( const FIndex*vert )
  :FCell(myGeoDescription,myIndices,myPositions,myPositionData,myTensors)
{
  for(int i=0;i<4;i++)
    vertexIndices[i] = vert[i];

  reuseInit = false;
  normalVectorsComputed = false;
}

//--------------------------------------------------------------------------- 

FTetrahedronCell::~FTetrahedronCell()
{
}

//--------------------------------------------------------------------------- 

FCell* FTetrahedronCell::getClone() const 
{
  try 
    {
      return new FTetrahedronCell( vertexIndices );
    }
  catch (FException& e) 
    {
      e.addTraceMessage("FCell* FTetrahedronCell::getClone() const");
      throw e;
    }
}

//--------------------------------------------------------------------------- 

positive FTetrahedronCell::sizeOfCellType() const
{
  return 4;
}

//--------------------------------------------------------------------------- 

void FTetrahedronCell::interpolate(FTensor& result,
				   const FPosition& position) const
{
  try {
    
    if(!tensors||!positions)
      throw FException("tensors or positions not set");

    double b[4];
    barycentricCoordinates(b, position);
    result = 
      b[0]*tensors[0] + 
      b[1]*tensors[1] + 
      b[2]*tensors[2] + 
      b[3]*tensors[3]; 
  }
  catch (FException& e) {
    e.addTraceMessage("void FTetrahedronCell::interpolate"
		      "(FTensor& result, const FPosition& position)");
    throw e;
  }
}
 
//--------------------------------------------------------------------------- 

void FTetrahedronCell::derivatives(FTensor& result, const FPosition&) const
{
  try{

    if(!interpolationOK||!normalVectorsComputed)      
      buildDeriv();    
    
    result = deriv;
    
  }
  catch (FException& e) {
    e.addTraceMessage("void FTetrahedronCell::derivatives"
		      "(FTensor& result,const FPosition&) const");
    throw e;
  }
}

//--------------------------------------------------------------------------- 
 
bool FTetrahedronCell::isInside(const FPosition& position) const
{
    ++FCell::ntests;

  try{
    if(!positions)
      throw FException("positions not set");
  
#ifndef NODEBUG
    if (position.getDimension() != 3)
      throw FInvalidDimensionException();
#endif

    if (!geometryOK)
      {
	buildBoundingBox();
	initBaryCoords();
	geometryOK=true;
      }

    double b[3];

    double B[3];
    
    //B=position-positions[3]
    //Fd3op3(B,=position,-positions[3],); 
	for(int i=0; i<3; i++)
		B[i] = position[i] -positions[3][i];
		
     //    cout<< normalVectors<<endl;
    //    cout<<"row0:"<<normalVectors.row(0)<<" row1:"<<normalVectors.row(1)
    //	<<" row2"<<normalVectors.row(2)<<endl;
    b[0] = Fd3prod(normalVectors[0],B);

    if (b[0] > -epsilon) 
    {
        b[1] = Fd3prod(normalVectors[1],B);

        if (b[1] > -epsilon) 
        {
            b[2] = Fd3prod(normalVectors[2],B);

            if( b[2] > -epsilon && 1 - b[0] - b[1] - b[2] > -epsilon ) 
            {

                // save results for reuse
                lastBs[0] = b[0];
                lastBs[1] = b[1];
                lastBs[2] = b[2];
                lastBs[3] = 1. - lastBs[0] - lastBs[1] - lastBs[2];
                lastPos = position;
                reuseInit = true;
                
                return true;
            }
        }
    }

    return false;
  }
  catch(FException e){
    e.addTraceMessage("bool FTetrahedronCell::isInside"
		      "(const FPosition& position) const");
    throw e;
  }
}

//--------------------------------------------------------------------------- 

bool FTetrahedronCell::isInside( const FPosition&  position,
				 const vector< FPosition >& vertices)
{
  // compute barycentric coordinates (no preliminary bounding box test)
  
  //removed the FMatrix and Fvector
  //to avoid function calls for access operators ([])
  //(makes approx. factor 2 in speed:
  //with 200000 calls 1.6 to 0.64 seconds
  //in the test with cellperf3D)

  //and I also introduced the criterium,
  // that not only the 4th baryz. coordinate >0,
  //but also the sum of the yet calculated barycentric coords smaller than 
  //the determinant (max)

  // The expressions of the inverse matrix and determinant are taken
  // from the optimized code of Max (cf. FMatrix::invert( void )).

  // matrix to be inverted
  double a11, a12, a13, a21, a22, a23, a31, a32, a33;
  a11 = vertices[0][0] - vertices[3][0];
  a21 = vertices[1][0] - vertices[3][0];
  a31 = vertices[2][0] - vertices[3][0];
  
  a12 = vertices[0][1] - vertices[3][1];
  a22 = vertices[1][1] - vertices[3][1];
  a32 = vertices[2][1] - vertices[3][1];
  
  a13 = vertices[0][2] - vertices[3][2];
  a23 = vertices[1][2] - vertices[3][2];
  a33 = vertices[2][2] - vertices[3][2];
  

  double 
    B0 = position[0]-vertices[3][0],
    B1 = position[1]-vertices[3][1], 
    B2 = position[2]-vertices[3][2];

  
  // 1st column of inverse matrix
  double 
    A00 = (a22 * a33 - a23 * a32),
    A01 = - (a21 * a33 - a23 * a31),
    A02 = (a21 * a32 - a22 * a31);

  double denom = A00*a11 + A01*a12 + A02*a13;

  double alpha = - epsilon * denom; 

  double bary = A00*B0+A01*B1+A02*B2;
  double restof1 = denom-bary;

  if ( denom > 0 
       ? (bary > alpha) & (restof1 > alpha)
       : (bary < alpha) & (restof1 < alpha)  ) {

    bary = 
      - (a12 * a33 - a13 * a32) * B0
      + (a11 * a33 - a13 * a31) * B1
      - (a11 * a32 - a12 * a31) * B2;
    //bary = b1 = A10*B0+A11*B1+A12*B2

    restof1-=bary;
    //restof1 is now denom-b0-b1
    
    if ( denom > 0 
	 ? (bary > alpha) & (restof1 > alpha) 
	 : (bary < alpha) & (restof1 < alpha)  ) {

      bary =
	+ (a12 * a23 - a13 * a22) * B0
	- (a11 * a23 - a13 * a21) * B1
	+ (a11 * a22 - a12 * a21) * B2;
      //bary = b2 = A20*B0+A21*B1+A22*B2
      restof1-=bary;
      //restof1 is now denom-b0-b1-b2 = b3

      if ( denom > 0 
	   ? (bary > alpha) & (restof1 > alpha) 
	   : (bary < alpha) & (restof1 < alpha)  ) {
	// 	cout << "bary = (" << b[0]/denom << ", " << b[1]/denom << ", "
	// 	     << b[2]/denom << ", " << 1. - (b[0] - b[1] - b[2])/denom
	// 	     << ")" << endl;
	
	return true;
      }
    }
  }

  return false;
}  

//--------------------------------------------------------------------------- 
  
void FTetrahedronCell::initBaryCoords() const
{  

  try {

    if( !positions[0].size() )
      throw FException("positions were not set");
    
    //  cout<<positions[0]<<" "<<positions[1] <<" "<<positions[2] <<" "<<positions[3]<<endl;
    

    double A[3][3];

    double *p0=positionData;
    double *p1=p0+3;
    double *p2=p1+3;
    double *p3=p2+3;

    //Fd3op3(A[0],=p0,-p3,);
    //Fd3op3(A[1],=p1,-p3,);
    //Fd3op3(A[2],=p2,-p3,);
	for(int i=0; i<3; i++){
		A[0][i] = p0[i] -p3[i];
		A[1][i] = p1[i] -p3[i];
		A[2][i] = p2[i] -p3[i];
	}

    //set normalvectors to inverted,transposed A
    Fd3kreuz(normalVectors[0],A[1],A[2]);
    Fd3kreuz(normalVectors[1],A[2],A[0]);
    Fd3kreuz(normalVectors[2],A[0],A[1]);

    denom = Fd3prod(normalVectors[0],A[0]);
    double invdenom=1/denom;

    //Fd3op1(normalVectors[0],*=invdenom);
    //Fd3op1(normalVectors[1],*=invdenom);
    //Fd3op1(normalVectors[2],*=invdenom);
	for(int i=0; i<3; i++){
		normalVectors[0][i] *=invdenom;
		normalVectors[1][i] *=invdenom;
		normalVectors[2][i] *=invdenom;
	}	
		

    normalVectorsComputed=true;
    reuseInit=false;

  }
  catch (FException& e) {
    e.addTraceMessage("void FTetrahedronCell::initBaryCoords()");
    throw e;
   }
}

//--------------------------------------------------------------------------- 

void FTetrahedronCell::barycentricCoordinates(double b[4], 
                                              const FPosition& position) const
{


  if(!geometryOK){
    buildBoundingBox();
    initBaryCoords();
    geometryOK = normalVectorsComputed = true;
  }

  if (reuseInit && position == lastPos) {
    b[0] = lastBs[0];
    b[1] = lastBs[1];
    b[2] = lastBs[2];
    b[3] = lastBs[3];

    return;
  }


  if(!normalVectorsComputed)
    initBaryCoords();

  //  cout << endl << "matrixmutl" << endl;
  double deltp[3];
  //Fd3op3(deltp,=position,-positions[3],);
  for(int i=0; i<3; i++)
	deltp[i]=position[i] -positions[3][i];

  b[0] = Fd3prod(normalVectors[0],deltp);
  b[1] = Fd3prod(normalVectors[1],deltp);
  b[2] = Fd3prod(normalVectors[2],deltp);


  b[3] = 1.-b[0]-b[1]-b[2];
}

//--------------------------------------------------------------------------- 

/**
 *\par Description:
 *\see FCell::neighborFaceForPos
 */
bool FTetrahedronCell::
neighborFaceForPos(const FArray& pos,FIndex&faceId ) const
{
  double b[4];

  barycentricCoordinates(b,pos);

  //if pos is inside, return false
  if((b[0]>0)&(b[1]>0)&(b[2]>0)&(b[3]>0)) 
    return false;

  //get number of smallest barycentric ccord
  int i=0;
  if(b[1]<b[0])i=1;
  if(b[2]<b[i])i=2;
  if(b[3]<b[i])i=3;
  
  //map barycentric coord number to face number

  static const FIndex map[4]={2,3,1,0};

  faceId=map[i];

  return true;
}
//--------------------------------------------------------------------------- 

void FTetrahedronCell::derivBaryCoord(double dbdx[4],double dbdy[4], 
				      double dbdz[4])
  const
{
  try{
    if(!positions)
      throw FException("positions not set");
    
    if(!geometryOK){
      buildBoundingBox();
      initBaryCoords();
      geometryOK = normalVectorsComputed = true;     
    }
    
    if(!normalVectorsComputed)
      initBaryCoords();
  
    for(int i=0;i<3;i++){
      dbdx[i]=normalVectors[i][0];
      dbdy[i]=normalVectors[i][1];
      dbdz[i]=normalVectors[i][2];
    }
   

    dbdx[3] = - dbdx[0] - dbdx[1] - dbdx[2];
    dbdy[3] = - dbdy[0] - dbdy[1] - dbdy[2];
    dbdz[3] = - dbdz[0] - dbdz[1] - dbdz[2];
  }

  catch(FException&e){
    e.addTraceMessage("void FTetrahedronCell::derivBaryCoord"
		      "(double dbdx[4],double dbdy[4], double dbdz[4]) const");
    throw e;
  }
}


//--------------------------------------------------------------------------- 

FCell::CellType FTetrahedronCell::getCellType(void) const
{
  return FCell::TETRAHEDRON;
}

//--------------------------------------------------------------------------- 

void FTetrahedronCell::getZeros( list<FAMSingularPoint>& result) const
{
  try {
    if( !positions[0].size() || !tensors[0].size() )
      throw FException("positions or tensors not set");

    if( !geometryOK ){
      buildBoundingBox();
      initBaryCoords();
      geometryOK = normalVectorsComputed = true;     
    }

#ifndef NODEBUG 
    if (tensors[0].getDimension() != 3) {
      FInvalidDimensionException e("ERROR: cannot locate a zero in a 3D cell if the vector field is not 3D");
      throw e;
    }

    if (tensors[0].getOrder() != 1) {
      FInvalidDimensionException e;
      throw e;
    }
#endif

    // --------------------------------------------------------------------
    //            CALCULATE THE SINGULARITY POSITION
    // --------------------------------------------------------------------

    FMatrix T(3,3);
    FVector C(3);

    for(int i=0;i!=3;i++)
      for(int j=0;j!=3;j++){
	T(j,i)=tensors[3](j)-tensors[i](j);
      }

    T.invert();

    C = T * tensors[3];

    // testing the barycentric coordinates one found
    for (positive i=0 ; i<3 ; i++) {
      if (C[i] < -epsilon)
	return;
    }
    if ( C[0] + C[1] + C[2] > 1. + epsilon)
      return;

    // a zero has been found inside the cell
    FVector singularityPos(3);
    singularityPos = ( positions[0]*C[0] + positions[1]*C[1] +
		       positions[2]*C[2] + positions[3]*(1.0-C[0]-C[1]-C[2]) );


    // setting the properties of the singularity
    FAMSingularPoint singPoint (singularityPos);

    singPoint.setOrder(tensors[0].getOrder());
    if(!interpolationOK||!normalVectorsComputed)
      buildDeriv();
    singPoint.setLinearNature(deriv);

    result.push_back( singPoint );
    
    // --------------------------------------------------------

  }

  catch( FException& e) {
    e.addTraceMessage("void FTetrahedronCell::getZeros(list<FAMSingularPoint>& result) const");
    throw e;
  }

}

//--------------------------------------------------------------------------- 

bool FTetrahedronCell::getZerosArbitrary( list<FAMSingularPoint>& result) const
{
  try {
    if( !positions[0].size() || !tensors[0].size() )
      throw FException("positions or tensors not set");

    if( !geometryOK ){
      buildBoundingBox();
      initBaryCoords();
      geometryOK = normalVectorsComputed = true;     
    }

#ifndef NODEBUG 
    if (tensors[0].getDimension() != 3) {
      FInvalidDimensionException e("ERROR: cannot locate a zero in a 3D cell if the vector field is not 3D");
      throw e;
    }

    if (tensors[0].getOrder() != 1) {
      FInvalidDimensionException e;
      throw e;
    }
#endif

    // --------------------------------------------------------------------
    //            CALCULATE THE SINGULARITY POSITION
    // --------------------------------------------------------------------

    FMatrix T(3,3);
    FVector C(3);

    for(int i=0;i!=3;i++)
      for(int j=0;j!=3;j++){
	T(j,i)=tensors[3](j)-tensors[i](j);
      }
    
    try
    {
      T.invert();
    }
    catch(FException)
    {
      return false;
    }

    C = T * tensors[3];

    // testing the barycentric coordinates one found
  //   for (positive i=0 ; i<3 ; i++) {
//       if (C[i] < -epsilon)
// 	return false;
//     }
    //  if ( C[0] + C[1] + C[2] > 1. + epsilon)
    //  return;

    // a zero has been found inside the cell
    FVector singularityPos(3);
    singularityPos = ( positions[0]*C[0] + positions[1]*C[1] +
		       positions[2]*C[2] + positions[3]*(1.0-C[0]-C[1]-C[2]) );


    // setting the properties of the singularity
    FAMSingularPoint singPoint (singularityPos);

    singPoint.setOrder(tensors[0].getOrder());
    if(!interpolationOK||!normalVectorsComputed)
      buildDeriv();
    // singPoint.setLinearNature(deriv);

    
      singPoint.setNature(FAMSingularPoint::CRITICAL); // will be changed if other ...
      singPoint.setOrder(1); // will be changed if other ...
      
      vector<FVector> eigenvectors(3);
      vector< complex<double> > l_eigenvalues(3);
      complex<double> tmp;
      FMatrix r_matrix=deriv;
    
      r_matrix.getEigenValues(l_eigenvalues);

      // --------------------------------------------------------------------
      //            MAKE VERY SMALL VALUES ZERO
      // --------------------------------------------------------------------

      for (char i = 0; i<3; i++) 
      {
	if (fabs(l_eigenvalues[i].imag()) < 1.0e-8 )
	  l_eigenvalues[i] = complex<double> (l_eigenvalues[i].real(), 0.0);
	
	if (fabs(l_eigenvalues[i].real()) < 1.0e-8 )
	  l_eigenvalues[i] = complex<double> (0.0, l_eigenvalues[i].imag());
      }

      // --------------------------------------------------------
      // --------------------------------------------------------
      // IMPORTANT NOTE :::::
      // Eigenvectors of saddle points are always stored so that
      // single signums are in front... - + + or + - - 
      // --------------------------------------------------------
      // --------------------------------------------------------

  
      // --------------------------------------------------------
      // NOW WE HAVE 3 real eigenvalues
      // --------------------------------------------------------
      if ( ( l_eigenvalues[0].imag() == 0.0 ) &&
	   ( l_eigenvalues[1].imag() == 0.0 ) &&
	   ( l_eigenvalues[2].imag() == 0.0 ) ) 
      {
	// sort the real values, small absolute values to front
	if ( fabs (l_eigenvalues[0].real()) > fabs (l_eigenvalues[2].real()) )
	{
	  tmp = l_eigenvalues[0]; l_eigenvalues[0]=l_eigenvalues[2]; l_eigenvalues[2]=tmp;
	}
	
	if ( fabs (l_eigenvalues[1].real()) > fabs (l_eigenvalues[2].real()) )
	{
	  tmp = l_eigenvalues[1]; l_eigenvalues[1]=l_eigenvalues[2]; l_eigenvalues[2]=tmp;
	}
	
	if ( fabs (l_eigenvalues[0].real()) > fabs (l_eigenvalues[1].real()) )
	{
	  tmp = l_eigenvalues[1]; l_eigenvalues[1]=l_eigenvalues[0]; l_eigenvalues[0]=tmp;
	}

	// --------------------------------------------------------
	// determine the eigenvectors
	// --------------------------------------------------------
    
	vector<FVector> nullspace, range;
	FVector X(3);
	FVector N(3);

	try {
	  for (char i=0; i<3; i++) {
	    FMatrix saveM=r_matrix;
	    saveM(0,0) -= l_eigenvalues[i].real();
	    saveM(1,1) -= l_eigenvalues[i].real();
	    saveM(2,2) -= l_eigenvalues[i].real();

	    //	    cout<<saveM;
	    nullspace.clear ();
        
	    saveM.solveSingularSystem(N, X, nullspace, range);

	    if(nullspace.size()<=0)
	    {
	      //return;
	      // THROW_EXCEPTION(FException,"Nullspace of tet has dimension zero!");
	    }

	    eigenvectors[i]=nullspace[0];
	    // if the nulspace is of dimension > 1 then 2 or more succeeding
	    // eigenvalues must have been the same -> take the results and skip
	    // the next ev's
	    if (nullspace.size() > 1) {
	      if (i==2) cout << "FATAL ERROR IN BARYTETTOOLS - NULSPACES.1" 
			     << endl;
	      cout << "BARYTETTOOLS - nulspace > dim 1" << endl;
	      eigenvectors[++i]=nullspace[1];
	    }
        
	    if (nullspace.size() > 2) {
	      if (i==2) cout << "FATAL ERROR IN BARYTETTOOLS - NULSPACES.2" 
			     << endl;
	      cout << "BARYTETTOOLS - nulspace > dim 2" << endl;
	      eigenvectors[++i]=nullspace[2];
	    }
	  }
	}
	catch (FException& e) {
	  cout << "EXCEPTION CAUGHT" << endl<< endl<< endl;
	  cout << e << endl;
	}

	//	cout << "Eigenvectors have been calculated:" << eigenvectors.size()
	//	     <<endl;
    
	if ( l_eigenvalues[0].real() == 0.0 )
	{
	  if ( l_eigenvalues[1].real() == 0.0 )
	  {
	    if ( l_eigenvalues[2].real() == 0.0 ) // the whole field is singular...
	    {
	      //		cout << "NO INTERESTING SINGULARITY FOUND" << endl;
	      singPoint.setNature(FAMSingularPoint::DEGENERATE);
	      singPoint.setType(FAMSingularPoint::NONE);
	    }
	    else // only one direction
	    {
	      //		cout << "NO INTERESTING SINGULARITY FOUND" << endl;
	      THROW_EXCEPTION( FNotImplementedException, "BARYTETTOOLS-singularities"); 
	    }
	  }
	  else // all happens in one plane
	  {
	    if ( (l_eigenvalues[1].real() > 0) + (l_eigenvalues[2].real() > 0) != 1 )
	    {
	      // same signum ...
	      if (l_eigenvalues[1].real() > 0) // repelling line in 3d
	      {
		//		  cout << "NO INTERESTING SINGULARITY FOUND" << endl;
		THROW_EXCEPTION( FNotImplementedException, "BARYTETTOOLS-singularities"); 
	      }
	      else // attracting line in 3d
	      {
		//		  cout << "NO INTERESTING SINGULARITY FOUND" << endl;
		THROW_EXCEPTION( FNotImplementedException, "BARYTETTOOLS-singularities"); 
	      }
	    }
	    else // different signum -> saddle line in 3d
	    {
	      //		cout << "NO INTERESTING SINGULARITY FOUND" << endl;
	      THROW_EXCEPTION( FNotImplementedException, "BARYTETTOOLS-singularities"); 
	    }
	  }
	}//end of  if ( l_eigenvalues[0].real() == 0.0 )
	else // 3 real eigenvalues
	  if (l_eigenvalues[0].real() > 0)
	    if (l_eigenvalues[1].real() > 0)
	      if (l_eigenvalues[2].real() > 0) // repelling node 3d
		{
		  singPoint.setType(FAMSingularPoint::REPELL_NODE_3D);
		}
	      else // saddle 3d
		{
		  singPoint.setType(FAMSingularPoint::SADDLE_3D);
		  FVector tmpVec(eigenvectors[2]);
		  complex < double > tmpClx = l_eigenvalues[2];
		  eigenvectors[2]=eigenvectors[0];
		  l_eigenvalues[2] = l_eigenvalues [0];
		  eigenvectors[0]=tmpVec;
		  l_eigenvalues[0] = tmpClx;
		}
	    else // saddle 3d
	      {
		singPoint.setType(FAMSingularPoint::SADDLE_3D);
		// maybe swap the single negative ev to front
		if (l_eigenvalues[2].real() > 0) 
		{
		  FVector tmpVec(eigenvectors[1]);
		  complex < double > tmpClx = l_eigenvalues[1];
		  eigenvectors[1]=eigenvectors[0];
		  l_eigenvalues[1] = l_eigenvalues [0];
		  eigenvectors[0]=tmpVec;
		  l_eigenvalues[0] = tmpClx;
		}
	      }
	  else
	    if (l_eigenvalues[1].real() < 0)
	      if (l_eigenvalues[2].real() < 0) // repelling node 3d
		{
		  singPoint.setType(FAMSingularPoint::ATTRACT_NODE_3D);
		}
	      else // saddle 3d
		{
		  singPoint.setType(FAMSingularPoint::SADDLE_3D);
		  FVector tmpVec(eigenvectors[2]);
		  complex < double > tmpClx = l_eigenvalues[2];
		  eigenvectors[2]=eigenvectors[0];
		  l_eigenvalues[2] = l_eigenvalues [0];
		  eigenvectors[0]=tmpVec;
		  l_eigenvalues[0] = tmpClx;
		}
	    else // saddle 3d
	      {
		singPoint.setType(FAMSingularPoint::SADDLE_3D);
		// maybe swap the single positive ev to front
		if (l_eigenvalues[2].real() < 0) 
		{
		  FVector tmpVec(eigenvectors[1]);
		  complex < double > tmpClx = l_eigenvalues[1];
		  eigenvectors[1]=eigenvectors[0];
		  l_eigenvalues[1] = l_eigenvalues [0];
		  eigenvectors[0]=tmpVec;
		  l_eigenvalues[0] = tmpClx;
		}
	      }
	singPoint.setEigenvectors(eigenvectors);
	singPoint.setEigenvalues(l_eigenvalues);

      }

      // --------------------------------------------------------
      // NOW WE HAVE 2 complex conjugate eigenvalues...
      // --------------------------------------------------------

      else {

	//	cout << "singularity with two complex evalues"<< endl;

	// take the real eigenvalue to top...
	if ( l_eigenvalues[0].imag() != 0.0 ) {
	  tmp = l_eigenvalues[0]; l_eigenvalues[0]=l_eigenvalues[1]; l_eigenvalues[1]=tmp; }
	if ( l_eigenvalues[0].imag() != 0.0 ) {
	  tmp = l_eigenvalues[0]; l_eigenvalues[0]=l_eigenvalues[2]; l_eigenvalues[2]=tmp; }
	if ( l_eigenvalues[0].imag() != 0.0 )
	  cout << "UIUIUI MORE THAN 2 COMPLEX EIGENVALUES IN BARYTETTOOLS"
	       << endl;

	FVector X(3);
	FVector N(3);

	// m - lamda * I
	r_matrix(0,0) -= l_eigenvalues[0].real();
	r_matrix(1,1) -= l_eigenvalues[0].real();
	r_matrix(2,2) -= l_eigenvalues[0].real();
    
	vector<FVector> nullspace, range;
	r_matrix.solveSingularSystem(N, X, nullspace, range);

	if (nullspace.size() != 1)
	  return false;

	eigenvectors[0] = nullspace[0];
    
	// in nullspace[0] we should have our eigenvector for lamda[0];

	// --------------------------------------------------------
	// test the 3 canonical base vectors for good results...
	// --------------------------------------------------------

	FVector E(3);
	FVector Y(3), Z(3);

	E[0]=1.0;
	r_matrix.mult(E, X);

	E[0]=0.0;
	E[1]=1.0;
	r_matrix.mult(E, Y);

	E[1]=0.0;
	E[2]=1.0;
	r_matrix.mult(E, Z);

	double a,b,c;

	a=(X*Y)/(X.norm()*Y.norm());
	b=(X*Z)/(X.norm()*Z.norm());
	c=(Y*Z)/(Y.norm()*Z.norm());

	// the smallest value is the biggest angle... take it...
	if (a<=b)
	  if (a<=c);
	  else 
	    if (b<=c)
	      Y=Z;
	    else
	      X=Z;

	eigenvectors[1] = X;
	eigenvectors[2] = Y;
    
	singPoint.setNature(FAMSingularPoint::CRITICAL);

	if ( l_eigenvalues[0].real() == 0.0 )
	  if ( l_eigenvalues[1].real() == 0.0 ) // center
	    {
	      singPoint.setType(FAMSingularPoint::CENTER_3D);
	    }
	  else if ( l_eigenvalues[1].real() > 0.0 ) // repelling spiral FOCUS
	    {
	      singPoint.setType(FAMSingularPoint::REPELL_FOCUS_3D);
	    }
	  else // attracting spiral FOCUS
	    {
	      singPoint.setType(FAMSingularPoint::ATTRACT_FOCUS_3D);
	    }
	else if ( l_eigenvalues[0].real() > 0.0 ) // axis repelling
	  if ( l_eigenvalues[1].real() > 0 ) // repelling focus
	    {
	      singPoint.setType(FAMSingularPoint::REPELL_FOCUS_3D);
	    }
	  else // spiral saddle repelling axis
	    {
	      singPoint.setType(FAMSingularPoint::SPIRAL_SADDLE_3D);
	    }
	else // axis attracting
	  if ( l_eigenvalues[1].real() < 0 ) // attracting focus
	    {
	      singPoint.setType(FAMSingularPoint::ATTRACT_FOCUS_3D);
	    }
	  else // spiral saddle attracting axis
	    {
	      singPoint.setType(FAMSingularPoint::SPIRAL_SADDLE_3D);
	    }

	singPoint.setEigenvectors(eigenvectors);
	singPoint.setEigenvalues(l_eigenvalues);
      }

    result.push_back( singPoint );
    
    // --------------------------------------------------------

  }

  catch( FException& e) {
    e.addTraceMessage("void FTetrahedronCell::getZeros(list<FAMSingularPoint>& result) const");
    throw e;
  }
  return true;

}

//--------------------------------------------------------------------------- 

void FTetrahedronCell::buildDeriv() const
{
  try{

#ifndef NODEBUG
    if(!tensors||!positions)
      throw FException("tensors or positions not set");
#endif    

    if(!geometryOK)
      {
	buildBoundingBox();
	initBaryCoords();
	geometryOK=true;
      }


    if(!normalVectorsComputed)
      initBaryCoords();
    
    
    deriv.resizeTensor(tensors[3].getDimension(),tensors[3].getOrder()+1);

    //temporary var. for calculation
    //r = t0*b0 + t1*b1 + t2*b2 + t3*(1-b0-b1-b2)
    //also r = t3 + (t0-t3) * nv0*pos + (t1-t3) * nv1*pos + (t2-t3)*nv2*pos;
    //also dr/dpos0 = (t0-t3) *nv0_0  + (t1-t3) * nv1_0 +(t2-t3)*nv2_0
    //(r = interpol.tensorwert, nv = normalVectors)


    FTensor 
      t0=tensors[0]-tensors[3],
      t1=tensors[1]-tensors[3],
      t2=tensors[2]-tensors[3];
    
    for(int i=0;i<3;i++){

      deriv[i] 
	= t0 * normalVectors[0][i] 
	+ t1 * normalVectors[1][i] 
 	+ t2 * normalVectors[2][i] ;
      // derivi is now a FTensor with a pointer to the same data array 
      // as deriv[i]
    }
    
    interpolationOK = true;

  }
  catch(FException e){
    e.addTraceMessage("void FTetrahedronCell::buildDeriv()");
    throw e;
    }
}


positive FTetrahedronCell::memSize() const
{
  return
    tensorData?tensors[0].size()*sizeof(double)*
    (sizeOfCellType()+tensors[0].getDimension() ):0
    +sizeof(*this);
}
