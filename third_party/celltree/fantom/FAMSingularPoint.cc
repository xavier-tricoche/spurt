//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMSingularPoint.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 13:16:30 $
// Author:    $Author: garth $
// Version:   $Revision: 1.35 $
//
//--------------------------------------------------------------------------- 

#include "FAMSingularPoint.hh"

#include "FMath.hh"
#include "FException.hh"

#include "FVector.hh"
#include "FMatrix.hh"

#include "FTensor.hh"
#include "FCell.hh"

#include <complex>
#include <iostream>

#include "FMathValues.hh" 

//---------------------------------------------------------------------------
#ifdef OUTLINE
#include "FAMSingularPoint.icc"
#endif
//--------------------------------------------------------------------------- 



const double TEST_DISTANCE = 1.0E-5;
const double ANGLE_TEST_TOLERANCE = 2.0E-1;

FAMSingularPoint::mainType FAMSingularPoint::whichType(FAMSingularPoint::singularityType type) 
{
  switch (type) {

  case FAMSingularPoint::SADDLE_2D:
  case  FAMSingularPoint::SADDLE_3D:
  case  FAMSingularPoint::SPIRAL_SADDLE_3D:
    return FAMSingularPoint::SADDLE;

  case FAMSingularPoint::ATTRACT_NODE:
  case FAMSingularPoint::ATTRACT_FOCUS:
  case FAMSingularPoint::ATTRACT_NODE_3D:
  case FAMSingularPoint::ATTRACT_FOCUS_3D:
    return FAMSingularPoint::SINK;

  case FAMSingularPoint::REPELL_NODE:
  case FAMSingularPoint::REPELL_FOCUS:
  case FAMSingularPoint::REPELL_NODE_3D:
  case FAMSingularPoint::REPELL_FOCUS_3D:
    return FAMSingularPoint::SOURCE;
    
  case FAMSingularPoint::CENTER_2D:
    return FAMSingularPoint::CENTER;

  case FAMSingularPoint::TRISECTOR_POINT:
    return FAMSingularPoint::TRISECTOR;

  case FAMSingularPoint::WEDGE_POINT:
    return FAMSingularPoint::WEDGE;
    
  case FAMSingularPoint::WEDGE_SINGLE_POINT:
    return FAMSingularPoint::WEDGE_SINGLE;

  default:
    return FAMSingularPoint::OTHER;
  }
}

//--------------------------------------------------------------------------- 

FAMSingularPoint::FAMSingularPoint() 
  : FObject(), pos(), nature(FAMSingularPoint::CRITICAL), type(FAMSingularPoint::NONE), order(0), eigenvector(0),
    eigenvalues(0), sepStart(), sepStop(),
    cellId(), enstrophy(-1.), max_norm(-1.)
{
}

//--------------------------------------------------------------------------- 

FAMSingularPoint::FAMSingularPoint(const FAMSingularPoint::singularityNature& nature) 
  : pos(), type(FAMSingularPoint::NONE), order(0), eigenvector(0), eigenvalues(0), sepAngles(0),
    sepStart(), sepStop(), cellId(), enstrophy(-1.), max_norm(-1.)
{
  this->nature = nature;
}

//--------------------------------------------------------------------------- 

FAMSingularPoint::FAMSingularPoint(const FAMSingularPoint& sing) 
  : FObject(),  pos(sing.pos), nature(sing.nature), type(sing.type), order(sing.order),
    eigenvector(0), sepStart(sing.sepStart), sepStop(sing.sepStop),
    cellId(sing.cellId), connections(sing.connections), 
    enstrophy(sing.enstrophy), max_norm(sing.max_norm)
{
  this->eigenvalues=sing.eigenvalues;
  
  eigenvector.resize(sing.eigenvector.size());
  for (positive i=0 ; i<sing.eigenvector.size() ; i++)
    eigenvector[i] = sing.eigenvector[i];
  sepAngles.resize(sing.sepAngles.size());
  for (positive i=0 ; i<sing.sepAngles.size() ; i++)
    sepAngles[i] = sing.sepAngles[i];
}

//--------------------------------------------------------------------------- 

FAMSingularPoint::FAMSingularPoint(const FPosition& pos)
  : pos(), nature(FAMSingularPoint::CRITICAL), type(FAMSingularPoint::NONE), order(0), eigenvector(0),
    eigenvalues(0),sepStart(), sepStop(), cellId(), connections(), 
    enstrophy(-1.), max_norm(-1.)
{
  this->pos = pos;
}

//--------------------------------------------------------------------------- 

FAMSingularPoint::FAMSingularPoint(const FPosition& pos, 
				   const FAMSingularPoint::singularityNature& nature)
  : pos(), type(FAMSingularPoint::NONE), order(0), eigenvector(0), eigenvalues(0),
    sepAngles(0), sepStart(), sepStop(), cellId(), connections(), 
    enstrophy(-1.), max_norm(-1.)
{
  this->pos = pos;
  this->nature = nature;
}

//--------------------------------------------------------------------------- 

FAMSingularPoint::FAMSingularPoint(const FPosition& pos, 
				   const FAMSingularPoint::singularityType& type,
				   positive order)
    : pos(), nature(FAMSingularPoint::CRITICAL), eigenvector(0), eigenvalues(0),
      sepStart(), sepStop(), cellId(), connections(), enstrophy(-1.),
      max_norm(-1.)
{
  this->pos = pos;
  this->type = type;
  this->order = order;
}

//--------------------------------------------------------------------------- 

FAMSingularPoint::FAMSingularPoint(const FPosition& pos,
                                   const FAMSingularPoint::singularityNature& nature,
                                   const FAMSingularPoint::singularityType& type,
                                   positive order)
    : pos(), eigenvector(0), eigenvalues(0), sepAngles(0), sepStart(),
      sepStop(), cellId(), connections(), enstrophy(-1.), max_norm(-1.)
{
  this->pos = pos;
  this->nature = nature;
  this->type = type;
  this->order = order;

}

//--------------------------------------------------------------------------- 

FAMSingularPoint::~FAMSingularPoint()
{
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::setIncludingCellIndex(const FIndex& cellId) {

  this->cellId = cellId;
}

//--------------------------------------------------------------------------- 

const FString& FAMSingularPoint::getClassName() const 
{
  static FString className = "FAMSingularPoint";

  return className;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::setLinearNature(const FMatrix& approx) {

  if (order != 1) {
    THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
  }
    
  // first order singularity
  order = FIndex(1);

  if (approx.getDimensionX() == 2) 
  {

    complex<double> eigenVal[2];
    vector< complex<double> > l(2);
    double eigenVec1[2], eigenVec2[2];
    
    FMath::Eigensystem(approx, eigenVal, eigenVec1 , eigenVec2);

    l[0] = eigenVal[0];
    l[1] = eigenVal[1];    

    setEigenvalues(l);

    if (!eigenVal[0].imag()) {
      // 2 real eigenvalues
      
      if (eigenVal[0].real() * eigenVal[1].real() == 0.) {
	cerr << "WARNING: degenerate or higher order singularity found!" 
	     << endl;
	type = FAMSingularPoint::NONE;
      }
      else if (eigenVal[0].real() * eigenVal[1].real() < 0.) {
	// saddle case
	
	type = FAMSingularPoint::SADDLE_2D;
	
	eigenvector.clear();
	eigenvector.resize(2);
	FVector vec(2);
	
	if (eigenVal[0].real() > 0.) {
	  // repelling direction
	  FVector vec(2);
	  vec[0] = eigenVec1[0];
	  vec[1] = eigenVec1[1];
	  vec.normalize();
	  eigenvector[0] = vec;
	  // attracting direction	
	  vec[0] = eigenVec2[0];
	  vec[1] = eigenVec2[1];
	  vec.normalize();
	  eigenvector[1] = vec;
	}
	else {
	  // repelling direction
	  vec[0] = eigenVec2[0];
	  vec[1] = eigenVec2[1];
	  vec.normalize();
	  eigenvector[0] = vec;
	  // attracting direction	
	  vec[0] = eigenVec1[0];
	  vec[1] = eigenVec1[1];
	  vec.normalize();
	  eigenvector[1] = vec;

	  l[1] = eigenVal[0];
	  l[0] = eigenVal[1];    
	  
	  setEigenvalues(l);
	}
      }
      else if (eigenVal[0].real() > 0.) {
	// repelling node case
	
	type = FAMSingularPoint::REPELL_NODE;
		
	eigenvector.clear();
	eigenvector.resize(2);

	if(eigenVal[0].real() > eigenVal[1].real() ){
	  FVector vec(2);
	  vec[0] = eigenVec1[0];
	  vec[1] = eigenVec1[1];
	  vec.normalize();
	  eigenvector[0] = vec;
	  
	  vec[0] = eigenVec2[0];
	  vec[1] = eigenVec2[1];
	  vec.normalize();
	  eigenvector[1] = vec;
	}
	else{
	  FVector vec(2);
	  vec[0] = eigenVec2[0];
	  vec[1] = eigenVec2[1];
	  vec.normalize();
	  eigenvector[0] = vec;

	  vec[0] = eigenVec1[0];
	  vec[1] = eigenVec1[1];
	  vec.normalize();
	  eigenvector[1] = vec;
	  
	  l[1] = eigenVal[0];
	  l[0] = eigenVal[1];    
	  
	  setEigenvalues(l);
	}

	//eigenvector.clear();

	enstrophy = .5 *fabs(approx(1,0)-approx(0,1));
      }
      else {
	// attracting node case
	
	type = FAMSingularPoint::ATTRACT_NODE;
		
	eigenvector.clear();
	eigenvector.resize(2);

	if( eigenVal[0].real() > eigenVal[1].real() ){
	  FVector vec(2);
	  vec[0] = eigenVec1[0];
	  vec[1] = eigenVec1[1];
	  vec.normalize();
	  eigenvector[0] = vec;
	  
	  vec[0] = eigenVec2[0];
	  vec[1] = eigenVec2[1];
	  vec.normalize();
	  eigenvector[1] = vec;
	}
	else{
	  FVector vec(2);
	  vec[0] = eigenVec2[0];
	  vec[1] = eigenVec2[1];
	  vec.normalize();
	  eigenvector[0] = vec;
	  
	  vec[0] = eigenVec1[0];
	  vec[1] = eigenVec1[1];
	  vec.normalize();
	  eigenvector[1] = vec;
	  
	  l[1] = eigenVal[0];
	  l[0] = eigenVal[1];    
	  
	  setEigenvalues(l);
	}
	
	// eigenvector.clear();

	enstrophy = .5 *fabs(approx(1,0)-approx(0,1));
      }
    }
    else {
      // 2 conjugate complex eigenvalues
      
      if (eigenVal[0].real() == 0.) {
	// center case
	
	type = FAMSingularPoint::CENTER_2D;
	eigenvector.clear();

	enstrophy = .5 *fabs(approx(1,0)-approx(0,1));
      }
      else if (eigenVal[0].real() > 0.) {
	// repelling focus case
	
	type = FAMSingularPoint::REPELL_FOCUS;
	eigenvector.clear();

	enstrophy = .5 *fabs(approx(1,0)-approx(0,1));
      }
      else {
	// attracting focus case
	
	type = FAMSingularPoint::ATTRACT_FOCUS;
	eigenvector.clear();

	enstrophy = .5 *fabs(approx(1,0)-approx(0,1));
      }
    }
  }

  //calculation for 3d taken from FTetrahedronCell
  else 
    if(approx.getDimensionX()==3)      
    {
      setNature(FAMSingularPoint::CRITICAL); // will be changed if other ...
      setOrder(1); // will be changed if other ...

      vector<FVector> eigenvectors(3);
      vector< complex<double> > l_eigenvalues(3);
      complex<double> tmp;
      FMatrix r_matrix=approx;
    
      r_matrix.getEigenValues(l_eigenvalues);

      // --------------------------------------------------------------------
      //            MAKE VERY SMALL VALUES ZERO
      // --------------------------------------------------------------------

      for (char i = 0; i<3; i++) {
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


	//	cout << "singularity with all real evalues"<< endl;
    
	// sort the real values, small absolute values to front
	if ( fabs (l_eigenvalues[0].real()) > fabs (l_eigenvalues[2].real()) ){
	  tmp = l_eigenvalues[0]; l_eigenvalues[0]=l_eigenvalues[2]; l_eigenvalues[2]=tmp; }
	if ( fabs (l_eigenvalues[1].real()) > fabs (l_eigenvalues[2].real()) ){
	  tmp = l_eigenvalues[1]; l_eigenvalues[1]=l_eigenvalues[2]; l_eigenvalues[2]=tmp; }
	if ( fabs (l_eigenvalues[0].real()) > fabs (l_eigenvalues[1].real()) ){
	  tmp = l_eigenvalues[1]; l_eigenvalues[1]=l_eigenvalues[0]; l_eigenvalues[0]=tmp; }

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
	      THROW_EXCEPTION(FException,"Nullspace of tet has dimension zero!");
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
	      setNature(FAMSingularPoint::DEGENERATE);
	      setType(FAMSingularPoint::NONE);
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
		  setType(FAMSingularPoint::REPELL_NODE_3D);
		}
	      else // saddle 3d
		{
		  setType(FAMSingularPoint::SADDLE_3D);
		  FVector tmpVec(eigenvectors[2]);
		  complex < double > tmpClx = l_eigenvalues[2];
		  eigenvectors[2]=eigenvectors[0];
		  l_eigenvalues[2] = l_eigenvalues [0];
		  eigenvectors[0]=tmpVec;
		  l_eigenvalues[0] = tmpClx;
		}
	    else // saddle 3d
	      {
		setType(FAMSingularPoint::SADDLE_3D);
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
		  setType(FAMSingularPoint::ATTRACT_NODE_3D);
		}
	      else // saddle 3d
		{
		  setType(FAMSingularPoint::SADDLE_3D);
		  FVector tmpVec(eigenvectors[2]);
		  complex < double > tmpClx = l_eigenvalues[2];
		  eigenvectors[2]=eigenvectors[0];
		  l_eigenvalues[2] = l_eigenvalues [0];
		  eigenvectors[0]=tmpVec;
		  l_eigenvalues[0] = tmpClx;
		}
	    else // saddle 3d
	      {
		setType(FAMSingularPoint::SADDLE_3D);
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
	setEigenvectors(eigenvectors);
	setEigenvalues(l_eigenvalues);

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
	  cout << "WHILE CALCULATING EIGENVECTOR OF TET-Vectorfield: nullspace"
	       << "does not match (!=1)" << endl;

	if(nullspace.size()<=0)
	{
	  //return;
	  THROW_EXCEPTION(FException,"Nullspace of tet has dimension zero!");
	}
	
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
    
	setNature(FAMSingularPoint::CRITICAL);

	if ( l_eigenvalues[0].real() == 0.0 )
	  if ( l_eigenvalues[1].real() == 0.0 ) // center
	    {
	      setType(FAMSingularPoint::CENTER_3D);
	    }
	  else if ( l_eigenvalues[1].real() > 0.0 ) // repelling spiral FOCUS
	    {
	      setType(FAMSingularPoint::REPELL_FOCUS_3D);
	    }
	  else // attracting spiral FOCUS
	    {
	      setType(FAMSingularPoint::ATTRACT_FOCUS_3D);
	    }
	else if ( l_eigenvalues[0].real() > 0.0 ) // axis repelling
	  if ( l_eigenvalues[1].real() > 0 ) // repelling focus
	    {
	      setType(FAMSingularPoint::REPELL_FOCUS_3D);
	    }
	  else // spiral saddle repelling axis
	    {
	      setType(FAMSingularPoint::SPIRAL_SADDLE_3D);
	    }
	else // axis attracting
	  if ( l_eigenvalues[1].real() < 0 ) // attracting focus
	    {
	      setType(FAMSingularPoint::ATTRACT_FOCUS_3D);
	    }
	  else // spiral saddle attracting axis
	    {
	      setType(FAMSingularPoint::SPIRAL_SADDLE_3D);
	    }

	setEigenvectors(eigenvectors);
	setEigenvalues(l_eigenvalues);
      }



    }
  else {
    THROW_DEFAULT_EXCEPTION( FNotImplementedException );
  }
}

//---------------------------------------------------------------------------

void FAMSingularPoint::setDegenerateType(double a,double
					 b,double c,double d,
					 const FCell *cell)
{
  if (nature==FAMSingularPoint::CRITICAL) {
    THROW_EXCEPTION( FException, "ERROR: function only for degenerate points!");
  }
  double delta;
  delta=a*d-b*c;
  
  if (delta==0) {        // merged degenerate points
    type = FAMSingularPoint::HIGHER_DEGENERATE;
    cerr << "WARNING: higher order degenerate point found!" <<endl;
  } else {
    if (delta<0) {
      type = FAMSingularPoint::TRISECTOR_POINT;
    } else {
      type = FAMSingularPoint::WEDGE_POINT;
    }
    complex<double> x[3];
    int nrSols, i;
    double tmpAngle;
    FTensor testRes(2,2);
    FPosition correction, corrStart;
    // needed for recomputation
    sepAngles.clear();
    
    // handle properly 2D+T cells
    unsigned int dim = cell->getDimension();
    if (dim == 2) {
      correction = FPosition(2);
      corrStart  = FPosition(2);
    }
    else {
      correction = FPosition(3);
      correction[2] = 0.;
      corrStart  = FPosition(3);
      corrStart[2] = pos[2];
    }

    nrSols=FMath::CubicEquation(d, (c+2.0*b), (2.0*a-d), (-1.0*c), x);
    for (i=0;i<nrSols;i++) {
      if (!x[i].imag()) {
        tmpAngle=atan(x[i].real());
        correction[0]=cos(tmpAngle);
        correction[1]=sin(tmpAngle);
        correction*=TEST_DISTANCE;
        corrStart=pos+correction;
        cell->interpolate(testRes, corrStart);
        
        FMatrix testMatr(testRes);
        complex<double> evalue[2];
        double ev0[2], ev1[2];
        double test;
        vector<FVector> testEVecs;
        FVector tmpVec(2), testDir(2);
        testDir[0]=correction[0];
        testDir[1]=correction[1];
        FMath::Eigensystem(testMatr, evalue, ev0, ev1);
        tmpVec[0]=ev0[0];
        tmpVec[1]=ev0[1];
        testEVecs.push_back(tmpVec);
        tmpVec[0]=ev1[0];
        tmpVec[1]=ev1[1];
        testEVecs.push_back(tmpVec);
        if (evalue[0].real() > evalue[1].real()) {
          test=(testEVecs[0]*testDir)/(testEVecs[0].norm()*testDir.norm());
        } else { // EV1 >=EV0
          test=(testEVecs[1]*testDir)/(testEVecs[1].norm()*testDir.norm());
        }
        if ( ( test <  1 - ANGLE_TEST_TOLERANCE ) &&
             ( test > -1 + ANGLE_TEST_TOLERANCE ) ) {
          correction=correction*-1.0;

          if (tmpAngle >= 0)
            tmpAngle -= M_PI;
          else tmpAngle += M_PI;
          
        } else {
        }
        if ( type == FAMSingularPoint::WEDGE_POINT && i == 2) {
          double diff01 = sepAngles[0] - sepAngles[1];
          if (diff01 >  M_PI) diff01 -= 2*M_PI;
          if (diff01 < -M_PI) diff01 += 2*M_PI;
          
          double diff02 = sepAngles[0] - tmpAngle;
          if (diff02 >  M_PI) diff02 -= 2*M_PI;
          if (diff02 < -M_PI) diff02 += 2*M_PI;
          
          double diff12 = sepAngles[1] - tmpAngle;
          if (diff12 >  M_PI) diff12 -= 2*M_PI;
          if (diff12 < -M_PI) diff12 += 2*M_PI;
          
          if (diff01 > 0.0 ) {
            if (diff02 > 0.0) {
              if (diff12 > 0.0) {
                sepAngles[1]=tmpAngle;
              }
              // else discard tmpAngle
            } else {  // diff01 > 0, diff12 < 0
              sepAngles[0]=tmpAngle;
            }
          } else {   // diff01 < 0
            if (diff12 >0.0) {
              if (diff02 >0.0) {
                sepAngles[0]=tmpAngle;
              }
              // else discard tmpAngle
            } else {  // diff01 < 0, diff12 < 0
              sepAngles[1]=tmpAngle;
            }
          }
        } else sepAngles.push_back(tmpAngle);
        
      }
    }
  }
}
      
//--------------------------------------------------------------------------- 

FAMSingularPoint& FAMSingularPoint::operator=(const FAMSingularPoint& sing) {

  if (&sing != this) {
    this->pos = sing.pos;
    
    this->type = sing.type;

    this->nature = sing.nature;

    this->order = sing.order;

    this->cellId = sing.cellId;

    this->eigenvector=sing.eigenvector;
    this->eigenvalues=sing.eigenvalues;

    this->sepAngles = sing.sepAngles;
    
    this->sepStart.clear();
    list<FIndex>::const_iterator iter = sing.sepStart.begin();
    while (iter != sing.sepStart.end()) {
      this->sepStart.push_back(*iter);
      iter++;
    }

    this->sepStop.clear();
    iter = sing.sepStop.begin();
    while (iter != sing.sepStop.end()) {
      this->sepStop.push_back(*iter);
      iter++;
    }

    this->enstrophy = sing.enstrophy;

    this->connections = connections;
  }

  return *this;
}

//---------------------------------------------------------------------------

ostream& operator<< (ostream &os, const FAMSingularPoint::singularityNature& nature) {
  switch (nature) {
  case FAMSingularPoint::CRITICAL:
    os << "FAMSingularPoint::CRITICAL POINT";
    break;
  case FAMSingularPoint::DEGENERATE:
    os <<"FAMSingularPoint::DEGENERATE POINT";
    break;
  case FAMSingularPoint::UNKNOWN:
    os <<"FAMSingularPoint::UNKNOWN NATURE";
    break;
  }
  return(os);
}

//---------------------------------------------------------------------------

ostream& operator<< (ostream &os, const FAMSingularPoint::singularityType& type) {
  switch (type) {
  case FAMSingularPoint::SADDLE_2D: 
    os << "FAMSingularPoint::SADDLE POINT";
    break;
  case FAMSingularPoint::SADDLE_3D: 
    os << "FAMSingularPoint::SADDLE POINT 3D";
    break;
  case FAMSingularPoint::ATTRACT_NODE:
    os << "ATTRACTING NODE";
    break;
  case FAMSingularPoint::REPELL_NODE:
    os << "REPELLING NODE";
    break;
  case FAMSingularPoint::ATTRACT_NODE_3D:
    os << "ATTRACTING NODE_3D";
    break;
  case FAMSingularPoint::REPELL_NODE_3D:
    os << "REPELLING NODE_3D";
    break;
  case FAMSingularPoint::ATTRACT_FOCUS:
    os << "ATTRACTING FOCUS";
    break;
  case FAMSingularPoint::REPELL_FOCUS:
    os << "REPELLING FOCUS";
    break;
  case FAMSingularPoint::ATTRACT_FOCUS_3D:
    os << "ATTRACTING FOCUS_3D";
    break;
  case FAMSingularPoint::REPELL_FOCUS_3D:
    os << "REPELLING FOCUS_3D";
    break;
  case FAMSingularPoint::CENTER_2D:
    os << "FAMSingularPoint::CENTER POINT";
    break;
  case FAMSingularPoint::CENTER_3D:
    os << "FAMSingularPoint::CENTER LINE 3D";
    break;
  case FAMSingularPoint::NONE:
    os << "FAMSingularPoint::UNKNOWN/UNDEFINED TYPE";
    break;
  case FAMSingularPoint::HIGHER:
    os << "FAMSingularPoint::HIGHER ORDER TYPE";
    break;
  case FAMSingularPoint::BOUNDARY:
    os << "FAMSingularPoint::BOUNDARY FAMSingularPoint::SOURCE/FAMSingularPoint::SINK";
    break;
  case FAMSingularPoint::TRISECTOR_POINT:
    os << "FAMSingularPoint::TRISECTOR POINT";
    break;
  case FAMSingularPoint::WEDGE_POINT:
    os << "FAMSingularPoint::WEDGE POINT";
    break;
  case FAMSingularPoint::WEDGE_SINGLE_POINT:
    os << "FAMSingularPoint::WEDGE SINGLE POINT";
    break;
  case FAMSingularPoint::HIGHER_DEGENERATE:
    os << "FAMSingularPoint::HIGHER FAMSingularPoint::DEGENERATE POINT";
    break;
  case FAMSingularPoint::SPIRAL_SADDLE_3D:
    os << "SPIRAL FAMSingularPoint::SADDLE_3D";
    break;
  }
  return os;
}
    
//---------------------------------------------------------------------------

ostream& operator<< (ostream &os, const FAMSingularPoint::mainType& type) {
  switch (type) {
  case FAMSingularPoint::SADDLE: 
    os << "FAMSingularPoint::SADDLE";
    break;
  case FAMSingularPoint::SOURCE:
    os << "FAMSingularPoint::SOURCE";
    break;
  case FAMSingularPoint::SINK:
    os << "FAMSingularPoint::SINK";
    break;
  case FAMSingularPoint::CENTER:
    os << "FAMSingularPoint::CENTER";
    break;
  case FAMSingularPoint::TRISECTOR:
    os << "FAMSingularPoint::TRISECTOR";
    break;
  case FAMSingularPoint::WEDGE:
    os << "FAMSingularPoint::WEDGE";
    break;
  case FAMSingularPoint::WEDGE_SINGLE:
    os << "FAMSingularPoint::WEDGE SINGLE";
    break;
  case FAMSingularPoint::OTHER:
    os << "UNDEFINED/FAMSingularPoint::UNKNOWN";
    break;
  }
  return os;
}

//---------------------------------------------------------------------------

ostream& operator<< (ostream &os, const FAMSingularPoint& singularity)
{
  os << "FAMSingularPoint:" << endl;
  os << "--position: " << singularity.pos << endl;
  os << "--type: " << singularity.type << endl;
  os << "--order: " << singularity.order << endl;
  os << "--enstrophy: " << singularity.enstrophy << endl;
  os << "--cellId: " << singularity.cellId << endl;
  for (positive i=0 ; i<singularity.eigenvector.size() ; i++)
    os << "--eigenvector " << i << ": " << singularity.eigenvector[i];
  for (positive i=0 ; i<singularity.sepAngles.size() ; i++)
    os << "--sepAngle " << i <<": "<< singularity.sepAngles[i];
  os << "connections : ";
  for (positive i=0 ; i<singularity.connections.size() ; i++)
    os << singularity.connections[i] << " ";
  os << endl;
  list<FIndex>::const_iterator iter;
  iter = singularity.sepStart.begin();
  while (iter != singularity.sepStart.end()) {
    cout << "--starting separatrices: " << (*iter) << endl;
    iter++;
  }
  iter = singularity.sepStop.begin();
  while (iter != singularity.sepStop.end()) {
    cout << "--stopping separatrices: " << (*iter) << endl;
    iter++;
  }

  return os;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::getEigenvectors(std::vector<FVector>& eigenVec) const{
  eigenVec = eigenvector;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::setEigenvectors(const std::vector<FVector>& eigenVec) {
  eigenvector = eigenVec;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::getEigenvalues(std::vector< std::complex<double> >& eigenVal) const{
  eigenVal = this->eigenvalues;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::setEigenvalues(const std::vector< std::complex<double> >& eigenVal) {
  this->eigenvalues = eigenVal;
}

//---------------------------------------------------------------------------

void FAMSingularPoint::getSeparatrixAngles(std::vector<double>& sAngles) 
{
  sAngles=sepAngles;
}

//---------------------------------------------------------------------------

void FAMSingularPoint::setSeparatrixAngles(const std::vector<double>& sAngles) 
{
  sepAngles.resize(sAngles.size());
  
  for (positive i=0; i<sAngles.size(); i++)
    sepAngles[i] = sAngles[i];
}

//---------------------------------------------------------------------------

void FAMSingularPoint::getPosition(FPosition& result) const {

  if (pos.getDimension() == 0) {
    THROW_DEFAULT_EXCEPTION( FEmptyObjectException );
  }
  result = pos;
}

//---------------------------------------------------------------------------

void FAMSingularPoint::setPosition(const FPosition& pos) {

  if ((pos.getDimension() < 2) ||(pos.getDimension() > 3)) {
    THROW_DEFAULT_EXCEPTION( FInvalidDimensionException );
  }
  this->pos = pos;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::getType(FAMSingularPoint::singularityType& result) const {

  if (pos.getDimension() == 0 || order == 0) {
    THROW_EXCEPTION( FException, "ERROR: no type has been set" );
  }
  result = type;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::setType(FAMSingularPoint::singularityType val) {
  this->type = val;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::getNature(FAMSingularPoint::singularityNature& result) const {

  if (pos.getDimension() == 0 || order == 0) {
    THROW_EXCEPTION( FException, "ERROR: no type has been set");
  }
  result = nature;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::setNature(FAMSingularPoint::singularityNature val) {
  this->nature = val;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::setOrder( positive order)
{
  this->order=order;
}

//--------------------------------------------------------------------------- 

positive FAMSingularPoint::getOrder() const {

  if (pos.getDimension() == 0 || order == 0) {
    THROW_EXCEPTION( FException, "ERROR: no order has been set");
  }
  return order;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::getIncludingCellIndex(FIndex& cellId) const {

  if (!this->cellId.isValid()) {
    THROW_EXCEPTION( FException, "ERROR: no cell index has been set");
  }

  cellId = this->cellId;
}

//---------------------------------------------------------------------------

void FAMSingularPoint::getStartingSeparatrices(std::list<FIndex>& starts) const
{
  starts.resize(0);

  list<FIndex>::const_iterator iter = sepStart.begin();

  while (iter != sepStart.end()) {
    
    starts.push_back(*iter);
    iter++;
  }
}

//---------------------------------------------------------------------------

void FAMSingularPoint::getStoppingSeparatrices(std::list<FIndex>& stops) const
{
  stops.resize(0);

  list<FIndex>::const_iterator iter = sepStop.begin();

  while (iter != sepStop.end()) {
    
    stops.push_back(*iter);
    iter++;
  }
}

//---------------------------------------------------------------------------

void FAMSingularPoint::addStartingSeparatrix(const FIndex& start)
{
  sepStart.push_back(start);
}

//---------------------------------------------------------------------------

void FAMSingularPoint::addStoppingSeparatrix(const FIndex& stop)
{
  sepStop.push_back(stop);
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::addConnection(const FIndex& connection)
{
  connections.push_back(connection);
  //cout << "connection added" << endl;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::setConnections(const std::vector< FIndex >& connections)
{
//   cout << "in setConnections: connections.size() = " 
//        << connections.size() << endl;

  try {
    switch (type) {
    case FAMSingularPoint::SADDLE_2D: {
      if (connections.size() != 4) {
	THROW_EXCEPTION( FInvalidDimensionException, "ERROR: invalid number of reached singularities");
      }
      else
	this->connections = connections;
      break;
    }
    case FAMSingularPoint::TRISECTOR_POINT: {
      if (connections.size() != 3) {
	THROW_EXCEPTION( FInvalidDimensionException, "ERROR: invalid number of reached singularities");
      }
      else
	this->connections = connections;
      break;
    }
    case FAMSingularPoint::WEDGE_POINT: {
      if (connections.size() != 2) {
	THROW_EXCEPTION( FInvalidDimensionException, "ERROR: invalid number of reached singularities");
      }
      else
	this->connections = connections;
      break;
    }
    case FAMSingularPoint::WEDGE_SINGLE_POINT: {
      if (connections.size() != 1) {
	THROW_EXCEPTION( FInvalidDimensionException, "ERROR: invalid number of reached singularities");
      }
      else
	this->connections = connections;
      break;
    }
    default: {
      this->connections = connections;
    }
    }
  }
  CATCH_N_RETHROW( FException );

//   cout << "leaving setConnections: " << this->connections.size() << endl;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::getConnections(std::vector< FIndex >& connections) const
{
  connections = this->connections;
}

//--------------------------------------------------------------------------- 

void FAMSingularPoint::setMaxNorm(double max_norm) {

  this->max_norm = max_norm;
} 

//---------------------------------------------------------------------------

