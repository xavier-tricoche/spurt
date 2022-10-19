//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSprMatrix.cc,v $
// Language:  C++
// Date:      $Date: 2003/03/27 08:29:53 $
// Author:    $Author: garth $
// Version:   $Revision: 1.7 $
//
//--------------------------------------------------------------------------- 

#include "FPosition.hh"
#include "FVector.hh"
#include "FMatrix.hh"
#include "FException.hh"

#include "FSprMatrix.hh"

#include <cmath>
#include <iostream>

using namespace std;

FSprMatrix::FSprMatrix( void )
{
  //  do nothing as this constructor hopefully never is called...
}


FSprMatrix::FSprMatrix(const FSprMatrix& m)
    : ija(), sa()
{

  ija=m.ija; // allocate memory
  sa=m.sa;

  fill_m=m.fill_m;
  fill_n=m.fill_n;
  
  dimension = m.dimension;
  _map = NULL;
  EPS=1.0e-14;
}

/*
  Converts a square matrix a[1..n][1..n] into row-indexed sparse storage mode.
  Only elements of a with magnitude  thresh are retained. Output is in two
  linear arrays with dimension nmax (an input parameter): sa[1..] contains
  array values, indexed by ija[1..]. The number of elements lled of sa and
  ija on output are both ija[ija[1]-1]-1 (see text).
*/

FSprMatrix::FSprMatrix (const FMatrix &m, double thresh)
{
  unsigned int i,j;
  unsigned int space;
  unsigned int k;

  dimension = m.getDimensionX();

#ifndef NODEBUG
  if (m.getDimensionY() != m.getDimensionX())
    {
      FInvalidDimensionException e ("Matrix must be square to fit into Sparse");
      e.addTraceMessage ("FSprMatrix::FSprMatrix (FMatrix &m, double thresh)");
      throw(e);
    }
#endif

  // do some guess how much space it may take...
  space = dimension +2; // just the minimum without any offdiagonal elements

  for (i=0;i<dimension;i++)
    for (j=0;j<dimension;j++)
      if (fabs(m(i,j)) >= thresh && i != j)
        space++; // look how many entries are needed...

  sa.resize(space);  // allocate the Memory
  ija.resize(space);
  

  for (j=0;j<dimension;j++) sa[j]=m(j,j); /* Store diagonal elements.*/
  ija[0]=dimension+1; /* Index to 1st row offdiagonal element, if any. */

  k=dimension;
  for (i=0;i<dimension;i++) { /* Loop over rows.*/
    for (j=0;j<dimension;j++) { /* Loop over columns.*/
      if (fabs(m(i,j)) >= thresh && i != j) {
        
#ifndef NODEBUG        
        if (++k > space) // end of array, resize !
          {
            FException e ("Error in Source Code (space not calculated correctly");
            e.addTraceMessage ("FSprMatrix::FSprMatrix (const FMatrix &m, double thresh)");
            throw e;            
          }
#endif        
        
        sa[k]=m(i,j);  /* Store offdiagonal elements and their columns. */
        ija[k]=j;
      }
    }
    ija[i+1]=k+1; /* As each row is completed, store index to next. */
  }
  _map = NULL;
  EPS=1.0e-14;
}

//--------------------------------------------------------------------------- 

FSprMatrix::FSprMatrix(unsigned int m, double fill_estimate)
{
  unsigned int i;

#ifndef NODEBUG        
  if ( fill_estimate > 1.0 )
    {
      FException e ("The estimate for the amount of nonzero elements must be 0<x<=1");
      e.addTraceMessage("FSprMatrix::FSprMatrix(unsigned int m, double fill_estimate)");
      throw e;
    }
#endif
  
  // points to 1st element that's not in list
  ija.resize (m+1 + (unsigned int) ((double)m*(double)m*fill_estimate) );
  sa.resize  (m+1 + (unsigned int) ((double)m*(double)m*fill_estimate) );

  fill_m=fill_n=0; // preset to see where fill-up begins...

  for(i=0; i<m+1;i++)
    ija[i]=m+1;
  //  ija [m+1] = m+2; // points to 1st element that's not in list

  dimension = m;
  _map = NULL;
  EPS=1.0e-14;
}

//--------------------------------------------------------------------------- 

FSprMatrix::~FSprMatrix()
{
  if(_map) delete _map;  // C++ doesn't initialize Pointers to NULL, therefore be CAREFUL
}

//--------------------------------------------------------------------------- 

void FSprMatrix::_reset_(unsigned int m)
{
  unsigned int i;

  // points to 1st element that's not in list

  ija.resize (m+1);
  sa.resize  (m+1);

  fill_m=fill_n=0; // preset to see where fill-up begins...

  for(i=0; i<m+1;i++)
    {
      ija[i]=m+1;
      sa[i]=0.0;
    }
  //  ija [m+1] = m+2; // points to 1st element that's not in list

  dimension = m;
  if (_map) delete _map;
  _map = NULL;
}


//--------------------------------------------------------------------------- 

FSprMatrix FSprMatrix::operator* (const double lambda)
{
  FSprMatrix result(*this);

  unsigned int i;

  for (i=0; i<ija[ija[0]-1]-1; i++)
    result.sa[i] *=  lambda;
  
  return (result);
}

//--------------------------------------------------------------------------- 

FSprMatrix& FSprMatrix::operator*= (const double lambda)
{
  unsigned int i;

  // when multiplying with a value geometry of the matrix doesn´t change...
  for (i=0; i<ija[ija[0]-1]-1; i++)
    sa[i] *=  lambda;
  
  return (*this);
}

//--------------------------------------------------------------------------- 

/*
  Multiply a matrix in row-index sparse storage arrays sa and ija by a
  vector x[1..n], giving a vector b[1..n].
*/

FVector FSprMatrix::operator*(const FVector &v)
{
  unsigned long i,k;
  unsigned long n = v.getDimension();

  FVector r (n);

#ifndef NODEBUG  
  if (ija[0] != n+1)
    {
      FInvalidDimensionException e ("matrix and vector don't match in dimensions");
      e.addTraceMessage("FSprMatrix::FVector operator*(const FVector &v)");
      throw e;
    }

#endif
  for (i=0;i<n;i++) {
    r[i]=sa[i]*v[i]; /* Start with diagonal term.*/
    for (k=ija[i];k<ija[i+1];k++) /* Loop over offdiagonal terms.*/
      r[i] += sa[k]*v[ija[k]];
  }
  return (r);
}

//--------------------------------------------------------------------------- 

/*
  Multiply the transpose of a matrix in row-index sparse storage arrays
  sa and ija by a vector x[1..n], giving a vector b[1..n].
*/

FVector FSprMatrix::tmult ( FVector &v)
{
  unsigned long i,k;
  unsigned long n = v.getDimension ();

  FVector r (n);
  
#ifndef NODEBUG  
  if (ija[0] != n+1)
    {
      FInvalidDimensionException e ("matrix and vector don't match in dimensions");
      e.addTraceMessage("FSprMatrix::FVector tmult(const FVector &v)");
      throw e;
    }

#endif

  for (i=0;i<n;i++) r[i]=sa[i]*v[i]; /* Start with diagonal terms. */

  for (i=0;i<n;i++) { /* Loop over o-diagonal terms.*/
    for (k=ija[i];k<ija[i+1];k++) {
      r[ija[k]] += sa[k]*v[i];
    }
  }
  return (r);
}


//--------------------------------------------------------------------------- 

/* IINDEXX */
#define SwaP(a,b) itemp=(a);(a)=(b);(b)=itemp;
#define M 7
#define NSTACK 50

void iindexx(unsigned long n, std::vector<unsigned int>& arr,unsigned int starta, std::vector<unsigned int>& indx ,unsigned int startb )
  /*  Indexes an array arr[starta+1..n], i.e., outputs the array indx[startb+1..n] such that arr[starta+indx[startb+j]] is
      in ascending order for j = 1;2;: : :;N. The input quantities n and arr are not changed.*/
{
  unsigned int i,indxt,ir=n,itemp,j,k,l=1,
    uplimit = NSTACK;
  unsigned int jstack=0;
  std::vector<unsigned int> istack( uplimit+1 ); // array of int, 50 is the initial value
  double a;

  for (j=1;j<=n;j++) indx[startb+j]=j;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
        indxt=indx[startb+j];
        a=arr[starta+indxt];
        for (i=j-1;i>=l;i--) {
          if (arr[starta+indx[startb+i]] <= a) break;
          indx[startb+i+1]=indx[startb+i];
        }
        indx[startb+i+1]=indxt;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SwaP(indx[startb+k],indx[startb+l+1]);
      if (arr[starta+indx[startb+l]] > arr[starta+indx[startb+ir]]) {
        SwaP(indx[startb+l],indx[startb+ir])
          }
      if (arr[starta+indx[startb+l+1]] > arr[starta+indx[startb+ir]]) {
        SwaP(indx[startb+l+1],indx[startb+ir])
          }
      if (arr[starta+indx[startb+l]] > arr[starta+indx[startb+l+1]]) {
        SwaP(indx[startb+l],indx[startb+l+1])
          }
      i=l+1;
      j=ir;
      indxt=indx[startb+l+1];
      a=arr[starta+indxt];
      for (;;) {
        do i++; while (arr[starta+indx[startb+i]] < a);
        do j--; while (arr[starta+indx[startb+j]] > a);
        if (j < i) break;
        SwaP(indx[startb+i],indx[startb+j])
          }
      indx[startb+l+1]=indx[startb+j];
      indx[startb+j]=indxt;
      jstack += 2;
      if (jstack > uplimit)
        {
          istack.resize (uplimit+=10);
        }
      if (ir-i+1 >= j-l) {
        istack[jstack]=ir;
        istack[jstack-1]=i;
        ir=j-1;
      } else {
        istack[jstack]=j-1;
        istack[jstack-1]=l;
        l=i;
      }
    }
  }
}

//--------------------------------------------------------------------------- 

/*
  Construct the transpose of a sparse square matrix, from row-index sparse
  storage arrays sa and ija into arrays sb and ijb.
*/

FSprMatrix& FSprMatrix::transpose( void )


  /*TRANSPONIEREN*/  
  // void sprstp(float sa[], unsigned long ija[], float sb[], unsigned long ijb[])
  
{
  //  void iindexx(unsigned long n, long arr[], unsigned long indx[]);
  //  Version of indexx with all float variables changed to long.

  unsigned int j,jl,jm,jp,ju,k,m,n2,noff,inc,iv;
  double v;

  // make a copy of self, then operate as if self is the destination.
  std::vector<unsigned int> ijb (ija); // copy 
  std::vector<double> sb (sa); // copy

  n2=ijb[0]; // Linear size of matrix plus 2.
  for (j=0;j<n2-1;j++) sa[j]=sb[j]; // Diagonal elements.
  
  iindexx(ijb[n2-1]-ijb[0],ijb, n2-1, ija, n2-1);
  
  // Index all off-diagonal elements by their columns.
  jp=0;
  jp--;
  for (k=ijb[0];k<ijb[n2-1];k++) { // Loop over output off-diagonal elements.
    m=ija[k]+n2-1; // Use index table to store by (former) columns.
    sa[k]=sb[m];
    for (j=jp+1;j<=ijb[m];j++) ija[j]=k;
    //Fill in the index to any omitted rows.
    jp=ijb[m];
    // Use bisection to nd which row element
    // m is in and put that into ija[k].
    jl=0;
    ju=n2-1;
    while (ju-jl > 1) {
      jm=(ju+jl)/2;
      if (ijb[jm] > m) ju=jm; else jl=jm;
    }
    ija[k]=jl;
  }
  for (j=jp+1;j<n2;j++) ija[j]=ijb[n2-1];
  for (j=0;j<n2-1;j++) {
    // Make a nal pass to sort each row by
    //Shell sort algorithm.
    jl=ija[j+1]-ija[j];
    noff=ija[j]-1;
    inc=1;
    do {
      inc *= 3;
      inc++;
    } while (inc <= jl);
    do {
      inc /= 3;
      for (k=noff+inc+1;k<=noff+jl;k++) {
        iv=ija[k];
        v=sa[k];
        m=k;
        while (ija[m-inc] > iv) {
          ija[m]=ija[m-inc];
          sa[m]=sa[m-inc];
          m -= inc;
          if (m-noff <= inc) break;
        }
        ija[m]=iv;
        sa[m]=v;
      }
    } while (inc > 1);
  }

  return ( *this );
  
}

//--------------------------------------------------------------------------- 

/*
  Solves A  x = b for x[1..n], givenb[1..n], by the iterative biconjugate
  gradient method. On input x[1..n] should be set to an initial guess of the
  solution (or all zeros); itol is 1,2,3, or 4, specifying which convergence
  test is applied (see text); itmax is the maximum number of allowed iterations;
  and tol is the desired convergence tolerance. On output, x[1..n] is reset to
  the improved solution, iter is the number of iterations actually taken, and err
  is the estimated error. The matrix A is referenced only through the user-supplied
  routines atimes, which computes the product of either A or its transpose on a
  vector; and asolve, which solves

  ( look into numerical recipes c++ page 86 bottom , 2.7.45 )
  
*/

// SOLVES for trivial diagonalpart matrix
void FSprMatrix::asolve(FVector& b,FVector& x)
{
  for (unsigned int i=0; i<dimension; i++)
    if (fabs(sa[i]) > EPS)
      x[i] = b[i] / sa[i];
    else
      x[i] = b[i];
}


//--------------------------------------------------------------------------- 

// Compute one of two norms for a vector sx[1..n], as signaled by itol. Usedbylinbcg.
double FSprMatrix::snrm(FVector& sx,unsigned int itol)
{
  unsigned int i,isamax;
  unsigned int n = dimension;
  double ans;
  if (itol <= 3) {
    ans = 0.0;
    for (i=0;i<n;i++) ans += sx[i]*sx[i]; // Vector magnitude norm.
    return sqrt(ans);
  } else {
    isamax=0;
    for (i=0;i<n;i++) { // Largest component norm.
      if (fabs(sx[i]) > fabs(sx[isamax])) isamax=i;
    }
    return fabs(sx[isamax]);
  }
}

FVector FSprMatrix::solve_bcg ( FVector &b, unsigned int itol, double tol, unsigned int itmax, unsigned int & iter, double & err)
{
  
  // void linbcg(TYPE* sa, int* ija, unsigned long n, TYPE b[], TYPE x[], int itol,
  //		TYPE tol, int itmax, int *iter, TYPE *err) {

  unsigned long j, n=dimension;
  double ak, akden, bk, 
    bkden=0., //initialized to avoid warning
    bknum, 
    bnrm=0., //initialized to avoid warning
    dxnrm, 
    xnrm, 
    zm1nrm, 
    znrm=0. ;//initialized to avoid warning
  FVector p(n), pp(n), r(n), rr(n), z(n), zz(n), x(n);

  for(unsigned int i=0; i<n; i++)
    p[i] = pp[i] = r[i] = rr[i] = z[i] = zz[i] = x[i] = 0.0;

  
  /* Calculate initial residual */
	iter = 0 ;
	r=(*this)*x ;
	for (j=0 ; j<n ; j++) {
		r[j] = b[j]-r[j] ;
		rr[j] = r[j] ;
	}

	//atimes(r,rr,0) ; 
	if (itol == 1) {
		bnrm = snrm(b, itol) ;
		asolve(r, z) ;
	}
	else if (itol == 2) {
		asolve(b, z  ) ;
		bnrm = snrm(z, itol) ;
		asolve(r, z) ;
	}
	else if (itol == 3 || itol == 4) {
		asolve(b, z) ;
		bnrm = snrm(z, itol) ;
		asolve(r, z) ;
		znrm = snrm(z, itol) ;
	} else
#ifndef NODEBUG
    {
      FException e ("illegal itol in linbcg");
      e.addTraceMessage ("FVector FSprMatrix::solve_bcg ( FVector & v, unsigned int itol, double tol, unsigned int itmax, unsigned int & iter, double & err)");
      throw e;
    }
#endif
  ;
  
	while (iter <= itmax) {
		++(iter) ;
		asolve(rr, zz) ;
		for (bknum=0.0, j=0 ; j<n ; j++) bknum += z[j]*rr[j] ;
		// calculate coefficient bk and direction vectors p and pp 
		if (iter == 1) {
			for (j=0 ; j<n ; j++) {
				p[j]=z[j] ;
				pp[j]=zz[j] ;
			}
		}
		else {
			bk = bknum/bkden ;
			for (j=0 ; j<n ; j++) {
				p[j] = bk*p[j]+z[j] ;
				pp[j] = bk*pp[j]+zz[j] ;
			}
		}
		bkden = bknum ;
		z=(*this)* p ;
		for (akden=0.0, j=0 ; j< n ; j++) akden += z[j]*pp[j] ;
		ak = bknum / akden ;
    zz= tmult( pp );
		for (j=0 ; j<n ; j++) {
			x[j] += ak*p[j] ;
			r[j] -= ak*z[j] ;
			rr[j] -= ak*zz[j] ;
		}
		asolve(r, z) ;
		if (itol == 1) 
			err=snrm( r, itol)/bnrm ;
		else if (itol == 2)
			err = snrm( z, itol)/bnrm ;
		else if (itol == 3 || itol == 4) {
			zm1nrm = znrm ;	
			znrm = snrm( z, itol) ;
			if (fabs(zm1nrm-znrm) > EPS*znrm) {
				dxnrm = fabs(ak)*snrm(p, itol) ;
				err = znrm/fabs(zm1nrm-znrm)*dxnrm ;
			}
			else {
				err = znrm/bnrm ;
				continue ;
			}
			xnrm = snrm(x, itol) ;
			if (err <= 0.5*xnrm) err /= xnrm ;
			else {
				err = znrm/bnrm ;
				continue ;
			}
		}
    if (err <= tol) 
			 break ;
  }

  return (x);
}

//--------------------------------------------------------------------------- 

unsigned int FSprMatrix::space( void ) const
{
  // this should be the amount of used entries,
  // since C++ runs arrays from 0 we should add 1
  return ( ija[ija[0]-1]);

  // otherwise, why not use "sizeof" ??
  //  return ( sizeof(*ija) );
}

//--------------------------------------------------------------------------- 

unsigned int FSprMatrix::getDimension(void)
{
  return (dimension);
}

//--------------------------------------------------------------------------- 

double FSprMatrix::getElement(unsigned int m, unsigned int n) const
{

#ifndef NODEBUG
  if ( (m>=dimension) || (n>=dimension) )
    {
      FIndexOutOfBoundsException e ("Index greater than Matrixdimension");
      e.addTraceMessage("double FSprMatrix::getElement(unsigned int m, unsigned int n) const");
      throw e;
    }
#endif
  
  if (n == m)
    return (sa[m]);
  else
    for (unsigned int i=ija[m];i<ija[m+1];i++)
      if (ija[i] == n)
        {
          return (sa[i]);
        }
  return (0.0); // kein Eintrag gefunden -> Nullelement
}

//--------------------------------------------------------------------------- 

#if 0
// buggy, see test results
double FSprMatrix::setElement(unsigned int m, unsigned int n, double val)
{


  //      cout << endl << "Halli :" << fill_m << " und " << fill_n;
  //    cout << endl << "Hallo :" << m << " und " << n;
  
  unsigned int endindex, i;
  if ( (m>dimension) || (n>dimension) )
    {
      FIndexOutOfBoundsException e ("Index greater than Matrixdimension");
      e.addTraceMessage("double& FSprMatrix::setElement(unsigned int m, unsigned int n)");
      throw e;
    }
  else
    if ( (m<fill_m) || ((m==fill_m)&&(n<fill_n))) // matrix iss filled up beyond the desired point already
      {
        FException e ("Matrix must be filled left - right and top - down");
        e.addTraceMessage("double& FSprMatrix::setElement(unsigned int m, unsigned int n)");
        throw e;
      }

  /* too small values shouldn't be taken */
  if (fabs(val) < EPS) return (0.0);

  if (n == m)
    return ( sa[n]=val ); // diagonal element should be no problem
  else
    {
      endindex = ija[ija[1]-1];
      ija[ija[1]-1]++; // push a bit
      if (ija.size() <= endindex) // vectorarrays too small
        {
          ija.resize(ija.size()+10);
          sa.resize(ija.size()+10);
        }

      // rest of array has to be updated
      for (i=m+1; i<ija[1]-1;i++)
        ija[i] = endindex+1;

      fill_n = n;
      fill_m = m;

      // return reference to first Element that's not in the list,
      // and set ija-elements accordingly
      ija[ endindex ] = n;

      sa[endindex] = val;
      return ( sa[endindex]);
    }

}
#endif

//--------------------------------------------------------------------------- 

void FSprMatrix::cueElement (unsigned int m, unsigned int n, double val)
{

#ifndef NODEBUG
  if ( (m>=dimension) || (n>=dimension) )
    {
      FIndexOutOfBoundsException e ("Index greater than Matrixdimension");
      e.addTraceMessage("void FSprMatrix::cueElement (unsigned int m, unsigned int n, double val)");
      throw e;
    }
#endif

  std::pair<unsigned int, unsigned int> _p;
  _p.first=m;
  _p.second=n;
  // new map 
  if (_map == NULL)
    {
      _map = new map<pair<unsigned int, unsigned int>, double> ();
    }
  
  
  if ( fabs(val) >= EPS )
    (*_map)[_p] = val; // insert pair...
}

//--------------------------------------------------------------------------- 

void FSprMatrix::commitCue ( void )
{
  unsigned int i,j;
  unsigned int space;
  unsigned int k;
  pair<unsigned int, unsigned int> _p;

  // do some guess how much space it may take...
  space = dimension +4 + (*_map).size(); // just the minimum without any offdiagonal elements

  sa.resize(space);  // allocate the Memory
  ija.resize(space);

  for (j=0;j<dimension;j++)
    {
      _p.first = _p.second = j;
      sa[j]=(*_map)[_p]; /* Store diagonal elements.*/
    }
  ija[0]=dimension+1; /* Index to 1st row offdiagonal element, if any. */

  k=dimension;
  for (i=0;i<dimension;i++) { /* Loop over rows.*/
    for (j=0;j<dimension;j++) { /* Loop over columns.*/
      _p.first=i;
      _p.second=j;
      if (i != j && ((*_map).find(_p) != (*_map).end())) {
        ++k;
#ifndef NODEBUG        
        if ( k >= space) // end of array, resize !
          {
            FException e ("Error in Source Code (space not calculated correctly");
            e.addTraceMessage ("void FSprMatrix::commitCue ( void )");
            throw e;
          }
#endif        

        sa[k]=(*_map)[_p];  /* Store offdiagonal elements and their columns. */
        ija[k]=j;
      }
    }
    ija[i+1]=k+1; /* As each row is completed, store index to next. */
  }
  // kill map..
  delete _map;
  _map = NULL;

  sa.resize(ija[ija[0]-1]);  // allocate the Memory
  ija.resize(ija[ija[0]-1]);

}

//--------------------------------------------------------------------------- 

void FSprMatrix::setTreshold ( double t )
{
  EPS = fabs(t);
}

//--------------------------------------------------------------------------- 

// be very cautious when using this operator, 100x100 matrices wouldn't be
// so nice on the screen
std::ostream& operator<<(std::ostream& os, const FSprMatrix& m)
{
  unsigned int i,j,k,l;

  for(i=0;i<m.dimension; i++)
    {
      os << endl;
      l = m.ija[i];
      k = m.ija[i+1];
      
      for (j=0;j<m.dimension;j++)
        {
          if (i==j)
            os << m.sa[i]; // diagonal element
          else if (m.ija[l]==j && l<k)
            os << m.sa[l++];
          else
            os << "0";
          os << "  "; // insert spaces for better readability
        }
    }

  os << endl;

  os << "INDXlooks like : ";
  for (i=0;i<m.space();i++)
    {
      os << i << "  ";
    }
  os << endl;
  os << "ija looks like : ";
  for (i=0;i<m.space();i++)
    {
      os << (m.ija[i]) << "  ";
    }
  os << endl;
  os << "sa  looks like : ";
  for (i=0;i<m.space();i++)
    {
      os << m.sa[i] << "  ";
    }
  os << endl;
  
  return os;
}
