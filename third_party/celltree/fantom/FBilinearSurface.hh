#ifndef FBILINEARSURFACE_HH
#define FBILINEARSURFACE_HH

#include <math.h>
#include "Fd3op.hh"




/** 
 * this class represents
 *  the bilinear surface 
 *  ( { (1-b)((1-a)F00 + a F10 ) 
 *       + b ((1-a)F01 + a F11 ) , a,b in R } )
 *  with Fij in R^3
 * by F00 + a*(df/da(0,0)) + b*(df/db)(0,0) + a*b*(d^2f/(dadb))  
 * and gives some basic geometric operations
 * for it
 */
class FBilinearSurface{

public:

  /**
     empty constructor
     needed for the use in FPyramidCell
   */
  FBilinearSurface()
  {
  }

  /**
    same as constructor, needed for FPyramidCell
   */
  template<class ArrayType>
  inline
  void init(const ArrayType&F00,const ArrayType&F10,
	    const ArrayType&F11,const ArrayType&F01);

  /**
   *\par
   * constructor. 
   *\rem ArrayType is a Type who has
   * an operator[] defined which returns a scalar.
   * \param F00,F01,F11,F10: 
   * Fij is the value of the bilinear function
   * with the parameters i=I, j=J
   * 
  */
  template<class ArrayType>
  inline
  FBilinearSurface(const ArrayType&F00,const ArrayType&F10,
		   const ArrayType&F11,const ArrayType&F01);




  /**
   *\par
   *  gives the points where the surface
   *  and the the line through start in direction dir
   *  intersect
   *  in local coordinates 
   *\par 
   * ( 1st and 2nd entries are the parameters
   *  of the bilinear surface,
   *  3rd entry is the parameter t of the line (x + t*s = intersection point ) )
   *\param x:
   *startpoint of the line
   *\param s: 
   * direction vector of the line
   *\param solutions: 2 vectors with the solutions 
   *in the format described above
   *\return number of solutions
   */
   template<class ArrayType>
  int cutWithLine(const ArrayType &x,const ArrayType &s,
		  ArrayType solutions[2],bool refine=true) const;

  /**
   *refines one solution of cutWithLine by newton iteration
   *\param start
   *startpoint of the line
   *\param dir: 
   * direction vector of the line
   *\param solutions: 
   * vector with solution
   * to be refined
   *\return number of refinemwnt steps used
   *-> =0 if no refinement was applied
   */
  template<class ArrayType>
  int refineCutWithLine(const ArrayType &x,const ArrayType &s,
			 ArrayType& solution) const;

  /**
   *interpolate position x at local coords a,b
   */
  template<class ArrayType>
  inline
  void interpolate(ArrayType &x,double a,double b) const;

  /**
   *\par
   *interpolate derivatives respective 1st local coord
   *\retval dxda: derivation after first local coord
   *\param b: 2nd local coord
   */
  template<class ArrayType>
  inline
  void deriva(ArrayType &dxda,double b) const;

  /**
   *\par
   *interpolate derivatives respective 2nd local coord
   *\retval dxdb: derivation after second local coord
   *\param a: 1st local coord
   */
  template<class ArrayType>
  inline
  void derivb(ArrayType &dxdb,double a) const;


  /**
   * get surface normal
   * at position a,b by crossproduct
   * dxda x dxdb
   *\param a 1st local coord
   *\param b 2nd local coord
   *\retval n
   * normal
   */
  template<class ArrayType>
  inline
  void normal(ArrayType &n,double a,double b) const;

 
protected:
  
  double f00[3],fab[3],fa0[3],f0b[3];

};

template<class ArrayType>
FBilinearSurface::
FBilinearSurface(const ArrayType&F00,const ArrayType&F10,
		 const ArrayType&F11,const ArrayType&F01)
{  
	for(int i=0; i<3; i++)
	{
		f00[i] = F00[i];
		fa0[i] = F10[i] - f00[i];
		f0b[i] = F01[i] - f00[i];
		fab[i] = F11[i] - F01[i] - fa0[i];
	}
}

//----------------------------------------------------------


template<class ArrayType>
void
FBilinearSurface::
init(const ArrayType&F00,const ArrayType&F10,
     const ArrayType&F11,const ArrayType&F01)
{  
	for(int i=0; i<3; i++)
	{
		f00[i] = F00[i];
		fa0[i] = F10[i] - f00[i];
		f0b[i] = F01[i] - f00[i];
		fab[i] = F11[i] - F01[i] - fa0[i];
	}
}


template<class ArrayType>
int 
FBilinearSurface::
cutWithLine(const ArrayType &x,const ArrayType &s,
		ArrayType solutions[2],bool refine) const
{

  /*

  first, some notation:
    < , > := scalar product, % := cross product    

    fAB := f(a,b) 
    # value of bilinear function with two arguments,
    # first set to a, 2nd set to b

    fA0 := f(a,0) 
    # value of bilinear function with two arguments,
    # first set to a, 2nd set to 0


    fAb := (df/db)(a) 
    # derivation after 2nd parameter
    # with first parameter set to a

    (big letters are parameters, 
     small ones derivations)

    in this function , we do the following:
    first  we solve the equation

    < x - fA0 , s % fAb > = 0 
    for a
    to obtain the first local coordinate a in the surface.
    ( this means that x lies on the surface that is
      parallel to s and the derivation of the surface after b 
      and goes through the point line f(a,0) )

    This expands to:
    < x - f00-a*fa0 , s % ( f0b + a * fab) > = 0 .

    This gives quadradic equation
    
    < x - f00 , s % fb0 > 
    + a * ( < x - f00 , s % fab > + < -fa0 , s % f0b > )
    + a^2 * < - fa0 , s % fab >

    for a.

    then we compute the other two local coordinates ( b and c)    
    by solving < x+ c*s -fA0, s % fAb %fAb > = 0 for c

    and b is equal to 
    | x+c*s-fA0 | / | fAb | = < x+c*s-fA0,fAb>/<fAb,fAb>
    because fAb is parallel to x+c*s-fA0

  */
    
  
  static double xf00[3],sXf0b[3],sXfab[3];
  
  Fd3kreuz(sXf0b,s,f0b);
  Fd3kreuz(sXfab,s,fab);
  
  for(int i=0; i<3; i++)
	xf00[i] = x[i] - f00[i];
	
  double sols[2];
  int nsols;

  //c,b,a = parameters for standard quadratic equation
  // a x^2 + b x + c = 0
  double 
    c= Fd3prod(xf00,sXf0b),
    b= Fd3prod(xf00,sXfab)-Fd3prod(fa0,sXf0b),
    a= -Fd3prod(fa0,sXfab) ;


  //solve quadratic(or linear) equation for first local coord

  if(a==0){
    sols[0]=-c/b;
    if(!finite(sols[0]))
      return 0;
    nsols=1;
  }
  else{
    double d=b*b-4*a*c;

    if(d<0)return 0;      

    d=sqrt(d);
    if((b+d)==b){
      sols[0]=-0.5*b/a;
      nsols=1;
    }
    else{
      if(b<0)d=-d;
      d = -.5*(b+d) ;
      int i = ( (b<0)^(a<0) );
      sols[ !i ] = c/d;
      sols[ i  ] = d/a;
      nsols=2;
    }

  }
      
  double fAb[3],u[3],v[3];
  
    
  //compute 2nd local coords + t-values for all solutions found
  for(int i=0;i<nsols;i++)
    {
      ArrayType &ret=solutions[i];
      //1st local coordinate
      ret[0] = a = sols[i];

	   for(int i=0; i<3; i++)
		fAb[i] = f0b[i] + a*fab[i];
      Fd3kreuz(u,s,fAb); 
      Fd3kreuz(v,u,fAb);   //v = s x fAb x fAb 

	   for(int i=0; i<3; i++)
		u[i] = xf00[i] - a*fa0[i]; //u= x-fA0
		
      //t value for line
      ret[2] = c = - Fd3prod(u,v)/Fd3prod(s,v);


	  for(int i=0; i<3; i++)
		u[i] += c*s[i]; //u=x+s*c-fA0

      //2nd local coordinate
      ret[1] = b = Fd3prod(u,fAb)/Fd3prod(fAb,fAb);	

      if(refine){

	//check for correctness and make newton-raphson if necessary
	
	for(int i=0; i<3; i++)
		u[i] -= b*fAb[i]; // u= x - ( fAB -s*c )
	
	static double 
	  mat[3][3], //jacobian
	  invMat[3][3]; //inverse,transposed jacobian


	double delt = Fd3prod(u,u);
	double odelt = delt;
	double eps = (Fd3prod(x,x)+Fd3prod(f00,f00)) * 1e-29;

	//num.prec must be better than 1e-15.5
	if(delt > eps)
	{
		
		for(int i=0; i<3; i++)
		{
			mat[0][i] = fa0[i] + b * fab[i]; //mat0=df/da(b)
			mat[1][i] = fAb[i];            //mat1=df/db(a)
			mat[2][i] = -s[i];             //mat2=-df/dc
		}
		
	  //compute transposed inverse
	  Fd3kreuz(invMat[0],mat[1],mat[2]);
	  Fd3kreuz(invMat[1],mat[2],mat[0]);
	  Fd3kreuz(invMat[2],mat[0],mat[1]);	

	  double invdenom = 1.0/Fd3prod(invMat[0],mat[0]);
	  
	  do{
	  
	    //ret+= u*invMat
	    ret[0]+=Fd3prod(u,invMat[0])*invdenom;
	    ret[1]+=Fd3prod(u,invMat[1])*invdenom;
	    ret[2]+=Fd3prod(u,invMat[2])*invdenom;
	  
	    //u = x - f(ret)
	    interpolate(v,ret[0],ret[1]);
	    for(int i=0; i<3; i++)
			u[i] = x[i] - ret[2]* mat[2][i] -v[i];
			
	    odelt=delt;
	    delt=Fd3prod(u,u);	  
	  
	  }while(odelt>delt && delt>eps);
	}
      }
    }

  return nsols;
    
}

template<class ArrayType>
int FBilinearSurface::
refineCutWithLine(const ArrayType &x,const ArrayType &s,
		  ArrayType & solution) const
{
  //make newton-iteration
  
  static double 
    u[3], //rest error 
    a[3], //actual solution
    mat[3][3], //jacobian
    invMat[3][3],invDet, //inverse,transposed jacobian
    delt,odelt; //fabs(u)

	for(int i=0; i<3; i++)  
		a[i] = solution[i];
		
  //compute rest
  interpolate(u,a[0],a[1]);
  for(int i=0; i<3; i++)  
	u[i] -=  x[i] + a[2] * s[i];
  
  delt = fabs(u[0])+fabs(u[1])+fabs(u[2]);
  odelt = delt;

  int iters=-1;

  //   printf("%2d:%.0a  a:(% a % a % a)\n",
  // 	 iters,delt,a[0],a[1],a[2]);

  // u= f(a[0],a[1])-(x+a[2]*s);
  //compute jacobian
  deriva(mat[0],a[1]);   //mat0=(du/da[0])(-,b,-)
  derivb(mat[1],a[0]);   //mat1=(du/da[1])(a,-,-)
	for(int i=0; i<3; i++)  
		mat[2][i] = -s[i];  //mat2=(du/da[2])(-,-,-)
	
  //compute transposed inverse
  Fd3kreuz(invMat[0],mat[1],mat[2]);
  Fd3kreuz(invMat[1],mat[2],mat[0]);
  Fd3kreuz(invMat[2],mat[0],mat[1]);	

  invDet = 1.0/Fd3prod(invMat[0],mat[0]);

  //inverse of derivation is computed only once,
  //because we assume we already 
  //have a good solution and the matrix does not
  //change much


  do
  {
     
    for(int i=0; i<3; i++)  
      solution[i] = a[i];
    
    //a -= (du/da)^(-1) u
    a[0] -= Fd3prod(u,invMat[0])*invDet;
    a[1] -= Fd3prod(u,invMat[1])*invDet;
    a[2] -= Fd3prod(u,invMat[2])*invDet;
	  
    //u = f(ret[0],ret[1])-(x+a[2]*s)
    interpolate(u,a[0],a[1]);
	for(int i=0; i<3; i++)  
		u[i] -= x[i] + a[2]* s[i];
    
    odelt = delt;
    delt = fabs(u[0])+fabs(u[1])+fabs(u[2]);	 
    ++iters;
       // printf("%2d:%.0a delta a:(% a % a % a)\n",
       //     iters,delt,a[0]-solution[0],a[1]-solution[1],a[2]-solution[2]);

    if( iters > 1000 )
        break;

  } while(odelt>delt);      

  return iters;
}

/**
 *interpolate position at local coords a,b
 */
template<class ArrayType>
void FBilinearSurface::interpolate(ArrayType &x,double a,double b) const
{
  double ab=a*b;
  for(int i=0; i<3; i++)  
	x[i] = f00[i] + a*fa0[i] + b*f0b[i] + ab*fab[i];
}


/**
 *interpolate derivatives respective 1st local coord
 *\param dxda: derivation after first local coord
 *\param b: 2nd local coord
 */


template<class ArrayType>
void FBilinearSurface::deriva(ArrayType &dxda,double b) const
{
  for(int i=0; i<3; i++) 
	dxda[i] = fa0[i] + b*fab[i];
}

/**
 *interpolate derivatives respective 2nd local coord
 *\param dxdb: derivation after second local coord
 *\param a: 1st local coord
 */
template<class ArrayType>
void FBilinearSurface::derivb(ArrayType &dxdb,double a) const
{
   for(int i=0; i<3; i++) 
	dxdb[i] = f0b[i] + a*fab[i];
}

/**
 *interpolate derivatives respective 2nd local coord
 *\param dxdb: derivation after second local coord
 *\param a: 1st local coord
 */
template<class ArrayType>
void FBilinearSurface::normal(ArrayType &n,double a,double b) const
{
  double dxda[3],dxdb[3];
  
  for(int i=0; i<3; i++)
  {
	dxda[i] = fa0[i] + b*fab[i];
	dxdb[i] = f0b[i] + a*fab[i];
  }
  Fd3kreuz(n,dxda,dxdb);
}

#endif //FBILINEARSURFACE_HH
