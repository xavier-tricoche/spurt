#ifndef FBILINEARHALFSPACE_HH
#define FBILINEARHALFSPACE_HH

#include "FBilinearSurface.hh"


/**
 * class which describes a halfspace 
 * bounded by a bilinear surface
 */

class FBilinearHalfspace:public FBilinearSurface{
  
public:
  /**
   * constructor. ArrayType is a Type who has
   * an operator[] defined which returns a scalar.
   * \param F00,F01,F11,F10: 
   * Fij is the value of the bilinear function
   * with the parameters i=I, j=J
   * 
  */
  template<class ArrayType>
  FBilinearHalfspace(const ArrayType&F00,const ArrayType&F10,
		     const ArrayType&F11,const ArrayType&F01);

  /**
   * describes if position is inside or outside the halfspace:
   * \param p
   * Position to check
   * \return
   *  if <0: point outside
   *  if >0: point inside
   *  if =0: point at surface
   */
  template<class ArrayType>
  double isInside(const ArrayType&p);

private:

  double fabXf0b[3],fabXfa0[3];
  double ca,cb,cab;
  
};

//----------------------------------------------------------

template<class ArrayType>
FBilinearHalfspace::
FBilinearHalfspace(const ArrayType&F00,const ArrayType&F10,
		   const ArrayType&F11,const ArrayType&F01)
  :FBilinearSurface(F00,F10,F11,F01)
{
  static double da,db;
    
  Fd3kreuz(fabXf0b,fab,f0b);
  Fd3kreuz(fabXfa0,fab,fa0);
  
  da = Fd3prod(fa0,fabXf0b);
  ca = Fd3prod(fab,fa0)/da;

  db = Fd3prod(f0b,fabXfa0);
  cb = Fd3prod(fab,f0b)/db;
  
  cab = Fd3prod(fab,fab)/(da*db);
  
}

//----------------------------------------------------------

template<class ArrayType>
double
FBilinearHalfspace::
isInside(const ArrayType&p)
{
  static double x[3],a,b;
  Fd3op3(x,=p,-f00,);
  a = Fd3prod(x,fabXf0b);
  b = Fd3prod(x,fabXfa0);
  
  return 
    Fd3prod(x,fab) - ca*a - cb*b - cab* a*b;     
}


#endif //FBILINEARHALFSPACE_HH
