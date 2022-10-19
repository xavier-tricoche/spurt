//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTrilinearInterpolation.hh,v $
// Language:  C++
// Date:      $Date: 2002/02/06 19:57:51 $
// Author:    $Author: mlangbe $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef FTRILINEAR_INTERPOLATION_H
#define FTRILINEAR_INTERPOLATION_H

#include<vector>
#include "FArray.hh"

using namespace std;

class FTrilinearInterpolation
{

public:

  /*
   *\par Description:
   * builds the parameters needed
   * for fast trilinear interpolation
   *\pre
   * p describes the vertices values
   * with counterclockwise-enumeration
   * where x describes the trilinear coordinates
   *\post
   * p describes the interpolant parameters
   *\param p: array of objects
   * with operators += , -= , *  defined
   */

  template<class T> 
  static void buildParameters(T p[8])
  {
  
    //old numering of hexeders (bitwise)
    //     p[7]-=p[6]; 
    
    //     p[6]-=p[4]; p[5]-=p[4]; p[3]-=p[2];
    
    //     p[1]-=p[0]; p[2]-=p[0]; p[4]-=p[0];
    
    //     p[3]-=p[1];
    
    //     p[7]-=p[5]; p[7]-=p[3];
    
    //     p[6]-=p[2]; p[5]-=p[1];
    
    //new numbering of hexaeders (counter clockwise)
    p[6]-=p[7]; 
    
    p[7]-=p[4]; p[5]-=p[4]; p[2]-=p[3];
    
    p[1]-=p[0]; p[3]-=p[0]; p[4]-=p[0];

    p[2]-=p[1];

    p[6]-=p[5]; p[6]-=p[2];

    p[7]-=p[3]; p[5]-=p[1];
  
  
  }

  //-------------------------------------------------------------------

  /*
   *\par Description:
   *calculate value at given local coordinates
   *\pre
   *bulidParameters(p) has been invoked
   *\post
   *none
   *\param p: interpolation parameters
   *\param x: local coordinates
   *\return: value at x
   */
  template<class T> static T interpolateT(const T p[8],const double x[3])
  {

    //old numering of hexeders (bitwise)
    //     return 
    //       p[0] 
    //       + x[0]*p[1] + x[1]*p[2] + x[2]*p[4]
    //       + (x[0]*x[1]) * p[3]  +  (x[0]*x[2]) * p[5]  +  (x[1]*x[2]) * p[6]
    //       + (x[0]*x[1]*x[2]) * p[7];
    //new numbering of hexaeders (counter clockwise)
    return 
      p[0] 
      + x[0]*p[1] + x[1]*p[3] + x[2]*p[4]
      + (x[0]*x[1]) * p[2]  +  (x[0]*x[2]) * p[5]  +  (x[1]*x[2]) * p[7]
      + (x[0]*x[1]*x[2]) * p[6];
  }

  //-------------------------------------------------------------------

  /*
   *\par Description:
   *calculate derivative value resp. first local coordinate
   *at given local coordinates
   *\pre
   *bulidParameters(p) has been invoked
   *\post
   *none
   *\param p: interpolation parameters
   *\param x: local coordinates
   *\return: value at x
   */
  template<class T> static T interpolatedTdx(const T p[8],const double x[3])
  {
    return  p[1] + x[1] * p[2] + x[2] * p[5] + (x[1]*x[2]) * p[6];
  }

  //-------------------------------------------------------------------

  /*
   *\par Description:
   *calculate derivative value resp. 2nd local coordinate
   *at given local coordinates
   *\pre
   *bulidParameters(p) has been invoked
   *\post
   *none
   *\param p: interpolation parameters
   *\param x: local coordinates
   *\return: value at x
   */
  template<class T> static T interpolatedTdy(const T p[8],const double x[3])
  {
    return  p[3] + x[0] * p[2] + x[2] * p[7] + (x[0]*x[2]) * p[6];
  }

  //-------------------------------------------------------------------

  /*
   *\par Description:
   *calculate derivative value resp. 3rd local coordinate
   *at given local coordinates
   *\pre
   *bulidParameters(p) has been invoked
   *\post
   *none
   *\param p: interpolation parameters
   *\param x: local coordinates
   *\return: value at x
   */
  template<class T> static T interpolatedTdz(const T p[8],const double x[3])
  {
    return  p[4] + x[0] * p[5] + x[1] * p[7] + (x[0]*x[1]) * p[6];
  }

  //--------------------------------------------------------------------
  
  //optimized function versions I(max) use in hexahedroncell:

  /*
   *\par Description:
   *interpolation functions like above
   *\param p:
   *pointer on array of interpolation parameters for each
   *interpolated value(that means of type double[][8])
   *\param end:
   *pointer after the end of the array
   *\param x:
   *pointer on trilinear coordinates
   *\param ret:
   *pointer on array of returned values
   */

  static void buildParameters(double * p, const double * end);

  static void interpolateT(const double * p, const double * end,
		    const double*x, double * ret);

  static void interpolatedTdx(const double * p, const double * end,
		       const double*x, double * ret);

  static void interpolatedTdy(const double * p, const double * end,
		       const double*x, double * ret);

  static void interpolatedTdz(const double * p, const double * end,
		       const double*x, double * ret);




  void buildParameters(vector<double>&p);

  void interpolateT(const FVector&x, FArray&ret) const;

  void interpolatedTdx(const FVector&x, FArray&ret) const;
  void interpolatedTdy(const FVector&x, FArray&ret) const;
  void interpolatedTdz(const FVector&x, FArray&ret) const;


private:

  vector<double> params;
  
};


#endif //FTrilinearInterpolation_HH
