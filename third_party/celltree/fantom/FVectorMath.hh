#ifndef FVectormath_hh
#define FVectormath_hh

#include "FVector.hh"
#include "FPosition.hh"

#include <cmath>

#ifndef NODEBUG
   #include "eassert.hh"
#endif

#include "FMathValues.hh"
#include <iostream>
#include <vector>


namespace FMath
{
	#define MY_EPSILON (1e-9)
	
	/**
	*	Distance between 2 Points
	*/
	inline double pointPointDistance(const FPosition& x, const FPosition& y)
	{
		return x.distance(y);
	}
	
	/**
	*	Minimal distance between point  and line segment
	*	p0,p1 : 	end - points of segment
	*	p:		point
	*/
	inline double pointSegmentDistance(const FPosition& p, const FPosition& p0, const FPosition& p1)
	{
		//segment vector
		FVector v = p1 - p0;
		FVector w = p - p0;

		//if ( c1 <= 0 ) then the p0 point is the closest to p
	    double c1 = w*v;
	    if ( c1 <= 0 )
			return p.distance(p0);

		//if ( c2 <= c1 ) then the p1 point is   closest to p
	    double c2 = v*v;
	    if ( c2 <= c1 )
	        return p.distance(p1);

	    double b = c1 / c2;
	    FPosition Pb = p0 + b * v;
	    return p.distance(Pb);
	}
   
   
	/**
	*	minimal distance between a point p and a plane given by a normal and a point
	*	distance = n * (p-e0)    where p the point, e0 is point of the plane, n is not unit normal
	*	
	*/
	inline double pointPlaneDistance(const FPosition& p, //point
									const FPosition& e0,  const FVector& n)   // point and normal vector of a plane 
	{
		FVector nn = n.normalized();
		FVector u = p-e0;
		return nn*u;
	}   
   
	/**
	*	minimal distance between a point and a plane
	*	d= n * (p-e0)    where n is unit normal vector, p the point, e0 is one of the tree given points of the plane
	*/
	inline double pointPlaneDistance(const FPosition& p, //point
									const FPosition& e0,  const FPosition& e1, const FPosition& e2)   // 3 points in space to define a plane
	{
		FVector n, v1, v2;
	    
	    v1 = e1-e0;
	    v2 = e2-e0;
	    
	    n = v1.crossProduct( v2);
		return pointPlaneDistance(p, e0, n);
	}
	
	
	
	/**
	*	minimal distance between a line segment and a plane
	*/
	inline double segmentPlaneDistance( const FPosition& e0,  const FPosition& e1, const FPosition& e2,   // 3 points in space to define a plane
									  const FPosition& p1, const FPosition& p2 )                        // 2 points in space to define a line
	{
	    FVector n, m, v1, v2;
	    
	    v1 = e1-e0;
	    v2 = e2-e0;
	    
	    n = v1.crossProduct( v2);
	    n.normalize();              // unit normal vector of the plane
	  
	    m = p1 - p2;  // calculate the direction vector of the line (this vector points from p2 to p1)
	    
	    if(fabs(n*m) < MY_EPSILON)         // plane and line are parallel
	    {   
			//calc distance here
			FVector u = p1-e0;
		
		    return  n*u;
	    }
		return 0.0;
	}

	/**
	*	minimal distance between 2 planes
	*/
	inline double planePlaneDistance( const FPosition& e0,  const FPosition& e1, const FPosition& e2,   // 3 points in space to define a plane
									const FPosition& p0,  const FPosition& p1, const FPosition& p2)   // 3 points in space to define a plane
	{
		FVector n, m, v1, v2, u1, u2;
	    
		//get vectors of the plane1 
	    v1 = e1-e0;
	    v2 = e2-e0;
		
		//get normal of the plane1
	    n = v1.crossProduct( v2);
	    n.normalize(); 
		 
		//get vectors of the plane2 
		u1 = p1-p0;
	    u2 = p2-p0;
		
	    //get normal of the plane2
	    m = u1.crossProduct( u2);
	    m.normalize();   
	 
		//not sure check this
	    if(fabs(n*m) < MY_EPSILON)         // plane1 and plane2 are parallel
	    {   
			//calc distance here
			FVector w = p0-e0;
		
			return n*w;
	    }
		return 0.0;
	}
	
	
  /** 
   * Minimal distance between two line segments defined by
   * (a0,a1) and (b0,b1). 
   */
  inline double segmentDistance(FPosition a0, FPosition a1,FPosition b0, FPosition b1)
  {
    double distance=0;
    
    FVector a = a1 - a0;
    FVector b = b1 - b0;
    FVector c = crossProduct(a,b);
    
#ifndef NODEBUG
    eassert( c != FVector(0.,0.,0.) );
#endif
    
    double param_a, param_b;
    FVector ab=( a0 - b0 );
    {
      FVector n = crossProduct(a,c);

      double l = n * ab;
      double b_projected = n * b;
      
      param_b = l/b_projected;
      param_b = param_b > 0 ? param_b : 0;
      param_b = param_b < 1 ? param_b : 1;
    }
    
    {
      FVector n = crossProduct(b,c);
      
      double l = -(n * ab);
      double a_projected = n * a;
      
      param_a = l/a_projected;
      param_a = param_a > 0 ? param_a : 0;
      param_a = param_a < 1 ? param_a : 1;
    }
    
    distance = ((a0+param_a*a)-(b0+param_b*b)).norm();
    
    return distance; 
  }

    
  /** 
   * This function tests for intersection of triangle (x,y,z)
   * and segment (p0,p1). It returns true if there is an intersection and 
   * sets t_parameter to the position of intersection on the segment.
   */
  inline bool intersectTriangleWithSegment( double& t_parameter, 
				            const FPosition& x,  const FPosition& y, const FPosition& z, 
				            const FPosition& p0, const FPosition& p1 )
  {
    FVector tvec, pvec, qvec;
    double det, inv_det, t;
    
    double u, v;
    
    FVector e1 = y - x;
    FVector e2 = z - x;
    pvec = (p1-p0).crossProduct(e2);
    det = e1 * pvec;
    
    if( fabs(det) < 10e-9 )      // singular case
      return false;
    
    inv_det = 1.0 / det;
    
    tvec = p0 - x;
    u = (tvec * pvec) * inv_det;
    
    if(u < 0.0 || u > 1.0)
      return false;
    
    qvec = tvec.crossProduct(e1);
    v = ((p1-p0) * qvec) * inv_det;
    
    if(v < 0.0 || u + v > 1.0)
      return false;
    
    t = (e2 * qvec) * inv_det;

    t_parameter = t;

    if((0 <= t) && (t < 1))
    {      
      return true;
    }
    
    return false;
  }

  
  /** 
   * This function is a conveniece function for intersectTriangleWithSegment( double& t_parameter, FPosition x, FPosition y, FPosition z, FPosition p0, FPosition p1 )
   * that allows only to test and not to get the t_parameter.
   */
  inline bool intersectTriangleWithSegment(const FPosition& x, const FPosition& y, const FPosition& z, const FPosition& p0, const FPosition& p1 )
  {
    double dummy=0;
    return intersectTriangleWithSegment( dummy,  x,  y, z,  p0,  p1 );
  }


  /**
   * \brief Calculates the point where line and plane intersect.
   * This method is doing what the name suggests:
   * It returns the point in space where the given line is intersecting the plane.
   * If the plane and the line are parallel the interecting point is some point in the plane (e0),
   * but the return value of the function will be "false" instead of true in the other case.
   */
  inline bool intersectPlaneWithLine( FPosition &cutpoint,                                              // the intersecting point
				      const FPosition& e0,  const FPosition& e1, const FPosition& e2,   // 3 points in space to define a plane
				      const FPosition& p1, const FPosition& p2 )                        // 2 points in space to define a line
  {
    FVector n, m, v1, v2;
    double t, d;
    
    v1 = e1-e0;
    v2 = e2-e0;
    
    n = v1.crossProduct( v2);
    n.normalize();              // normal vector of the plane
    
    d = n * e0;
  
    m = p1 - p2;  // calculate the direction vector of the line (this vector points from p2 to p1)
    
    if(fabs(n*m) < MY_EPSILON)         // plane and line are parallel
    {   cutpoint = e0;   // otherwise it would be undefined
        return false;
    }
    
    t = (d - n*p2) / (n*m);
    
    cutpoint = p2 + t*m;
    
    return true;
  }  

  
/**
 * \brief Method to project a bunch of coordniates to a plane.
 * Method to project a bunch of coordinates to a plane specified by its normal vector "normal". The normal vector is
 * considered to be at the position p2.
 * "coordinates" containes the coordinates to be projected and "projcoords" are the projected coordinates
 */

inline void projectCoordsToPlaneBatch( const std::vector<FPosition> &coordinates, std::vector<FPosition> &projcoords, FVector normal, FPosition p2)
{
    projcoords.clear();
    projcoords.resize(coordinates.size());
    
    FVector p;
    normal.normalize(); 
    
    for(unsigned int i=0; i<projcoords.size(); i++)
    {   p = coordinates[i] - p2;
        projcoords[i] = coordinates[i] - ( normal*p) * normal;  // actually simple but cool stuff
    }   
}


/**
 * \brief Method to project one coordniate to a plane.
 * Method to project one coordinate to a plane specified by its normal vector "normal". The normal vector is
 * considered to be at the position p2.
 * "coordinate" containes the coordinate to be projected and "projcoord" is the projected coordinate
 */
inline void projectCoordsToPlane( const FPosition &coordinate, FPosition &projcoord, FVector normal, FPosition p2)
{
    projcoord.clear();
    projcoord.resize(coordinate.size());
    
    FVector p;
    normal.normalize(); 
    
    p = coordinate - p2;
    projcoord = coordinate - ( normal*p) * normal;  // actually simple but cool stuff
       
}

    
  /**
   *  \author Dominic Schneider
   *  \brief Returns barycentric coordinates for a given point and a triangle.
   *  This method returns the barycentric coordinates given three points to define the triangle
   *  and one point the barycentric coordinates should be calculated for.
   */  
  inline bool barycentricCoordinates(double &lambda1, double &lambda2, double &lambda3, FPosition p, FPosition p1, FPosition p2, FPosition p3 )
  {
      double A = p1[0] - p3[0];
      double B = p2[0] - p3[0];
      double C = p3[0] -  p[0];
      double D = p1[1] - p3[1];
      double E = p2[1] - p3[1];
      double F = p3[1] -  p[1];
      double G = p1[2] - p3[2];
      double H = p2[2] - p3[2];
      double I = p3[2] -  p[2];

      double N = A*(E+H) - B*(D+G);  // Nenner
      double Q = 0.0;
      double M = 0.0;
      double R = 0.0;
      
      if( N*N < 0.0001)   // bad case, triangle is perpendicular to the X-axis
      {    N = D*(B+H) - E*(A+G);   // alternative definition switch A with D, B with E and C with F
           Q = E*(C+I) - F*(B+H);   // taken from wikipedia.org english version
	    
	   M = E*(A+G) - D*(B+H);
	   R = D*(C+I) - F*(A+G);
      }
      else  
      {    Q = B*(F+I) - C*(E+H);
           
	   M = B*(D+G) - A*(E+H);
	   R = A*(F+I) - C*(D+G);      
      }

      if((N==0.0) ||  (M==0.0)) return false;  // bad stuff happend
                  
      lambda1 = Q / N;     
      lambda2 = R / M;
      lambda3 = 1.0 - lambda1 - lambda2;
     
      return true;
  }      

  
/**
 * \author Dominic Schneider
 * \brief Method to determine whether a point is located within a triangle or not.
 * Algorithm to determine whether a point is located within the specified triangle or not using barycentric coordinates
 * The provided coordinates need to be in the plane otherwise one will get wrong results
 */
inline bool insideTriangle3D(FPosition p, FPosition p1, FPosition p2, FPosition p3)
{
    double lambda1=0.0, lambda2=0.0, lambda3=0.0;
    
  
    if(!barycentricCoordinates( lambda1, lambda2, lambda3, p, p1, p2, p3)) 
    {   
        //cout << " either barycentric rubbish or wiered coordinates" << endl;
        return false;
    }

    FVector v1, v2, n;
    double d;
    
    v1 = p2 - p1;
    v2 = p3 - p1;
    
    n = v2.crossProduct(v1);
    n.normalize();
    d = n*p1;        
    
    if( (p*n -d)*(p*n -d) > 0.01)     // the point p is NOT located within the plane of the triangle
    {   //cout << "Arrrg, point is not in plane" << endl;  // due to numeric errors it cant be exactly 0.0
        return false;
    }
        
    double threshold = -0.00001;  // due to numeric errors it has been prooven necessary to impose a threshold
    if((lambda1 < threshold) || (lambda2 < threshold) || (lambda3 < threshold)) return false;
    
    return true;    
}  
    
};  
#endif
