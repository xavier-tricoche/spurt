#ifndef FD3OP_HH
#define FD3OP_HH

/*
 *checks if p is in the bbox define by max, and min
 */

#define Fd3IsInside(max,min,p)\
 (min[0]<=p[0]&&p[0]<=max[0]  \
  &&min[1]<=p[1]&&p[1]<=max[1] \
  &&min[2]<=p[2]&&p[2]<=max[2])

/*
 * a is set to the crossproduct of b and c
 */
#define Fd3kreuz(a,b,c)\
a[0]=b[1]*c[2]-b[2]*c[1];\
a[1]=b[2]*c[0]-b[0]*c[2];\
a[2]=b[0]*c[1]-b[1]*c[0];

/*
  scalar product
 */

#define Fd3prod(a,b)\
(a[0]*b[0]+a[1]*b[1]+a[2]*b[2])



#define Fd3out(p)\
'('<<p[0]<<' '<<p[1]<<' '<<p[2]<<')'

/* 
 *returns norm of x
*/
#define Fd3norm(x)\
sqrt(Fd3prod(x,x))

/*
  Fd3opI
  makes 3 lines out of the expression with
  the commas replaced by [0] in the first line
  the commas replaced by [1] in the 2nd line
  the commas replaced by [2] in the 3rd line
  with exactly I commas
*/

#endif
