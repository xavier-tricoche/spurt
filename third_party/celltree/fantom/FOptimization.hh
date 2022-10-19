#ifndef FOptimization_hh
#define FOptimization_hh

#include "FMathValues.hh"
#include "FArray.hh"
#include "FException.hh"
#include "FMatrix.hh"
#include "eassert.hh"

//#include <iostream>
namespace FMath
{
  /**
   * Methods for Optimization (i.e. Maximization or minimization of function values)
   */
  namespace FOptimization
  {
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

    /**
     * given function func
     * initial points ax,bx. Searches in downhill direction (defined by the function as evaluated at the initial points) and returns new points ax, bx, cs that bracket a minimum of the function. Also returned are the function values at the three points fa, fb and fc
     *
     *
     * NR pp 400f
     */
    template<typename FUNCTION> void mnbrak(double &ax, double &bx, double&cx, double&fa, double&fb, double&fc, FUNCTION &func)
    {
#ifndef NODEBUG
      try{
#endif
      double ulim, u,r,q,fu,dum;
      fa=(func)(ax);
      fb=(func)(bx);
      if(fb > fa)
      {
        SHFT(dum, ax, bx, dum);
        SHFT(dum,fb,fa,dum);
      }
      cx=bx+GOLD*(bx-ax);
      fc=(func)(cx);
      while(fb>fc)
      {
        r=(bx-ax)*(fb-fc);
        q=(bx-cx)*(fb-fa);
        u=bx-(bx-cx)*q-(bx-ax)*r/(2.0*F_SIGN2(F_MAX_DBL(fabs(q-r),TINY),q-r));
        ulim=bx+GLIMIT*(cx-bx);
        // text various possibilities
        if((bx-u)*(u-cx) >0.0) // parabolic u is between b and c, try it
        {
          fu=(func)(u);
          if(fu < fc)
          {
            ax = bx;
            bx = u;
            fa = fb;
            fb = fu;
            return;
          }
          else if (fu > fb ) // got a minimum between a and u
          {
            cx=u;
            fc=fu;
            return;
          }
          u=cx+GOLD*(cx-bx); // parabolic fit was not used. Use default magnification
          fu=(func)(u);
        }
        else if((cx-u)*(u-ulim) > 0.0) // parabolic fit is bewween c and its allowed limit
        {
          fu = (func)(u);
          if(fu<fc)
          {
            SHFT(bx,cx,u,cx+GOLD*(cx-bx));
            SHFT(fb,fc,fu,(func)(u));
          }
        }
        else if((u-ulim)*(ulim-cx) >= 0.0)
        {
          u=ulim;
          fu=(func)(u);
        }
        else
        {
          u=cx+GOLD*(cx-bx);
          fu=(func)(u);
        }
        SHFT(ax,bx,cx,u);
        SHFT(fa,fb,fc,fu);
      }
#ifndef NODEBUG
      } CATCH_N_RETHROW(FException);
#endif
    }


#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10

    /**
     *  Given function f and a given bracketing triplet of abscissas ax,bx,cx (such that bx is between ax and cx and f(bx) is less than both f(ax and f(bx). This routine isolates the minimum to a fractional precision of about tol using Brent's method. The abscissa of the minimum is returned as xmin and the minimum function value is returned.
     *  \param ax ...
     *  \param bx ...
     *  \param cx ...
     *  \param xmin (returned) abszissa value of the minimum (position)
     *  \returns minimum function value
     *
     *  NR pp. 404f
     */
    template<typename FUNCTION>double brent(double ax, double bx, double cx, FUNCTION f, double tol, double &xmin)
    {
#ifndef NODEBUG
      try{
#endif
      // HELP:
      // d is uninitialized here, I just set it to e, because it is associated with it, but this may not work
      unsigned int iter;
      double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
      double e=0.0;
      d=e; // FIXME!

      a=(ax<cx? ax:cx);
      b=(ax>cx? ax:cx);
      x=w=v=bx;
      fw=fv=fx=f(x);
      for(iter =1; iter <= ITMAX; ++iter)
      {
//        std::cout << "brend iteration " << iter << std::endl;
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
//        std::cout << "tol2=" << tol2 << " fabs: " << fabs(x-xm) << " ?? " << tol2-0.5*(b-a) << std::endl;
        if(fabs(x-xm) <=(tol2-0.5*(b-a)))
        {
          xmin=x;
          return fx;
        }
        if(fabs(e) > tol1)
        {
          r =(x-w)*(fx-fv);
          q=(x-v)*(fx-fw);
          p=(x-v)*q-(x-w)*r;
          q=2.0*(q-r);
          if(q>0.0) p =-p;
          q = fabs(q);
          etemp=e;
          e=d;
          if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
            d=CGOLD*(e=(x >= xm? a-x : b-x));
          else
          {
            d=p/q;
            u=x+d;
            if(u-a < tol2 || b-u < tol2)
              d=F_SIGN2(tol1,xm-x);
          }
        }
        else
        {
          d=CGOLD*(e=(x>= xm? a-x:b-x));
        }
        u=(fabs(d) >= tol1 ? x+d : x+F_SIGN2(tol1,d));
        fu=(f)(u); // single function evaluation per iteration
        if(fu <=fx)
        {
          if(u >= x ) a=x; else b=x;
          SHFT(v,w,x,u);
          SHFT(fv,fw,fx,fu);
        }
        else
        {
          if(u<x) a=u; else b=u;
          if(fu <=fw || w==x )
          {
            v=w;
            w=u;
            fv=fw;
            fw=fu;
          }
          else if(fu <= fv || v == x || v == w )
          {
            v=u;
            fv=fu;
          }
        }
      }
      THROW_DEFAULT_EXCEPTION( FTooManyIterationsException );
      //xmin =x;
      //return fx;
#ifndef NODEBUG
      }CATCH_N_RETHROW(FException);
#endif
    }
#define TOL 2.0e-4
#define EPS 1.0e-10
#define MOV3(a,b,c,   d,e,f) (a)=(d); (b)=(e); (c)=(f);
    /**
     * to a given function f and its derivative df and a given bracketing triplet of abscissas ax, bx, cx ( such that bx is between ax and cx and f(bx) is less than both f(ax) and f(cx) this routine isolates the minimum to a fractional precision of about tol using a modification of Brent's method that uses derivatives. Th abscissa of the minimum is returned as xmin. The function value at its minimum is returned.
     * \param ax
     * \param bx
     * \param cx interval
     * \param f function
     * \param df derivative of function
     * \param tol tolerance
     * \param xmin value at minimum
     */
    template<typename FUNCTION, typename DERIV>
      double dbrent(double ax, double bx, double cx, FUNCTION f, DERIV df, double tol, double *xmin)
      {
#ifndef NODEBUG
        try{
#endif
        unsigned int iter, ok1, ok2;
        double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
        double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

        a=(ax<cx?ax:cx);
        b=(ax>cx?ax:cx);
        x=w=v=bx;
        fw=fv=fx=f(x);
        dw=dv=dx=df(x);
        for(iter=1;iter<=ITMAX;++iter)
        {
          xm=0.5*(a+b);
          tol1=tol*fabs(x)+ZEPS;
          tol2=2.0*tol1;
          if(fabs(x-xm) <=(tol2-0.5*(b-a)))
          {
            *xmin=x;
            return fx;
          }
          if( fabs(e) > tol1)
          {
            d1=2.0*(b-a);
            d2=d1;
            if(dw != dx) d1=(w-x)*dx/(dx-dw); // Warning: check for epsilon?
            if(dv != dx) d2=(v-x)*dx/(dx-dv); // Warning: check for epsilon?
            u1=x+d1;
            u2=x+d2;
            ok1=(a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
            ok2=(a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
            olde=e;
            e=d;
            if(ok1 || ok2) {
              if(ok1 && ok2)
                d=(fabs(d1) < fabs(d2) ? d1:d2);
              else if(ok1)
                d=d1;
              else
                d=d2;
              if(fabs(d) <= fabs(0.5*olde))
              {
                u=x+d;
                if(u-a<tol2 || b-u < tol2)
                  d=F_SIGN2(tol1,xm-x);
              }
              else
              {
                d=0.5*(e=(dx >= 0.0? a-x: b-x));
              }
            }
            else
            {
              d=0.5*(e=(dx >= 0.0 ? a-x: b-x));
            }
          }
          else
          {
            d=0.5*(e=(dx >= 0.0 ? a-x: b-x));
          }

          if(fabs(d) >= tol1)
          {
            u=x+d;
            fu=f(u);
          }
          else
          {
            u=x+F_SIGN2(tol1,d);
            fu=f(u);
            if(fu>fx) // minimum step downhill goes uphill, we are done
            {
              *xmin=x;
              return fx;
            }
          }
          du=df(u);
          if(fu <= fx)
          {
            if(u>= x) a=x; else b=x;
            MOV3(v,fv,dv, w,fw,dw);
            MOV3(w,fw,dw, x,fx,dx);
            MOV3(x,fx,dx, u,fu,du);
          }
          else
          {
            if(u<x) a=u; else b=u;
            if(fu<=fw || w==x)
            {
              MOV3(v,fv,dv, w,fw,dw);
              MOV3(w,fw,dw, u,fu,du);
            }
            else if( fu < fv || v==x || v==w)
            {
              MOV3(v,fv,dv, u,fu,du);
            }
          }
        }
        THROW_DEFAULT_EXCEPTION( FTooManyIterationsException );
#ifndef NODEBUG
        } CATCH_N_RETHROW(FException);
#endif
      }



    // this does not seem to work, maybe a buggy helper function
    // use the one without derivatives below :(
    template< typename FUNCTION>
    class FDerivativesLinMin
    {
      private:
        unsigned int dim;
        const FUNCTION& func;
        const FUNCTION& dfunc;
        FArray pcom, xicom;

        double f1dim(double x) const
        {
          double f;
          FArray xt(dim);

          for(unsigned j=0; j<dim; ++j) xt[j]=pcom[j]+x*xicom[j];
          f=func(xt);
          return f;
        }

        double df1dim(double x) const
        {
          //unsigned int j;
          double df1=0.0;
          //FArray xt(dim);
          FArray df(dim);
          //for(j=0; j<dim; ++j) xt[j]=pcom[j]+x*xicom[j];
          FArray xt = pcom+xicom*x;
          dfunc(df, xt);
          //for(j=0; j<dim; ++j) df1 += df[j]*xicom[j];
          df1 = df*xicom; // scalar product
          return df1;
        }

        struct Callf1dim{
          const FDerivativesLinMin &linmin;
          Callf1dim(FDerivativesLinMin& linmin) : linmin(linmin){}
          double operator()(double x) const { return linmin.f1dim(x); }
        };
        struct Calldf1dim{
          const FDerivativesLinMin &linmin;
          Calldf1dim(FDerivativesLinMin& linmin) : linmin(linmin){}
          double operator()(double x) const { return linmin.df1dim(x); }
        };


      public:
        FDerivativesLinMin(unsigned int dim, FUNCTION &func, FUNCTION &dfunc) : dim(dim), func(func), dfunc(dfunc)
        {}
        void linmin(FArray &p, FArray xi, double *fret)
        {
          Callf1dim callf1dim(*this);
          Calldf1dim calldf1dim(*this);

          unsigned int j;
          double xx,xmin,fx,fb,fa,bx,ax;

          pcom=FArray(dim);
          xicom=FArray(dim);
          for(j=0; j<dim; ++j)
          {
            pcom[j]=p[j];
            xicom[j]=xi[j];
          }
          ax=0.0;
          xx=1.0;
          mnbrak(ax,xx,bx,fa,fx,fb,callf1dim);

          *fret=dbrent(ax,xx,bx,callf1dim,calldf1dim, TOL, &xmin);
          for(j=0; j<dim; ++j)
          {
            xi[j] *= xmin;
            p[j] += xi[j];
          }
        }
    };

    template< typename FUNCTION >
    class FLinMin
    {
      private:
        unsigned int dim;
        FUNCTION &func;

        // tmp values needed by f1dim
        FArray pcom, xicom;



        double f1dim(double x) const
        {
          //FArray xt(dim);
          //for(unsigned j=0; j<dim; ++j) xt[j]=pcom[j]+x*xicom[j];
          FArray xt( pcom+xicom*x );
          return func(xt);
        }

        struct Callf1dim{
          const FLinMin &linmin;
          Callf1dim(FLinMin& linmin) : linmin(linmin){}
          double operator()(double x) const { return linmin.f1dim(x); }
        };

      public:
        FLinMin(unsigned int dim, FUNCTION &func) : dim(dim), func(func)
        {
        }

        void linmin( FArray&p, FArray &xi, double *fret)
        {
#ifndef NODEBUG
          try{
#endif
          Callf1dim callf1dim(*this);
          //unsigned int j;
          double xx, xmin, fx, fb, fa, bx, ax;

//          pcom = FArray(dim);
//          xicom= FArray(dim);
/*
          for(j=0; j<dim;++j)
          {
            pcom[j]=p[j];
            xicom[j]=xi[j];
          }
          */
          pcom = p;
          xicom = xi;
          ax=0.0;
          xx=1.0;
          mnbrak(ax,xx,bx,fa,fx,fb,callf1dim);
//          std::cout << "ax" << ax << "xx" << xx << "bx " << bx << " fa " << fa << " fx " << fx << " fb " << fb << std::endl;
          *fret=brent(ax,xx,bx,callf1dim,TOL,xmin);
          /*
          for(j=0;j<dim;++j)
          {
            xi[j]*=xmin;
            p[j]+=xi[j];
          }
          */
          xi *= xmin;
          p  += xi; // is xi ok? or should this be xicom?
#ifndef NODEBUG
        } CATCH_N_RETHROW(FException);
#endif
        }
   };

    /**
     * conjugate gradient method for multidimensions
     */
    template < typename FUNCTION, typename LINEMIN >
    class FFletcherReeves
    {
      private:
        const unsigned int dim;
      public:

        FFletcherReeves( unsigned int dim ) : dim(dim)
        { }
        
        /**
         * NR pp. 423f (frprmn)
         * Given a starting point p, Fletcher-Reeves minimization is performed on a function func using its gradient dfunc.
         * \param p position of minimum
         * \param ftol convergence tolerance
         * \param iter returns number of iterations performed
         * \param fret returns minimum value of the function
         * \param itmax maximum number of iterations
         */
        void minimize( FArray &p, double ftol, int *iter, double *fret, FUNCTION func, FUNCTION dfunc, LINEMIN &linemin, unsigned int itmax)
        {
#ifndef NODEBUG
          try{
#endif
          eassert(dim == p.getDimension());
//          FLinMin<FUNCTION> linemin( dim, func );
          FArray g(dim);
          FArray h(dim);
          FArray xi(dim);
          double dgg;
          double gg;
          double gam;
          double fp=(func)(xi,p);
          /*
          for( unsigned int j=0; j<dim; ++j)
          {
            g[j] = -xi[j];
            xi[j]=h[j]=g[j];
          }
          */
          g = -xi;
          xi = h = g;
          for(unsigned int its=1; its<=itmax; ++its)
          {
            *iter = its;
            linemin.linmin(p,xi,fret);
            if(2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) return;
            
            fp = *fret;
            (dfunc)(xi,p);
            dgg=gg=0.0;
            for(unsigned int j=0; j<dim;++j)
            {
              gg+=g[j]*g[j];
              dgg+=xi[j]*xi[j];         // Fletcher-Reeves
              //dgg+=(xi[j]+g[j])*xi[j];  // Polak-Ribiere
            }
            if(gg == 0.0) return;

            gam=dgg/gg;
            for(unsigned int j=0; j<dim;++j)
            {
              g[j] = -xi[j];
              xi[j]=h[j]=g[j]+gam*h[j];
            }
          }
#ifndef NODEBUG
          }CATCH_N_RETHROW(FException);
#endif
        }
    };

    /**
     * conjugate gradient method for multidimensions
     */
    template < typename FUNCTION, typename LINEMIN>
    class FPolakRiviere
    {
      private:
        const unsigned int dim;
      public:

        FPolakRiviere( unsigned int dim ) : dim(dim)
        {
        }

        /**
         * NR pp. 423f (frprmn)
         * Given a starting point p, Fletcher-Reeves minimization is performed on a function func using its gradient dfunc.
         * \param p position of minimum
         * \param ftol convergence tolerance
         * \param iter returns number of iterations performed
         * \param fret returns minimum value of the function
         * \param itmax maximum number of iterations
         */
        void minimize( FArray &p, double ftol, int *iter, double *fret, FUNCTION func, FUNCTION dfunc, LINEMIN &linemin, const unsigned int itmax)
        {
#ifndef NODEBUG
          try{
#endif
          eassert(dim == p.getDimension());
          //FLinMin<FUNCTION> linemin( dim, func );
          FArray g(dim);
          FArray h(dim);
          FArray xi(dim);
          double dgg;
          double gg;
          double gam;
          double fp=(func)(xi,p);
          /*
          for( unsigned int j=0; j<dim; ++j)
          {
            g[j] = -xi[j];
            xi[j]=h[j]=g[j];
          }
          */
          g = -xi;
          xi = h = g;
          for(unsigned int its=1; its<=itmax; ++its)
          {
            *iter = its;
            linemin.linmin(p,xi,fret);
            if(2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) return;

            fp = *fret;
            (dfunc)(xi,p);
            dgg=gg=0.0;
            for(unsigned int j=0; j<dim;++j)
            {
              gg+=g[j]*g[j];
              //dgg+=xi[j]*xi[j];         // Fletcher-Reeves
              dgg+=(xi[j]+g[j])*xi[j];  // Polak-Ribiere
            }
            if(gg == 0.0) return;
            
            gam=dgg/gg;
            for(unsigned int j=0; j<dim;++j)
            {
              g[j] = -xi[j];
              xi[j]=h[j]=g[j]+gam*h[j];
            }
          }
#ifndef NODEBUG
          }CATCH_N_RETHROW(FException);
#endif
        }
    };

#define SQR(a) ((a)*(a))
    /**
     * minimization of function func.
     * Similar to those above
     */
    template< typename FUNCTION, typename LINEMIN >
      class FPowell
      {
        unsigned int dim;
        FMatrix xi; //< column vectors contain the initial set of directions

        public:
        FPowell(unsigned int dim) : dim(dim) 
        {
          // initial Directions
          xi.makeDiagonal(1.0,dim);
        }

        // only needed if initial directions should be modified
        FMatrix& getMatrix()
        {
          return xi;
        }

        void minimize( FArray &p, double ftol, int *iter, double *fret, FUNCTION func, LINEMIN &linemin/*, const unsigned int itmax*/)
        {
          unsigned int i, ibig, j;
          double del, fp, fptt, t;
          FArray pt(dim), ptt(dim), xit(dim);

          *fret = func(p);
          pt = p;
          for(*iter=1;;++(*iter))
          {
            fp=(*fret);
            ibig=0;
            del = 0.0;
            for(i=0;i<dim;++i)
            {
              for(j=0;j<dim;++j)
                xit[j]=xi(j,i);
              fptt=(*fret);
              linemin.linmin(p,xit,fret);
              if(fptt-(*fret) > del)
              {
                del=fptt-(*fret);
                ibig=i;
              }
            }
            if( 2.0*(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))+TINY)
            {
              return;
            }
            if(*iter == ITMAX)
              THROW_DEFAULT_EXCEPTION( FTooManyIterationsException );
            for(j=0; j<dim;++j)
            {
              ptt[j]=2.0*p[j]-pt[j];
              xit[j]=p[j]-pt[j];
              pt[j]=p[j];
            }
            fptt=func(ptt);
            if(fptt < fp)
            {
              t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
              if(t<0.0)
              {
                linemin.linmin(p,xit,fret);
                for(j=0;j<dim;++j)
                {
                  xi(j,ibig)=xi(j,dim-1);
                  xi(j,dim-1)=xit[j];
                }
              }
            }
          }
        }
        
      };
  } // FOptimiziation
} // FMath
#endif
