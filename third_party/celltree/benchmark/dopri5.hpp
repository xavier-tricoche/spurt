#ifndef __dopri5_hpp
#define __dopri5_hpp

#include <cmath>
#include <cassert>
#include <iostream>
#include <limits>

namespace 
{

const double safe = 0.9;
const double epsilon = std::numeric_limits<double>::epsilon();
const double facl = 0.2;
const double facr = 10.0;
const double beta = 0.04;
const unsigned int nstiff = 100;

const double a21=0.2, a31=3.0/40.0, a32=9.0/40.0, a41=44.0/45.0,
    a42=-56.0/15.0, a43=32.0/9.0, a51=19372.0/6561.0, a52=-25360.0/2187.0,
    a53=64448.0/6561.0, a54=-212.0/729.0, a61=9017.0/3168.0, 
    a62=-355.0/33.0, a63=46732.0/5247.0, a64=49.0/176.0, a65=-5103.0/18656.0,
    a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0, a75=-2187.0/6784.0, a76=11.0/84.0;

const double c2=0.2, c3=0.3, c4=0.8, c5=8.0/9.0;

const double d1=-12715105075.0/11282082432.0, d3=87487479700.0/32700410799.0,
    d4=-10690763975.0/1880347072.0, d5=701980252875.0/199316789632.0,
    d6=-1453857185.0/822651844.0, d7=69997945.0/29380423.0;

const double e1=71.0/57600.0, e3=-71.0/16695.0, e4=71.0/1920.0, 
    e5=-17253.0/339200.0, e6=22.0/525.0, e7=-1.0/40.0;

}

// -------------------------------------------------------------------------

template<typename V, unsigned int D>
struct dopri5
{
    typedef V value_type;
    
    enum result
    {
        OK,
        T_MAX_REACHED,
        OUTSIDE_DOMAIN,
        STEPSIZE_UNDERFLOW,
        STIFFNESS_DETECTED,
    };
    
    value_type reltol;
    value_type abstol;
    
    value_type h;
    value_type h_max;
    value_type h_fixed;
    // value_type h_max;
    
    unsigned int n_accepted;
    unsigned int n_rejected;
    unsigned int n_steps;
    unsigned int n_eval;
    
    value_type     t, t_max;
    value_type y[D], k1[D];

    // stepsize control stabilization
    value_type     facold;

    // stiffness detection
    value_type     hlamb;
    int            iasti;
    int            nonsti;
    bool           adaptive;
    
public:

    // struct step
    // {
    //     value_type p[5];
    //     value_type y0, y1;
    //     double     t0, t1;

    //     value_type y( const double& t ) const
    //     {
    //         double a = (t - t0) / (t1 - t0);
    //         double b = 1.0 - a;
        
    //         return p[0]+a*(p[1]+b*(p[2]+a*(p[3]+b*p[4])));
    //     }

    //     double y( const double& t, unsigned int n ) const
    //     {
    //         double a = (t - t0) / (t1 - t0);
    //         double b = 1.0 - a;
        
    //         return p[0][n]+a*(p[1][n]+b*(p[2][n]+a*(p[3][n]+b*p[4][n])));
    //     }

    //     value_type dy( const double& t ) const
    //     {
    //         double a = (t - t0) / (t1 - t0);
    //         double b = 1.0 - a;

    //         value_type pc = p[3]+b*p[4];

    //         return ( p[1] + b*(p[2]+a*pc) -
    //                  a*(p[2] + a*pc -  b*(p[3]+(b-a)*p[4]) ) ) / (t1-t0);
    //     }

    //     double dy( const double& t, unsigned int n ) const
    //     {
    //         double a = (t - t0) / (t1 - t0);
    //         double b = 1.0 - a;

    //         double pc = p[3][n]+b*p[4][n];

    //         return ( p[1][n] + b*(p[2][n]+a*pc) -
    //                  a*(p[2][n] + a*pc -  b*(p[3][n]+(b-a)*p[4][n]) ) ) / (t1-t0);
    //     }
    // };

    dopri5()
    {
        reltol = 1e-6;
        abstol = 1e-6;
        
        h       = 0.0;
        h_max   = 0.0;
        
        n_accepted = 0;
        n_rejected = 0;
        n_steps    = 0;
        n_eval     = 0;
        
        facold = 1e-5;
        hlamb  = 0.0;
        iasti  = 0;
    }
    
    template<typename RHS> result step( const RHS& rhs )
    {
        const double direction = t_max < t ? -1 : 1;

        // compute k1 for FSA and h_init()
        if( n_steps == 0 )
        {
            ++n_eval;
            if( !rhs( t, y, k1 ) )
                return OUTSIDE_DOMAIN;
        }
        
        // if h is zero, estimate a suitable initial stepsize
        if( h == 0.0 )
        {
            value_type dnf = 0, dny = 0, tmp, sy;

            for( unsigned int d=0; d<D; ++d )
            {
                sy = abstol + reltol * std::abs(y[d]);
        
                tmp = k1[d] / sy;
                dnf += tmp*tmp;
                tmp = y[d] / sy;
                dny += tmp*tmp;
            }

            if( (dnf <= 1.0e-10) || (dny <= 1.0e-10) ) 
                h = 1.0e-6;
            else 
                h = std::sqrt( dny/dnf ) * 0.01;
              
            h *= direction;
          
            // explicit Euler step
            value_type y1[D];

            for( unsigned int d=0; d<D; ++d )
                y1[d] = y[d] + h * k1[d];

            value_type k2[D];

            ++n_eval;
            if( !rhs( t+h, y1, k2 ) )
                return OUTSIDE_DOMAIN;

            // estimate the second derivative of the solution
            // step size is computed such that h^(1/5) * max( norm(k1), norm(der2) ) = 0.01
            value_type der2 = 0.0; 
 
            for( unsigned int d=0; d<D; ++d )
            {
                sy = abstol + reltol * std::abs(y[d]);

                tmp = (k2[d] - k1[d]) / sy;
                der2 += tmp*tmp;
            }

            der2 = std::sqrt( der2 ) / h;

            value_type der12 = std::max( std::abs( der2 ), std::sqrt( dnf ) );

            value_type h1;

            if( der12 <= 1.0e-15 ) 
                h1 = std::max( 1.0e-6, std::abs(h) * 1.0e-3 );
            else 
                h1 = pow( 0.01/der12, 0.2 );

            h = direction * std::min( value_type(100.0)*std::abs(h), h1 );

            // fprintf( stderr, "dopri5 hinit: h = %f\n", h );
        }

        // --- loop until an acceptable step is taken
        
        bool reject = false;
        
        while( true )
        {
            assert( h * direction >= 0 );

            bool last = false;
            value_type y_new[D], y_stiff[D];

            // check for stepsize underflow
            if( 0.1*std::abs(h) <= std::abs(t)*epsilon ) 
                return STEPSIZE_UNDERFLOW;

            // if the following step ends up very close to t_max, 
            // mark it as the last step
            if( (t + 1.01*h - t_max) * direction > 0.0 ) 
            {
                last = true;
                h = t_max - t;
            }

            // --- perform stages
            
            // k1 is implicitly set from FSAL
            for( unsigned int d=0; d<D; ++d )
                y_new[d] = y[d] + h * a21*k1[d];

            value_type k2[D], k3[D], k4[D], k5[D], k6[D], k7[D]; 

            ++n_eval;
            if( !rhs( t+c2*h, y_new, k2 ) )
                return OUTSIDE_DOMAIN;

            for( unsigned int d=0; d<D; ++d )
                y_new[d] = y[d] + h * ( a31*k1[d] + a32*k2[d] );

            ++n_eval;
            if( !rhs( t+c3*h, y_new, k3 ) )
                return OUTSIDE_DOMAIN;

            for( unsigned int d=0; d<D; ++d )
                y_new[d] = y[d] + h * ( a41*k1[d] + a42*k2[d] + a43*k3[d] );

            ++n_eval;
            if( !rhs( t+c4*h, y_new, k4 ) )
                return OUTSIDE_DOMAIN;

            for( unsigned int d=0; d<D; ++d )
                y_new[d] = y[d] + h * ( a51*k1[d] + a52*k2[d] + a53*k3[d] + a54*k4[d] );

            ++n_eval;
            if( !rhs( t+c5*h, y_new, k5 ) )
                return OUTSIDE_DOMAIN;

            for( unsigned int d=0; d<D; ++d )
                y_new[d] = y_stiff[d] = y[d] + h * ( a61*k1[d] + a62*k2[d] + a63*k3[d] + a64*k4[d] + a65*k5[d] );         

            ++n_eval;
            if( !rhs( t+h, y_new, k6 ) )
                return OUTSIDE_DOMAIN;

            for( unsigned int d=0; d<D; ++d )
                y_new[d] = y[d] + h * ( a71*k1[d] + a73*k3[d] + a74*k4[d] + a75*k5[d] + a76*k6[d] );

            ++n_eval;
            if( !rhs( t+h, y_new, k7 ) )
                return OUTSIDE_DOMAIN;

            // fprintf( stderr, "dopri5 step: t = %f, h = %f, y = %f %f %f, yn = %f %f %f\n", 
            //          t, h, y[0], y[1], y[2], y_new[0], y_new[1], y_new[2] );


            // --- compute error estimate and perform stepsize control

            value_type h_new = h;

            // error estimation

            value_type ee[D], tmp, err = 0.0;

            for( unsigned int d=0; d<D; ++d )
            {
                ee[d] = h * ( e1*k1[d] + e3*k3[d] + e4*k4[d] + e5*k5[d] + e6*k6[d] + e7*k7[d] );

                tmp = ee[d] / (abstol + reltol * std::max( std::abs(y[d]), std::abs(y_new[d]) ));
                err += tmp*tmp;
            }

            err = sqrt( err / D );

            // compute next (potential) stepsize
            double fac11 = pow( err, 0.2 - beta*0.75 );
            
            // Lund-stabilization (require facl <= h_new/h <= facr)
            double fac = fac11 / pow( facold, beta );
            fac = std::max( 1.0/facr, std::min( 1.0/facl, fac/safe ) );
                        
            h_new = h / fac;

            if( err > 1.0 ) 
            {
                // step rejected
                h /= std::min( 1.0/facl, fac11/safe );
                reject = true;
                
                if( n_accepted >= 1 ) 
                    n_rejected++;

                // printf( "dopri5 reject: h = %f\n", h );

                // restart the step
                continue;
            }

            // stiffness detection
            if( (n_accepted % nstiff) == 0 || iasti > 0 ) 
            {
                double stnum, stden;

                for( unsigned int d=0; d<D; ++d )
                {
                    tmp = k7[d]-k6[d];
                    stnum += tmp*tmp;
                    tmp = y_new[d] - y_stiff[d];
                    stden += tmp*tmp;
                }
                
                if( stden > 0.0 ) 
                    hlamb = std::abs(h) * sqrt( stnum/stden );
                
                if( hlamb > 3.25 ) 
                {
                    nonsti = 0;
                    iasti++;
                    
                    if( iasti == 15 ) 
                        return STIFFNESS_DETECTED;
                }
                else 
                {
                    nonsti++;
                    if( nonsti == 6 ) 
                        iasti = 0;
                }
            }

            // step accepted
            facold = std::max( err, value_type(1.0e-4) );
            
            if( reject ) 
            {
                h_new = direction * std::min( std::abs(h_new), std::abs(h) );
                reject = false;
            }
                
            for( unsigned int d=0; d<D; ++d )
                y[d] = y_new[d];

            t = t + h;
            h = h_new;

            // FSAL
            for( unsigned int d=0; d<D; ++d )
                k1[d] = k7[d];
            
            return last ? T_MAX_REACHED : OK;
        }
    }
};

#endif // __dopri5_hpp
