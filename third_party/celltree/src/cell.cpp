#include "cell.hpp"
#include <cmath>
#include <cstdio>

inline void solve3( const float m[4][3], float* r, float& det )
{
    float tmp[3];

    tmp[0] = m[1][1]*m[2][2] - m[1][2]*m[2][1];
    tmp[1] = m[1][2]*m[2][0] - m[1][0]*m[2][2];
    tmp[2] = m[1][0]*m[2][1] - m[1][1]*m[2][0];

    det = m[0][0]*tmp[0] + m[0][1]*tmp[1] + m[0][2]*tmp[2];

    r[0] = m[3][0]*tmp[0] + m[3][1]*tmp[1] + m[3][2]*tmp[2];
    r[0] = -r[0];

    tmp[0] = m[0][1]*m[2][2] - m[0][2]*m[2][1];
    tmp[1] = m[0][2]*m[2][0] - m[0][0]*m[2][2];
    tmp[2] = m[0][0]*m[2][1] - m[0][1]*m[2][0];

    r[1] = m[3][0]*tmp[0] + m[3][1]*tmp[1] + m[3][2]*tmp[2];

    tmp[0] = m[0][1]*m[1][2] - m[0][2]*m[1][1];
    tmp[1] = m[0][2]*m[1][0] - m[0][0]*m[1][2];
    tmp[2] = m[0][0]*m[1][1] - m[0][1]*m[1][0];

    r[2] = m[3][0]*tmp[0] + m[3][1]*tmp[1] + m[3][2]*tmp[2];
    r[2] = -r[2];
}


inline float det3( const float* v0, const float* v1, const float* v2 )
{
    return v0[0]*(v1[1]*v2[2] - v2[1]*v1[2]) + 
           v1[0]*(v2[1]*v0[2] - v0[1]*v2[2]) +
           v2[0]*(v0[1]*v1[2] - v1[1]*v0[2]);
}

inline float cross( float* r, const float* v0, const float* v1 )
{
    float l = 0.0;

    r[0] = v0[1]*v1[2] - v0[2]*v1[1];
    l += r[0]*r[0];
    r[1] = v0[2]*v1[0] - v0[0]*v1[2];
    l += r[1]*r[1];
    r[2] = v0[0]*v1[1] - v0[1]*v1[0];
    l += r[2]*r[2];
    
    return l;
}

inline float dot( const float* v0, const float* v1 )
{
    return v0[0]*v1[0] + v0[1]*v1[1] + v0[2]*v1[2];
}

inline bool intersect_triangle( float& _t, const float* p0, const float* p1, const float* p2, 
                                const float* origin, const float* direction )
{
    
    float u[3] = { p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2] };
    float v[3] = { p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2] };
    float n[3];

    if( cross( n, u, v ) == 0 )
        return false;
    
    float b = dot( n, direction );
    
    if(fabs(b) < 1e-6) 
        return false;

    float w[3] = { origin[0]-p0[0], origin[1]-p0[1], origin[2]-p0[2] };

    _t = -dot( n, w ) / b;

    // if( _t < 0.0 )
    //     return false;

    for( unsigned int d=0; d<3; ++d )
        w[d] = w[d] + _t*direction[d];
    
    // is I inside T?
    const float uu = dot( u, u );
    const float uv = dot( u, v );
    const float vv = dot( v, v );
    const float wu = dot( w, u );
    const float wv = dot( w, v );
    const float D = uv * uv - uu * vv;

    // get and test parametric coords
    float s, t, r;
    s = (uv * wv - vv * wu) / D;
    t = (uv * wu - uu * wv) / D;
    r = 1.0-s-t;

    return s >= 0 && s <= 1 && t >= 0 && t <= 1 && r >= 0 && r <= 1;
}

inline bool intersect_quad( float& t, const float* p0, const float* p1, const float* p2, const float* p3, 
                            const float* o, const float* d )
{
    const float diag1[3] = { p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2] };
    const float diag2[3] = { p3[0]-p1[0], p3[1]-p1[1], p3[2]-p1[2] };
    
    const float d1 = diag1[0]*diag1[0] + diag1[1]*diag1[1] + diag1[2]*diag1[2];
    const float d2 = diag2[0]*diag2[0] + diag2[1]*diag2[1] + diag2[2]*diag2[2];
  
    if( d1 < d2 )
    {
        if( intersect_triangle( t, p0, p1, p2, o, d ) )
        {
            return true;
            //pcoords[0] = pcoords[0] + pcoords[1];
        }
        else if( intersect_triangle( t, p2, p3, p0, o, d ) )
        {
            // pcoords[0] = 1.0 - (pcoords[0]+pcoords[1]);
            // pcoords[1] = 1.0 - pcoords[1];
            return true;
        }
    }
    else
    {
        if( intersect_triangle( t, p0, p1, p3, o, d ) )
        {
            return true;
        }
        else if( intersect_triangle( t, p2, p3, p1, o, d ) )
        {
            return true;
            // pcoords[0] = 1.0 - pcoords[0];
            // pcoords[1] = 1.0 - pcoords[1];
        }
    }

    return false;
}

// -------------------------------------------------------------------------

template<unsigned int KIND>
class cell_traits
{
};

// --- tetrahedron ------------------------------------------------------------

template<>
struct cell_traits<celltree::TETRAHEDRON>
{
    static void interpolant( float* result, const float* values, const float* param, 
                             const unsigned int dim )
    {
        const float w[4] = { param[0], param[1], param[2], 1.0f-param[0]-param[1]-param[2] };
        
        for( int i=0; i<dim; ++i )
        {
            result[i] = 0;
    
            for( int j=0; j<4; ++j )
                result[i] += w[j] * values[dim*j+i];
        }
    }

    static bool invert( float* coeff, const float* pts, const float* pos )
    {
        const float epsilon = 1e-5;
        
        float v[4][3], c[4], vol;
    
        for( int i=0; i<3; ++i )
        {
            v[3][i] = pos[i] - pts[9+i];
        
            for( int j=0; j<3; ++j )
                v[j][i] = pts[9+i] - pts[3*j+i];
        }

        solve3( v, c, vol );
        c[3] = vol - c[0] - c[1] - c[2];

        for( int i=0; i<4; ++i )
            c[i] /= vol;

        if( c[0] < -epsilon || c[1] < -epsilon || 
            c[2] < -epsilon || c[3] < -epsilon )
            return false;
    
        for( int i=0; i<3; ++i )    
            coeff[i] = c[i];

        return true;
    }    

    static unsigned int intersect( float* t, const float* pts, const float* origin, const float* direction )
    {
        unsigned int k = 0;
        
        if( k<2 && intersect_triangle( t[k], pts+0, pts+3, pts+9, origin, direction ) )
            ++k;
        if( k<2 && intersect_triangle( t[k], pts+3, pts+6, pts+9, origin, direction ) )
            ++k;
        if( k<2 && intersect_triangle( t[k], pts+6, pts+0, pts+9, origin, direction ) )
            ++k;
        if( k<2 && intersect_triangle( t[k], pts+0, pts+6, pts+3, origin, direction ) )
            ++k;
            
        return k;
    }
};

// ----------------------------------------------------------------------------

template<>
struct cell_traits<celltree::HEXAHEDRON>
{
    static void interpolant( float* result, const float* values, const float* c,
                             const unsigned int dim )    
    {
        const float d[3] = { 1.0f-c[0], 1.0f-c[1], 1.0f-c[2] };
        
        for( int i=0; i<dim; ++i )
        {
            result[i] = values[3*0+i]*d[0]*d[1]*d[2] + 
                        values[3*1+i]*c[0]*d[1]*d[2] + 
                        values[3*2+i]*c[0]*c[1]*d[2] + 
                        values[3*3+i]*d[0]*c[1]*d[2] + 
                        values[3*4+i]*d[0]*d[1]*c[2] + 
                        values[3*5+i]*c[0]*d[1]*c[2] + 
                        values[3*6+i]*c[0]*c[1]*c[2] + 
                        values[3*7+i]*d[0]*c[1]*c[2];  
        }
    }
    
    static void derivative( float* result, const float* values, const float* c,
                            const unsigned int dim )
    {
        const float d[3] = { 1.0f-c[0], 1.0f-c[1], 1.0f-c[2] };

        for( int i=0; i<dim; ++i )
        {
            result[dim*0+i] = -values[dim*0+i]*d[1]*d[2] + values[dim*1+i]*d[1]*d[2] 
                              +values[dim*2+i]*c[1]*d[2] - values[dim*3+i]*c[1]*d[2]
                              -values[dim*4+i]*d[1]*c[2] + values[dim*5+i]*d[1]*c[2] 
                              +values[dim*6+i]*c[1]*c[2] - values[dim*7+i]*c[1]*c[2];

            result[dim*1+i] = -values[dim*0+i]*d[0]*d[2] - values[dim*1+i]*c[0]*d[2]
                              +values[dim*2+i]*c[0]*d[2] + values[dim*3+i]*d[0]*d[2]
                              -values[dim*4+i]*d[0]*c[2] - values[dim*5+i]*c[0]*c[2] 
                              +values[dim*6+i]*c[0]*c[2] + values[dim*7+i]*d[0]*c[2];

            result[dim*2+i] = -values[dim*0+i]*d[0]*d[1] - values[dim*1+i]*c[0]*d[1]
                              -values[dim*2+i]*c[0]*c[1] - values[dim*3+i]*d[0]*c[1]
                              +values[dim*4+i]*d[0]*d[1] + values[dim*5+i]*c[0]*d[1] 
                              +values[dim*6+i]*c[0]*c[1] + values[dim*7+i]*d[0]*c[1];
        }
    }

    static bool invert( float* c, const float* pts, const float* pos )
    {
        const int   maxiter = 8;
        const float epsilon = 1e-5;

        float h[3], d[9], p[3], denom;

        for( int i=0; i<3; ++i )
            c[i] = 0.5;
   
        for( int iter=0; iter<maxiter; ++iter )
        {
            interpolant( p, pts, c, 3 );
            derivative( d, pts, c, 3 );
            denom = det3( d, d+3, d+6 );

            for( int i=0; i<3; ++i )
                p[i] -= pos[i];

            c[0] -= (h[0] = det3( p, d+3, d+6 ) / denom);
            c[1] -= (h[1] = det3( d+0, p, d+6 ) / denom);
            c[2] -= (h[2] = det3( d+0, d+3, p ) / denom);

            if( std::abs(h[0])<epsilon && std::abs(h[1])<epsilon && std::abs(h[2])<epsilon )
                break;
        }

        return c[0] > -0.001 && c[1] > -0.001 && c[2] > -0.001 && 
               c[0] <  1.001 && c[1] <  1.001 && c[2] <  1.001;
    }
    
    static unsigned int intersect( float* t, const float* pts, const float* origin, const float* direction )
    {
        unsigned int k = 0;
        
        if( k<2 && intersect_quad( t[k], pts+0,  pts+12, pts+21, pts+9,  origin, direction ) )
            ++k;
        if( k<2 && intersect_quad( t[k], pts+3,  pts+6,  pts+18, pts+15, origin, direction ) )
            ++k;
        if( k<2 && intersect_quad( t[k], pts+0,  pts+3,  pts+15, pts+12, origin, direction ) )
            ++k;
        if( k<2 && intersect_quad( t[k], pts+9,  pts+21, pts+18, pts+6,  origin, direction ) )
            ++k;
        if( k<2 && intersect_quad( t[k], pts+0,  pts+9,  pts+6,  pts+3,  origin, direction ) )
            ++k;
        if( k<2 && intersect_quad( t[k], pts+12, pts+15, pts+18, pts+21, origin, direction ) )
            ++k;
            
        return k;
    }
};
   
// --- prism ------------------------------------------------------------------

template<>
struct cell_traits<celltree::PRISM>
{
    static void interpolant( float* result, const float* values, const float* c,
                             const unsigned int dim )    
    {
        const float d[2] = { 1.0f-c[0]-c[1], 1.0f-c[2] };

        for( int i=0; i<dim; ++i )
        {
            result[i] = values[dim*0+i]*d[0]*d[1] +
                        values[dim*1+i]*c[0]*d[1] +
                        values[dim*2+i]*c[1]*d[1] + 
                        values[dim*3+i]*d[0]*c[2] + 
                        values[dim*4+i]*c[0]*c[2] + 
                        values[dim*5+i]*c[1]*c[2];
        }
    }

    static void derivative( float *result, const float* values, const float* c,
                            const unsigned int dim )
    {
        const float d[2] = { 1.0f-c[0]-c[1], 1.0f-c[2] };

        for( int i=0; i<dim; ++i )
        {
            result[dim*0+i] = d[1]*( values[dim*1+i] - values[dim*0+i] ) + 
                              c[2]*( values[dim*4+i] - values[dim*3+i] );
                
            result[dim*1+i] = d[1]*( values[dim*2+i] - values[dim*0+i] ) +
                              c[2]*( values[dim*5+i] - values[dim*3+i] );

            result[dim*2+i] = d[0]*( values[dim*3+i] - values[dim*0+i] ) +
                              c[0]*( values[dim*4+i] - values[dim*1+i] ) +
                              c[1]*( values[dim*5+i] - values[dim*2+i] );
        }
    }

    static bool invert( float* c, const float* pts, const float* pos )
    {
        const int   maxiter = 8;
        const float epsilon = 1e-5;

        float h[3], d[9], p[3], denom;

        // for( int i=0; i<3; ++i )
        c[0] = 0.33;
        c[1] = 0.33;
        c[2] = 0.5;
   
        // Newton iteration

        int iter = 0;
        for( iter=0; iter<maxiter; ++iter )
        {
            interpolant( p, pts, c, 3 );

            for( int i=0; i<3; ++i )
                p[i] -= pos[i];
            
            derivative( d, pts, c, 3 );

            float denom = det3( d+0, d+3, d+6 );
                
            if( fabs(denom) < 1e-20 )
                return false;
        
            c[0] -= (h[0] = det3( p, d+3, d+6 ) / denom);
            c[1] -= (h[1] = det3( d+0, p, d+6 ) / denom);
            c[2] -= (h[2] = det3( d+0, d+3, p ) / denom);

            if( std::abs(h[0])<epsilon && std::abs(h[1])<epsilon && std::abs(h[2])<epsilon )
                break;
        }

        return c[0] > -0.001 && c[1] > -0.001 && (1.0-c[0]-c[1]) > -0.001 && 
               c[2] > -0.001 && c[2] <  1.001;
    }
    
    static unsigned int intersect( float* t, const float* pts, const float* origin, const float* direction )
    {
        unsigned int k = 0;
        
        if( k<2 && intersect_triangle( t[k], pts+0,  pts+3,  pts+6,  origin, direction ) )
            ++k;
        if( k<2 && intersect_triangle( t[k], pts+9,  pts+15, pts+12, origin, direction ) )
            ++k;
        if( k<2 && intersect_quad( t[k], pts+0,  pts+9,  pts+12,  pts+3,  origin, direction ) )
            ++k;
        if( k<2 && intersect_quad( t[k], pts+3,  pts+12,  pts+15,  pts+6,  origin, direction ) )
            ++k;
        if( k<2 && intersect_quad( t[k], pts+6,  pts+15,  pts+12,  pts+0,  origin, direction ) )
            ++k;

        return k;
    }
};

// --- pyramid -------------------------------------------------------------

template<>
struct cell_traits<celltree::PYRAMID>
{
    static void interpolant( float* result, const float* values, const float* c,
                             const unsigned int dim )    
    {
        const float m[3] = { 1.0f-c[0], 1.0f-c[1], 1.0f-c[2] };
        const float d[5] = { m[0]*m[1]*m[2], c[0]*m[1]*m[2], c[0]*c[1]*m[2], m[0]*c[1]*m[2], c[2] };

        for( int i=0; i<dim; ++i )
        {
            result[i] = values[dim*0+i]*d[0] +
                        values[dim*1+i]*d[1] +
                        values[dim*2+i]*d[2] + 
                        values[dim*3+i]*d[3] + 
                        values[dim*4+i]*d[4];
        }
    }

    static void derivative( float *result, const float* values, const float* c,
                            const unsigned int dim )
    {
        const float rm = 1.0 - c[0];
        const float sm = 1.0 - c[1];
        const float tm = 1.0 - c[2];

        float d[15];

        // r-derivatives
        d[0] = -sm*tm;
        d[1] = sm*tm;
        d[2] = c[1]*tm;
        d[3] = -c[1]*tm;
        d[4] = 0.0;

        // s-derivatives
        d[5] = -rm*tm;
        d[6] = -c[0]*tm;
        d[7] = c[0]*tm;
        d[8] = rm*tm;
        d[9] = 0.0;

        // t-derivatives
        d[10] = -rm*sm;
        d[11] = -c[0]*sm;
        d[12] = -c[0]*c[1];
        d[13] = -rm*c[1];
        d[14] = 1.0;

        for( int i=0; i<dim; ++i )
        {
            result[dim*0+i] = d[0]*values[dim*0+i] + d[1]*values[dim*1+i] + 
                              d[2]*values[dim*2+i] + d[3]*values[dim*3+i] + 
                              d[4]*values[dim*4+i];
                
            result[dim*1+i] = d[5]*values[dim*0+i] + d[6]*values[dim*1+i] + 
                              d[7]*values[dim*2+i] + d[8]*values[dim*3+i] + 
                              d[9]*values[dim*4+i];

            result[dim*2+i] = d[10]*values[dim*0+i] + d[11]*values[dim*1+i] + 
                              d[12]*values[dim*2+i] + d[12]*values[dim*3+i] + 
                              d[14]*values[dim*4+i];
        }
    }

    static bool invert( float* c, const float* pts, const float* pos )
    {
        const int   maxiter = 8;
        const float epsilon = 1e-5;

        float h[3], d[9], p[3], denom;

        // for( int i=0; i<3; ++i )
        c[0] = 0.5;
        c[1] = 0.5;
        c[2] = 0.5;
   
        // Newton iteration
        int iter = 0;
        for( iter=0; iter<maxiter; ++iter )
        {
            interpolant( p, pts, c, 3 );

            for( int i=0; i<3; ++i )
                p[i] -= pos[i];
            
            derivative( d, pts, c, 3 );

            float denom = det3( d+0, d+3, d+6 );
                
            if( fabs(denom) < 1e-20 )
                return false;
        
            c[0] -= (h[0] = det3( p, d+3, d+6 ) / denom);
            c[1] -= (h[1] = det3( d+0, p, d+6 ) / denom);
            c[2] -= (h[2] = det3( d+0, d+3, p ) / denom);

            if( std::abs(h[0])<epsilon && std::abs(h[1])<epsilon && std::abs(h[2])<epsilon )
                break;
        }

        return c[0] > -0.001 && c[1] > -0.001 && c[2] > -0.001 && 
               c[0] <  1.001 && c[1] <  1.001 && c[2] <  1.001;
    }
    
    static unsigned int intersect( float* t, const float* pts, const float* origin, const float* direction )
    {
        unsigned int k = 0;
        
        if( k<2 && intersect_quad( t[k], pts+0, pts+9, pts+6, pts+3,  origin, direction ) )
            ++k;
        if( k<2 && intersect_triangle( t[k], pts+0,  pts+3, pts+12, origin, direction ) )
            ++k;
        if( k<2 && intersect_triangle( t[k], pts+3,  pts+6, pts+12, origin, direction ) )
            ++k;
        if( k<2 && intersect_triangle( t[k], pts+6,  pts+9, pts+12, origin, direction ) )
            ++k;
        if( k<2 && intersect_triangle( t[k], pts+9,  pts+0, pts+12, origin, direction ) )
            ++k;

        return k;
    }
};

// -------------------------------------------------------------------------

bool tet_invert( float* c, const float* pts, const float* pos )
{
    return cell_traits<celltree::TETRAHEDRON>::invert( c, pts, pos );
}

// -------------------------------------------------------------------------
bool celltree::invert_cell( float* c, const float* pts, const float* pos, cell_kind kind )
{
    typedef bool (*invert_func)( float* c, const float* pts, const float* pos );
    
    const invert_func funcs[] = {
        cell_traits<TETRAHEDRON>::invert,
        cell_traits<HEXAHEDRON>::invert,
        cell_traits<PYRAMID>::invert,
        cell_traits<PRISM>::invert,
    };
    
    return funcs[kind]( c, pts, pos );
}

void celltree::interpolate_cell( float* result, const float* values, const float* c,
                       const unsigned int dim, cell_kind kind )
{
    typedef void (*interpolate_func)( float* result, const float* values, 
                                      const float* c, const unsigned int dim );
    
    const interpolate_func funcs[] = {
        cell_traits<TETRAHEDRON>::interpolant,
        cell_traits<HEXAHEDRON>::interpolant,
        cell_traits<PYRAMID>::interpolant,
        cell_traits<PRISM>::interpolant,
    };
    
    return funcs[kind]( result, values, c, dim );
}


unsigned int celltree::intersect_cell( float* t, const float* pts, const float* origin, const float* direction,
                             cell_kind kind )
{
    typedef unsigned int (*intersect_func)( float* t, 
                                            const float* pts, const float* origin, const float* direction );
    
    const intersect_func funcs[] = {
        cell_traits<TETRAHEDRON>::intersect,
        cell_traits<HEXAHEDRON>::intersect,
        cell_traits<PYRAMID>::intersect,
        cell_traits<PRISM>::intersect,
    };
    
    return funcs[kind]( t, pts, origin, direction );
}
