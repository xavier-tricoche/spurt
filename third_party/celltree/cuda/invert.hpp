#ifndef __invert_hpp
#define __invert_hpp

#define CELL_INVERT_USE_NEWTON 0

#include <cutil_math.h>

__host__ __device__
float det3( float3 v0, float3 v1, float3 v2 )
{
    return v0.x*(v1.y*v2.z - v2.y*v1.z) + 
           v1.x*(v2.y*v0.z - v0.y*v2.z) +
           v2.x*(v0.y*v1.z - v1.y*v0.z);
}

__host__ __device__
inline void solve3( float3 m0, float3 m1, float3 m2, float3 r, float3& c )
{
    float3 tmp;

    tmp.x = m1.y*m2.z - m1.z*m2.y;
    tmp.y = m1.z*m2.x - m1.x*m2.z;
    tmp.z = m1.x*m2.y - m1.y*m2.x;

    float det = 1.0f / (m0.x*tmp.x + m0.y*tmp.y + m0.z*tmp.z);

    c.x = (r.x*tmp.x + r.y*tmp.y + r.z*tmp.z) * det;

    tmp.x = m2.y*m0.z - m2.z*m0.y;
    tmp.y = m2.z*m0.x - m2.x*m0.z;
    tmp.z = m2.x*m0.y - m2.y*m0.x;

    c.y = (r.x*tmp.x + r.y*tmp.y + r.z*tmp.z) * det;

    tmp.x = m0.y*m1.z - m0.z*m1.y;
    tmp.y = m0.z*m1.x - m0.x*m1.z;
    tmp.z = m0.x*m1.y - m0.y*m1.x;

    c.z = (r.x*tmp.x + r.y*tmp.y + r.z*tmp.z) * det;
}

__host__ __device__
bool invert_tet( float3& c, float3 p0, float3 p1, float3 p2, float3 p3, float3 pos )
{
    solve3( p0-p3, p1-p3, p2-p3, pos-p3, c );
    return c.x > -0.001f && c.y > -0.001f && c.z > -0.001f && c.x+c.y+c.z < 1.001f;
}

// --------------------------------------------------------------------------

template<typename T>
__host__ __device__ 
T interpolate_tet( const T v[4], const float3& c )
{
    return c.x * v[0] + c.y * v[1] +
           c.z * v[2] + (1.0f-c.x-c.y-c.z) * v[3];
}

// --------------------------------------------------------------------------

__host__ __device__
bool invert_tet( const float3 pts[4], float3& c, float3 pos )
{
    if( pos.x < fminf( fminf( pts[0].x, pts[1].x ), fminf( pts[2].x, pts[3].x ) ) ||
        pos.x > fmaxf( fmaxf( pts[0].x, pts[1].x ), fmaxf( pts[2].x, pts[3].x ) ) || 
        pos.y < fminf( fminf( pts[0].y, pts[1].y ), fminf( pts[2].y, pts[3].y ) ) ||
        pos.y > fmaxf( fmaxf( pts[0].y, pts[1].y ), fmaxf( pts[2].y, pts[3].y ) ) || 
        pos.z < fminf( fminf( pts[0].z, pts[1].z ), fminf( pts[2].z, pts[3].z ) ) ||
        pos.z > fmaxf( fmaxf( pts[0].z, pts[1].z ), fmaxf( pts[2].z, pts[3].z ) ) )
        return false;

    return invert_tet( c, pts[0], pts[1], pts[2], pts[3], pos );
}

// --------------------------------------------------------------------------

template<typename T>
__host__ __device__ 
T interpolate_hex( const T v[8], float3& c )
{
    float3 v0 = (1.0f-c.x) * v[0] + c.x * v[1];
    float3 v1 = (1.0f-c.x) * v[3] + c.x * v[2];
    float3 v2 = (1.0f-c.x) * v[4] + c.x * v[5];
    float3 v3 = (1.0f-c.x) * v[7] + c.x * v[6];
    
    v0 = (1.0f-c.y)*v0 + c.y*v1;
    v1 = (1.0f-c.y)*v2 + c.y*v3;
    
    return (1.0f-c.z)*v0 + c.z*v1;
}

__host__ __device__ 
bool invert_hex( const float3 pts[8], float3& c, float3 pos )
{
#if CELL_INVERT_USE_NEWTON

    const int maxiter = 8;

    c = make_float3( 0.5f, 0.5f, 0.5f );
            
    for( int iter=0; iter<maxiter; ++iter )
    {
        float3 d = make_float3( 1.0f-c.x, 1.0f-c.y, 1.0f-c.z );

        float3 rp = pts[0]*d.x*d.y*d.z + pts[1]*c.x*d.y*d.z + 
                    pts[2]*c.x*c.y*d.z + pts[3]*d.x*c.y*d.z + 
                    pts[4]*d.x*d.y*c.z + pts[5]*c.x*d.y*c.z + 
                    pts[6]*c.x*c.y*c.z + pts[7]*d.x*c.y*c.z - pos;  

        float3 d0 =  pts[1]*d.y*d.z - pts[0]*d.y*d.z + pts[2]*c.y*d.z - pts[3]*c.y*d.z
                    -pts[4]*d.y*c.z + pts[5]*d.y*c.z + pts[6]*c.y*c.z - pts[7]*c.y*c.z;
        
        float3 d1 =  pts[2]*c.x*d.z + pts[3]*d.x*d.z - pts[0]*d.x*d.z - pts[1]*c.x*d.z
                    -pts[4]*d.x*c.z - pts[5]*c.x*c.z + pts[6]*c.x*c.z + pts[7]*d.x*c.z;
        
        float3 d2 =  pts[4]*d.x*d.y + pts[5]*c.x*d.y  +pts[6]*c.x*c.y + pts[7]*d.x*c.y
                    -pts[0]*d.x*d.y - pts[1]*c.x*d.y - pts[2]*c.x*c.y - pts[3]*d.x*c.y;

        float denom = det3( d0, d1, d2 );

        c.x -= det3( rp, d1, d2 ) / denom;
        c.y -= det3( rp, d2, d0 ) / denom;
        c.z -= det3( rp, d0, d1 ) / denom;
    }

    return c.x > -0.001f && c.y > -0.001f && c.z > -0.001f && 
           c.x <  1.001f && c.y <  1.001f && c.z <  1.001f;

#else

    if( invert_tet( c, pts[1], pts[4], pts[6], pts[3], pos ) )                                                                                   
    {
        c = make_float3( c.x+c.z, 1.0f-c.x-c.y, c.y+c.z );
        return true;
    }

    if( invert_tet( c, pts[3], pts[0], pts[1], pts[4], pos ) )                                                                                   
    {
        c = make_float3( c.z, c.x, 1.0f-c.x-c.y-c.z );                                              
        return true;
    }

    if( invert_tet( c, pts[1], pts[4], pts[5], pts[6], pos ) ) 
    {
        c = make_float3( 1.0f-c.y, 1.0f-c.x-c.y-c.z, 1.0f-c.x );
        return true;
    }                                                                                                                                            
                                                                                                                      
    if( invert_tet( c, pts[3], pts[6], pts[2], pts[1], pos ) )
    {
        c = make_float3( 1.0f-c.x, c.x+c.y+c.z, c.y );
        return true;
    }                                                                                                                                            

    if( invert_tet( c, pts[3], pts[6], pts[7], pts[4], pos ) )
    {
        c = make_float3( c.y, c.x+c.y+c.z, 1.0f-c.x );
        return true;
    }

    return false;

#endif
}

// --------------------------------------------------------------------------

template<typename T>
__host__ __device__ 
T interpolate_prism( const T v[6], float3& c )
{
    float3 v0 = (1.0f-c.x-c.y) * v[0] + c.x*v[1] + c.y*v[2];
    float3 v1 = (1.0f-c.x-c.y) * v[3] + c.x*v[4] + c.y*v[5];
    return (1.0f-c.z)*v0 + c.z*v1;
}

__host__ __device__
bool invert_prism( const float3 pts[6], float3& c, float3 pos )
{
#if CELL_INVERT_USE_NEWTON

    const int maxiter = 8;

    c = make_float3( 0.5f, 0.5f, 0.5f );
            
    for( int iter=0; iter<maxiter; ++iter )
    {
        float2 d = make_float2( 1.0f-c.x-c.y, 1.0f-c.z );

        float3 rp = pts[0]*d.x*d.y + pts[1]*c.x*d.y + 
                    pts[2]*c.y*d.y + pts[3]*d.x*c.z + 
                    pts[4]*c.x*c.z + pts[5]*c.y*c.z;

        float3 d0 = d.y*(pts[1]-pts[0]) + c.z*(pts[4]-pts[3]);
        float3 d1 = d.y*(pts[2]-pts[0]) + c.z*(pts[5]-pts[3]);
        float3 d2 = d.x*(pts[3]-pts[0]) + c.x*(pts[4]-pts[1]) + 
                    c.y*(pts[5]-pts[2]);

        float denom = det3( d0, d1, d2 );

        c.x -= det3( rp, d1, d2 ) / denom;
        c.y -= det3( rp, d2, d0 ) / denom;
        c.z -= det3( rp, d0, d1 ) / denom;
    }

    return c.x > -0.001f && c.y > -0.001f && c.z > -0.001f && 
           c.x <  1.001f && c.y <  1.001f && c.z <  1.001f;

#else 

    // FIXME: figure out correct prism coordinates
    //        from tet barycentric coordinates
    #warning "FIXME: using prisms with tet decomposition yields wrong coordinates"

    if( invert_tet( c, pts[0], pts[2], pts[1], pts[3], pos ) )                                                                                   
    {
        // FIXME: here
        c = make_float3( 1.0f, 0.0f, 0.0f );
        return true;
    }

    if( invert_tet( c, pts[1], pts[2], pts[5], pts[3], pos ) )                                                                                   
    {
        // FIXME: here
        c = make_float3( 2.0f, 0.0f, 0.0f );
        return true;
    }

    if( invert_tet( c, pts[1], pts[3], pts[4], pts[5], pos ) ) 
    {
        // FIXME: here
        c = make_float3( 3.0f, 0.0f, 0.0f );
        return true;
    }                                                                                                                                            

    return false;

#endif
}

// -------------------------------------------------------------------------

template<typename T>
__host__ __device__
T interpolate_pyramid( const T v[5], float3& c )
{
    const float3 d = make_float3( 1.0f-c.x, 1.0f-c.y, 1.0f-c.z );

    return d.z * (d.x*d.y*v[0] + c.x*d.y*v[1] + 
                  c.x*c.y*v[2] + d.x*c.y*v[3] ) + c.z * v[4];
}

__host__ __device__
bool invert_pyramid( const float3 pts[5], float3& c, float3 pos )
{
#if CELL_INVERT_USE_NEWTON

    const int maxiter = 8;

    c = make_float3( 0.5f, 0.5f, 0.5f );
            
    for( int iter=0; iter<maxiter; ++iter )
    {
        float3 d = make_float3( 1.0f-c.x, 1.0f-c.y, 1.0f-c.z );

        float3 rp = pts[0]*d.x*d.y*d.z + pts[1]*c.x*d.y*d.z + 
                    pts[2]*c.x*c.y*d.z + pts[3]*d.x*c.y*d.z + 
                    pts[4]*c.z;

        float3 d0 = pts[1]*d.y*d.z - pts[0]*d.y*d.z + pts[2]*c.y*d.z - pts[3]*c.y*d.z;
        float3 d1 = pts[2]*c.x*d.z - pts[3]*c.x*d.z - pts[0]*d.x*d.z - pts[1]*c.x*d.z;
        float3 d2 = make_float3( 0, 0, 0 ) - pts[0]*d.x*d.y - pts[1]*c.x*d.y - pts[2]*c.x*c.y - pts[3]*d.x*c.y;

        float denom = det3( d0, d1, d2 );

        c.x -= det3( rp, d1, d2 ) / denom;
        c.y -= det3( rp, d2, d0 ) / denom;
        c.z -= det3( rp, d0, d1 ) / denom;
    }

    return c.x > -0.001f && c.y > -0.001f && c.z > -0.001f && 
           c.x <  1.001f && c.y <  1.001f && c.z <  1.001f;

#else

    // FIXME: figure out correct pyramid coordinates
    //        from tet barycentric coordinates
    #warning "FIXME: using pyramids with tet decomposition yields wrong coordinates"

    if( invert_tet( c, pts[0], pts[1], pts[2], pts[4], pos ) )                                                                                   
    {
        // FIXME: here
        c = make_float3( 0.5f, 0.5f, 0.5f );
        return true;
    }

    if( invert_tet( c, pts[0], pts[2], pts[3], pts[4], pos ) )
    {
        // FIXME: here
        c = make_float3( 0.5f, 0.5f, 0.5f );
        return true;
    }

    return false;

#endif
}


#endif // __invert_hpp
