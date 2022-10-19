#include <cmath>
#include <cstdio>

#include "trackball.hpp"

trackball::trackball()
{
    m_cur.r[0] = 0.0;
    m_cur.r[1] = 0.0;
    m_cur.r[2] = 0.0;
    m_cur.r[3] = 1.0;

    m_cur.c[0] = 0;
    m_cur.c[1] = 0;
    m_cur.c[2] = 0;

    m_cur.d = 10.0;

    m_home = m_cur;
    m_mode = NONE;

    update();
}

// --------------------------------------------------------------------------

void trackball::rotate_vq( double v[3], const double q[4] )
{
    double xx = 2.0 * q[0] * q[0];
    double xy = 2.0 * q[0] * q[1];
    double xz = 2.0 * q[0] * q[2];
    double yy = 2.0 * q[1] * q[1];
    double yz = 2.0 * q[1] * q[2];
    double zz = 2.0 * q[2] * q[2];
    double wx = 2.0 * q[3] * q[0];
    double wy = 2.0 * q[3] * q[1];
    double wz = 2.0 * q[3] * q[2];

    double c[3] = { v[0], v[1], v[2] };

    v[0] = (1.0-yy-zz)*c[0] + (xy-wz)*c[1] + (xz+wy)*c[2];
    v[1] = (xy+wz)*c[0] + (1.0-xx-zz)*c[1] + (yz-wx)*c[2];
    v[2] = (xz-wy)*c[0] + (yz+wx)*c[1] + (1.0-xx-yy)*c[2];
}

// --------------------------------------------------------------------------

void trackball::rotate_qq( double q[4], const double r[4] )
{
    double x = r[3]*q[0] + r[0]*q[3] + r[1]*q[2] - r[2]*q[1];
    double y = r[3]*q[1] - r[0]*q[2] + r[1]*q[3] + r[2]*q[0];
    double z = r[3]*q[2] + r[0]*q[1] - r[1]*q[0] + r[2]*q[3];
    double w = r[3]*q[3] - r[0]*q[0] - r[1]*q[1] - r[2]*q[2];
    
    q[0] = x;
    q[1] = y;
    q[2] = z;
    q[3] = w;
}

// --------------------------------------------------------------------------

void trackball::normalize( double* v )
{
    double l = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];

    if( l != 0 )
    {
        l = sqrt( l );

        v[0] /= l;
        v[1] /= l;
        v[2] /= l;
    }
}

// --------------------------------------------------------------------------

void trackball::apply_trackball( const double* p1, const double* p2 )
{
    const double tbsize = 0.8;
    const double r2 = tbsize*tbsize / 2.0;

    double d1 = p1[0]*p1[0] + p1[1]*p1[1];
    double d2 = p2[0]*p2[0] + p2[1]*p2[1];

    double sp1[3] = { p1[0], p1[1], d1 < tbsize ? sqrt( 2.0*tbsize-d1 ) : tbsize/sqrt(d1) };
    double sp2[3] = { p2[0], p2[1], d2 < tbsize ? sqrt( 2.0*tbsize-d2 ) : tbsize/sqrt(d2) };

    double axis[3] = { sp2[1]*sp1[2] - sp2[2]*sp1[1],
                       sp2[2]*sp1[0] - sp2[0]*sp1[2],
                       sp2[0]*sp1[1] - sp2[1]*sp1[0] };

    double saxis[3] = { axis[0], axis[1], axis[2] };

    rotate_vq( axis, m_cur.r );
    normalize( axis );

    sp2[0] -= sp1[0];
    sp2[1] -= sp1[1];
    sp2[2] -= sp1[2];

    double angle = sqrt( sp2[0]*sp2[0] + sp2[1]*sp2[1] + sp2[2]*sp2[2] ) / tbsize;

    if( angle >  1.0 )  angle =  1.0;
    if( angle < -1.0 )  angle = -1.0;

    angle = asin( angle );

    const double ch = cos( 0.5*angle );
    const double sh = sin( 0.5*angle );

    double qr[4] = { axis[0]*sh, axis[1]*sh, axis[2]*sh, ch };
    rotate_qq( m_cur.r, qr );
}

// --------------------------------------------------------------------------

void trackball::begin( Mode mode, float mx, float my )
{
    m_m0[0] = m_m1[0] = mx;
    m_m0[1] = m_m1[1] = my;

    m_mode = mode;
}

// --------------------------------------------------------------------------

void trackball::end()
{
    m_mode = NONE;
}

// --------------------------------------------------------------------------

void trackball::move( float mx, float my )
{
    m_m1[0] = mx;
    m_m1[1] = my;

    switch( m_mode )
    {
    case ROTATE:
        apply_trackball( m_m0, m_m1 );
        break;
    case ZOOM:
    {
        double s = 1.0 + (m_m1[1] - m_m0[1]);

        if( m_cur.d * s > 0.05 )
        {
            m_cur.d *= s;
        }
        else
        {
            double d[3] = { 0, 0, -(m_m1[1]-m_m0[1])*s };

            rotate_vq( d, m_cur.r );
            m_cur.c[0] += d[0];
            m_cur.c[1] += d[1];
            m_cur.c[2] += d[2];
        }
        break;
    }
    case PUSH:
    {
        double d[3] = { 
            0, 
            0, 
            (m_m1[1] - m_m0[1])*m_cur.d 
        };

        rotate_vq( d, m_cur.r );

        m_cur.c[0] += d[0];
        m_cur.c[1] += d[1];
        m_cur.c[2] += d[2];
        break;
    }
    case PAN:
    {
        double d[3] = { 
            -0.3 * m_cur.d * (m_m1[0] - m_m0[0]),
            -0.3 * m_cur.d * (m_m1[1] - m_m0[1]),
            0
        };

        rotate_vq( d, m_cur.r );

        m_cur.c[0] += d[0];
        m_cur.c[1] += d[1];
        m_cur.c[2] += d[2];
        break;
    }
    default:
        break;
    }
    
    update();

    m_m0[0] = m_m1[0];
    m_m0[1] = m_m1[1];
}

// --------------------------------------------------------------------------

void trackball::update()
{
    m_mat[0 ] = 1.0 - 2.0*m_cur.r[1]*m_cur.r[1] - 2.0*m_cur.r[2]*m_cur.r[2];
    m_mat[1 ] = 2.0*m_cur.r[0]*m_cur.r[1]-2.0*m_cur.r[2]*m_cur.r[3];
    m_mat[2 ] = 2.0*m_cur.r[0]*m_cur.r[2]+2.0*m_cur.r[1]*m_cur.r[3];
    m_mat[3 ] = 0;

    m_mat[4 ] = 2.0*m_cur.r[0]*m_cur.r[1]+2.0*m_cur.r[2]*m_cur.r[3];
    m_mat[5 ] = 1.0 - 2.0*m_cur.r[0]*m_cur.r[0] - 2.0*m_cur.r[2]*m_cur.r[2];
    m_mat[6 ] = 2.0*m_cur.r[1]*m_cur.r[2]-2.0*m_cur.r[0]*m_cur.r[3];
    m_mat[7 ] = 0;

    m_mat[8 ] = 2.0*m_cur.r[0]*m_cur.r[2]-2.0*m_cur.r[1]*m_cur.r[3];
    m_mat[9 ] = 2.0*m_cur.r[1]*m_cur.r[2]+2.0*m_cur.r[0]*m_cur.r[3];
    m_mat[10] = 1.0 - 2.0*m_cur.r[0]*m_cur.r[0]- 2.0*m_cur.r[1]*m_cur.r[1];
    m_mat[11] = 0;

    m_mat[12] = -(m_mat[0]*m_cur.c[0] + m_mat[4]*m_cur.c[1] + m_mat[8]*m_cur.c[2]);
    m_mat[13] = -(m_mat[1]*m_cur.c[0] + m_mat[5]*m_cur.c[1] + m_mat[9]*m_cur.c[2]);
    m_mat[14] = -(m_mat[2]*m_cur.c[0] + m_mat[6]*m_cur.c[1] + m_mat[10]*m_cur.c[2] + m_cur.d);
    m_mat[15] = 1;
}

// --------------------------------------------------------------------------

const double* trackball::matrix() const
{
    return m_mat;
}

// --------------------------------------------------------------------------

void trackball::home( double eyex, double eyey, double eyez,
                      double cenx, double ceny, double cenz,
                      double upx, double upy, double upz )
{
    double look[3], side[3], up[3];

    look[0] = cenx - eyex;
    look[1] = ceny - eyey;
    look[2] = cenz - eyez;

    up[0] = upx;
    up[1] = upy;
    up[2] = upz;

    m_cur.c[0] = cenx;
    m_cur.c[1] = ceny;
    m_cur.c[2] = cenz;

    m_cur.d = sqrt( look[0]*look[0] + look[1]*look[1] + look[2]*look[2] );

    if( m_cur.d )
    {
        look[0] /= m_cur.d;
        look[1] /= m_cur.d;
        look[2] /= m_cur.d;
    }

    side[0] = look[1]*up[2] - look[2]*up[1];
    side[1] = look[2]*up[0] - look[0]*up[2];
    side[2] = look[0]*up[1] - look[1]*up[0];

    double ls = sqrt( side[0]*side[0] + side[1]*side[1] + side[2]*side[2] );

    side[0] /= ls;
    side[1] /= ls;
    side[2] /= ls;

    double lu = sqrt( up[0]*up[0] + up[1]*up[1] + up[2]*up[2] );

    up[0] = side[1]*look[2] - side[2]*look[1];
    up[1] = side[2]*look[0] - side[0]*look[2];
    up[2] = side[0]*look[1] - side[1]*look[0];

    double m[3][3];

    m[0][0] = side[0]; m[1][0] = side[1]; m[2][0] = side[2];
    m[0][1] = up[0];   m[1][1] = up[1];   m[2][1] = up[2];
    m[0][2] = -look[0]; m[1][2] = -look[1]; m[2][2] = -look[2];

    double tq[4];
    tq[0] = 1 + m[0][0] + m[1][1] - look[2];
    tq[1] = 1 + m[0][0] - m[1][1] + look[2];
    tq[2] = 1 - m[0][0] + m[1][1] + look[2];
    tq[3] = 1 - m[0][0] - m[1][1] - look[2];

    unsigned int j = 0;

    for( unsigned int i=1; i<4; ++i )
        if( tq[i] > tq[j] )
            j = i;

    double s = 0.5*sqrt( 1.0/tq[j] );

    if( j == 0 )
    {
        m_cur.r[0] = s*(m[1][2]-m[2][1]);
        m_cur.r[1] = s*(m[2][0]-m[0][2]);
        m_cur.r[2] = s*(m[0][1]-m[1][0]);
        m_cur.r[3] = s*tq[0];
    }
    else if( j == 1 )
    {
        m_cur.r[0] = s*tq[1];
        m_cur.r[1] = s*(m[0][1]+m[1][0]);
        m_cur.r[2] = s*(m[2][0]+m[0][2]);
        m_cur.r[3] = s*(m[1][2]-m[2][1]);
    }
    else if( j == 2 )
    {
        m_cur.r[0] = s*(m[0][1]+m[1][0]);
        m_cur.r[1] = s*tq[2];
        m_cur.r[2] = s*(m[1][2]+m[2][1]);
        m_cur.r[3] = s*(m[2][0]-m[0][2]);
    }
    else 
    {
        m_cur.r[0] = s*(m[2][0]+m[0][2]);
        m_cur.r[1] = s*(m[1][2]+m[2][1]);
        m_cur.r[2] = s*tq[3];
        m_cur.r[3] = s*(m[0][1]-m[1][0]);
    }

    m_home = m_cur;
    update();
}

// --------------------------------------------------------------------------

void trackball::home()
{
    m_cur = m_home;
    update();
}



