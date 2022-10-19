#ifndef __trackball_hpp
#define __trackball_hpp

#include <cmath>

class trackball
{
public:
    
    trackball();

    const double* matrix() const;

    enum Mode
    {
        NONE,
        ROTATE,
        ROTATE_Z,
        PAN,
        ZOOM,
        PUSH
    };

    void begin( Mode mode, float mx, float my );
    void move( float mx, float my );
    void end();

    void home( double ex, double ey, double ez,
               double cx, double cy, double cz,
               double ux, double uy, double uz );

    void home( double min[3], double max[3] );

    void home();

protected:

    void apply_trackball( const double* m0, const double* m1 );

    void rotate_vq( double* v, const double* q );
    void rotate_qq( double* q1, const double* q2 );
    void normalize( double* v );

    void update();

    struct position
    {
        double      r[4];
        double      c[3];
        double      d;
    };

    position     m_cur, m_home;
    double       m_m0[2], m_m1[2];
    double       m_mat[16];
    Mode         m_mode;

};

#endif // __trackball_hpp
