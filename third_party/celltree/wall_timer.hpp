/*************************************************************************
* celltree: boundary interval hierarchy tree-based solution for point 
*           location and interpolation in unstructured meshes. This
*           technique was first proposed in:
*           "Fast, Memory-Efficient Cell location in Unstructured Grids
*            for Visualization"
*                by Christop Garth and Kenneth I Joy, VisWeek 2011.
* Author: Christoph Garth
* Copyright (C) 2011 Christoph Garth
* All rights reserved
**************************************************************************/

#ifndef __wall_timer_hpp
#define __wall_timer_hpp

#include <sys/time.h>
#include <ctime>

// A wall clock timer, modeled loosely on boost::timer.

namespace CellTree {

class wall_timer
{
public:
    
    wall_timer()
    {
        restart();
    }
    
    double restart() 
    { 
        timeval current;
        gettimeofday( &current, 0 );
        
        double elapsed = delta( m_start, current );
        
        m_start = current;
        return elapsed;
    }
    
    double elapsed() const
    { 
        timeval current;
        gettimeofday( &current, 0 );
        
        return delta( m_start, current );
    }

private:

    static double delta( const timeval& t1, const timeval& t2 )
    {
        timeval r;

        r.tv_sec  = t2.tv_sec  - t1.tv_sec;
        r.tv_usec = t2.tv_usec - t1.tv_usec;

        if( r.tv_usec < 0 )
        {
            --r.tv_sec;
            r.tv_usec += 1000000;
        }

        return r.tv_sec + r.tv_usec / 1000000.0;
    }
    
    timeval m_start;
};

} // namespace CellTree

#endif // __wall_timer_hpp
