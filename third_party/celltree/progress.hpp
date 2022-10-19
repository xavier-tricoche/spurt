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

#ifndef __progress_hpp
#define __progress_hpp

#include <cstdio>
#include <cmath>

namespace CellTree {

class progress
{
    size_t       m_goal;
    size_t       m_current;
    size_t       m_check;
    float        m_last;
    const char*  m_msg;
    
public:
    
    progress() : m_goal(0), m_msg(0)
    {
    }
    
    progress( const char* msg, size_t count ) : m_goal(0)
    {
        restart( msg, count );
    }

    ~progress()
    {
        finished();
    }    
    
        
    void finished()
    {
        if( !m_goal )
            return;
        
        fprintf( stderr, "\r%50s\r", " " );
        fflush( stderr );
        
        m_goal = 0; 
    }    
    
    void restart( const char* msg, size_t count )
    {
        finished();
        
        m_msg = msg;
        m_goal = count;
        m_check = count / 10101 + 1;
        m_current = 0;
        m_last = 0;
        
        output();
    }

    void output()
    {
        if( m_current == m_goal )
        {
            finished();
            return;
        }
        
        if( (m_current % m_check) )
            return;
        
        fprintf( stderr, "\r%s: %.2f%%", m_msg, (100.0f * m_current) / m_goal );
        fflush( stderr );
    }
    
    double elapsed() const
    {
		return 0.0;
    }
    
    void operator++() 
    {
        #pragma omp atomic
        ++m_current;
        output();
    }
    
    void operator++( int )
    {
        ++(*this);
    }
    
    void operator+=( size_t inc )
    {
        #pragma omp atomic
        m_current += inc;
        output();
    }
};

} // namespace CellTree

#endif // __progress_hpp
