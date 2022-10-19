#ifndef __progress_hpp
#define __progress_hpp

#include <cstdio>
#include <cmath>

#include <wall_timer.hpp>

namespace celltree {

class progress
{
    size_t       m_goal;
    size_t       m_current;
    size_t       m_check;
    float        m_last;
    const char*  m_msg;
    wall_timer   m_timer;
    
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
        m_timer.restart();
        
        output();
    }

    void output()
    {
        if( m_current == m_goal )
        {
            finished();
            return;
        }
        
        if( (m_current % m_check) || m_timer.elapsed() - m_last < 0.1 )
            return;
        
        fprintf( stderr, "\r%s: %.2f%%", m_msg, (100.0f * m_current) / m_goal );
        fflush( stderr );
        
        m_last = m_timer.elapsed();
    }
    
    double elapsed() const
    {
        return m_timer.elapsed();
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

} // namespace celltree

#endif // __progress_hpp
