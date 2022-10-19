#ifndef __fps_timer_hpp
#define __fps_timer_hpp

#include <wall_timer.hpp>

class fps_timer
{
public:

    fps_timer()
    {
        m_frames = 0;
        m_rate   = 0;
        m_time   = 0;
    }

    void frame()
    {
        ++m_frames;

        if( m_timer.elapsed() < 0.2 )
            return;

        m_time = m_timer.restart();
        m_rate = double(m_frames)/m_time;
        m_frames = 0;
    }

    double rate() const
    {
        return m_rate;
    }

private:

    wall_timer   m_timer;
    unsigned int m_frames;
    double       m_rate;
    double       m_time;
};

#endif // __fps_timer_hpp
