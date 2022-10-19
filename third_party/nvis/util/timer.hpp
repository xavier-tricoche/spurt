#ifndef __timer_hpp
#define __timer_hpp

#include <sys/time.h>
#include <ctime>

namespace nvis
{
    // --- a wall clock timer modeled after boost::timer ---

    class timer
    {
    public:

	timer() 
	{
	    gettimeofday( &_start_time, 0 );
	}
	
	void restart() 
	{ 
	    gettimeofday( &_start_time, 0 );
	}
	
	double elapsed() const
	{ 
	    timeval cur;
	    gettimeofday( &cur, 0 );
	    
	    return timeval_diff( _start_time, cur );
	}

    private:
	
	static double timeval_diff( const timeval& t1,
				    const timeval& t2 )
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
	
	timeval _start_time;
    };

} // namespace nvis

#endif //__timer_hpp
