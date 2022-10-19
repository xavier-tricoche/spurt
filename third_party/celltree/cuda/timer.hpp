#ifndef __cuda_timer_hpp
#define __cuda_timer_hpp

#include <cuda.h>

class event_timer
{
    cudaEvent_t start, stop;
    
public:
    
    event_timer()
    {
        cudaEventCreate( &start );
        cudaEventCreate( &stop );
    }
    
    ~event_timer()
    {
        cudaEventDestroy( start );
        cudaEventDestroy( stop );
    }
    
    void tic()
    {
        cudaEventRecord( start, 0 );
    }
    
    void toc()
    {
        cudaEventRecord( stop, 0 );
    }
    
    float elapsed()
    {
        cudaEventSynchronize( stop );
        
        float interval;
        cudaEventElapsedTime( &interval, start, stop );
        
        return interval;
    }
};

#endif // __cuda_timer_hpp