#include <cstdlib>
#include <cstdio>
#include <new>

using namespace std;

// -------------------------------------------------------------------------

std::size_t __memstat_allocated_current = 0;
std::size_t __memstat_allocated_max     = 0;

void memstat_reset()
{
    #pragma omp single 
    {
        __memstat_allocated_current = 0;
        __memstat_allocated_max     = 0;
    }
}

std::size_t memstat_max()
{
    return __memstat_allocated_max;
}

std::size_t memstat_cur()
{
    return __memstat_allocated_current;
}
    
// -------------------------------------------------------------------------

void *operator new(std::size_t sz) throw (std::bad_alloc)
{
    unsigned char *p = (unsigned char*)malloc( sz + sizeof(size_t) );

    if( p == 0 )
        throw bad_alloc();

    *(size_t*)p = sz;

    #pragma omp atomic
    __memstat_allocated_current += sz;

    if( __memstat_allocated_current > __memstat_allocated_max )
    {
        #pragma omp atomic
        __memstat_allocated_max = __memstat_allocated_current;
    }

    return p + sizeof(size_t);
}

// -------------------------------------------------------------------------

void operator delete(void* p) throw ()
{
    if( p )
    {
        size_t* ps = (size_t*)p;
        --ps;

        #pragma omp atomic
        __memstat_allocated_current -= *ps;

        std::free( ps );
    }
}

// -------------------------------------------------------------------------
