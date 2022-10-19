#ifndef FSynchronizedQueue_hh
#define FSynchronizedQueue_hh

#include "FMutex.hh"
#include "FThread.hh"

#include <queue>

// synced queue, will go in its own header file if everything is ready for it
template< typename T >
class FSynchronizedQueue
{
    public:
        FSynchronizedQueue(){}
        ~FSynchronizedQueue(){};

        T pop()
        {
            FAutoLock l(mutex);

            while( queue.empty() )
                condition.wait( l.mutex );

            T t = queue.front();
            queue.pop();
            return t;
        }
        
        bool try_pop( T& t )
        {
            FAutoLock l(mutex);
            if( queue.empty() )
                return false;
            t= queue.front();
            queue.pop();
            return true;
        }

        void push( const T& t)
        {
            FAutoLock l(mutex);
            queue.push( t );

            condition.signal();
        }

    private:
        FMutex mutex;
        FCondition condition;
        std::queue<T> queue;
};

#endif
