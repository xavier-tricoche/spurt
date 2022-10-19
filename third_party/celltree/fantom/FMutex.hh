#ifndef FMutex_hh
#define FMutex_hh
//#ifndef USE_QT_THREADS
//#define USE_QT_THREADS
//#endif

#if defined(USE_QT_THREADS) && defined( USE_PTHREADS)
#error do not define both USE_QT_THREADS and USE_PTHREADS
#endif

#if !defined(USE_QT_THREADS) && !defined(USE_PTHREADS)
#error do define one of USE_QT_THREADS or USE_PTHREADS
#endif


#ifdef USE_PTHREADS
#include <pthread.h>
#endif

#ifdef USE_QT_THREADS
#ifdef QT3
#include <qmutex.h>
#include <qthread.h>
#endif
#ifdef QT4
#include <QMutex>
#include <QThread>
#include <QWaitCondition>
#endif
#endif

/**
 * FMutex encapsulates a general mutex.
 * Currently it uses the pthread mutex mechanism,
 */
class FMutex
{
	public:
	FMutex();		//< Constructor, initializes mutex with 0
	~FMutex();		//< Destructor, destroys mutex
	
	void lock();	//< Locks mutex
	void unlock();	//< Unlocks mutex
	
	private:
	void operator=(const FMutex&){}; // not implemented

	protected:
#ifdef USE_PTHREADS
	pthread_mutex_t theMutex;
#endif
#ifdef USE_QT_THREADS
    QMutex theMutex;
#endif
    friend class FCondition;
};

inline void FMutex::lock()
{
#ifdef USE_PTHREADS
	pthread_mutex_lock(&(this->theMutex));
#endif
#ifdef USE_QT_THREADS
    theMutex.lock();
#endif
}

inline void FMutex::unlock()
{
#ifdef USE_PTHREADS
	pthread_mutex_unlock(&(this->theMutex));
#endif
#ifdef USE_QT_THREADS
    theMutex.unlock();
#endif
}


class FAutoLock 
{
    public:
        FAutoLock(FMutex &mutex): mutex(mutex)
        {
            mutex.lock();
        }

        ~FAutoLock()
        {
            mutex.unlock();
        }
        FMutex &mutex;
};


/**
 * The FAnToM Condition Variable object
 *
 * Semantics currently are those of pthreads.
 */
class FCondition
{
  public:
    FCondition();
    ~FCondition();

    void wait(FMutex& lock);
    void signal();
    void broadcast();

  private:
    void operator=(const FCondition&){}; // not implemented

  protected:
#ifdef USE_PTHREADS
    pthread_cond_t condition;
#endif
#ifdef USE_QT_THREADS
    QWaitCondition condition;
#endif
};

inline void FCondition::wait(FMutex &lock)
{
#ifdef USE_PTHREADS
    pthread_cond_wait( &condition, &lock.theMutex );
#endif
#ifdef USE_QT_THREADS
    condition.wait( &(lock.theMutex) );
#endif
}

inline void FCondition::signal()
{
#ifdef USE_PTHREADS
    pthread_cond_signal( &condition );
#endif
#ifdef USE_QT_THREADS
    condition.wakeOne();
#endif
}

inline void FCondition::broadcast()
{
#ifdef USE_PTHREADS
    pthread_cond_broadcast( &condition );
#endif
#ifdef USE_QT_THREADS
    condition.wakeAll();
#endif
}

#endif
