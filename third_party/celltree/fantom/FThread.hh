#ifndef FThread_hh
#define FThread_hh

/**
 * An environment independant thread class for FAnToM.
 *
 * All threads should be derived from this class and
 * have to overload the main() function which is 
 * executed by the newly created thread.
 *
 * To create a thread:
 * class MyThread : public FThread {...};
 *
 * MyThread t;
 * t.run(); //< we return immediately here
 * do_something_else();
 * t.join();//< block until thread has finished work
 * 
 */
class FThread
{
    public:
        FThread();
        virtual ~FThread();


        static const char* getThreadLibName();

        /** 
         * Spawns the new thread and returns immediately
         *
         * WARNING: Because of unknown scheduling, execution
         * of the newly created thread may have finished before 
         * we return!
         */
        void run();

        /**
         * Block until the thread has finished its execution
         */
        void join();

        /**
         * cancel the current thread.
         * WARNING: untested!
         */
        /*
        inline void cancel()
        {
#ifdef USE_PTHREADS
            pthread_cancel(tid);
#endif
        }
        */
    protected:
        /**
         * Main function.
         * This is executed by the new thread and does all the
         * work.
         */
        virtual void *main() = 0;

    public:
        struct FThreadPrivateData;
        friend struct FThreadPrivateData;
    private:
        FThreadPrivateData *pd;
};

#endif
