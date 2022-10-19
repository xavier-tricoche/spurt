#ifndef FObserver_hh
#define FObserver_hh


class FObservable;
class FObserver
{
  public:
  FObserver(){}
  virtual ~FObserver(){}


  /**
   * called while holding a mutex so do not call
   * addObserver or removeObsererver here.
   * 
   * thread context is the thread that has triggered the updateObservers 
   * call.
   */
  virtual void update( FObservable* ) = 0;
};

#endif

