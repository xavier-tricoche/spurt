#ifndef FObservable_hh
#define FObservable_hh

#include "FObserver.hh"

class FObservable
{
  struct FObservablePrivateData;
  
  public:
  FObservable();
  virtual ~FObservable();
  void addObserver( FObserver* observer );
  void removeObserver( FObserver* observer );

  unsigned int nbObservers() const;
  protected:
  void updateObservers();

  private:
  FObservablePrivateData * pd;
};

#endif

