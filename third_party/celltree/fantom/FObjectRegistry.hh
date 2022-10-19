#ifndef FObjectRegistry_hh
#define FObjectRegistry_hh

#include "FString.hh"
#include "FMutex.hh"
#include <iosfwd>
#include <list>

class FObject;
class FObjectRegistry
{
  public:
    ~FObjectRegistry(){}
    void dump(std::ostream &os) const;
    FString getDump() const;

    FObject* getObject( int id ) const;


    static FObjectRegistry* getRegistry()
    {
      return &registry;
    }

    size_t size() const { return objects.size();}

  protected:
    friend class FObject;

    void registerObject(FObject* object);
    void unregisterObject( FObject* object);

  private:
    FObjectRegistry(){}

    mutable FMutex mutex;
    std::list<FObject*> objects; 


    static FObjectRegistry registry;
};

#endif
