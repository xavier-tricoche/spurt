#ifndef FStaticDeleter_hh
#define FStaticDeleter_hh

/**
 * This class is used to delete global pointer-type objects that are
 * created dynamically and last as long as the application is running.
 *
 * This makes it easy to delete them in a clean way using their destructor.
 * 
 * use example:
 * \code
 * class MyThing;
 * FStaticDeleter<MyThing> thingDeleter;
 * 
 * class MyThing
 * {
 *   public:
 *   static MyThing *getThing(){ 
 *      if(!thing) 
 *      {
 *          // we are not allowed to break between
 *          // these lines, otherwise thing won't
 *          // be deleted.
 *          
 *          // critical
 *          thing = new MyThing;
 *          thingDeleter.setObject(obj);
 *          // end critical
 *      }
 *      return thing;
 *   }
 *   ...
 *
 *   private:
 *   Mything(){...};
 * };
 * \endcode
 */
template< class Type > class FStaticDeleter
{
    public:
        FStaticDeleter() { deleteit = 0; }

        Type *setObject( Type *obj, bool isArray = false )
        {
            this->deleteit = obj;
            this->isArray  = isArray;
            return obj;
        }

        virtual void destructObject()
        {
            if(isArray)
                delete [] deleteit;
            else
                delete deleteit;
            deleteit = 0;
        }

        virtual ~FStaticDeleter()
        {
            destructObject();
        }

    private:
        Type *deleteit;
        bool isArray;
};      

#endif
