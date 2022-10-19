#ifndef FDynamicPlugin_hh
#define FDynamicPlugin_hh

#include "FObject.hh"

#ifdef WIN32
#include <windows.h>
#endif
/**
 * Class as an object oriented wrapper to
 * dlfcn and similar code relating dynamic
 * plugin handling
 */
class FDynamicPlugin : public FObject
{

  public:
    // define basic types to handlers
#ifdef WIN32
  typedef HINSTANCE Handle;
#else
  typedef void* Handle;
#endif
  typedef void* FuncAddr;
  static const FuncAddr InvalidFuncAddr;
 
  /**
   * destructor, closes the library
   */
  virtual ~FDynamicPlugin();

  /**
   * return the address of symbol in the loaded library
   * returns InvalidFuncAddr on error
   */
  FuncAddr getFuncAddress( const FString& symbol ) const;
  
  /**
   * Load the library with pathname libraryName
   * returns the reference to FDynamicPlugin created
   * returns 0 on error.
   */
  static FDynamicPlugin* loadLibrary( const std::string& libraryName );

  //! overloaded FObject
  virtual const FString& getClassName() const;

  virtual bool isValid() const;
 
  FString getMessage() const;
  
  protected:
 /**
   *
   */
  FDynamicPlugin( const std::string& libraryName );

  private:
  /**
   * constructor
   */
  FDynamicPlugin( Handle h );
  
   /**
   * handle pointing to the library
   */
  Handle handle;

  private:
  mutable std::string msg;
};

#endif

