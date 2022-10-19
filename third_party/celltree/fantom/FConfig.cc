#include "FConfig.hh"
#include "FStaticDeleter.hh"

#include <iostream>

FStaticDeleter<FConfig> staticConfigDeleter;

FConfig *FConfig::_theConfig=NULL;

FConfig *FConfig::theConfig()
{
    if(!_theConfig)
    {
        _theConfig = new FConfig;
        staticConfigDeleter.setObject(_theConfig);
    }
    return _theConfig;
}

bool FConfig::setFromEnvironment( const std::string& name ) const
{
   if ( getenv( name.c_str() ) == 0 )
   {
     std::cout << "WARNING: config entry " << name << " not found in environment." << std::endl;
     string_entries[ name ] == "";
     return false;
   }
   else
   {
     const char* str = getenv( name.c_str() );
     string_entries[ name ] = str;
     return true;
   }
}


