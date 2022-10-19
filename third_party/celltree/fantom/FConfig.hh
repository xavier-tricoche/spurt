#ifndef FConfig_hh
#define FConfig_hh

#include <string>
#include <map>
#include <ostream>

/** 
 * Global configuration object for FAnToM
 * 
 * Currently only string data can be "stored" here :(
 */
struct FConfig
{
    //! get (and create if needed) the global instance
    static FConfig *theConfig();

    //! add a (new) entry or overwrite an existing one
    void addEntry(const std::string &name, const std::string &entry)
    {
        string_entries[name] = entry;
    }

    //! add a (new) bool entry or overwrite an existing one
    void addBoolEntry(const std::string &name, const bool &entry)
    {
        bool_entries[name] = entry;
    }
    
    
    //! get an entry, returns "" if no entry is stored and appends
    //! an empty entry to the configuration
    const std::string &getEntry(const std::string &name) const
    {
      if ( string_entries.find( name ) == string_entries.end() )
      {
        setFromEnvironment( name );
      }
      return string_entries[name];
    }

    bool setFromEnvironment( const std::string &name ) const; //< thingy is mutable so we can do this

    //! get an entry, returns "" if no entry is stored and appends
    //! an empty entry to the configuration
    const bool &getBoolEntry(const std::string &name) const
    {
        return bool_entries[name];
    }

    //! get a reference to the internal map. use with care because
    //! there is no synchronization!
    const std::map<std::string, std::string>& getMap() const
    {
        return string_entries;
    }

    //! get a reference to the internal map. use with care because
    //! there is no synchronization!
    const std::map<std::string, bool>& getBoolMap() const
    {
        return bool_entries;
    }

    ~FConfig()
    {
        _theConfig = 0;
    };

    void dump(std::ostream& o)
    {
      for(std::map<std::string, std::string>::const_iterator it
          = string_entries.begin(); it != string_entries.end(); ++it)
      {
        o << it->first << " -> " << it->second << std::endl;
      }
      for(std::map<std::string, bool>::const_iterator it
          = bool_entries.begin(); it != bool_entries.end(); ++it)
      {
        o << it->first << " -> " << it->second << std::endl;
      }
    }
    
    private:
    mutable std::map<std::string, std::string> string_entries;
    mutable std::map<std::string, bool> bool_entries;
    static FConfig *_theConfig;

    FConfig(){}; // private constructor
};

#endif
