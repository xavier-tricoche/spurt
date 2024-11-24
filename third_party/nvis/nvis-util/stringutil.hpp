#ifndef __stringutil_hpp
#define __stringutil_hpp

#include <string>
#include <sstream>
#include <map>

// ------------------------------------------------------------------------

template<typename T> std::string to_string( const T& t )
{
    std::ostringstream out;
    out << t;
    return out.str();
}

// ------------------------------------------------------------------------

template<typename T> T from_string( const std::string& s )
{
    T t;
    std::istringstream(s) >> t;
    return t;
}

template<typename T> bool from_string( const std::string& s, T& t )
{
    std::istringstream in(s);
    in >> t;
    return !in.fail();
}

template<typename T> bool from_string_n( const std::string& s, T t[], std::size_t n )
{
    std::istringstream in(s);

    for( std::size_t i=0; i<n && in.good(); ++i )
        in >> t[i];

    return !in.fail();
}

// ------------------------------------------------------------------------

template<typename T> class string_map
{
public:
    typedef std::map< std::string, T >  map_type;

    string_map& add( const std::string& name, const T& m )
    {
        map_.insert( typename map_type::value_type( name, m ) );
        return *this;
    }

    bool map( const std::string& name, T& t ) const
    {
        typename map_type::const_iterator i = map_.find( name );
        
        if( i == map_.end() )
            return false;

        t = i->second;
        return true;
    }

private:
    map_type map_;
};

// ------------------------------------------------------------------------

#endif // __stringutil_hpp
