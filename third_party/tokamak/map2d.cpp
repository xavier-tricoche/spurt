#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

#include "tokamak_nimrod.hpp"
#include "poincare_map.hpp"
#include "standard_map.hpp"
#include "map2d.hpp"

map2d* map2d::load( const std::string& mapdesc )
{
    using namespace boost;
    using namespace boost::algorithm;
    
    std::vector<std::string> parts;
    split( parts, mapdesc, is_any_of(":") );

    if( parts.size() < 2 )
        throw std::runtime_error( "invalid map description" );

    if( parts[0] == "std" )
        return new standard_map( lexical_cast<double>(parts[1]) );
    else if( parts[0] == "nimrod" || parts[0] == "nim" )
    {
        if( parts.size() < 3 )
            throw std::runtime_error( "invalid map description" );
        
        tokamak_field* field = new tokamak_nimrod( parts[1], parts[2], parts.size()<4 );
        return new poincare_map( field );
    }
    else
        throw std::runtime_error( "unknown map type" );
}
