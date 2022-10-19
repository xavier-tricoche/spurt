//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSerializableObject.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/30 07:31:38 $
// Author:    $Author: garth $
// Version:   $Revision: 1.4 $
//
//--------------------------------------------------------------------------- 
//#define DEBUG_SERIALIZABLE
#include "FSerializableObject.hh"
#include "binio.hh"

#include <map>
#include <cassert>
#include <iostream>

typedef std::map< std::string, _factory > _factory_map_t;

static _factory_map_t *the_map = 0;

// ---------------------------------------------------------------

void _register_factory( const std::string& type, _factory f )
{
#ifdef DEBUG_SERIALIZABLE
  std::cout << "registering: " << type << std::endl;
#endif
    if( !the_map ) 
	the_map = new _factory_map_t;

    if( the_map->end() != the_map->find(type) )
    {
// 	std::cout << "already have factory for type " << type << "\n";
	return;
    }

//     std::cout << "registering factory for type " << type << "\n";
    the_map->operator[](type) = f;
}

// ---------------------------------------------------------------

FSerializableBase *_create_type( const std::string& type, std::istream& in )
{
    if( !the_map )
	return 0;

    _factory_map_t::iterator i;

    if( the_map->end() == ( i = the_map->find(type) ) )
	return 0;

    return (i->second)( in );
}

// ---------------------------------------------------------------

FSerializableBase* rebuildObject( std::istream& in )
{
    std::string type;
    binread( in, type );
#ifdef DEBUG_SERIALIZABLE
    std::cout << "rebuild: \"" << type << "\" (" << type.size() << ")" << std::endl;
#endif
    return  _create_type( type, in );
}

// ---------------------------------------------------------------

void serializeObject( std::ostream& out, FSerializableBase *obj )
{
#ifdef DEBUG_SERIALIZABLE
  std::cout << "serialize: \"" << obj->_type_id() << "\"." << std::endl;
#endif
    binwrite( out, obj->_type_id() );
    obj->serialize( out );
}
