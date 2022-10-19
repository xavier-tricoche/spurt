//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FSerializableObject.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/30 13:01:56 $
// Author:    $Author: garth $
// Version:   $Revision: 1.5 $
//
//--------------------------------------------------------------------------- 

#ifndef __FSerializableObject_hh
#define __FSerializableObject_hh

#include <typeinfo>
#include <iosfwd>
#include <string>
#include <stdexcept>

// ---------------------------------------------------------------

struct FSerializableBase
{
    virtual void serialize( std::ostream& out ) const = 0;

    virtual const char* _type_id() const = 0;

    virtual ~FSerializableBase() {};
};

// ---------------------------------------------------------------

typedef FSerializableBase* (*_factory)( std::istream& );

void _register_factory( const std::string&, _factory );

FSerializableBase *_create_type( const std::string&, std::istream& );

// ---------------------------------------------------------------

template<typename S> struct _generic_factory
{
    static FSerializableBase *create( std::istream& in )
    {
	return S::rebuild( in );
    }

    _generic_factory()
    {
	_register_factory( S::_static_type_id(), create );
    }

    void _noop(){} // to make the compiler happy

    static _generic_factory<S> _this_factory;

};

template<typename S> _generic_factory<S> _generic_factory<S>::_this_factory;

// ---------------------------------------------------------------

// the "_this_factory" reference has been there to guarantee, that
// the object is referenced, thus constructed. Some compiler complained,
// that that statement has no effect (even though it did, what we wanted
// it to do), so I introduced the noop function that is called here, only
// for the purpose, that _this_factory is created and _register_factory(...)
// will be called.
#define MAKE_SERIALIZABLE( class )                                           \
struct _##class##_serializer:                                                \
    public _generic_factory<class>                                           \
{                                                                            \
    _##class##_serializer() { _this_factory._noop(); }                       \
};                                                                           \
public:                                                                      \
static const char* _static_type_id() { return #class; }                      \
virtual const char* _type_id() const { return _static_type_id(); }           \
private:

// ---------------------------------------------------------------

// high level methods 

FSerializableBase* rebuildObject( std::istream& );
void serializeObject( std::ostream&, FSerializableBase* );

template<typename T> T* rebuildObject( std::istream& in )
{
    return dynamic_cast<T*>( rebuildObject(in) );
}

#endif //__FSerializableObject_hh
