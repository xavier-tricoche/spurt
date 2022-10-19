//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: parseutils.hh,v $
// Language:  C++
// Date:      $Date: 2003/09/11 09:00:35 $
// Author:    $Author: garth $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef __parseutils_hh
#define __parseutils_hh

#include "FException.hh"
#include <boost/spirit/core.hpp>
#include "mmapped_file.hh"

#include <sstream>
#include <string>
#include <stdexcept>

// this causes name clashes
//using namespace boost::spirit;

using boost::spirit::parser;
using boost::spirit::parse_info;
using boost::spirit::space_p;


//--------------------------------------------------------------------------- 

template<typename T> class set_actor
{
public:
    
    explicit set_actor( T& ref_, T val_ )
        : ref(ref_), val(val_) {}
    
    template<typename T2> void operator()(T2 const&) const 
    {
	ref = val; 
    }

    template<typename IterT> void operator()( IterT const&, 
					      IterT const& ) const
    {
	ref = val;
    }
private:
   
    T& ref;
    T  val;
};


template<typename T> const set_actor<T> set( T& t, T v )
{
    return set_actor<T>(t, v);
}


//--------------------------------------------------------------------------- 

// moved to exception class
//typedef std::runtime_error parse_failure;

class parse_helper
{
    mmapped_file file;
    const char *file_begin_;
    const char *begin_, *end_;

public:

    const char* get_pos() const
    {
      return begin_;
    }

    const char* end() const
    {
      return end_;
    }
    
    const char* begin() const
    {
      return file_begin_;
    }
    
    void goto_pos( const char* pos )
    {
      // this is an hard assert because
      // you should check before calling
      // this function for validy to be
      // able to show an appropriate error
      // message
      assert( pos >= file_begin_);
      assert( pos <= end_ );
      begin_ = pos;
    }      
    
    void goto_sol()
    {
	while( begin_<end_ && (*begin_ != '\n') )
	    begin_++;

	if( begin_ == end_ )
	    throw parse_failure( "parse failure: input exhausted (goto_sol)" );

	begin_++;
    }

    template<typename T> void parse( parser<T> const& r )
    {
	if( begin_ == end_ )
	    throw parse_failure( "parse failure: input exhausted (parse)" );

	parse_info<> info = boost::spirit::parse( begin_, end_, r, space_p );

	if( !info.hit )
	{
	    begin_ = info.stop - 1;
        std::ostringstream os;
        os << "parse failure at position " << (begin_ - file_begin_) << " before '" << std::string( begin_, std::min( begin_+20, end_ )) << ( begin_+20 >= end_ ? "(EOF)" : "" ) << "'.";
   	    throw parse_failure( os.str() );
        //"parse failure before '" +
		//                 std::string( begin_, begin_+10 ) + '\',' + " offset: ");
	}

	begin_ = info.stop-1;
    }

    template<typename T> bool try_parse( parser<T> const& r )
    {
	if( begin_ == end_ )
	    return false;

	parse_info<> info = boost::spirit::parse( begin_, end_, r, space_p );

	if( !info.hit )
	    return false;

	begin_ = info.stop-1;
	return true;
    }

    parse_helper( const std::string& path ) : 
	file(path), file_begin_( file.begin() ), begin_(file.begin()), end_(file.end())
    {
    }
};

#endif //__parseutils_hh
