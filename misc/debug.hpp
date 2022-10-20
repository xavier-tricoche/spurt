#ifndef __debug_hpp__
#define __debug_hpp__

#include <string>
#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/stream.hpp>

namespace spurt {
    
typedef boost::iostreams::basic_null_sink<char> null_sink;
typedef boost::iostreams::stream<nullsink>      null_stream;

}

#endif // __debug_hpp__
