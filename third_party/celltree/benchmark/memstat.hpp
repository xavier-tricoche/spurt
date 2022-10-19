#ifndef __memstat_hpp
#define __memstat_hpp

#include <cstdlib>

// A somewhat half-assed attempt at a memory allocation tracker...
// works correctly on at least one machine, if that means anything :-/

// Reset memory allocation counters.
void   memstat_reset();

// Return the number of bytes allocated since the 
// last reset.
std::size_t memstat_cur();

// Return the maximum number of bytes allocated 
// at any time between now and the last reset.
std::size_t memstat_max(); 

#endif // __memstat_hpp