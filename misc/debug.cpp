#include "debug.hpp"
#include <iostream>

#if defined NDEBUG
_null_stream dout;
#else
std::ostream& dout = std::cout;
#endif
