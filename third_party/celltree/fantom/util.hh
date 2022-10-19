#ifndef util_hh 
#define util_hh

#ifdef WIN32
char *basename (const char* name);

//! I strongly advise NOT to use this function as it is a bad wrapper around windows rand() which is said to be a bad random number generator for most cases
#ifndef drand48
double drand48();
#else
#error drand48 already defined
#endif
#endif

#endif

