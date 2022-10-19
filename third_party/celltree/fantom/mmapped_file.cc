//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: mmapped_file.cc,v $
// Language:  C++
// Date:      $Date: 2003/09/10 13:12:10 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#include "mmapped_file.hh"
#include <cstring>

#ifdef NO_MMAPPED
#include <stdexcept>
#include <unistd.h>
#include <fcntl.h>
#include <errno.h>
#include <iostream>

mmapped_file::mmapped_file(const std::string& path )
{
  std::cout << "WARNING: Using non mmapped version of mmapped_file, which may be slow" << std::endl;

  fd = open( path.c_str(), O_RDONLY );
  if( fd < 0 )
	throw std::runtime_error( "unable to open " + path +
				  strerror(errno) );

  length = (size_t)lseek( fd, 0, SEEK_END );
  lseek(fd, 0, SEEK_SET);
  
  ptr = (const char*)malloc( length );

#ifndef WIN32
  // TODO: implement in loop, if size is larger?
  // Windows has no SSIZE_MAX defined here
  assert( length < SSIZE_MAX );
#endif
 
  // we have to cast the const away for writing of the data
  size_t s = read( fd, const_cast<char*>( ptr ), length );
  assert( s > 0 );
}

mmapped_file::~mmapped_file()
{
  // we have to cast the const away for freeing the data
  free(const_cast<char*>(ptr));
  if(fd >=0)
    close( fd );
  fd = -1;
}

bool mmapped_file::isSupported()
{ return false; }
#else

#include <stdexcept>

#include <unistd.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>

mmapped_file::mmapped_file( const std::string& path )
{
    fd = open( path.c_str(), O_RDONLY );
    
    if( fd < 0 )
	throw std::runtime_error( "unable to open " + path +
				  strerror(errno) );
    
    length = (size_t)lseek( fd, 0, SEEK_END );
    
    ptr = (const char*)mmap( 0, length, PROT_READ, 
			     MAP_PRIVATE, fd, 0 );
    
    if( !ptr )
    {
	close( fd );
	
	throw std::runtime_error( "unable to mmap " + path +
				  strerror(errno) );
    }
}

mmapped_file::~mmapped_file()
{
  if(ptr && fd >=0)
    close( fd );
}

bool mmapped_file::isSupported()
{ return true; }

#endif

const char *mmapped_file::begin() const
{
    return ptr;
}

const char *mmapped_file::end() const
{
    return ptr + length + 1;
}
