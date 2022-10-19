//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: sockutil.hh,v $
// Language:  C++
// Date:      $Date: 2003/07/01 16:19:38 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __sockutil_hh
#define __sockutil_hh

#ifdef WIN32
#error Sockutil not yet portet to Win32
#else
#include <iostream>
#include <string>
#include <fcntl.h>
#include <vector>


int listenTcpStream( unsigned int port );

void clearPendingConnections(int socket_fd);

int acceptSocket( int socket_fd );

int connectTcpStream(const char*hostname, unsigned int port );



class tcp_socket
{
    int      _fd;
    fd_set   _fd_set;

    void set_nonblocking();
    void set_blocking();

public:

    // construct with a given filedescriptor, preferably an already open socket
    tcp_socket( int fd );

    // construct by connecting a socket to hostname/port
    tcp_socket( const char* hostname, unsigned int port );

    ~tcp_socket();

    bool     write( void* ptr, size_t size );
    size_t   read( void* ptr, size_t size );

    size_t   try_read( void* ptr, size_t size );
    bool     wait_for_input();

    // listen for an incoming connection on port
    static tcp_socket* listen( unsigned int port );

    static unsigned int wait_for_input( const std::vector<tcp_socket*>& sockets );

    int fd() { return _fd; }
};

class socketstreambuf: public std::streambuf
{
    typedef char         char_type;
    typedef int          int_type;
    typedef tcp_socket   io_type;

    char_type        *_get_buffer, *_put_buffer;
    size_t            _buffer_size;
    io_type          *_io;

public:

    explicit socketstreambuf( io_type *io );
    ~socketstreambuf();

protected:

    // reimplemented from std::streambuf
    virtual int_type sync();
    virtual int_type underflow();
    virtual int_type overflow( int_type ch = EOF );

//     virtual streamsize xsgetn( char_type * s, streamsize n );
};
#endif
#endif // __sockutil_h
