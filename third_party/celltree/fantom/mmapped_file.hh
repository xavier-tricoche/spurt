//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: mmapped_file.hh,v $
// Language:  C++
// Date:      $Date: 2003/09/10 13:12:10 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __mmapped_file_hh
#define __mmapped_file_hh

#include <string>

// --------------------------------------------------------------------------

class mmapped_file
{
    int fd;

    const char *ptr;
    size_t      length;

    mmapped_file( const mmapped_file& );


    /** Retruns true, iff mmap is supported on this platform
      otherwise returns false. mmap will then be simulated
      by loading the file into the main memory.
      */
    static bool isSupported();
    
public:

    mmapped_file( const std::string& path );
    ~mmapped_file();

    const char *begin() const;
    const char *end() const;
};

#endif //__mmapped_file_hh
