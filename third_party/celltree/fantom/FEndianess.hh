#ifndef FEndianess_hh
#define FEndianess_hh

#ifndef NODEBUG
#include "eassert.hh"
#endif


#if defined (__DARWIN_X11__)
// BYTE_ORDER magically defined somewhere else
#elif defined (__APPLE__)
#include <machine/endian.h>
#elif defined (WIN32)
#include <windows.h>
#define	LITTLE_ENDIAN	1234
#define	BIG_ENDIAN	4321
#define	PDP_ENDIAN	3412
#define BYTE_ORDER LITTLE_ENDIAN
#else // linux
#include <endian.h>
#endif

#if BYTE_ORDER == PDP_ENDIAN
#error "PDP Endian code not implemeted in FAnToM"
#endif

#if BYTE_ORDER == BIG_ENDIAN
#define FBIG_ENDIAN
#define BE(a) a
#define LE(a)
#endif
#if BYTE_ORDER == LITTLE_ENDIAN
#define FLITTLE_ENDIAN
#define BE(a)
#define LE(a) a
#endif


namespace FEndianess
{
#define SWAP(a,b) \
  {char tmp = (a); (a)=(b); (b)=tmp;}
  
  inline void Swap( char* data, unsigned int size )
  {
#ifndef NODEBUG
    eassert( (size & 0x01) == 0 );
#endif
    for(unsigned int i=0;i<(size>>1); ++i)
      SWAP(data[i], data[size-i-1]);
  }

  inline void SwapRange( char* data, unsigned int dataSize, unsigned int nbChars )
  {
#ifndef NODEBUG
    eassert( (dataSize & 0x01) == 0);
    eassert( nbChars % dataSize == 0);
#endif
    unsigned int loops = nbChars / dataSize;
    for(unsigned int i=0; i< loops; ++i)
      Swap(&data[i*dataSize], dataSize);
  }
#undef SWAP


  inline bool isBigEndian()
  {
    return BYTE_ORDER == BIG_ENDIAN;
  }

  inline bool isLittleEndian()
  {
    return BYTE_ORDER == LITTLE_ENDIAN;
  }

  /** 
   * when saving data, convert host to big endian format
   */
  inline void HostToBigEndian(char*& LE(data), int LE(size))
  {
#if BYTE_ORDER == BIG_ENDIAN
    return;
#endif
#if BYTE_ORDER == LITTLE_ENDIAN
    Swap( data, size );
#endif
  }

  /**
   * when saving data, convert host to little endian format
   */
  inline void HostToLittleEndian(char* & BE(data), int BE(size) )
  {
#if BYTE_ORDER == LITTLE_ENDIAN
    return;
#endif
#if BYTE_ORDER == BIG_ENDIAN
    Swap( data, size );
#endif
  }

  /**
   * when loading data, convert host to little endian format
   */
  inline void LittleEndianToHost(char* &BE(data), int BE(size) )
  {
#if BYTE_ORDER == LITTLE_ENDIAN
    return;
#endif
#if BYTE_ORDER == BIG_ENDIAN
    Swap(data, size);
#endif
  }
  
  /**
   * when loading data, convert host to little endian format
   */
  inline void BigEndianToHost(char* &LE(data), int LE(size) )
  {
#if BYTE_ORDER == LITTLE_ENDIAN
    Swap(data, size);
#endif
#if BYTE_ORDER == BIG_ENDIAN
    return;
#endif
  }













  /** 
   * when saving data, convert host to big endian format
   */
  inline void HostToBigEndianRange(char* &LE(data), int LE(size), int LE(length))
  {
#if BYTE_ORDER == BIG_ENDIAN
    return;
#endif
#if BYTE_ORDER == LITTLE_ENDIAN
    SwapRange( data, size, length );
#endif
  }

  /**
   * when saving data, convert host to little endian format
   */
  inline void HostToLittleEndianRange(char* &BE(data), int BE(size), int BE(length) )
  {
#if BYTE_ORDER == LITTLE_ENDIAN
    return;
#endif
#if BYTE_ORDER == BIG_ENDIAN
    SwapRange( data, size, length );
#endif
  }

  /**
   * when loading data, convert host to little endian format
   */
  inline void LittleEndianToHostRange(char* & BE(data), int BE(size), int BE(length) )
  {
#if BYTE_ORDER == LITTLE_ENDIAN
    return;
#endif
#if BYTE_ORDER == BIG_ENDIAN
    SwapRange(data, size, length);
#endif
  }
  
  /**
   * when loading data, convert host to little endian format
   */
  inline void BigEndianToHostRange(char* &LE(data), int LE(size), int LE(length) )
  {
#if BYTE_ORDER == LITTLE_ENDIAN
    SwapRange(data, size, length);
#endif
#if BYTE_ORDER == BIG_ENDIAN
    return;
#endif
  }









  
  struct EndianessHelper
  {
    enum Endianess
    {
      BigEndian,
      LittleEndian
    } endianess;
   
    /**
     * \param e endianess of file
     */
    EndianessHelper(Endianess e ): endianess(e)
    { }

    inline void toHost( char* data, int size) const
    {
      if(endianess==BigEndian) BigEndianToHost( data, size );
      else LittleEndianToHost( data, size );
    }

    inline void toHostRange( char* data, int size, int length ) const
    {
      if(endianess==BigEndian) BigEndianToHostRange( data, size, length );
      else LittleEndianToHostRange( data, size, length );
    }

    inline void toFile( char* data, int size) const
    {
      if(endianess==BigEndian) HostToBigEndian( data, size );
      else HostToLittleEndian( data, size );
    }

    inline void toFile( char* data, int size, int length) const
    {
      if(endianess==BigEndian) HostToBigEndianRange( data, size, length );
      else HostToLittleEndianRange( data, size, length );
    }

    void setFileEndianess( Endianess endianess )
    {
      this->endianess = endianess;
    }
  };
}

// free the defines that are used here
#undef BE
#undef LE
#endif

