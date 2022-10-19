#ifndef F_MMAPPED_ARRAY_HH
#define F_MMAPPED_ARRAY_HH

#include <string>
#include <vector>
#include <list>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include "stdAliases.hh"
#include "FIntervalArray.hh"

using namespace std;
using namespace boost;


class FIntervalArray;

//10 blocks around to protect
//#define PROT_SIZE (10<<12)

/**
 *class which represents a block-structured
 *array distributed over multiple files
 *and then mapped into memory via mmap
 */
class FmmappedArray
{

//-----------------------------------------------------------------------------



public:


  /**
   * constructor
   * builds up an array of blocks with size blocksize
   * distributed over different files
   * described by the vector of parts
   * 
   *\param fileNames
   * filenames of the different files constituting the array
   *
   *\param blocksPerFile
   * (desired)number of blocks to be present in each file
   *
   *\param blocksize
   *size of a block in bytes
   *
   *\param preferredMmapSize
   * preferred size for the mmapped memory regions;
   * the real size will be bigger
   *
   *\param blocksOverlap
   * percentage of overlap of blocks in one mmapped region
   *
   *\param readonly
   * if true,indicates that only read access is allowed
   */
  FmmappedArray( const vector<string> & fileNames, 
		 const vector<positive> & blocksPerFile,
		 positive blockSize, 
		 positive preferredMmapSize=128*1024, //128KB
		 positive blocksOverlap=10,
		 bool readonly=false
		 );
  ~FmmappedArray();

  /// 
  positive getBlockSize() const{return blockSize;}

  ///get total number of blocks in array
  positive getNbBlocks() const{return numBlocks;}
  
  ///
  void getFileNames( vector<string> &names) const;
  
  ///get file number of file in which block number i lies
  positive getFileForBlock(positive i) const;

  ///
  void getBlocksPerFile( vector<positive> &sizes) const;


  ///reads and writes a datum in part id i
  void syncNFS(positive partId);

  /**
   * Gives a pointer on an mmapped block with index i.
   *\param i
   * blocknumber
   *\retval numconsecutive
   * number of consecutively accessible blocks
   *(including the actual one)
   * at the returned address
   *\return
   * pointer on the block with index i
   */
  void * getBlock( positive & numconsecutive , positive i );


  static void setDeleteAtEnd(bool d=true)
  {deleteAtEnd=d;}

  ///helper classes defined in FmmappedArray.cc
  class FFileHandle;

  struct FmmappedRegion{
  
    //to hold the file handle 
    shared_ptr < FFileHandle > fileHandle;
    char * pStart;
    int size;


#ifdef PROT_SIZE
    void * pProtectUpper, * pProtectLower;
#endif

    //pointer on entry in regionMaps 
    //used to set it to 0 in destructor
    list<FmmappedRegion>::iterator * pIter;



    /// block interval in array which is mapped
    /// start index of blocks in this
    positive startBlockId;

    //number of blocks in this mmapped mem region
    positive numBlocks;

    // pointer on block with id startBlockId,
    // can be different from pstart
    // because of rounding
    char * startBlock; 

    void unmap() ;


  };

private:



  ///blocksize in which array will be accessed
  positive blockSize;

  ///total number of blocks in array
  positive numBlocks;

  ///memory page size of system
  positive pageSize;

  /// ~(pageSize-1)
  /// (has bits set to 0 which are lower than ln2(pagesize))
  positive notPageSizeM1;


  ///number of blocks in one mmapped region
  positive mmappedBlocks;
  
  ///log2(mmappedBlocks)
  int pMmappedBlocks;

  ///number of overlapping blocks
  positive mmapOverlap;



  typedef  list< FmmappedRegion > mapRegions;
  static mapRegions EmptyList;

  mapRegions::iterator lastMapped;


  //...........................................................................
  //  globally defined variables:
  //...........................................................................

  friend class FmmappedRegion;
  friend class FFileHandle;


  //list of mmapped memory regions
  static mapRegions actuallyMapped;

  //this is for efficiency:
  //here, unused elements of the list are stored for potential re-use
  static mapRegions unusedElems;


  static positive usedAdressSpace;
  static const positive maxAdressSpaceToUse = 1LU<<30; //1GB

  static positive nbOpenedFiles;
  static const positive maxNbOpenedFiles = 500;


  ///should files be deleted in destructor ?
  static bool deleteAtEnd;


  //...........................................................................



  ///gets the index of a mapped region in a file from
  ///the block index relative to the first block in the file
  positive getRegionIndex(positive blockIndex);


  struct part
  {
    part(){}
    ~part(){}

    string fileName;

    ///id returned by an open() command
    ///(in unistd.h)
    weak_ptr<FFileHandle> handle;

    ///returns handle if not empty,
    ///otherwise constructs new shared ptr
    void  getFileHandle(shared_ptr<FFileHandle> & sharedHandle);

    ///index of first block in file
    positive firstIndex;
 
    ///number of blocks in file
    positive numBlocks;


    //is this to be opened in rdonly mode ?
    bool readonly;
    

    /**
     *gives for every region
     *(a region has size mmappedBlocks)
     *an iterator into the actuallyMapped list
     */
    vector< mapRegions::iterator > regionMaps;
    
  };

  vector<part> parts;

  FIntervalArray partsMap;

  /// helper function for getBlock
  void mapNewBlock(  mapRegions::iterator & aIter ,part & actPart, positive regionId ) ;


};


positive greatestCommonDivisor(positive a,positive b);

#ifndef OUTLINE
#include "FmmappedArray.icc"
#endif


#endif // F_MMAPPED_ARRAY_HH
