#ifndef FDirectory_hh
#define FDirectory_hh

#include "FString.hh"

/**
 * class that helps to manage directories and scan for files
 */
class FDirectory
{
  public:
    /**
     * constructor
     * \param dirname name of the directory to open
     */
    FDirectory( const char* dirname );
    
    ~FDirectory();
    
    /**
     * returns true if the current structure is valid
     */
    bool isValid();

    /**
     * get the number of files present in the directory
     */
    unsigned int getNbFiles() const;

    /**
     * returns the name of the i-th file in the directory
     */
    FString getFile(unsigned int i) const;

    /**
     * returns the full path of the i-th file in the directory
     */
    FString getFullFilePath(unsigned int i) const;

    /**
     * might fail on some platforms or filesystems, not tested yet.
     */
    bool isDirectory(unsigned int i) const;

    struct FDirectoryPrivateData;
  protected:
    //int idx;
    FDirectoryPrivateData *pd;
};

#endif

