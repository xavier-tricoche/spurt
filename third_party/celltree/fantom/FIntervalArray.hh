#ifndef F_INTERVALARRAY_HH
#define F_INTERVALARRAY_HH

#include <vector>
#include "stdAliases.hh"

/**
 *class which implements fast access to interval indices.
 *(helper class for FmmappedArray)
 *it is required, that the minimum interval size
 *is not too small (except the size od the last one )
 *( because an array of size (sum of sizes) / (min size) will be built) )
*/



class FIntervalArray
{


public:
  /** 
   * constructor.
   * contructs a number of consecutive 
   * intervals with size sizes[i]
   * beginning at 0
   *
   *\param borders 
   *start indices of the arrays
   */
  FIntervalArray(const std::vector<positive> & sizes);


  /** returns the number of the interval in which
   * i lies
   *(the intervals are separated by the entries in 
   * the array borders given in the constructor)
   */
  positive getInterval(positive i) const;

  inline positive operator [](positive i) const
  {return getInterval(i);}

private:

  struct hashIt{
    positive interval;
    positive border;
    /** operator used in getInterval for std::upper_bound */
    friend bool operator < (positive i,const hashIt&h)
    {return i < h.border;}
  };


  std::vector< hashIt > table;
  
  positive sumSize;
  int pMinSize;

};


#ifndef OUTLINE
#include "FIntervalArray.icc"
#endif

#endif // F_INTERVALARRAY_HH




