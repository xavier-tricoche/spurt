//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMCycle.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:33:33 $
// Author:    $Author: garth $
// Version:   $Revision: 1.21 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMCycle_hh
#define __FAMCycle_hh

#include "FObject.hh"
#include "FVector.hh"

#include "FAMSingularPoint.hh"

#include <vector>
#include <list>

/** 
 * Structure to handle Cycles. There is a singular point that is contained within the cycle
 */
class FAMCycle : public FObject
{
public:

  /** 
   *\par Description:
   *Constructor: returns an empty FAMCycle.
   */
  FAMCycle();

  /** 
   *\par Description:
   * copy constructor
   *\param cyc
   * Cycle to copy
   */
  FAMCycle(const FAMCycle& c);

  /** 
   *\par Description:
   * Destructor: cleans up.
   */
  virtual ~FAMCycle();

  /**
   *\par Description:
   * Get (one of) the singular point(s) which leads into this cycle.
   */
  void getSingularity (const FAMSingularPoint *s) const;

  /** 
   *\par Description:
   * Set (one of) the singular point(s) which leads into this cycle.
   */
  void setSingularity (const FAMSingularPoint &s);

  /** 
   *\par Description:
   * Get a pointer to the approximation of the cycle.
   *\param ap
   * The pointer to the objects internal std::vector. You are not allowed
   * to modify or delete the memory that ap points to.
   */
  void getApproximation (std::vector<FVector>* &ap) const;

  /** 
   *\par Description:
   * Set a pointer to the approximation of the cycle.
   *\param ap
   * Pointer to a std::vector containing the information. After calling, the object
   * ap is property of FAMCycle, it will be deleted here.
   */
  void setApproximation (std::vector<FVector> *ap);

  /** 
   *\par Description:
   * Get a pointer to the std::list of cells which the cycle passes through.
   */
  void getCycleCells (std::list<FIndex>* &cells) const;

  /** 
   *\par Description:
   * Set a pointer to the std::list of cells which the cycle passes through.
   *\param cells
   * Pointer to a std::vector containing the information. After calling, the object
   * ap is property of FAMCycle, it will be deleted here.
   */
  void setCycleCells (std::list<FIndex> *cells);

  /**
   *\par Description:
   * Returns the class name.
   */
  virtual const FString& getClassName() const;

  /**
   *\par Description:
   * operator =.
   */
  FAMCycle& operator= (const FAMCycle& c);

private:

  FAMSingularPoint sing; //! (One of) the singular point(s) that lead(s)
                         //! into the cycle.
  std::vector<FVector> *aprox;
  std::list<FIndex> *cycleCells;
};

//===========================================================================

#endif // __FAMCycle_hh
