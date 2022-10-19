//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMElement.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:41 $
// Author:    $Author: garth $
// Version:   $Revision: 1.11 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMElement_hh
#define __FAMElement_hh
 
#include "FObject.hh"

/** 
 * Abstract base class for objects held by the analysis module data class.
 * The introduced functionality covers flags for emtpy and changed objects
 * and methods to save the internal state.
 */
class FAMElement : public FObject
{
public:
  /**
   * \par Description:
   *   Default Constructor.
   */
  FAMElement ();

  /**
   * Destructor
   */
  virtual ~FAMElement();

  /**
   * \par Description:
   *   Returns wether the object contains any data. This is done based
   *   on a flag that has to be set by the subclass accordingly.
   */
  bool isEmpty() const;

  /**
   * \par Description:
   *   Stuff needed by the FObject inheritance.
   */
  virtual const FString& getClassName() const = 0;

  /**
   * \par Description:
   * Saves the information held within to the file specified by fileNameBase
   * the FAMElement may use different files by adding chars to fileNameBase
   * \pre
   * the Element exists...
   * \post
   * the file has been created/changed to current values.
   * \exception
   * filesystem errors maybe
   */
  virtual void save ( const FString& fileNameBase ) = 0;

  /**
   * \par Description:
   * Loads the information held within the file so it has not to be calculated
   * again.
   * The FAMElement may use different files by adding chars to fileNameBase
   * \pre
   * the Element exists...
   * \post
   * the file has been created/changed to current values.
   * \exception
   * filesystem errors maybe
   */
  virtual void load ( const FString& fileNameBase ) = 0;
  
  /**
   * \par Description:
   * Set the changed bit of the data. This method is to be called whenever
   * the internal state of the object is modified.
   * \pre
   * none
   * \post
   * The change flag has been set.
   * \exception
   * none
   */
  void change (bool flag = true);

  /**
   * \par Description:
   * Get the changed bit of the data.
   * \return
   *    bool to determine iff data has changed
   */ 
  bool getChanged () const;
  
protected:
  bool empty;     //! Flag indicating if the object is in use.
  bool changed;   //! Flag indicating wether the object's internal state was changed.
};

//===========================================================================

#endif // __FAMElement_hh
