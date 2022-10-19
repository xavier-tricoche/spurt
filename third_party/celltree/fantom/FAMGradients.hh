//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAMGradients.hh,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:41 $
// Author:    $Author: garth $
// Version:   $Revision: 1.11 $
//
//--------------------------------------------------------------------------- 

#ifndef __FAMGradients_hh
#define __FAMGradients_hh
 
#include "FAMElement.hh"

#include <vector>
#include "FVector.hh"
/** 
 * Object that holds directions and strength of gradient pulsed
 * applied to DTI datasets. They are additional information needed
 * for raw MRI/DTI data to be interpretet.
 *
 * 
 */
class FAMGradients : public FAMElement
{
public:
  /**
   * \par Description:
   *   Default Constructor for empty gradients object
   */
  FAMGradients ();
  virtual ~FAMGradients();

  bool isEmpty() const;

  /**
   * \par Description:
   *   Stuff needed by the FObject inheritance.
   */
  virtual const FString& getClassName() const
  {
      static FString classname("FAMGradients");
      return classname;
  }

  /**
   * \par Description:
   * Saves the information held within to the file specified by fileNameBase
   * the FAMGradients may use different files by adding chars to fileNameBase
   * \pre
   * the Element exists...
   * \post
   * the file has been created/changed to current values.
   * \exception
   * filesystem errors maybe
   */
  virtual void save ( const FString& fileNameBase );

  /**
   * \par Description:
   * Loads the information held within the file so it has not to be calculated
   * again.
   * The FAMGradients may use different files by adding chars to fileNameBase
   * \pre
   * the Element exists...
   * \post
   * the file has been created/changed to current values.
   * \exception
   * filesystem errors maybe
   */
  virtual void load ( const FString& fileNameBase );
  
  void change (bool flag = true);

//  bool getChanged () const;

  const std::vector<FVector>& getGradients() const { return gradients;}
  void setGradients(const std::vector<FVector>& g) { gradients = g;}
  
protected:
//  bool empty;     //! Flag indicating if the object is in use.
//  bool changed;   //! Flag indicating wether the object's internal state was changed.

private:
  std::vector<double> bvalues;
  std::vector<FVector> gradients;
};

//===========================================================================

#endif // __FAMGradients_hh
