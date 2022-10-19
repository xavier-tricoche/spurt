#ifndef __FQUEUEELEMENT_HH
#define __FQUEUEELEMENT_HH
//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile
// Language:  C++
// Date:      $Date: 2001/01/17 10:08:31 $
// Author:    $Author: jfrey $
// Version:   $Revision: 1.4 $
//
//---------------------------------------------------------------------------

/**
 * Undocumented ...
 */
class FQueueElement
{
public:
  // declare destructor virtual
  virtual ~FQueueElement(){}
  /**
   * \brief
   * executes the queued command using its internal data.
   */
  virtual void execute () = 0;
};

#endif
