//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FNewHandler.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/13 18:35:48 $
// Author:    $Author: wiebel $
// Version:   $Revision: 1.5 $
//
//--------------------------------------------------------------------------- 

#ifndef __FNewHandler_hh
#define __FNewHandler_hh

//===========================================================================

/** 
 * The FNewHandler class provides a way to throw an exceptions when an 
 * operator new() fails. Therefore the method
 * FNewHander::new_handler must be installed in main:
 * void (*old_handler() =
 * set_new_handler(&FNewHandler::new_handler);} 
*/
class FNewHandler
{
public:

  /** This method simply throws a FOutOfMemoryException when called.
   * \exception
   * FOutOfMemoryException
   */
  static void new_handler();

};

//===========================================================================
#endif // __FNewHandler_hh
 
