#ifndef __FCOMMANDQUEUE_HH
#define __FCOMMANDQUEUE_HH
//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile
// Language:  C++
// Date:      $Date: 2004/07/14 06:35:20 $
// Author:    $Author: garth $
// Version:   $Revision: 1.10 $
//
//---------------------------------------------------------------------------


class FQueueElement;

/**
 * A queue for commands. Each class containing a command must be
 * derived from FQueueElement and allocated by the user. After
 * adding to the queue the element is owned by the queue and must not
 * be deleted by the user (contrasting the programming rules)!
 */
class FCommandQueue
{
public:
  /**
   * \brief
   * Constructor.
   */
  FCommandQueue ();

  /**
   * \brief
   * Destructor.
   */
  ~FCommandQueue ();

  /**
   * \brief
   * stop the queue.
   */
  void stop ();

  /**
   * \brief
   * append a command to the queue.
   */
  void append (FQueueElement *element);

  /**
   * \brief
   * executes the next queued command using its internal data.
   */
  bool executeNext ();

  /**
   * \brief
   * Is this queue empty?
   */
  bool isEmpty();


private:
  struct FCommandQueuePrivateData;
  FCommandQueuePrivateData *pd;
};

#endif
