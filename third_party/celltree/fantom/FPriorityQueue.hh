
//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FPriorityQueue.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:10 $
// Author:    $Author: garth $
// Version:   $Revision: 1.12 $
//
//--------------------------------------------------------------------------- 

// forward-declaration
template <class _Tp, class _Compare, class _PositionAccess>
class  FPriorityQueue;

#ifndef __FPriorityQueue_hh
#define __FPriorityQueue_hh

#include <vector>
#include <functional>
#include <iosfwd>

using namespace std;

/**
 * This code is partly "stolen" from the class priority_queue
 * as defined in stl_queue.h.
 * It extends the code so that every element knows its
 * position, i.e. the element has a member
 *     int FPriorityQueuePosition
 * which does what it looks like.
 * Otherwise it has to be specified as second template parameter.
 */

template<class T>
class FPriorityQueueDefaultPositionAccess {
public:
     long int & operator () ( T & x );
};

template<class T>
class FPriorityQueueDefaultPositionAccess< T * > {
public:
     long int & operator () ( T * & x );
};

// forward for the friends
template <class __Tp, class __Compare, class __PositionAccess>
ostream & operator<<(ostream & os, const FPriorityQueue<__Tp,__Compare,__PositionAccess> & pq);


/**
 * FPriorityQueue is an extension of the classical priority queue implemented
 * on a heap.
 * In addition to the push, pop, and top operations an update operation is 
 * introduced, that keeps the priority queue consistent if the priority of
 * one of its elements changes.
 *
 * The code is partly "stolen" from the class priority_queue
 * as defined in stl_queue.h.
 * Every element has to know its position, i.e. the element has a member
 *     int FPriorityQueuePosition
 * which does what it looks like.
 * Otherwise it has to be specified as second template parameter.
 *
 * If the elements are not compared by operator<, the comparison
 * operator has to be given as second template parameter.
 */
template <class _Tp, 
     class _Compare = less<typename vector<_Tp>::value_type>,
     class _PositionAccess = FPriorityQueueDefaultPositionAccess<_Tp> >
class  FPriorityQueue {
public:
	 /**
	  * \par Description:
	  * Constructor. Provides an empty priority queue.
	  * \pre
	  * none
	  * \post
	  * none
	  * \exception
	  * none
	  */
     FPriorityQueue();
     
	 /**
	  * \par Description:
	  * Constructor. Provides an empty priority queue.
	  * \pre
	  * none
	  * \post
	  * none
	  * \exception
	  * none
	  * \param
	  * __x: initializer for the _PositionAccess class
	  */
     explicit FPriorityQueue(const _PositionAccess& __x);
     
	 /**
	  * \par Description:
	  * Constructor. Provides an empty priority queue.
	  * \pre
	  * none
	  * \post
	  * none
	  * \exception
	  * none
	  * \param
	  * __x: initializer for the _Compare class
	  */
     explicit FPriorityQueue(const _Compare& __x);
     
	 /**
	  * \par Description:
	  * Constructor. Provides an empty priority queue.
	  * \pre
	  * none
	  * \post
	  * none
	  * \exception
	  * none
	  * \param
	  * __x: initializer for the _PositionAccess class
	  * __y: initializer for the _Compare class
	  */
     explicit FPriorityQueue(const _PositionAccess & __x, const _Compare& __y );

	 /**
	  * \par Description:
	  * Returns true iff the priority queue is empty
	  * \pre
	  * none
	  * \post
	  * none
	  * \exception
	  * none
	  * \return
	  * True iff the priority queue is empty
	  */
     bool empty() const;

	 /**
	  * \par Description:
	  * Returns the size, i.e. the number of entries of the priority queue
	  * \pre
	  * none
	  * \post
	  * none
	  * \exception
	  * none
	  * \return
	  * The number of entries in *this
	  */
     vector<_Tp>::size_type size() const;

	 /**
	  * \par Description:
	  * Returns the element with the highest priority
	  * \pre
	  * The priority queue is not empty
	  * \post
	  * none
	  * \exception
	  * none
	  * \param
	  * \return
	  * The element with the highest priority
	  */
     const _Tp & top() const;
	 
	 /**
	  * \par Description:
	  * Inserts the element __x into the priority queue
	  * \pre
	  * none
	  * \post
	  * none
	  * \exception
	  * none
	  * \param
	  * __x: the element that will be inserted
	  */
     void push(const _Tp& __x);

	 /**
	  * \par Description:
	  * Removes the element with the highest priority
	  * \pre
	  * The priority queue is not empty
	  * \post
	  * none
	  * \exception
	  * none
	  */
     void pop();

	 /**
	  * \par Description:
	  * Returns the element at position i
	  * \pre
	  * i < this->size()
	  * \post
	  * none
	  * \exception
	  * none
	  * \param
	  * i: the position which should be retrieved
	  * \return
	  * The element at position i
	  */
     const _Tp & operator[] (long int i) const;

	 /**
	  * \par Description:
	  * Erases the element at position i
	  * \pre
	  * i < this->size()
	  * \post
	  * none
	  * \exception
	  * none
	  * \param
	  * i: the position to erase
	  */
     void erase(long int i);

	 /**
	  * \par Description:
	  * Signals to the priority queue that the priority of the element at 
	  * position i has changed.
	  * The priority queue now updates itself
	  * \pre
	  * i < this->size()
	  * \post
	  * none
	  * \exception
	  * none
	  * \param
	  * i: the position where the priority has changed
	  */
	 void update(long int i);

protected:
	 /**
	  * \par Description:
	  * Method used internally by the priority queue that moves the element
	  * at position i upwards until the priority is consistent again
	  * \pre
	  * i < this->size()
	  * \post
	  * none
	  * \exception
	  * none
	  * \param
	  * i: the position that will has to be moved up
	  */
	 void move_up(long int i);
	 
     vector<_Tp> c;
     _PositionAccess posAcc;
     _Compare comp;

	 friend ostream & operator<< <_Tp,_Compare,_PositionAccess> (ostream & os, const FPriorityQueue<_Tp,_Compare,_PositionAccess> & pq);
};


#include "FPriorityQueue.icc"

#endif // __FPriorityQueue_hh
