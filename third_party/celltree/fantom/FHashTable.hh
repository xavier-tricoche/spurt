//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FHashTable.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:06 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FHashTable_hh
#define __FHashTable_hh

#include <boost/shared_ptr.hpp>
#include "FCell.hh"

using namespace boost;

//===========================================================================

/**
 * Undocumented.
 */
class FHashElement
{
public:
  /// The hashed object itself.
  shared_ptr<FCell> object;
  /// The object's id.
  unsigned int id;
};

/**
 * A hash table which implements a hash-rehash cache.
 */
class FHashTable
{
public:
  /**
   * default constructor
   */
  FHashTable(unsigned int size = 100);

  /**
   * destructor
   */
  ~FHashTable();

  /**
   * initialize hash table
   */
  void init() const;

  /**
   * return last element
   * \return requested element
   */
  shared_ptr<FCell> back() const;

  /**
   * return specified element
   * \param c position of element which will be returned
   * \return requested element
   */
  shared_ptr<FCell> operator[](unsigned int c) const;

  /**
   * add given element
   * \param element element to add
   * \param id of element
   */
  void add( shared_ptr<FCell> pos, unsigned int id);

  /**
   * clear the hash table
   */
  void flood ();

  positive memSize() const;

private:
  /** copy constructor
   */
  FHashTable(const FHashTable &object);

  const unsigned int sizeOfTable; /// size of hash table
  unsigned int counter; /// number of elements in hash table

  mutable FHashElement* table;
};

//===========================================================================

#endif  
