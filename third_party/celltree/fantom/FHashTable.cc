//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FHashTable.cc,v $
// Language:  C++
// Date:      $Date: 2004/09/02 14:51:09 $
// Author:    $Author: hlawit $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#include "FCell.hh"

#include "FHashTable.hh"
#include "FString.hh"

#include "FTetrahedronCell.hh"

#include <iostream>
//---------------------------------------------------------------------------

FHashTable::FHashTable(unsigned int size)
  :sizeOfTable( size )
{
    counter = 0;

    init();
}

//---------------------------------------------------------------------------

FHashTable::~FHashTable()
{
    flood();
    delete [] table;
}

//---------------------------------------------------------------------------

void FHashTable::init () const
{
    table = new FHashElement[sizeOfTable];

    for (unsigned int i=0; i<sizeOfTable; i++)
    {
	table[i].object = shared_ptr<FCell>();
	table[i].id = 0;
    }
}

//---------------------------------------------------------------------------

shared_ptr<FCell> FHashTable::back () const
{
    return operator[](counter - 1);
}

//---------------------------------------------------------------------------

shared_ptr<FCell> FHashTable::operator[]( unsigned int c ) const
{
    unsigned int pos1 = c % sizeOfTable;

    if( table[pos1].id == c && table[pos1].object )
	return table[pos1].object;

    unsigned int pos2 =
	(pos1 + (c % (sizeOfTable - 2) + 1)) % sizeOfTable;

    if( table[pos2].id == c && table[pos2].object )
	return table[pos2].object;

    return shared_ptr<FCell>();
}

//---------------------------------------------------------------------------

void FHashTable::add( shared_ptr<FCell> element, unsigned int id )
{
    unsigned int pos1 = id % sizeOfTable;

    if( table[pos1].object )
    {
	unsigned int pos2 =
	    (pos1 + (id % (sizeOfTable - 2) + 1)) % sizeOfTable;

	table[pos2].object = element;
	table[pos2].id = id;
    }
    else
    {
	table[pos1].object = element;
	table[pos1].id = id;
    }
}

//---------------------------------------------------------------------------

void FHashTable::flood ()
{
    for( unsigned int i=0; i<sizeOfTable; i++ )
	table[i].object = shared_ptr<FCell>();
}

//---------------------------------------------------------------------------

positive FHashTable::memSize() const 
{
    //cout << "Size of FHashTable (" << sizeOfTable << ") is ";

    unsigned int tmp = sizeOfTable * sizeof( FHashElement );
    for(positive i=0;i<sizeOfTable;i++)
      if(table[i].id)
	tmp+=table[i].object->memSize();

    //cout << tmp << endl;
    return tmp;
}

//===========================================================================
