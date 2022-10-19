//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTensorFieldReader.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/16 06:34:12 $
// Author:    $Author: wiebel $
// Version:   $Revision: 1.5 $
//
//---------------------------------------------------------------------------

#ifndef __FDataSetReader_hh
#define __FDataSetReader_hh

#include "FObject.hh"

#include "stdAliases.hh"

#include <cstdio>
#include <vector>
#include <ctime>
#include <string>
#include <list>
#include "FCell.hh"

#include <boost/shared_ptr.hpp>

using namespace boost;

class FTensorField;
class FPositionSet;
class FTensorSet;
class FCellDefinitions;
class FCell;

namespace FTensorFieldReader
{
    shared_ptr<FPositionSet>      loadPositionSetDLR( const std::string& file );
    shared_ptr<FCellDefinitions>  loadCellDefsDLR( const std::string& file );
    shared_ptr<FTensorSet>        loadTensorSetDLR( const std::string& file );

    shared_ptr<FTensorField>      loadDLR( const std::string& file );
}

#endif // __FDataSetReader_hh
