//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FFeature.hh,v $
// Language:  C++
// Date:      $Date: 2003/11/19 09:21:01 $
// Author:    $Author: tricoche $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef __FFeature_hh
#define __FFeature_hh

#include <string>
#include <iostream>

struct FFeature
{
    FFeature() : name() {};

    virtual std::string getType() const 
	{
	    cout << "not defined!!!" << endl;
	    throw FException();
	};
    
    std::string name;
    std::string type;

    virtual ~FFeature() {};
};

#endif // __FFeature_hh
