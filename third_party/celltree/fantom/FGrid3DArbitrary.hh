//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FGrid3DArbitrary.hh,v $
// Language:  C++
// Date:      $Date: 2004/07/14 06:37:00 $
// Author:    $Author: garth $
// Version:   $Revision: 1.2 $
//
//--------------------------------------------------------------------------- 

#ifndef __FGrid3DArbitrary_hh
#define __FGrid3DArbitrary_hh

#include "FGrid.hh"

class FGrid3DArbitrary : public FGrid
{
public:
    
    FGrid3DArbitrary( shared_ptr<FPositionSet> posSetPtr,
		      shared_ptr<FCellDefinitions> cellDef,
		      const std::string& newname );


    virtual const FString& getClassName() const;
    
    virtual ~FGrid3DArbitrary();
    
    void buildKdTree();
};

#endif
