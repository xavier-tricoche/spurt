//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTimeDependentTensorField.hh,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:13 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 

#ifndef __FTimeDependentTensorField_hh
#define __FTimeDependentTensorField_hh

//===========================================================================

#include "FObject.hh"
#include "FAnalysisModuleData.hh"
#include "FTensorField.hh"
#include "stdAliases.hh"

#include <string>
#include <vector>
#include <map>

#include "eassert.hh"

#include <boost/shared_ptr.hpp>

// needs new documentation badly :)

using namespace boost;


/*!  The FTimeDependentTensorField class represents a time dependent
  tensor field consisting of several FTensorFields where a time value
  for each tensor field is given. Using this time value one can
  interpolate in temporal direction using two tensor fields and fetch
  values inbetween.
  */
class FTimeDependentTensorField : public FObject
{
    typedef shared_ptr<FTensorField>    field_ptr;
    typedef std::map<double, field_ptr> time_map;

public:

    FTimeDependentTensorField( bool is_periodic = false ) : 
	name( "empty" ), periodic( is_periodic )
    {
    }

    ~FTimeDependentTensorField()
    {
    }

    // --- 

    const std::string& getName() const
    {
	return name;
    }

    void setName( const std::string& _name )
    {
        name = _name;
    }

    // ---

    const double& getMinT() const
    {
	return tmap.begin()->first;
    }

    const double& getMaxT() const
    {
	return tmap.rbegin()->first;
    }

    // --- access time steps by index

    positive getNbTimeSteps() const
    {
	return tmap.size();
    }

    field_ptr getTimeStep( const positive& id ) const;
    void      setTimeStep( const positive& id, field_ptr field );

    const double& getTimeValue( const positive &id ) const;
    void          setTimeValue( const positive &id, const double &t );

    // --- 

    void addTimeStep( field_ptr field, const double& t );
    void removeTimeStep( const positive& id );

    // --- interpolation / derivatives

    bool interpolate( FTensor& result, double t, const FArray& pos ) const;
    bool derivatives( FTensor& result, double t, const FArray& pos ) const;

    // --- misc.

    FAnalysisModuleData* getAnalysisModuleData() const
    {
	return &analysisModuleData;
    }

private:

    // --- access time steps by time

    bool getTimeSteps( double& t,
		       std::pair<time_map::const_iterator, 
		                 time_map::const_iterator>& steps ) const;

    time_map::const_iterator field_by_id( const positive& id ) const;
    time_map::iterator field_by_id( const positive& id );

    // strangely, map::lower_bound() is not const
    // therefore, need to declare tmap mutable
    time_map            tmap;
    std::string         name;
    bool                periodic;

    mutable FAnalysisModuleData analysisModuleData;
};

std::ostream& operator<<( std::ostream&, const FTimeDependentTensorField& );

#endif // __FTimeDependentTensorField_hh
