//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FTimeDependentTensorField.cc,v $
// Language:  C++
// Date:      $Date: 2003/08/18 11:58:13 $
// Author:    $Author: garth $
// Version:   $Revision: 1.1 $
//
//--------------------------------------------------------------------------- 


#include "FTimeDependentTensorField.hh"
#include "FException.hh"
#include "FTensor.hh"

#include "eassert.hh"

FTimeDependentTensorField::time_map::const_iterator FTimeDependentTensorField::field_by_id( const positive& id ) const
{
    eassert( id < tmap.size() );
    
    time_map::const_iterator ti = tmap.begin();
    std::advance( ti, id );
    
    return ti;
}

FTimeDependentTensorField::time_map::iterator FTimeDependentTensorField::field_by_id( const positive& id )
{
    eassert( id < tmap.size() );
    
    time_map::iterator ti = tmap.begin();
    std::advance( ti, id );
    
    return ti;
}

// --------------------------------------------------------------------------

FTimeDependentTensorField::field_ptr 
FTimeDependentTensorField::getTimeStep( const positive& id ) const
{
    return field_by_id( id )->second;
}

// --------------------------------------------------------------------------

void FTimeDependentTensorField::setTimeStep( const positive& id, field_ptr field )
{
    field_by_id( id )->second = field;
}

// --------------------------------------------------------------------------

const double& FTimeDependentTensorField::getTimeValue( const positive &id ) const
{
    return field_by_id( id )->first;
}

// --------------------------------------------------------------------------

void FTimeDependentTensorField::setTimeValue( const positive &id, const double &t )
{
    time_map::iterator ti = field_by_id( id );
    
    field_ptr field = ti->second;
    tmap.erase( ti );
    
    addTimeStep( field, t );
}

// --------------------------------------------------------------------------

void FTimeDependentTensorField::addTimeStep( field_ptr field, const double& t )
{
    // can't insert existing time, check
    time_map::iterator ti = tmap.find( t );
    
    if( ti != tmap.end() )
    {
	std::ostringstream out;
	out << "timestep t = " << t << " already in FTimeDependentTensorField";
	
	throw std::runtime_error( out.str() );
    }
    
    tmap.insert( std::make_pair( t, field ) );
}


// --------------------------------------------------------------------------

void FTimeDependentTensorField::removeTimeStep( const positive& id )
{
    eassert( id < tmap.size() );
    
    time_map::iterator ti = tmap.begin();
    std::advance( ti, id );
    
    tmap.erase( ti );
}

// --------------------------------------------------------------------------

bool FTimeDependentTensorField::getTimeSteps( double& t,
					      std::pair<time_map::const_iterator, 
					                time_map::const_iterator>& steps ) const
{
    if( tmap.empty() )
	return false;
    
    if( periodic )
    {
	double period = tmap.rbegin()->first - tmap.begin()->first;
	t -= period * floor( (t - tmap.begin()->first)/period );
    }
    else if( t > tmap.rbegin()->first || t < tmap.begin()->first )
	return false;

    time_map::const_iterator ti = tmap.upper_bound( t );
    
    if( ti == tmap.begin() )
	return false;

    time_map::const_iterator tl = ti, tu = ti;
    --tl;
    
    if( tu == tmap.end() )
    {
	--tu;
	--tl;
    }
    
    steps = std::make_pair( tl, tu );
    return true;
}

// --------------------------------------------------------------------------

bool FTimeDependentTensorField::interpolate( FTensor& result, double t, 
					     const FArray& pos ) const
{
    std::pair<time_map::const_iterator, time_map::const_iterator> ti;
    
    if( !getTimeSteps( t, ti ) )
	return false;
    
    const double& t0 = ti.first->first;
    const double& t1 = ti.second->first;
    
    FTensor r0, r1;
    
    bool okay = 
	ti.first->second->interpolate( r0, pos ) &&
	ti.second->second->interpolate( r1, pos );
    
    if( !okay )
	return false;
    
    result = (t-t0)/(t1-t0)*r0 + (t1-t)/(t1-t0)*r1;
    return true;
}

// --------------------------------------------------------------------------

bool FTimeDependentTensorField::derivatives( FTensor& result, 
					     double t, 
					     const FArray& pos ) const
{
    std::pair<time_map::const_iterator, time_map::const_iterator> ti;
    
    if( !getTimeSteps( t, ti ) )
	return false;
    
    const double& t0 = ti.first->first;
    const double& t1 = ti.second->first;
    
    FTensor r0, r1;
    
    bool okay = 
	ti.first->second ->derivatives( r0, pos ) &&
	ti.second->second->derivatives( r1, pos );
    
    if( !okay )
	return false;
    
    result = (t-t0)/(t1-t0)*r0 + (t1-t)/(t1-t0)*r1;
    return true;
}

// --------------------------------------------------------------------------

std::ostream& operator<<( ostream& os, const FTimeDependentTensorField& tdtf ) 
{
    os << "=====================================\n"
       << "= FTimeDependentTensorField Report ==\n"
       << "=====================================\n"
       << "Name: " << tdtf.getName() << "\n"
       << tdtf.getNbTimeSteps() << " timesteps stored.\n";
	
    for( unsigned int id=0; id<tdtf.getNbTimeSteps(); ++id )
	os << tdtf.getTimeValue(id) << "\t" << tdtf.getTimeStep(id)->getName() << "\n";
  
    return os;
}
