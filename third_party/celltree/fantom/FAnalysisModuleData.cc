//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FAnalysisModuleData.cc,v $
// Language:  C++
// Date:      $Date: 2003/06/25 07:22:43 $
// Author:    $Author: garth $
// Version:   $Revision: 1.9 $
//
//---------------------------------------------------------------------------

#include "FAnalysisModuleData.hh"
#include <iostream>

#include <cassert>

using namespace std;

  /**
   * \brief
   *generates an empty ModuleData
   * \param 
   * none
   * \return 
   * none
   */
FAnalysisModuleData::FAnalysisModuleData()
{
}   
   
  /**
   * \brief
   * destructor
   */
FAnalysisModuleData::~FAnalysisModuleData()
{
} 

  /**
   * \brief
   * Returns its class name
   * \return 
   * classname
   */
const FString& FAnalysisModuleData::getClassName() const
{
  static FString className("FAnalysisModuleData");

  return className;
}


  /**
   * \brief
   * Add an element derived from FAMElement.
   * \pre
   * 
   * \post
   * If another element of same class had been stored, it is now overwritten.
   * \exception
   * none
   * \param 
   * add_object Object to be added to the analysis module.
   * \return 
   * none
   */
void FAnalysisModuleData::addAMObject (boost::shared_ptr<FAMElement> add_object) 
{	
	deleteAMObject (add_object->getClassName ());
	
	// store new object
  	storedObjects.push_back (add_object);
}


  /**
   * \brief
   * Paste corresponding stored object into given object.
   * \pre
   * Object of needed class had been stored.
   * \post
   * If needed object had not been stored: get_object == NULL
   * \exception
   * none
   * \param 
   * get_object Object to paste the needed one into.
   * \return 
   * none
   */
void FAnalysisModuleData::getAMObject  (boost::shared_ptr<FAMElement> get_object) const
{
  assert( get_object ); // ensure the uses gives us a valid object (avoids nasty boost errors)

	for (list< boost::shared_ptr<FAMElement> >::const_iterator storedObjects_itr = storedObjects.begin();
	     storedObjects_itr != storedObjects.end ();
             storedObjects_itr++ )
	{
		if ( (*storedObjects_itr)->getClassName() == get_object->getClassName ()) {
			get_object = *storedObjects_itr;
			break;
		}
	}
}
 

  /**
   * \brief
   * Return object of given class.
   * \pre
   * none
   * \post
   * If needed object had not been stored: <returnedObject> == NULL
   * \exception
   * none
   * \param 
   * objName Class name of needed object.
   * \return 
   * Stored object or empty pointer.
   */
boost::shared_ptr<FAMElement> FAnalysisModuleData::getAMObject  (std::string objName) const
{

	for (list< boost::shared_ptr<FAMElement> >::const_iterator storedObjects_itr = storedObjects.begin();
	     storedObjects_itr != storedObjects.end ();
             storedObjects_itr++ )
	{
		if ( (*storedObjects_itr)->getClassName() == objName ) {
			return *storedObjects_itr;
		}
	}

	return boost::shared_ptr<FAMElement>();
}


  /**
   * \brief
   * Delete object with given class name from analysis module
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \param 
   * none
   * \return 
   * none
   */
void FAnalysisModuleData::deleteAMObject (std::string className)
{	
	for (list< boost::shared_ptr<FAMElement> >::iterator storedObjects_itr = storedObjects.begin();
	     storedObjects_itr != storedObjects.end ();
             storedObjects_itr++)
	{
		if ( (*storedObjects_itr)->getClassName() == className ) {
			storedObjects.erase (storedObjects_itr);
			return;
		}
	}
}


  /**
   * \brief
   * Clear list of stored objects.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \param 
   * none
   * \return 
   * none
   */
void FAnalysisModuleData::deleteAllAMObjects ()
{
	storedObjects.clear();	
}


  /**
   * \brief
   * Return list of all stored objects.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \param 
   * none
   * \return 
   * list of elements
   */
const std::list< boost::shared_ptr<FAMElement> >& FAnalysisModuleData::getAllAMObjects () const
{
	return storedObjects;
}

  /**
   * \brief
   * Return list of classnames of all stored objects.
   * \pre
   * none
   * \post
   * none
   * \exception
   * none
   * \param 
   * none
   * \return 
   * vector of classnames
   */
std::vector< std::string > FAnalysisModuleData::getAllClassNames () const
{
    std::vector< std::string > classnames;
    std::list< boost::shared_ptr<FAMElement> >::const_iterator it = storedObjects.begin();
    while(it != storedObjects.end())
    {
        classnames.push_back( (*it++)->getClassName() );
    }
    return classnames;
}
	
  /**
   * \brief
   * Method for debugging: Get number of stored objects.
   * \pre
   *none
   * \post
   * none
   * \exception
   * none
   * \param 
   * none
   * \return 
   * Number of stored objects
   */
positive FAnalysisModuleData::size () const
{
  	return storedObjects.size();
}


  /**
   * \brief
   * Method for debugging: Print list of all stored objects
   * \pre
   *none
   * \post
   * none
   * \exception
   * none
   * \param 
   * none
   * \return 
   * none
   */
void FAnalysisModuleData::print () const
{
	std::cout << "\nObjects stored in Analysis-Module:" << endl;
	std::cout << "\n\t" <<"Index\tClassName" << endl;
	std::cout <<   "\t" <<"------------------------------" << endl;

	int i=0;
	for (list< boost::shared_ptr<FAMElement> >::const_iterator storedObjects_itr = storedObjects.begin();
	     storedObjects_itr != storedObjects.end ();
         storedObjects_itr++, i++ )
	{
		std::cout << "\t" << i << "\t" << (*storedObjects_itr)->getClassName() << endl;
	}
	std::cout << endl;
}
