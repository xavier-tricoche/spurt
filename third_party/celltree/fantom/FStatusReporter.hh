//---------------------------------------------------------------------------
//
// Project:   FAnToM
// Module:    $RCSfile: FStatusReporter.hh,v $
// Language:  C++
// Date:      $Date: 2003/09/08 12:18:59 $
// Author:    $Author: garth $
// Version:   $Revision: 1.3 $
//
//--------------------------------------------------------------------------- 

#ifndef __FStatusReporter_hh
#define __FStatusReporter_hh

#include <string>
#include "stdAliases.hh"

/* an environment-independent status reporter base class
 * 
 * The Status Reporter can only be used in the algorithm
 * environment (thread), otherwise the main thread can
 * block esp. on ask and notify calls!
 */
class FStatusReporter
{
public:
    virtual ~FStatusReporter() {};

    /** set the status display text
     */ 
    virtual void say( const std::string& what ) const = 0;

    /** start a new progress indication
     */
    virtual void progressBegin( const std::string& what ) const = 0;

    /** end a progress indication
     */
    virtual void progressEnd() const = 0;

    /** indicate progress in percent from the fraction iter/max
     */
    virtual void progress( positive iter, positive max ) const = 0;

    /**
     * indicate progress in percent
     */
    virtual void progress( double percent ) const = 0;

	/** notify the user of an event and wait for acknowledgement
	 */
	virtual void notify( const std::string& what ) const = 0;

	/** ask the user a question
	 */
	virtual bool ask( const std::string& what, 
					  const std::string& positiveAnswer = "",
					  const std::string& negativeAnswer = "" ) const = 0;
};


/** the status reporter singleton, 
 *  (has to be) setup by the execution environment
 */
extern FStatusReporter *theStatusRep;

#endif
 
