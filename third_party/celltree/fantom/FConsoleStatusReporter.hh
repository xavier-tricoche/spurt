#ifndef __FConsoleStatusReporter_hh
#define __FConsoleStatusReporter_hh

#include <string>
#include "stdAliases.hh"
#include "FStatusReporter.hh"

/**
 * A status reporter dumping its output to console
 */

class FConsoleStatusReporter : public FStatusReporter
{
  mutable positive old_percent;

public:
  FConsoleStatusReporter();
  virtual ~FConsoleStatusReporter();

  /** set the status display text
   */ 
  void say( const std::string& what ) const;

  /** start a new progress indication
   */
  void progressBegin( const std::string& what ) const;

  /** end a progress indication
   */
  void progressEnd() const;

  /** indicate progress in percent from the fraction iter/max
   */
  void progress( positive iter, positive max ) const;

  /** indicate progress in percent
   */
  void progress( double percent ) const;

  /** notify the user of an event and wait for acknowledgement
   */
  void notify( const std::string& /*what*/ ) const {}

  /** ask the user a question
   */
  bool ask( const std::string& /*what*/, 
	        const std::string& /*positiveAnswer*/ = "",
			const std::string& /*negativeAnswer*/ = "" ) const 
  {
	return false;
  }
};

#endif
 
