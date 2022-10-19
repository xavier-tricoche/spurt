#ifndef __FStatusReporterAdapter_hh
#define __FStatusReporterAdapter_hh

#include "FStatusReporter.hh"
#include <sstream>


class FStatusReporterAdapter
{
public:
  FStatusReporterAdapter(FStatusReporter* statusRep )
    :myStatusRep(statusRep)
    {
    }
  
  /*
   * This function extends the functionality of FStatusReporter's progressBegin(...).
   * It allows to pass two arguments that indicate an iteration of iterations.
   * The metaIteration-th iteration of iterations is diplayed as "metaIteration / nbMetaIterations " 
   * at the end of the message passed by the argument "what". 
   * nbMetaIterations here represents the total number of iterations of iterations. 
   * If only the metaIteration is passed as argument, it is just appended the message passed by "what".
   * If neither metaIteration nor nbMetaIterations is passed this function behaves like 
   * FStatusReporter::progressBegin
   */
  void progressBegin( const std::string& what, const positive& metaIteration=0, const positive& nbMetaIterations=0 ) const
    {
      
      std::ostringstream out;
      out<<what;
      if(metaIteration==0&&nbMetaIterations==0)
      {
      }
      else if(nbMetaIterations==0)
      {
	out<<' '<<metaIteration;
      }
      else
      {
	out<<' '<<metaIteration<<'/'<<nbMetaIterations;
      }
      myStatusRep->progressBegin(out.str());
      
      
    }
  
  /*
   * The function just calls the function of the same name of the FStatusReporter passed to the constructor
   */
  void say( const std::string& what ) const 
    {
      myStatusRep->say(what);
    }
  
  /*
   * The function just calls the function of the same name of the FStatusReporter passed to the constructor
   */
  void progressEnd() const 
    {
      myStatusRep->progressEnd();
    }

  /*
   * The function just calls the function of the same name of the FStatusReporter passed to the constructor
   */
  void progress( positive iter, positive max ) const
    {
      myStatusRep->progress( iter, max );
    }

  /*
   * The function just calls the function of the same name of the FStatusReporter passed to the constructor
   */
  void progress( double percent ) const
    {
      myStatusRep->progress( percent );
    }

  /*
   * The function just calls the function of the same name of the FStatusReporter passed to the constructor
   */
  void notify( const std::string& what ) const 
    {
      myStatusRep->notify( what );
    }

  /*
   * The function just calls the function of the same name of the FStatusReporter passed to the constructor
   */
  bool ask( const std::string& what, 
	    const std::string& positiveAnswer = "",
	    const std::string& negativeAnswer = "" ) const
    {
      return myStatusRep->ask( what, positiveAnswer, negativeAnswer );
    }

  FStatusReporter* myStatusRep;
};

#endif
