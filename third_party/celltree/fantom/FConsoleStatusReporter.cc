#include "FConsoleStatusReporter.hh"
#include <iostream>
#include <cmath>


FConsoleStatusReporter::FConsoleStatusReporter() : old_percent(101)
{

}

FConsoleStatusReporter::~FConsoleStatusReporter()
{

}

void FConsoleStatusReporter::say( const std::string& what ) const
{
    std::cout << what << std::endl;
}

void FConsoleStatusReporter::progressBegin( const std::string& what ) const
{
    say( what );
    old_percent = 101;
}

void FConsoleStatusReporter::progressEnd() const
{
    say( "\rdone" );
}

void FConsoleStatusReporter::progress( positive iter, positive max ) const
{
    // indicate some progress...
    positive percent = (positive)rint( (double)iter * 100.0 / (double)max );
    
    if( percent != old_percent ) 
    {
	std::cout << '\r' << percent << '%' << std::flush;
	old_percent = percent;
    }
}

void FConsoleStatusReporter::progress( double p ) const
{
    // indicate some progress...
    positive percent = (positive)rint( (double)p );
    
    if( percent != old_percent ) 
    {
	std::cout << '\r' << percent << '%' << std::flush;
	old_percent = percent;
    }
}
