#include <vector>
#include <iostream>
#include <fstream>
#include <teem/nrrd.h>
#include <sstream>
#include "new_creaseline.hpp"

using namespace spurt;
using namespace gage_interface;

int main( int argc, char* argv[] )
{
    if ( argc != 4 )
    {
        std::cout << "USAGE: " << std::flush;
        std::cout << argv[0] << " <nrrd data> <input> <output>" 
            << std::endl;
        return -1;
    }

    Nrrd *nin = nrrdNew();
    if ( nrrdLoad( nin, argv[1], NULL ) )
    {
        std::cerr << "unable to open " << argv[1] << std::endl;
        return -1;
    }
    
    crease::MeasureWrapper *the_wrapper = new crease::MeasureWrapper( nin, 1 );
    
    std::fstream input( argv[2], std::ios::in );
    std::fstream output( argv[3], std::ios::out );

    vec3 p;
    while( !input.eof() )
    {
        double px, py, pz;
        char ch;

        input >> ch;

        if( input.eof() )
            break;

        switch( ch )
        {
            case 'p':
            input >> p[0] >> p[1] >> p[2];
            double val = the_wrapper->eigenvalue( p, 1 );
            output << "p " << p[0] << " " << p[1] << " " << p[2] << " " << -val << std::endl;
            break;
            case 'n':
            output << "n" << std::endl;
            break;
        }

        std::string tmp;
        std::getline( input, tmp );
    }
    
    input.close();
    output.close();

    return 0;
}
