#include <vector>
#include <iostream>
#include <fstream>
#include <teem/nrrd.h>
#include <sstream>
#include "creaseline3d.hpp"
#include <string>
#include <string.h>
#include <stdio.h>

using namespace xavier;
using namespace gage_interface;

std::string in_name, out_name;
unsigned int up;
unsigned int aniso;
bool do_ridge;

void parse( int argc, char* argv[] )
{
    if ( argc < 5 )
    {
        std::cerr << "USAGE: " << argv[0] << " -i <input> -o <output> [-u <upsampling>] [-a <aniso>] [-ridge / -valley]"
            << std::endl;
        exit( -1 );
    }

    up = 1;
    aniso = 0;
    do_ridge = true;

    bool set_in, set_out;
    set_in = set_out = false;

    bool failed = false;

    for ( unsigned int i=1 ; i<argc && !failed ; i++ )
    {
        if ( argv[i][0] != '-' || strlen( argv[i] )<2 ) 
        {
            failed = true;
            break;
        }
        switch( argv[i][1] )
        {
            case 'i':
            {
                if ( i==argc-1 )
                {
                    failed = true;
                }
                else
                {
                    set_in = true;
                    in_name = argv[i+1];
                    ++i;
                }
                break;
            }
            case 'o':
            {
                if ( i==argc-1 )
                {
                    failed = true;
                }
                else
                {
                    set_out = true;
                    out_name = argv[i+1];
                    ++i;
                }
                break;
            }
            case 'u':
            {
                if ( i==argc-1 )
                {
                    failed = true;
                }
                else
                {
                    up = atoi( argv[i+1] );
                    ++i;
                }
                break;
            }
            case 'a':
            {
                if ( i==argc-1 )
                {
                    failed = true;
                }
                else
                {
                    aniso = atoi( argv[i+1] );
                    ++i;
                }
                break;
            }
            case 'r':
            {
                if ( strcmp( argv[i], "-ridge" ) )
                {
                    failed = true;
                }
                else
                {
                    do_ridge = true;
                }
                break;
            }
            case 'v':
            {
                if ( strcmp( argv[i], "-valley" ) )
                {
                    failed = true;
                }
                else
                {
                    do_ridge = false;
                }
                break;
            }
            default:
            {
                std::cout << "default: " << argv[i] << std::endl;
                failed = true;
                break;
            }
        }
    }

    if ( failed )
    {
        std::cerr << "USAGE: " << argv[0] << " -i <input> -o <output> [-u <upsampling>] [-a <aniso>] [-ridge / -valley]"
            << std::endl;
        exit( -1 );
    }
}

int main( int argc, char* argv[] )
{
    parse( argc, argv );

    Nrrd *nin = nrrdNew();
    if ( nrrdLoad( nin, in_name.c_str(), NULL ) )
    {
        std::cerr << "unable to open " << argv[1] << std::endl;
        return -1;
    }

    bool is_tensor = ( nin->dim == 4 && nin->axis[0].size == 7 );
    if ( is_tensor )
    {
        if ( aniso == 0 ) aniso == 1;
    }

    gageShape *shape = gageShapeNew();
    if ( gageShapeSet( shape, nin, nin->dim-3 ) )
    {  
        std::cerr << "ERROR in FDTICoherence: " 
            << biffGetDone( GAGE ) << std::endl;
    }

    nvis::vec3 p;
    unsigned int id;

    xavier::crease::threshold = 0.80;
    xavier::crease::subdiv = 1;
    xavier::crease::eps = 0.05;
    xavier::crease::upsample = up;
    xavier::crease::extract_lines( nin, do_ridge, aniso );

    std::ostringstream os;
    os << out_name << "-all.vcl";

    std::fstream output( os.str().c_str(), std::ios::out );
    for ( unsigned int i=0 ; i<xavier::crease::all_edges.size() ; i++ )
    {
        id = xavier::crease::all_edges[i].first;
        p = xavier::crease::all_face_points[id];
        double v = xavier::crease::crease_strength[id];
        output << "p " << p[0] << " " << p[1] << " " << p[2] << " " << v << std::endl;

        id = xavier::crease::all_edges[i].second;
        p = xavier::crease::all_face_points[id];
        v = xavier::crease::crease_strength[id];
        output << "p " << p[0] << " " << p[1] << " " << p[2] << " " << v << std::endl;
        output << "n" << std::endl;
    }
    output.close();

    os.clear();
    os.str( "" );
    os << out_name << "-wc-connected.vcl";
    output.open( os.str().c_str(), std::ios::out );
    for ( unsigned int i=0 ; i<xavier::crease::components.size() ; i++ )
    {
        double wc[3];
        unsigned int n=0;
        for ( unsigned int j=0 ; j<xavier::crease::components[i].size() ; j++ ) 
        {
            std::cout << n++ << "-" << std::flush;
            p = xavier::crease::all_face_points[xavier::crease::components[i][j]];
            double ic[3] = { p[0], p[1], p[2] };
            gageShapeItoW( shape, wc, ic );
            double v = xavier::crease::crease_strength[xavier::crease::components[i][j]];
            output << "p " << wc[0] << " " << wc[1] << " " << wc[2] << " " << v << std::endl;
        }
        output << "n" << std::endl;
        std::cout << "done" << std::endl;
    }
    output.close();

    os.clear();
    os.str( "" );
    os << out_name << "-connected.vcl";
    output.open( os.str().c_str(), std::ios::out );
    for ( unsigned int i=0 ; i<xavier::crease::components.size() ; i++ )
    {
        std::cout << "component #" << i << ", length = " 
            << xavier::crease::components[i].size() << std::endl;
        unsigned int n=0;
        for ( unsigned int j=0 ; j<xavier::crease::components[i].size() ; j++ ) 
        {
            std::cout << n++ << "-" << std::flush;
            p = xavier::crease::all_face_points[xavier::crease::components[i][j]];
            double v = xavier::crease::crease_strength[xavier::crease::components[i][j]];
            output << "p " << p[0] << " " << p[1] << " " << p[2] << " " << v << std::endl;
        }
        output << "n" << std::endl;
        std::cout << "done" << std::endl;
    }
    output.close();

    os.clear();
    os.str( "" );
    os << out_name << "-connected_dot_vals.vcl";
    output.open( os.str().c_str(), std::ios::out );
    for ( unsigned int i=0 ; i<xavier::crease::components.size() ; i++ )
    {
        for ( unsigned int j=0 ; j<xavier::crease::components[i].size() ; j++ ) 
        {
            p = xavier::crease::all_face_points[xavier::crease::components[i][j]];
            double u = xavier::crease::crease_strength[xavier::crease::components[i][j]];
            double v = xavier::crease::grad_dot_evec[xavier::crease::components[i][j]];
            double w = xavier::crease::measure_value[xavier::crease::components[i][j]];
            output << "p " << p[0] << " \t" << p[1] << " \t" << p[2] << " \t dot: " 
                << v << " \t strength: " << u << " \t value: " << w
                << std::endl;
        }
    }
    output.close();
    std::cout << std::endl;

    return 0;
}
