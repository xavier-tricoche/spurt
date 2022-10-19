#include <fstream>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <map>
#include <math/fixed_vector.hpp>

std::vector< std::vector< unsigned int > > ueber_components;
std::vector< std::vector< unsigned int > > components;
std::vector< vec3 > points;

double eps;

struct EpsilonSpatialOrdering
{
    bool operator()( const vec3& p0, const vec3& p1 )
    {
        return 
            ( p0[0]<p1[0]-eps || 
              ( fabs( p1[0]-p0[0] )<=eps && 
                ( p0[1]<p1[1]-eps || 
                  ( fabs( p1[1]-p0[1] )<=eps && 
                    ( p0[2]<p1[2]-eps ) ) ) ) );
    }
};

int main( int argc, char* argv[] )
{
    if ( argc != 4 )
    {
        std::cout << "USAGE: " << argv[0] 
            << " <components> <ueber_components> <eps>"
            << std::endl;

        return -1;
    }

    std::fstream input( argv[1], std::ios::in );
    components.clear();
    unsigned int nb_pos = 0;

    components.push_back( std::vector< unsigned int >() );
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
            input >> px >> py >> pz;
            points.push_back( vec3( px, py, pz ) );
            components.back().push_back( nb_pos++ );
            break;
            case 'n':
            components.push_back( std::vector< unsigned int >() );
            break;
        }

        std::string tmp;
        std::getline( input, tmp );
    }
    components.pop_back();
    
    unsigned int length;
    unsigned int max = 0;
    std::cout << "INPUT: " << std::endl
        << "\t " << components.size() << " components" << std::endl;
    for ( unsigned int i=0 ; i<components.size() ; i++ )
    {
        unsigned int l = components[i].size();
        if ( l>max ) max = l;
        length += components[i].size();
    }
    std::cout << "\t average length = " << ( double )length/( double )components.size() << std::endl;
    std::cout << "\t max length = " << max << std::endl;

    ueber_components.clear();    
    
    eps = atof( argv[3] );

    // uniquify found locations
    std::map< vec3, unsigned int, EpsilonSpatialOrdering > unique_points;
    std::map< vec3, unsigned int, EpsilonSpatialOrdering >::iterator it;
    
    // test ordering
    if ( false )
    {
        std::cout << "epsilon = " << eps << std::endl;
        
        unique_points.clear();
        vec3 q0( 0, 0, 0 ), q;
        unique_points[q0] = 10;
        
        unsigned int wrong = 0;
        for ( unsigned int i=0 ; i<100000 ; i++ )
        {
            q[0] = eps*( -1.5 + 3.0*drand48() );
            q[1] = eps*( -1.5 + 3.0*drand48() );
            q[2] = eps*( -1.5 + 3.0*drand48() );
            
            bool inside = 
                ( q[0] >= -eps ) && ( q[0] <= eps ) && 
                ( q[1] >= -eps ) && ( q[1] <= eps ) && 
                ( q[2] >= -eps ) && ( q[2] <= eps );
                
            it = unique_points.find( q );
            if ( inside && it == unique_points.end() )
            {
                ++wrong;
            }
            else if ( !inside && it != unique_points.end() )
            {
                ++wrong;
                std::cout << "things are seriously fucked up: " << it->second << std::endl;
                unique_points.clear();
                unique_points[q0] = 10;
            }
        }
        
        std::cout << wrong << " wrong answers" << std::endl;
    }    
    
    unique_points.clear();
    
    std::vector< unsigned int > old2new( points.size() );
    unsigned int sz = 0;
    for ( unsigned int i=0 ; i<points.size() ; i++ )
    {
        std::cout << "point #" << i << ", nb inserted: " << unique_points.size() 
            << ", nb unique: " << sz << std::endl;
        it = unique_points.find( points[i] );
        if ( it != unique_points.end() )
        {
            std::cout << " * existing * " << std::endl;
            old2new[i] = it->second;
        }
        else
        {
            //unique_points.insert( std::pair< vec3, unsigned int >( points[i], sz ) );
            unique_points[points[i]] = sz;
            it = unique_points.find( points[i] );
            std::cout << "after insertion: element in map = "
                << ( ( it != unique_points.end() ) ? "TRUE" : "FALSE" )
                << std::endl;
            if ( it != unique_points.end() )
            {
                std::cout << "corresponding position: " << it->first 
                    << ", corresponding index: " << it->second << std::endl;
            }
            else
            {
                std::cout << "unique_points[points[i]] = " << unique_points[points[i]] << std::endl;
                std::cout << ( unique_points.find( points[i] ) == unique_points.end()  ?
                    "still not inside though" : "finally added to map" ) << std::endl;
            } 
        
            old2new[i] = sz;
            ++sz;
        }
    }
    std::cout << "last entered index is " << sz-1 << std::endl;
    sz = unique_points.size();
    std::cout << std::endl
        << "uniquification completed: " << sz << " unique points (" 
        << points.size() << ")" << std::endl;

    // connect loose ends of connected components
    std::vector< std::pair< unsigned int, unsigned int > > ends;
    for ( unsigned int i=0 ; i<components.size() ; i++ )
    {
        ends.push_back( std::pair< unsigned int, unsigned int >( old2new[components[i][0]], i ) );
        ends.push_back( std::pair< unsigned int, unsigned int >( old2new[components[i].back()], i ) );
    }

    std::vector< std::list< unsigned int > > p2c( sz );
    for ( unsigned int i=0 ; i<ends.size() ; i++ )
        p2c[ends[i].first].push_back( ends[i].second );
        
    std::cout << "point to curves connections:" << std::endl;
    unsigned int mc = 0;
    for ( unsigned int i=0 ; i<p2c.size() ; i++ )
    {
        unsigned int s = p2c[i].size();
        if ( s > mc ) mc = s;
        std::list< unsigned int >& links = p2c[i];
        std::list< unsigned int >::iterator it;
        std::cout << i << ": " << std::flush;
        for ( it=links.begin() ; it!=links.end() ; it++ )
            std::cout << *it << " " << std::flush;
        std::cout << std::endl; 
    }
    std::cout << "max number of connections at a point: " << mc << std::endl;

    std::vector< std::list< unsigned int > > chain_of_curves;

    // ***
    std::vector< bool > included( components.size(), false );
    for ( unsigned int i=0 ; i<p2c.size() ; i++ )
    {
        std::list< unsigned int >& links = p2c[i];
        if ( links.size() < 2 ) continue;

        if ( links.size() > 2  )
        {
            std::cout << links.size() << "-way branch detected" << std::endl;
            continue; 
        }
        unsigned int curve1 = links.front();
        unsigned int curve2 = links.back();

        // follow both links
        chain_of_curves.push_back( std::list< unsigned int >() );
        std::list< unsigned int >& mylist = chain_of_curves.back();

        // 1st component
        unsigned int cur_p = i;
        unsigned int cur_c = curve1;
        while ( true )
        {
            if ( included[cur_c] ) break;
            mylist.push_back( cur_c );

            unsigned int end1 = components[cur_c].front();
            unsigned int end2 = components[cur_c].back();

            // identify next position
            if ( old2new[end1] == cur_p )
                cur_p = old2new[end2];
            else
                cur_p = old2new[end1];

            // identify next curve
            if ( p2c[cur_p].size()==2 )
            {
                if ( p2c[cur_p].front() == cur_c )
                    cur_c = p2c[cur_p].back();
                else
                    cur_c = p2c[cur_p].front();
            }
            else
                break;
        }

        // 2nd component
        cur_p = i;
        cur_c = curve2;
        while ( true )
        {
            if ( included[cur_c] ) break;
            mylist.push_front( cur_c );

            included[cur_c] = true;

            unsigned int end1 = components[cur_c].front();
            unsigned int end2 = components[cur_c].back();

            if ( old2new[end1] == cur_p )
                cur_p = old2new[end2];
            else
                cur_p = old2new[end1];

            if ( p2c[cur_p].size()==2 )
            {
                if ( p2c[cur_p].front() == cur_c )
                    cur_c = p2c[cur_p].back();
                else
                    cur_c = p2c[cur_p].front();
            }
            else
                break;
        }

        if ( chain_of_curves.back().size() )
        {
            std::vector< unsigned int > curves( links.begin(), links.end() );
            unsigned int curve1, curve2;
            unsigned int end1, end2, n_end1, n_end2;

            ueber_components.push_back( std::vector< unsigned int >() );

            curve1 = curves[0];
            for ( unsigned int i=1 ; i<curves.size() ; i++ )
            {
                curve2 = curves[i];

                // both ends of current curve
                end1 = old2new[components[curve1].front()];
                end2 = old2new[components[curve1].back()];

                // both ends of next curve
                n_end1 = old2new[components[curve2].front()];
                n_end2 = old2new[components[curve2].back()];

                if ( end1==n_end1 || end1==n_end2 )
                {
                    std::vector< unsigned int >::reverse_iterator it;
                    for ( it=components[curve1].rbegin() ; it!=components[curve1].rend() ; it++ )
                        ueber_components.back().push_back( *it );
                }
                else
                {
                    std::vector< unsigned int >::iterator it;
                    for ( it=components[curve1].begin() ; it!=components[curve1].end() ; it++ )
                        ueber_components.back().push_back( *it );
                }

                curve1 = curve2;
            }
            
            // append last curve
            {    
                // both ends of current curve
                end1 = old2new[components[curve1].front()];
                end2 = old2new[components[curve1].back()];

                // end of last curve
                unsigned int n_end = ueber_components.back().back();

                if ( end1==n_end )
                {
                    std::vector< unsigned int >::iterator it;
                    for ( it=components[curve1].begin() ; it!=components[curve1].end() ; it++ )
                        ueber_components.back().push_back( *it );
                }
                else
                {
                    std::vector< unsigned int >::reverse_iterator it;
                    for ( it=components[curve1].rbegin() ; it!=components[curve1].rend() ; it++ )
                        ueber_components.back().push_back( *it );
                }
            }
            
        }
        else
            chain_of_curves.pop_back();
    }

    for ( unsigned int i=0 ; i<components.size() ; i++ )
    {
        if ( included[i] ) continue;

        ueber_components.push_back( std::vector< unsigned int >() );
        for ( unsigned int j=0 ; j<components[i].size() ; j++ )
            ueber_components.back().push_back( components[i][j] );
    }    
    
    max = 0;
    length = 0;
    std::cout << "OUTPUT: " << std::endl
        << "\t " << ueber_components.size() << " components" << std::endl;
    for ( unsigned int i=0 ; i<ueber_components.size() ; i++ )
    {
        unsigned int l = ueber_components[i].size();
        if ( l>max ) max = l;
        length += components[i].size();
    }
    std::cout << "\t average length = " << ( double )length/( double )ueber_components.size() << std::endl;
    std::cout << "\t max length = " << max << std::endl;
    
    // export
    std::fstream output( argv[2], std::ios::out );
    for ( unsigned int i=0 ; i<ueber_components.size() ; i++ )
    {
        for ( unsigned int j=0 ; j<ueber_components[i].size() ; j++ )
        {
            unsigned int id = ueber_components[i][j];
            output << "p " << points[id][0] << " " << points[id][1] << " " << points[id][2]
                << std::endl;
        }
        output << "n" << std::endl;
    }

    output.close();
    
    return 0;
}
