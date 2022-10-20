#include "new_creaseline.hpp"
#include <algorithm>
#include <thread>

using namespace std;
using namespace nvis;
using namespace spurt;
using namespace gage_interface;


double spurt::crease::threshold;
double spurt::crease::eps;
unsigned int spurt::crease::subdiv;

vector< crease::Edge > crease::problematic_edges;
vector< crease::Point > crease::all_face_points;
vector< pair< unsigned int, unsigned int > > crease::all_edges;
map< crease::FaceId, unsigned int > face_found;
vector< list< unsigned int > > crease::components;
crease::MeasureWrapper* the_wrapper;


MeasureWrapper::MeasureWrapper( const scalar_wrapper& scal )
    : _scal_wrap( scal ), _tens_wrap(0), _evecs(3), _evals(3), _measure(0)
    {}

MeasureWrapper::MeasureWrapper( const tensor_wrapper& tens, int aniso )
    : _scal_wrap(0), _tens_wrap( tens ), _evecs(3), _evals(3), _measure(aniso+1)
    {}

Vector MeasureWrapper::eigenvector( const Point& p, unsigned int idx ) const
{
    switch ( _measure )
    {
        case 0:
        {
            _scal_wrap.hess_evecs( p, _evecs );
            return _evecs[idx];
        }
        case 1:
        {
            double u[] = { p[0], p[1], p[2] };
            _tens_wrap.fa_hess_evec( u, _evecs_array );
            return Vector( _evecs_array[idx], _evecs_array[idx+1], _evecs_array[idx+2] );
        }
        case 2:
        {
            double u[] = { p[0], p[1], p[2] };
            _tens_wrap.mode_hess_evec( u, _evecs_array );
            return Vector( _evecs_array[idx], _evecs_array[idx+1], _evecs_array[idx+2] );
        }
        default:
        {
            assert( false );
        }
    }
}

double MeasureWrapper::eigenvalue( const Point& p, unsigned int idx )
{
    switch ( _measure )
    {
        case 0:
        {
            _scal_wrap.hess_evals( p, _evals );
            return _evals[idx];
        }
        case 1:
        {
            double u[] = { p[0], p[1], p[2] };
            _tens_wrap.fa_hess_evals( u, _evals_array );
            return _evals_array[idx];
        }
        case 2:
        {
            double u[] = { p[0], p[1], p[2] };
            _tens_wrap.mode_hess_evals( u, _evals_array );
            return _evals_array[idx];
        }
        default:
        {
            assert( false );
        }
    }
}

Vector MeasureWrapper::gradient( const Point& p )
{
    switch ( _measure )
    {
        case 0:
        {
            _scal_wrap.gradient( p, _grad );
            return _grad;
        }
        case 1:
        {
            double u[] = { p[0], p[1], p[2] };
            _tens_wrap.fa_grad( u, _grad_array );
            return Vector( _grad_array[0], _grad_array[1], _grad_array[2] );
        }
        case 2:
        {
            double u[] = { p[0], p[1], p[2] };
            _tens_wrap.mode_grad( u, _grad_array );
            return Vector( _grad_array[0], _grad_array[1], _grad_array[2] );
        }
        default:
        {
            assert( false );
        }
    }
}

Vector gradient( const Point& p ) const;
double eigenvalue( const Point& p, unsigned int idx ) const;

private:
    const scalar_wrapper& _scal_wrap;
    const tensor_wrapper& _tens_wrap;
    int _measure;
};

crease::Vector crease::eigenvector( const Point& p, const scalar_wrapper& gH,
    unsigned int idx )
{
    std::vector< Vector > evecs(3);
    gH.hess_evecs( p, evecs );
    return evecs[idx];
}

crease::Vector crease::gradient( const Point& p, const scalar_wrapper& gH )
{
    Vector g;
    gH.gradient( p, g );
    return g;
}

crease::Vector crease::basis_vector( const Vector& e )
{
    Vector ref( 0, 0, 0 );

    unsigned int minid = 0;
    double d = fabs(e[0]);
    if ( fabs(e[1])<d ) { d = fabs(e[1]); minid = 1; }
    if ( fabs(e[2])<d ) { d = fabs(e[2]); minid = 2; }

    ref[minid] = 1;
    Vector out = cross( e, ref );
    out *= 1/norm( out);

    return out;
}

double crease::zero_crossing( const Vector& ev0, const Vector& ev1,
    const Vector& g0, const Vector& g1 )
{
    double dot0 = inner( g0, ev0 );
    double dot1 = inner( g1, ev1 );
    if ( dot0*dot1 < 0 )
    {
        return -dot0/( dot1-dot0 );
    }

    return -1; // invalid
}

void crease::next_evec( Vector& out, const Vector& in,
    const Edge& edge, const scalar_wrapper& gH,
    unsigned int idx )
{
    out = eigenvector( edge.second, gH, idx );
    if ( nvis::inner(in,out)<0 )
        out *= -1;
}

bool crease::
track_eigenvector_orientation( Vector& out, const Vector& in,
    const Edge& edge, const scalar_wrapper& gH,
    unsigned int idx )
{
    if ( !subdiv )
    {
        next_evec( out, in, edge, gH, idx );
        return true;
    }

    vector< Vector > evecs(3);
    Vector prev, next;
    Edge current;
    prev = in;
    current.first = edge.first;

    double u = 0;
    double v = 0.5;
    if ( eps<0 || eps>=1 ) eps = 0.05;

    while ( u<v && v<=1 )
    {
        current.second = (1-v)*edge.first + v*edge.second;
        next_evec( next, prev, current, gH, idx );
        if ( inner( prev, next )<threshold )
        {
            // angular distance is too large: halve step size
            v = 0.5*(u+v);
            // check for step size underflow
            if ( v-u<eps ) break;

            continue;
        }
        u = v;
        v = std::min( 1., u+0.5 );
        prev = next;
    }

    if ( u==1 )
    {
        out = prev;
        return true;
    }

    return false;
}

bool crease::
orient_face( vector< Vector >& evecs, const vector< Point >& face,
    const scalar_wrapper& gH, unsigned int idx )
{
    evecs.resize(4);

    Vector cur, next;
    vector< Vector > vecs(3);
    Edge edge;

    // initialize orientation
    evecs[0] = cur = eigenvector( face[0], gH, idx );

    // loop over face edges
    for ( unsigned int i=0 ; i<3 ; i++ )
    {
        edge.first = face[i];
        edge.second = face[i+1];
        if ( !track_eigenvector_orientation( next, cur, edge, gH, idx ) )
        {
            return false;
        }

        evecs[i+1] = next;
        cur = next;
    }

    // check that we can consistently close the loop
    edge.first = face[3];
    edge.second = face[0];
    if ( !track_eigenvector_orientation( next, evecs[3], edge, gH, idx ) ||
        inner( next, evecs[0] )<0 )
        return false;

    return true;
}

bool crease::
extract_point_with_valid_emedium( Point& out, const vector< Vector >& evecs,
    const vector< Point >& face,
    const scalar_wrapper& gH, unsigned int idx )
{
    vector< Point > zero_xing;
    Edge edge;
    Vector g0, g1;
    for ( unsigned int i=0 ; i<4 ; i++ )
    {
        edge.first = face[i];
        edge.second = face[(i+1)%4];
        g0 = gradient( edge.first, gH );
        g1 = gradient( edge.second, gH );
        double u = zero_crossing( evecs[i], evecs[(i+1)%4], g0, g1 );
        if ( u>=0 )
        {
            zero_xing.push_back( (1-u)*edge.first + u*edge.second );
        }
    }

    if ( zero_xing.size()==2 )
    {
        Vector ev0, ev1;
        edge.first = zero_xing[0];
        edge.second = zero_xing[1];
        ev0 = eigenvector( edge.first, gH, idx );
        next_evec( ev1, ev0, edge, gH, idx );
        g0 = gradient( edge.first, gH );
        g1 = gradient( edge.second, gH );
        double u = zero_crossing( ev0, ev1, g0, g1 );
        if ( u>=0 )
        {
            out = (1-u)*edge.first + u*edge.second;
            return true;
        }
    }
    else if ( zero_xing.size()!=0 )
    {
        std::cerr << "wrong number of zero crossings on edges: " << zero_xing.size()
            << std::endl;
    }

    return false;
}

bool crease::extract_point_without_valid_emedium( Point& out,
    const vector< Point >& face,
    const scalar_wrapper& gH, unsigned int idx )
{
    vector< Vector > evecs_ref(4);
    vector< Vector > evecs_med(4);
    if ( orient_face( evecs_ref, face, gH, idx ) )
    {
        // we have a valid reference normal vector

        // pick an eigenvector in Span{orth(e_ref)} at each face vertex
        for ( unsigned int i=0 ; i<4 ; i++ )
        {
            evecs_med[i] = basis_vector( evecs_ref[i] );
        }

        vector< pair< double, unsigned int  > > zero_xing;
        Edge edge;
        Vector g0, g1;
        for ( unsigned int i=0 ; i<4 ; i++ )
        {
            edge.first = face[i];
            edge.second = face[(i+1)%4];
            g0 = gradient( edge.first, gH );
            g1 = gradient( edge.second, gH );
            double u = zero_crossing( evecs_med[i], evecs_med[(i+1)%4],
                g0, g1 );
            if ( u>=0 )
            {
                // save local coordinate and edge id
                zero_xing.push_back( pair< double, unsigned int >
                    ( u, i ) );
            }
        }

        if ( zero_xing.size()==2 )
        {
            Vector ev0, ev1, ev0_med, ev1_med, ev0_ref, ev1_ref;
            Vector evec0, evec1;

            double u0 = zero_xing[0].first;
            double u1 = zero_xing[1].first;
            unsigned int i0 = zero_xing[0].second;
            unsigned int i1 = zero_xing[1].second;

            // construct segment connecting both zero crossing positions
            edge.first = (1-u0)*face[i0] + u0*face[(i0+1)%4];
            edge.second = (1-u1)*face[i1] + u1*face[(i1+1)%4];

            // compute the reference and medium eigenvectors associated with
            // both segment vertices. Reference eigenvector is obtained by
            // linear interpolation (consistent with Marching Cubes) while the
            // second one is obtained by linear interpolation + projection
            // onto eigenspace normal to reference vector and subsequent
            // normalization. Finally, the eigenvector to be checked against
            // the gradient is obtained by cross product between reference
            // and medium eigenvector.

            // 1st vertex
            ev0_ref = (1-u0)*evecs_ref[i0] + u0*evecs_ref[(i0+1)%4];
            ev0_med = (1-u0)*evecs_med[i0] + u0*evecs_med[(i0+1)%4];
            ev0_med -= inner( ev0_med, ev0_ref )*ev0_ref;
            ev0_med *= 1/norm( ev0_med );
            evec0 = cross( ev0_ref, ev0_med );

            // 2nd vertex
            ev1_ref = (1-u1)*evecs_ref[i1] + u1*evecs_ref[(i1+1)%4];
            ev1_med = (1-u1)*evecs_med[i1] + u1*evecs_med[(i1+1)%4];
            ev1_med -= inner( ev1_med, ev1_ref )*ev1_ref;
            ev1_med *= 1/norm( ev1_med );
            evec1 = cross( ev1_ref, ev1_med );

            // check for zero crossing of dot product with gradient
            g0 = gradient( edge.first, gH );
            g1 = gradient( edge.second, gH );
            double u = zero_crossing( evec0, evec1, g0, g1 );

            if ( u>=0 )
            {
                out = (1-u)*edge.first + u*edge.second;
                return true;
            }
        }
        else if ( zero_xing.size()!=0 )
        {
            std::cerr << "wrong number of zero crossings on edges: " << zero_xing.size()
                << std::endl;
        }

        return false;
    }
    else
    {
        std::cerr << "unable to orient eigenvector provided as reference"
            << std::endl;
    }

    return false;
}

void crease::connect_segments()
{
    // compute point <- point -> point connections
    std::vector< std::pair< int, int > > connected_to( face_found.size(),
        std::pair< int, int >( -1, -1 ) );
    for ( unsigned int i=0 ; i<all_edges.size() ; i++ )
    {
        // face indices for current edge
        unsigned int f0, f1;
        f0 = all_edges[i].first;
        f1 = all_edges[i].second;

        int a, b;
        a = connected_to[f0].first;
        b = connected_to[f1].first;

        if ( a == -1 )
            connected_to[f0].first = f1;
        else
            connected_to[f0].second = f1;

        if ( b == -1 )
            connected_to[f1].first = f0;
        else
            connected_to[f1].second = f0;
    }

    components.clear();
    std::vector< bool > inserted( all_edges.size(), false );
    for ( unsigned int i=0 ; i<all_edges.size() ; i++ )
    {
        if ( inserted[i] ) continue;

        unsigned int cur, prev;
        int link0, link1;

        // start a new connected component
        components.push_back( std::list< unsigned int >() );
        std::list< unsigned int >& my_list = components.back();

        // edge end points
        unsigned int i0, i1;
        i0 = all_edges[i].first;
        i1 = all_edges[i].second;

        // initialize connected component with these two points
        my_list.push_back( i0 );
        my_list.push_back( i1 );

        // append forward
        prev = i0;
        cur = i1;
        inserted[prev] = true;
        for ( ; ; )
        {
            inserted[cur] = true;
            link0 = connected_to[cur].first;
            link1 = connected_to[cur].second;
            if ( ( unsigned int )link0 != prev && !inserted[link0] )
                my_list.push_back( link0 );
            else if ( link1 >= 0 && ( unsigned int )link1 != prev && !inserted[link1] )
                my_list.push_back( link1 );
            else break;
            prev = cur;
            cur = my_list.back();
        }

        // append backward
        cur = i0;
        prev = i1;
        for ( ; ; )
        {
            inserted[cur] = true;
            link0 = connected_to[cur].first;
            link1 = connected_to[cur].second;
            if ( ( unsigned int )link0 != prev && !inserted[link0] )
                my_list.push_front( link0 );
            else if ( link1 >= 0 && ( unsigned int )link1 != prev && !inserted[link1] )
                my_list.push_front( link1 );
            else break;
            prev = cur;
            cur = my_list.front();
        }
    }
}

void crease::extract_lines( const Nrrd* nrrd, bool ridge )
{
    // check if this is a tensor field
    bool is_tensor = false;
    if ( nrrd->dim == 4 )
        is_tensor = true;

    // collect information about size and bounding box
    unsigned int size[3];
    double min[3], max[3];
    unsigned int shift = ( is_tensor ? 1 : 0 );
    for ( unsigned int i=0 ; i<3 ; i++ )
    {
        min[i] = nrrd->axis[i+shift].min;
        max[i] = nrrd->axis[i+shift].max;
        size[i] = nrrd->axis[i+shift].size;
    }
    Grid grid( size, min, max );

    scalar_wrapper gH( nrrd );

    set< FaceId > face_done;
    all_face_points.clear();
    all_edges.clear();
    face_found.clear();

    // loop over mesh cells
    unsigned int faces[6][4] =
        { { 0, 1, 2, 3 },
        { 4, 5, 6, 7 },
        { 0, 1, 5, 4 },
        { 1, 2, 6, 5 },
        { 2, 3, 7, 6 },
        { 3, 0, 4, 7 } };

    unsigned int nbok = 0;
    unsigned int nb_segments = 0;
    unsigned int nb_fakebasis = 0;
    unsigned int nb_orientable = 0;

    // compute some statistics
    vec3 evals;
    vector< double > means;
    double mean = 0;

    //unsigned int counter=0;
    for ( unsigned int i=0 ; i<size[0] ; i++ )
    {
        for ( unsigned int j=0 ; j<size[1] ; j++ )
        {
            for ( unsigned int k=0 ; k<size[2] ; k++ )
            {
                gH.hess_evals( grid(i,j,k), evals );
                if ( ( ridge && evals[1]<0 ) ||
                    ( !ridge && evals[1]>0 ) )
                {
                    mean += evals[1];
                    means.push_back( evals[1] );
                }
            }
        }
    }

    std::cout << means.size() << " ("
        << ( float )means.size()/( float )( size[0]*size[1]*size[2] )*100.
        << "%) candidate vertices" << std::endl;

    sort( means.begin(), means.end() );

    mean /= ( float )means.size();
    std::cout << "mean medium valid eigenvalue = " << mean << std::endl;

    std::cout << "median valid eigenvalue = " << means[( int )(0.5*means.size())] << std::endl;

    if ( !ridge )
        std::cout << "max eigenvalue = " << means.back() << std::endl;
    else
        std::cout << "min eigenvalue = " << means[0] << std::endl;

    float answer;
    std::cout << "threshold (in %): ";
    std::cin >> answer;

    threshold = answer;

    unsigned int _id = ( !ridge
        ? ( unsigned int )( means.size()*threshold )
        : ( unsigned int )( means.size()*( 1-threshold ) ) );
    double _threshold = means[_id];

    std::cout << "corresponding threshold = " << _threshold << std::endl;
    std::this_thread::sleep_for(5000ms);
    std::cout << std::endl;

    cout << "threshold = " << _threshold << " (" << threshold << ")" << endl;
    cout << "median valid value =  " << means[means.size()/2] << endl;

    for ( unsigned int k=0 ; k<size[2]-1 ; k++ )
    {
        for ( unsigned int j=0 ; j<size[1]-1 ; j++ )
        {
            for ( unsigned int i=0 ; i<size[0]-1 ; i++ )
            {
                // voxel vertices
                vector< vec3 > v(8);
                v[0] = grid( i  , j  , k   );
                v[1] = grid( i+1, j  , k   );
                v[2] = grid( i+1, j+1, k   );
                v[3] = grid( i  , j+1, k   );
                v[4] = grid( i  , j  , k+1 );
                v[5] = grid( i+1, j  , k+1 );
                v[6] = grid( i+1, j+1, k+1 );
                v[7] = grid( i  , j+1, k+1 );

                vector< unsigned int > ids(8);
                ids[0] = grid.id( i  , j  , k   );
                ids[1] = grid.id( i+1, j  , k   );
                ids[2] = grid.id( i+1, j+1, k   );
                ids[3] = grid.id( i  , j+1, k   );
                ids[4] = grid.id( i  , j  , k+1 );
                ids[5] = grid.id( i+1, j  , k+1 );
                ids[6] = grid.id( i+1, j+1, k+1 );
                ids[7] = grid.id( i  , j+1, k+1 );

                // check sign of eigenvalues
                vec3 evals;
                bool candidate = true;
                for ( unsigned int l=0 ; l<8 ; l++ )
                {
                    gH.hess_evals( v[l], evals );
                    if ( ( ridge && evals[1]>=_threshold ) ||
                        ( !ridge && evals[1]<=_threshold ) )
                    {
                        candidate = false;
                        break;
                    }
                }
                if ( !candidate ) continue;

                //       cout << "valid voxel" << endl;

                nbok++;

                // loop over voxel faces
                vector< unsigned int > face_point;
                for ( unsigned int f=0 ; f<6 ; f++ )
                {
                    // face vertices
                    vector< vec3 > p(4);
                    for ( unsigned int i=0 ; i<4 ; i++ )
                        p[i] = v[faces[f][i]];

                    // check if we have processed that face already
                    FaceId fid( ids[faces[f][0]],
                        ids[faces[f][1]],
                        ids[faces[f][2]],
                        ids[faces[f][3]] );
                    if ( face_done.find(fid) != face_done.end() )
                    {
                        // add corresponding position, if any
                        map< FaceId, unsigned int >::const_iterator it =
                            face_found.find(fid);
                        if ( it != face_found.end() )
                        {
                            face_point.push_back( it->second );
                        }

                        // move to next voxel face
                        continue;
                    }

                    // mark current face as processed
                    face_done.insert( fid );

                    vector< Vector > evecs(4);
                    if ( false && orient_face( evecs, p, gH, 1 ) )
                    {
                        ++nb_orientable;
                        Point xing;
                        if ( extract_point_with_valid_emedium( xing, evecs, p, gH,
                            ( ridge ? 2 : 0 ) ) )
                        {
                            all_face_points.push_back( xing );
                            face_found[fid] = all_face_points.size()-1;
                            face_point.push_back( all_face_points.size()-1 );
                        }
                    }
                    else
                    {
                        ++nb_fakebasis;
                        Point xing;
                        // here we need to resort to plan B which consists in
                        // creating an artificial but consistent (e2,e3) basis
                        // of Span(orth(e1)) over the whole face and use it instead
                        // of the eigenvectors computed by Gage.
                        if ( extract_point_without_valid_emedium( xing, p, gH,
                            ( ridge ? 0 : 2 ) ) )
                        {
                            all_face_points.push_back( xing );
                            face_found[fid] = all_face_points.size()-1;
                            face_point.push_back( all_face_points.size()-1 );
                        }
                    }
                }

                if ( face_point.size() &&
                    face_point.size() != 2 )
                {
                    cout << "found " << face_point.size()
                        << " crease point(s) on voxel faces"
                        << endl;
                }
                else if ( face_point.size() == 2 )
                {
                    ++nb_segments;
                    cout << "found a crease line!!" << endl;
                    all_edges.push_back( pair< unsigned int, unsigned int >
                        ( face_point[0], face_point[1] ) );
                }
            }
        }
    }

    connect_segments();

    unsigned int nb_cells = ( size[0]-1 )*( size[1]-1 )*( size[2]-1 );
    cout << "fake eigenspace basis: " << nb_fakebasis
        << " / orientanble eigenvectors: " << nb_orientable
        << " ("
        << ( double )nb_fakebasis/( double )nb_orientable*100
        << "%)" << endl;
    cout << "percentage of candidate voxels: "
        << ( double )nbok/( double )nb_cells
        << endl;

    cout << "number of segments = " << all_edges.size()
        << ", number of connected components: " << components.size()
        << std::endl;
}
