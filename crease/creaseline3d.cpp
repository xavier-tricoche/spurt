#include "creaseline3d.hpp"
#include <algorithm>

using namespace std;
using namespace nvis;
using namespace spurt;
using namespace gage_interface;

/*
comments:

- fake eigenbasis is not guaranteed to be handled consistently throughout the computation
    - the basic logic remains to determine point-wise values at vertices and to combine interpolation
    + projection to determine intermediate values along edges
    - computing zero crossings may require to track involved vectors along the way. relying on linear
    interpolation seems unsatisfactory

*/


double spurt::crease::threshold;
double spurt::crease::eps;
unsigned int spurt::crease::subdiv;
unsigned int spurt::crease::upsample=1;

vector< crease::Edge > crease::problematic_edges;
vector< vec3 > crease::all_face_points;
vector< double > crease::crease_strength;
vector< double > crease::grad_dot_evec;
vector< double > crease::measure_value;

vector< vector< vector< spurt::crease::value > > > crease::nonorientable_faces;
vector< vector< vec3 > > crease::nonorientable_face_points;

vector< pair< unsigned int, unsigned int > > crease::all_edges;
map< crease::FaceId, unsigned int > face_found;
vector< vector< unsigned int > > crease::components;
crease::MeasureWrapper* the_wrapper;

bool verbose_mode = false;
bool _ridge;

unsigned int nb_evaluations = 0;
unsigned int nb_subdivided = 0;

// following values record required sampling resolution to capture
// consistent orientation of eigenvectors
double min_dus[4];
double min_du;

inline bool val_ok( double val, double threshold, bool ridge )
{
    if ( ridge ) return val>threshold;
    else return val<threshold;
}

inline bool eval_ok( double eval, double threshold, bool ridge )
{
    if ( ridge ) return eval<threshold;
    else return eval>threshold;
}

template< typename T >
void spurt::display_stats( const std::vector< T >& data, const std::string& msg )
{
    // examine the value distribution of the considered scalar quantity
    std::vector< T > vals( data.begin(), data.end() );
    T mean = 0;
    unsigned int sz = data.size();

    if ( sz<20 )
    {
        for ( unsigned int i=0 ; i<sz ; i++ )
            std::cout << "#" << i << ": " << data[i] << std::endl;
        return;
    }

    // compute mean
    for ( unsigned int i=0 ; i<sz ; i++ )
    {
        mean += data[i];
    }
    mean /= ( T )data.size();

    std::sort( vals.begin(), vals.end() );
    // compute variance
    T variance = 0;
    for ( unsigned int i=0 ; i<vals.size() ; i++ )
    {
        variance += ( vals[i]-mean )*( vals[i]-mean );
    }
    variance /= ( T )vals.size();

    unsigned int incr = vals.size()/100;
    std::cout << "Information on " << msg << ": " << std::endl
        << "mean: " << mean << ", variance: " << variance
        << ", median: " << vals[vals.size()/2] << std::endl
        << "distribution: "
        << vals[0] << " (0%), "
        << vals[1*incr] << " (1%), "
        << vals[2*incr] << " (2%), "
        << vals[5*incr] << " (5%), " << std::endl
        << vals[10*incr] << " (10%), "
        << vals[15*incr] << " (15%), "
        << vals[25*incr] << " (25%), " << std::endl
        << vals[50*incr] << " (50%), "
        << vals[75*incr] << " (75%), "
        << vals[85*incr] << " (85%), " << std::endl
        << vals[90*incr] << " (90%), "
        << vals[95*incr] << " (95%), "
        << vals[98*incr] << " (98%), "
        << vals[99*incr] << " (99%), "
        << vals.back() << " (100%)"
        << std::endl;
}

// --------------------------------------------------------------------------

crease::MeasureWrapper::MeasureWrapper( const Nrrd* nrrd, int aniso )
    : _scal_wrap(0), _tens_wrap(0), _evecs(3), _evals(3), _measure(aniso)
{
        // check if this is a tensor field
    bool is_tensor = ( nrrd->dim == 4 && nrrd->axis[0].size == 7 );

    std::cout << "measure is " << _measure << std::endl;

    if ( is_tensor )
    {
        assert( _measure == 1 || _measure == 2 );
        _tens_wrap = new tensor_wrapper( nrrd );
    }
    else
    {
        assert( _measure == 0 );
        _scal_wrap = new scalar_wrapper( nrrd );
    }
}

int crease::MeasureWrapper::what_measure() const
{
    return _measure;
}

double crease::MeasureWrapper::confidence( const vec3& p ) const
{
    double v;
    switch ( _measure )
    {
        case 0:
        {
            v = 1;
            break;
        }
        case 1:
        case 2:
        {
            _tens_wrap->confidence( p, v );
            break;
        }
        default:
        {
            assert( false );
        }
    }

    return v;
}

double crease::MeasureWrapper::value( const vec3& p ) const
{
    double v;
    switch ( _measure )
    {
        case 0:
        {
            _scal_wrap->value( p, v );
            break;
        }
        case 1:
        {
            _tens_wrap->fa( p, v );
            break;
        }
        case 2:
        {
            _tens_wrap->mode( p, v );
            break;
        }
        default:
        {
            assert( false );
        }
    }

    return v;
}

vec3 crease::MeasureWrapper::eigenvector( const vec3& p, unsigned int idx ) const
{
    switch ( _measure )
    {
        case 0:
        {
            _scal_wrap->hess_evecs( p, _evecs );
            break;
        }
        case 1:
        {
            _tens_wrap->fa_hess_evecs( p, _evecs );
            break;
        }
        case 2:
        {
            _tens_wrap->mode_hess_evecs( p, _evecs );
            break;
        }
        default:
        {
            assert( false );
        }
    }

    return _evecs[idx];
}

double crease::MeasureWrapper::eigenvalue( const vec3& p, unsigned int idx ) const
{
    ++nb_evaluations;
    nvis::vec3 _tmp;
    switch ( _measure )
    {
        case 0:
        {
            _scal_wrap->hess_evals( p, _tmp );
            break;
        }
        case 1:
        {
            _tens_wrap->fa_hess_evals( p, _tmp );
            break;
        }
        case 2:
        {
            _tens_wrap->mode_hess_evals( p, _tmp );
            break;
        }
        default:
        {
            assert( false );
        }
    }
    return _tmp[idx];
}

vec3 crease::MeasureWrapper::gradient( const vec3& p ) const
{
    switch ( _measure )
    {
        case 0:
        {
            _scal_wrap->gradient( p, _grad );
            break;
        }
        case 1:
        {
            _tens_wrap->fa_grad( p, _grad );
            break;
        }
        case 2:
        {
            _tens_wrap->mode_grad( p, _grad );
            break;
        }
        default:
        {
            assert( false );
        }
    }
    return _grad;
}

vec6 crease::MeasureWrapper::hessian( const vec3& p ) const
{
    ++nb_evaluations;
    nvis::vec6 _tmp;
    switch ( _measure )
    {
        case 0:
        {
            _scal_wrap->hessian( p, _tmp );
            break;
        }
        case 1:
        {
            _tens_wrap->fa_hess( p, _tmp );
            break;
        }
        case 2:
        {
            _tens_wrap->mode_hess( p, _tmp );
            break;
        }
        default:
        {
            assert( false );
        }
    }
    return _tmp;
}

// --------------------------------------------------------------------------

bool basis_vectors( std::vector< vec3 >& basis,
    const std::vector< vec3 >& evecs )
{
    // assumptions:
    // 1) eigenvector field varies rougly linearly over the face
    // 2) solid angle spanned by eigenvector covers (much) less than half of S^2

    // identify the solid angle spanned by the 4 vectors
    bool valid[3] = { true, true, true };

    // simple triangulation of gauss map
    unsigned int tris[2][3] = { { 0, 1, 2 }, { 0, 2, 3 } };
    std::vector< vec3 > P(3);
    for ( unsigned int i=0 ; i<2 ; i++ )
    {
        for ( unsigned int j=0 ; j<3 ; j++ )
            P[j] = evecs[tris[i][j]];

        // ( t*e_i - P0 ).normal = 0 => t*e_i.normal = P0.normal
        vec3 normal = nvis::cross( P[1]-P[0], P[2]-P[0] );
        double a = nvis::inner( P[0], normal );

        // loop over 3 axes
        for ( unsigned int l=0 ; l<3 ; l++ )
        {
            if ( !valid[l] ) continue;

            /*
            double t = a/normal[l];
            vec3 Q( 0, 0, 0 );
            Q[l] = t; // e_l[i] = \delta(i,l)
            */

            // project triangle onto plane orthogonal to e_l
            std::vector< vec2 > p(3);
            for ( unsigned int m=0 ; m<3 ; m++ )
            {
                unsigned int c=0;
                for ( unsigned int n=0 ; n<3 ; n++ )
                {
                    if ( n == l ) continue; // skip l
                    p[m][c++] = P[m][n];
                }
            }

            // compute barycentric coordinates
            double delta = nvis::cross( p[1]-p[0], p[2]-p[0] );
            if ( delta == 0 ) continue;

            // b1*(p1-p0) + b2*( p2-p0 ) = q-p0
            vec2 q(0,0);

            double b1 = nvis::cross( q-p[0], p[2]-p[0] )/delta;
            if ( b1<0 || b1>1 ) continue;
            double b2 = nvis::cross( p[1]-p[0], q-p[0] )/delta;
            if ( b2<0 || b2>1 ) continue;

            valid[l] = false;
        }
    }

    bool done = false;
    for ( unsigned int i=0 ; i<3 ; i++ )
    {
        if ( !valid[i] ) continue;

        vec3 e( 0, 0, 0);
        e[i] = 1;
        for ( unsigned int j=0 ; j<4 ; j++ )
        {
            basis[j] = nvis::cross( e, evecs[j] );
            basis[j] *= 1/nvis::norm( basis[j] );
        }

        return true;
    }

    return false;
}

vec3 crease::basis_vector( const vec3& e )
{
    vec3 ref( 0, 0, 0 );

    unsigned int minid = 0;
    double d = fabs(e[0]);
    if ( fabs(e[1])<d ) { d = fabs(e[1]); minid = 1; }
    if ( fabs(e[2])<d ) { d = fabs(e[2]); minid = 2; }

    ref[minid] = 1;
    nvis::vec3 out = cross( e, ref );
    out *= 1/norm( out);

    return out;
}

double crease::zero_crossing( const vec3& ev0, const vec3& ev1,
    const vec3& g0, const vec3& g1 )
{
    double dot0 = inner( g0, ev0 );
    double dot1 = inner( g1, ev1 );
    if ( dot0*dot1 < 0 )
    {
        return -dot0/( dot1-dot0 );
    }

    return -1; // invalid
}


inline void __next_evec( nvis::vec3& out, const nvis::vec3& in,
    const crease::Edge& edge, unsigned int idx, double u )
{
    nvis::vec3 p = (1.0-u)*edge.first + u*edge.second;
    out = the_wrapper->eigenvector( p, idx );
    if ( nvis::inner(in,out)<0 )
        out *= -1;
}

void crease::next_evec( vec3& out, const vec3& in,
    const Edge& edge, unsigned int idx )
{
    out = the_wrapper->eigenvector( edge.second, idx );
    if ( nvis::inner(in,out)<0 )
        out *= -1;
}

bool track_evec_on_edge( std::vector< spurt::crease::value >& evecs,
    const crease::Edge& edge, unsigned int idx )
{
    using namespace crease;

    evecs.clear();
    vec3 ev0, ev1;

    ev0 = the_wrapper->eigenvector( edge.first, idx );
    evecs.push_back( value( ev0, 0 ) );

    // if no subdivision is requested, assume trivial orientation consistency
    // and quasi-linear behavior of both vector quantities
    if ( !subdiv )
    {
        next_evec( ev1, ev0, edge, idx );
        evecs.push_back( value( ev1, 1 ) );

//        std::cout << evecs.size() << " vecs" << std::endl;
        return true;
    }

    double u0, u1;
    u0 = 0;
    u1 = 0.5;
    if ( eps < 0 || eps >= 1 ) eps = 0.05;

    while ( u0<u1 )
    {
        __next_evec( ev1, ev0, edge, idx, u1 );

        // can we trust the new found
        if ( inner( ev0, ev1 ) < 0.9 )
        {
            ++nb_subdivided;
            // angular distance is too large: halve step size
            u1 = 0.5*( u0+u1 );
            if ( u1-u0<eps )
            {
//                std::cout << evecs.size() << " vecs #" << std::endl;
                return false; // step size underflow
            }
        }
        evecs.push_back( value( ev1, u1 ) );

        u0 = u1;
        u1 = std::min( 1., u0+0.5 );
        ev0 = ev1;
    }

//    std::cout << evecs.size() << " vecs" << std::endl;
    return true;
}

bool crease::
track_eigenvector_orientation( vec3& out, const vec3& in,
    const Edge& edge, unsigned int idx )
{
    min_du = 1;
    if ( !subdiv )
    {
        next_evec( out, in, edge, idx );
        return true;
    }

    vector< vec3 > evecs(3);
    vec3 prev, next;
    Edge current;
    prev = in;
    current.first = edge.first;

    double u = 0;
    double v = 0.5;
    if ( eps<0 || eps>=1 ) eps = 0.05;

    while ( u<v )
    {
        if ( v-u < min_du ) min_du = v-u;

        current.second = (1-v)*edge.first + v*edge.second;
        next_evec( next, prev, current, idx );
        if ( inner( prev, next )<threshold )
        {
            ++nb_subdivided;
            // angular distance is too large: halve step size
            v = 0.5*(u+v);
            // check for step size underflow
            if ( v-u<eps ) break;

            continue;
        }
        u = v;
        v = std::min( 1., u+0.5 );
        prev = next;
        current.first = current.second;
    }

    if ( u==1 )
    {
        out = prev;
        return true;
    }

    return false;
}

bool crease::
orient_face( vector< vec3 >& evecs, const vector< vec3 >& face,
    unsigned int idx )
{
    evecs.resize(4);

    vec3 cur, next;
    vector< vec3 > vecs(3);
    Edge edge;

    // initialize orientation
    evecs[0] = cur = the_wrapper->eigenvector( face[0], idx );

    // loop over face edges
    for ( unsigned int i=0 ; i<3 ; i++ )
    {
        edge.first = face[i];
        edge.second = face[i+1];
        if ( !track_eigenvector_orientation( next, cur, edge, idx ) )
        {
            return false;
        }

        min_dus[i] = min_du;

        evecs[i+1] = next;
        cur = next;
    }

    // check that we can consistently close the loop
    edge.first = face[3];
    edge.second = face[0];
    if ( !track_eigenvector_orientation( next, evecs[3], edge, idx ) ||
        inner( next, evecs[0] )<0 )
        return false;

    min_dus[3] = min_du;

    return true;
}

void zero_crossing_on_edge( std::vector< nvis::vec3 >& zeros, const crease::Edge& edge,
    const std::vector< crease::value >& evecs )
{
    using namespace crease;

    zeros.clear();

    double dot0, dot1, u0, u1;
    nvis::vec3 ev0, ev1, g0, g1;

    ev0 = evecs[0].first;
    g0 = the_wrapper->gradient( edge.first );
    dot0 = nvis::inner( ev0, g0 );
    u0 = 0;
    for ( unsigned int i=1 ; i<evecs.size() ; i++ )
    {
        u1 = evecs[i].second;
        assert( u1>=0 && u1<=1 );
        nvis::vec3 p = ( 1.-u1 )*edge.first + u1*edge.second;
        g1 = the_wrapper->gradient( p );
        dot1 = nvis::inner( evecs[i].first, g1 );

        if ( dot0*dot1<0 )
        {
            double v = -dot0/( dot1-dot0 );
            double u = ( 1.-v )*u0 + v*u1;
            assert( u>=0 && u<=1 );
            zeros.push_back( (1.-u)*edge.first + u*edge.second );
        }
        dot0 = dot1;
        u0 = u1;
    }
}

int
    all_in_one_face( nvis::vec3& crossing, const vector< nvis::vec3 >& face, bool ridge )
{
    return -2;

    using namespace crease;

    Edge edge;
    std::vector< std::vector< value > > edge_values(4);

    unsigned int idx1, idx2;
    idx1 = ( ridge ? 2 : 0 ); // looking for crease surface first
    idx2 = 1;

    // loop over face edges
    for ( unsigned int i=0 ; i<4 ; i++ )
    {
        edge.first = face[i];
        edge.second = face[(i+1)%4];
        bool worked = track_evec_on_edge( edge_values[i], edge, idx1 );
        if ( !worked )
        {
            if ( verbose_mode )
                std::cout << "unable to orient minor eigenvector field on face (edge)"
                << std::endl;
            return -2;
        }
    }

    // check that we have a globally consistent orientation
    unsigned int nb_switch = 0;
    for ( unsigned int i=0 ; i<4 ; i++ )
    {
        if ( nvis::inner( edge_values[i].back().first, edge_values[(i+1)%4].front().first ) < 0 )
            ++nb_switch;
    }
    if ( nb_switch % 2 )
    {
        if ( verbose_mode )
            std::cout << "unable to orient minor eigenvector field on face (loop)"
            << std::endl;

        nonorientable_faces.push_back( std::vector< std::vector< value > >() );
        std::vector< std::vector< value > >& ref = nonorientable_faces.back();
        for ( unsigned int i=0 ; i<4 ; i++ )
        {
            ref.push_back( std::vector< value >() );
//            std::cout << "adding " << edge_values[i].size() << " vectors for this edge" << std::endl;
            for ( unsigned int j=0 ; j<edge_values[i].size() ; j++ )
            {
                ref.back().push_back( value( edge_values[i][j].first, edge_values[i][j].second ) );
            }
        }
        nonorientable_face_points.push_back( face );
        return -2;
    }

    // compute zero crossings along egdes
    std::vector< std::vector< nvis::vec3 > > zeros(4);
    unsigned int nb_zeros = 0;

    for ( unsigned int i=0 ; i<4 ; i++ )
    {
        edge.first = face[i];
        edge.second = face[(i+1)%4];
        zero_crossing_on_edge( zeros[i], edge, edge_values[i] );
        nb_zeros += zeros[i].size();
    }

    if ( nb_zeros != 2 ) return -1;

    // identify all possible edges of crease surface
    std::vector< Edge > all_edges;
    for ( unsigned int i=0 ; i<3 ; i++ )
    {
        if ( !zeros[i].size() ) continue;

        for ( unsigned int ii=0 ; ii<zeros[i].size() ; ii++ )
        {
            edge.first = zeros[i][ii];
            for ( unsigned int j=i+1 ; j<4 ; j++ )
            {
                if ( !zeros[j].size() ) continue;

                for ( unsigned int jj=0 ; jj<zeros[j].size() ; jj++ )
                {
                    edge.second = zeros[j][jj];
                    all_edges.push_back( edge );
                }
            }
        }
    }
    if ( !all_edges.size() ) return 0;

    assert( all_edges.size() == 1 );

    if ( verbose_mode && all_edges.size() > 1 )
        std::cout << all_edges.size() << " edges to examine on face!" << std::endl;

    // look for zero crossings along all edges
    std::vector< nvis::vec3 > __found;
    std::vector< value > vals;
    unsigned int nb_found = 0;
    for ( unsigned int i=0 ; i<all_edges.size() ; i++ )
    {
        edge.first = all_edges[i].first;
        edge.second = all_edges[i].second;
        if ( !track_evec_on_edge( vals, edge, idx2 ) )
        {
            return -1;
        }
        zero_crossing_on_edge( __found, edge, vals );
        if ( nb_found + __found.size() > 1 )
        {
            if ( verbose_mode )
                std::cout << "multiple crease point on face" << std::endl;
            return -1;
        }
        if ( __found.size() == 1 )
        {
            ++nb_found;
            crossing = __found[0];
            if ( verbose_mode )
                std::cout << "crossing = " << crossing << std::endl;
        }
    }

    if ( nb_found ) std::cout << "found a point on orientable face" << std::endl;
    return nb_found; // either 0 or 1
}

/*
bool newton( nvis::vec3& p3d, const vector< nvis::vec3 >& face,
    const vec3& b0, const vec3& b1 )
{
    // compute 2D basis
    vec3 e0 = face[1]-face[0];
    vec3 e1 = face[3]-face[0];
    vec3 n = nvis::cross( e0, e1 );
    double _norm = norm( n );
    unsigned int dim=0;
    for ( unsigned int i=1 ; i<3 ; i++ )
    {
        if ( fabs(n[i])>0.99*_norm )
        {
            dim = i;
            break;
        }
    }
    vec3 b0, b1; // assuming they are set to 0;
    b0[(dim+1)%3] = 1;
    b1[(dim+2)%3] = 1;

    p3d = 0.5*( face[0]+face[2] ); // middle of face is first guess
    vec2 p( inner( p3d-face[0], b0 ), inner( p3d-face[0], b1 ) );

    vec6 _hess6d;
    vec3 _grad3d;

    while ( true )
    {
        _grad3d = the_wrapper->gradient( p3d );
        _hess6d = the_wrapper->hessian( p3d );
        mat3 _hess3d;
        unsigned int c=0;
        for ( unsigned int i=0 ; i<3 ; i++ )
        {
            for ( unsigned int j=i ; j<3 ; j++ )
            {
                _hess3d[i][j] = _hess3d[j][i] = _hess6d[c++];
            }
        }

        vec2 grad( inner( _grad3d, b0 ), inner( _grad3d, b1 ) );
        mat2 hess;
        vec3 _tmp = _hess3d*b0;
        hess[0][0] = inner( _tmp, b0 );
        hess[0][1] = hess[1][0] = inner( _tmp, b1 );
        _tmp = _hess3d*b1;
        hess[1][1] = inner( _tmp, b1 );

        // now we have 2D problem

        // invert 2x2 hessian
        double delta = inner::det( hess );
        if ( !delta ) retun false;

        hess *= 1/delta;
        mat2 inv_hess;
        inv_hess[0][0] = hess[1][1];
        inv_hess[0][1] = inv_hess[1][0] = -hess[0][1];
        inv_hess[1][1] = hess[0][0];

        // compute step
        vec2 step = inv_hess*grad;
    }
}

int
    all_in_one_face_extreme_mode( nvis::vec3& crossing, const vector< nvis::vec3 >& face, bool ridge )
{
    // look for mode = +/- 1 on face

    // solve as an optimization problem

    using namespace nvis;


}
*/

bool crease::
extract_point_with_valid_emedium( vec3& out, const vector< vec3 >& evecs,
    const vector< vec3 >& face, unsigned int idx )
{
    vector< vec3 > zero_xing;
    Edge edge;
    vec3 g0, g1;
    for ( unsigned int i=0 ; i<4 ; i++ )
    {
        edge.first = face[i];
        edge.second = face[(i+1)%4];
        g0 = the_wrapper->gradient( edge.first );
        g1 = the_wrapper->gradient( edge.second );
        double u = zero_crossing( evecs[i], evecs[(i+1)%4], g0, g1 );
        if ( u>=0 && u<=1 )
        {
            zero_xing.push_back( (1-u)*edge.first + u*edge.second );
        }
    }

    if ( zero_xing.size()==2 )
    {
        vec3 ev0, ev1;
        edge.first = zero_xing[0];
        edge.second = zero_xing[1];
        ev0 = the_wrapper->eigenvector( edge.first, idx );
        if ( !track_eigenvector_orientation( ev1, ev0, edge, idx ) )
            return false;
        g0 = the_wrapper->gradient( edge.first );
        g1 = the_wrapper->gradient( edge.second );
        double u = zero_crossing( ev0, ev1, g0, g1 );
        if ( u>=0 && u<=1 )
        {
            out = (1-u)*edge.first + u*edge.second;
            return true;
        }
    }
    else if ( zero_xing.size()!=0 )
    {
        if ( verbose_mode )
            std::cerr << "wrong number of zero crossings on edges: " << zero_xing.size()
            << std::endl;
    }

    return false;
}

bool crease::extract_point_without_valid_emedium( vec3& out,
    const vector< vec3 >& face, unsigned int idx )
{
    assert( ( _ridge && idx==0 ) || ( !_ridge && idx==2 ) );

    vector< vec3 > evecs_ref(4);
    vector< vec3 > evecs_med(4);
    if ( orient_face( evecs_ref, face, idx ) )
    {
        // we have a valid reference normal vector

/*
        // pick an eigenvector in Span{orth(e_ref)} at each face vertex
        for ( unsigned int i=0 ; i<4 ; i++ )
        {
            evecs_med[i] = basis_vector( evecs_ref[i] );
        }
*/
        if ( !basis_vectors( evecs_med, evecs_ref ) )
        {
            if ( verbose_mode )
            {
                std::cout << "unable to assign fake basis!" << std::endl;
                double val=0, aval;
                for ( unsigned int i=0 ; i<4 ; i++ )
                {
                    aval = the_wrapper->value( face[i] );
                    val += aval;
                }
                std::cout << "average value of considered measure is: " << 0.25*val << std::endl;
            }
            return false;
        }

        vector< pair< double, unsigned int  > > zero_xing;
        Edge edge;
        vec3 g0, g1;
        for ( unsigned int i=0 ; i<4 ; i++ )
        {
            edge.first = face[i];
            edge.second = face[(i+1)%4];
            g0 = the_wrapper->gradient( edge.first );
            g1 = the_wrapper->gradient( edge.second );

            // we need a smart zero crossing
            double u = zero_crossing( evecs_med[i], evecs_med[(i+1)%4],
                g0, g1 );
            if ( u>=0 )
            {
                // checking
                vec3 p = ( 1.-u)*edge.first + u*edge.second;
                vec3 g = the_wrapper->gradient( p );
                g *= 1./( double )nvis::norm( g );
                vec3 e1 = the_wrapper->eigenvector( p, idx );
                vec3 e2 = ( 1.-u )*evecs_med[i] + u*evecs_med[(i+1)%4];
                e2 -= nvis::inner( e2, e1 )*e1;
                e2 *= 1./( double )nvis::norm( e2 );

                if ( verbose_mode && fabs( nvis::inner( g, e2 ) )>0.05 )
                    std::cout << "WARNING: < g, e_medium > = "
                    << fabs( nvis::inner( g, e2 ) ) << std::endl;

                // save local coordinate and edge id
                zero_xing.push_back( pair< double, unsigned int >
                    ( u, i ) );
            }
        }

        if ( zero_xing.size()==2 )
        {
            vec3 ev0, ev1, ev0_med, ev1_med, ev0_ref, ev1_ref;
            vec3 evec0, evec1;

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

            if ( false )
            {
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
            }
            else
            {
                // reference (major/minor) eigenvector
                ev0_ref = the_wrapper->eigenvector( edge.first, idx );
                if ( !track_eigenvector_orientation( ev1_ref, ev0_ref, edge, idx ) )
                    return false;

                // medium eigenvectors
                ev0_med = (1-u0)*evecs_med[i0] + u0*evecs_med[(i0+1)%4];
                ev0_med -= inner( ev0_med, ev0_ref )*ev0_ref;
                ev0_med *= 1/norm( ev0_med );
                evec0 = cross( ev0_ref, ev0_med );

                ev1_med = (1-u1)*evecs_med[i1] + u1*evecs_med[(i1+1)%4];
                ev1_med -= inner( ev1_med, ev1_ref )*ev1_ref;
                ev1_med *= 1/norm( ev1_med );
                evec1 = cross( ev1_ref, ev1_med );
            }

            // check for zero crossing of dot product with gradient
            g0 = the_wrapper->gradient( edge.first );
            g1 = the_wrapper->gradient( edge.second );
            double u = zero_crossing( evec0, evec1, g0, g1 );

            if ( u>=0 && u<=1 )
            {
                out = (1-u)*edge.first + u*edge.second;
                if ( verbose_mode )
                    std::cout << "found point on non-orientable face" << std::endl;
                return true;
            }
        }
        else if ( zero_xing.size()!=0 )
        {
            if ( verbose_mode )
                std::cerr << "wrong number of zero crossings on edges: " << zero_xing.size()
                << std::endl;
        }

        return false;
    }
    else
    {
        if ( verbose_mode )
            std::cerr << "unable to orient eigenvector provided as reference"
            << std::endl;
    }

    return false;
}


struct PairSort
{
    bool operator()( const std::pair< double, unsigned int >& p1,
        const std::pair< double, unsigned int >& p2 ) const
    {
        return ( p1.first < p2.first );
    }
};

void crease::connect_segments()
{
    // compute point <- point -> point connections
    std::vector< std::pair< int, int > > connected_to( all_face_points.size(),
        std::pair< int, int >( -1, -1 ) );
/*    std::cout << "all edges: " << std::endl;
*/    for ( unsigned int i=0 ; i<all_edges.size() ; i++ )
    {
        // face indices for current edge
        unsigned int f0, f1;
        f0 = all_edges[i].first;
        f1 = all_edges[i].second;
        /*
                std::cout << f0 << "-" << f1 << std::endl;
                assert( f0<all_face_points.size() && f1<all_face_points.size() );*/

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
    /*
        std::cout << "all connected to: " << std::endl;
        for ( unsigned int i=0 ; i<all_face_points.size() ; i++ )
        {
            std::cout << connected_to[i].first << " <- " << i << " -> " << connected_to[i].second << std::endl;
        }*/

    components.clear();
    std::vector< bool > inserted( all_face_points.size(), false );
    for ( unsigned int i=0 ; i<all_edges.size() ; i++ )
    {
        // edge end points
        unsigned int i0, i1;
        i0 = all_edges[i].first;
        i1 = all_edges[i].second;

        assert( i0<all_face_points.size() && i1<all_face_points.size() );

        if ( inserted[i0] || inserted[i1] ) continue;

        unsigned int cur, prev;
        int link0, link1;

        // start a new connected component
        std::list< unsigned int > my_list;
        if ( verbose_mode )
            std::cout << "starting connected component #" << components.size() << std::endl;

        // initialize connected component with these two points
        my_list.push_back( i0 );
        my_list.push_back( i1 );

        // append forward
        prev = i0;
        cur = i1;

        if ( verbose_mode )
            std::cout << prev << "-" << cur << "-" << std::flush;

        inserted[prev] = true;
        for ( ; ; )
        {
            inserted[cur] = true;
            link0 = connected_to[cur].first; // always >= 0
            link1 = connected_to[cur].second;
            assert( link0>=0 );
            if ( ( unsigned int )link0 != prev && !inserted[link0] )
                my_list.push_back( link0 );
            else if ( link1 >= 0 && ( unsigned int )link1 != prev && !inserted[link1] )
                my_list.push_back( link1 );
            else break;

            if ( verbose_mode )
                std::cout << my_list.back() << "-" << std::flush;

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
            assert( link1<0 || link1<all_face_points.size() );
            if ( ( unsigned int )link0 != prev && !inserted[link0] )
                my_list.push_front( link0 );
            else if ( link1 >= 0 && ( unsigned int )link1 != prev && !inserted[link1] )
                my_list.push_front( link1 );
            else break;

            if ( verbose_mode )
                std::cout << my_list.front() << "-" << std::flush;

            prev = cur;
            cur = my_list.front();
        }
        std::cout << std::endl;

        components.push_back( std::vector< unsigned int >( my_list.begin(), my_list.end() ) );
    }

    std::cout << "verifying results:" << std::endl;
    for ( unsigned int i=0 ; i<components.size() ; i++ )
    {
        std::cout << "component #" << i << ": (" << components[i].size() << ")" << std::flush;
        for ( unsigned int j=0 ; j<components[i].size() ; j++ )
        {
            std::cout << components[i][j] << " " << std::flush;
        }
        std::cout << std::endl;
    }
}

void crease::extract_lines( const Nrrd* nrrd, bool ridge, unsigned int aniso )
{
    _ridge = ridge;

    std::cout << "subdiv = " << subdiv << std::endl;
    std::cout << "eps = " << eps << std::endl;
    std::cout << "upsample = " << upsample << std::endl;

    // check if this is a tensor field
    bool is_tensor = ( nrrd->dim == 4 && nrrd->axis[0].size == 7 );

    gage_interface::grid ref_grid( nrrd );
    gage_interface::grid sample_grid( nrrd, upsample );
    // use FA for now in the tensor case
    the_wrapper = new MeasureWrapper( nrrd, aniso );

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

    // examine the value distribution of the considered scalar quantity
    unsigned int M, N, P;
    M = ref_grid.size[0];
    N = ref_grid.size[1];
    P = ref_grid.size[2];
    std::vector< double > vals;
    for ( unsigned int i=0 ; i<M ; i++ )
    {
        for ( unsigned int j=0 ; j<N ; j++ )
        {
            for ( unsigned int k=0 ; k<P ; k++ )
            {
                if ( the_wrapper->confidence( ref_grid(i,j,k) ) < 0.5 )
                    continue;
                double val = the_wrapper->value( ref_grid(i,j,k) );
                vals.push_back( val );
            }
        }
    }
    display_stats( vals, "measured scalar" );

    float answer;
    std::cout << "Value threshold: " << std::flush;
    std::cin >> answer;
    double thresh_val = answer;

    // determine ridge/valley strength threshold
    std::vector< double > means;
    for ( unsigned int i=0 ; i<M ; i++ )
    {
        for ( unsigned int j=0 ; j<N ; j++ )
        {
            for ( unsigned int k=0 ; k<P ; k++ )
            {
                if ( the_wrapper->confidence( ref_grid(i,j,k) ) < 0.5 )
                    continue;
                double val = the_wrapper->value( ref_grid(i,j,k) );
                if ( !val_ok( val, thresh_val, ridge ) )
                    continue;
                double eval = the_wrapper->eigenvalue( ref_grid(i,j,k), 1 );
                if ( eval_ok( eval, 0, ridge ) )
                {
                    means.push_back( eval );
                }
            }
        }
    }
    std::cout << "there are " << means.size() << " candidate vertices"
        << std::endl;
    display_stats( means, "crease strength" );

    std::cout << "Crease threshold: " << std::flush;
    std::cin >> answer;

    double thresh_strength = answer;

    std::this_thread::sleep_for(5000ms);
    std::cout << std::endl;

    M = sample_grid.size[0];
    N = sample_grid.size[1];
    P = sample_grid.size[2];

    unsigned int nb_wrong_val = 0;
    unsigned int nb_wrong_eval = 0;
    unsigned int nb_wrong_both = 0;
    for ( unsigned int k=0 ; k<P-1 ; k++ )
    {
        std::cout << std::endl
            << " ** k = " << k
            << ", ( " << 100.*( double )k/( double )( P-1 )
            << "% ), " << nb_evaluations << " eval., "
            << nb_subdivided << " subdiv., "
            << nbok << " inspected cells ("
            << 100.*( double )nbok/( double )( k*( M-1 )*( N-1 ) ) << "%), "
            << nb_segments << " segments ("
            << 100.*( double )nb_segments/( double )( k*( M-1 )*( N-1 ) )
            << "%) ** "
            << std::endl;
        for ( unsigned int j=0 ; j<N-1 ; j++ )
        {
            for ( unsigned int i=0 ; i<M-1 ; i++ )
            {
                // voxel vertices
                vector< vec3 > v(8);
                v[0] = sample_grid( i  , j  , k   );
                v[1] = sample_grid( i+1, j  , k   );
                v[2] = sample_grid( i+1, j+1, k   );
                v[3] = sample_grid( i  , j+1, k   );
                v[4] = sample_grid( i  , j  , k+1 );
                v[5] = sample_grid( i+1, j  , k+1 );
                v[6] = sample_grid( i+1, j+1, k+1 );
                v[7] = sample_grid( i  , j+1, k+1 );

                bool skip = false;
                for ( unsigned int l=0 ; l<8 ; l++ )
                {
                    if ( the_wrapper->confidence( v[l] ) < 0.5 )
                    {
                        skip = true;
                        break;
                    }
                }
                if ( skip ) continue;

                vector< unsigned int > ids(8);
                ids[0] = sample_grid.id( i  , j  , k   );
                ids[1] = sample_grid.id( i+1, j  , k   );
                ids[2] = sample_grid.id( i+1, j+1, k   );
                ids[3] = sample_grid.id( i  , j+1, k   );
                ids[4] = sample_grid.id( i  , j  , k+1 );
                ids[5] = sample_grid.id( i+1, j  , k+1 );
                ids[6] = sample_grid.id( i+1, j+1, k+1 );
                ids[7] = sample_grid.id( i  , j+1, k+1 );

                // check eigenvalue at each vertex
                /*
                                bool candidate = true;
                                for ( unsigned int l=0 ; l<8 ; l++ )
                                {
                                    double eval = the_wrapper->eigenvalue( v[l], 1 );
                                    double val = the_wrapper->value( v[l] );
                                    bool ok_value = val_ok( val, thresh_val, ridge );
                                    bool ok_strength = eval_ok( eval, thresh_strength, ridge );
                                    if ( !ok_strength || !ok_value )
                                    {
                                        if ( !ok_value ) ++nb_wrong_val;
                                        if ( !ok_strength ) ++nb_wrong_eval;
                                        if ( !ok_value && !ok_strength ) ++nb_wrong_both;
                                        candidate = false;
                                        break;
                                    }
                                }*/

                bool candidate = false;
                for ( unsigned int l=0 ; l<8 ; l++ )
                {
                    double eval = the_wrapper->eigenvalue( v[l], 1 );
                    double val = the_wrapper->value( v[l] );
                    bool ok_value = val_ok( val, thresh_val, ridge );
                    bool ok_strength = eval_ok( eval, thresh_strength, ridge );
                    if ( ok_strength && ok_value )
                    {
                        candidate = true;
                        break;
                    }
                }
                if ( !candidate ) continue;

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

                    vector< vec3 > evecs(4);

                    /*
                    if ( false && orient_face( evecs, p, 1 ) )
                    {
                        ++nb_orientable;
                        vec3 xing;
                        if ( extract_point_with_valid_emedium( xing, evecs, p,
                            ( ridge ? 2 : 0 ) ) )
                        {
                            all_face_points.push_back( xing );
                            face_found[fid] = all_face_points.size()-1;
                            face_point.push_back( all_face_points.size()-1 );
                        }
                    }
                    else
                    {*/
                        vec3 xing;
                    int nb = all_in_one_face( xing, p, ridge );
                    if ( nb > 0 )
                    {
                        ++nb_orientable;
                        all_face_points.push_back( xing );
                        face_found[fid] = all_face_points.size()-1;
                        face_point.push_back( all_face_points.size()-1 );
                    }
                    else if ( nb == -2 )
                    {
                        ++nb_fakebasis;
                        vec3 xing;
                            // here we need to resort to plan B which consists in
                            // creating an artificial but consistent (e2,e3) basis
                            // of Span(orth(e1)) over the whole face and use it instead
                            // of the eigenvectors computed by Gage.
                        if ( extract_point_without_valid_emedium( xing, p,
                            ( ridge ? 0 : 2 ) ) )
                        {
                            all_face_points.push_back( xing );
                            face_found[fid] = all_face_points.size()-1;
                            face_point.push_back( all_face_points.size()-1 );
                        }
                    }
                }

                if ( verbose_mode && face_point.size() == 1 )
                {
                    cout << "found " << face_point.size()
                        << " crease point(s) on voxel faces"
                        << endl;
                }
                else if ( face_point.size() >= 2 )
                {
                    if ( face_point.size()>2 )
                    {
                        continue;

                        // if connection is ambiguous, select best two vertices
                        unsigned int nzeros = face_point.size();
                        std::vector< std::pair< double, unsigned int > > vals( nzeros );
                        for ( unsigned int i=0 ; i<nzeros ; i++ )
                        {
                            double eval = the_wrapper->eigenvalue( all_face_points[face_point[i]], 1 );
                            vals[i].first = fabs( eval );
                            vals[i].second = i;
                        }

                        // select two highest values of crease strength
                        std::sort( vals.begin(), vals.end(), PairSort() );
                        unsigned int i0 = face_point[vals.back().second];
                        unsigned int i1 = face_point[vals[nzeros-2].second];
                        face_point[0] = i0;
                        face_point[1] = i1;

                        std::cout << "points at " << all_face_points[face_point[0]]
                            << " and " << all_face_points[face_point[1]]
                            << " were obtained as max of " << face_point.size() << " points"
                            << std::endl;
                    }

                    ++nb_segments;
                    all_edges.push_back( pair< unsigned int, unsigned int >
                        ( face_point[0], face_point[1] ) );
                }
            }
        }
    }

    // assess ridge strength along found ridge lines
    crease_strength.resize( all_face_points.size() );
    grad_dot_evec.resize( all_face_points.size() );
    measure_value.resize( all_face_points.size() );
    for ( unsigned int i=0 ; i<all_face_points.size() ; i++ )
    {
        double val;
        val = the_wrapper->eigenvalue( all_face_points[i], 1 );
        if ( ridge && val>=0 ) val=0;
        else if ( !ridge && val<=0 ) val=0;
        crease_strength[i] = val;
        val = the_wrapper->value( all_face_points[i] );
        measure_value[i] = val;

        nvis::vec3 grad, evec;
        grad = the_wrapper->gradient( all_face_points[i] );
        evec = the_wrapper->eigenvector( all_face_points[i], ( ridge ? 0 : 2 ) );
        grad *= 1/nvis::norm( grad );
        evec *= 1/nvis::norm( evec );
        grad_dot_evec[i] = fabs( nvis::inner( grad, evec ) );
    }

    // identify connected components
    connect_segments();

    unsigned int nb_cells = ( M-1 )*( N-1 )*( P-1 );
    cout << "fake eigenspace basis: " << nb_fakebasis
        << " / orientable eigenvectors: " << nb_orientable
        << " ("
        << ( nb_orientable ? ( double )nb_fakebasis/( double )nb_orientable*100 : 100 )
        << "%)" << endl;

    cout << "percentage of candidate voxels: "
        << ( double )nbok/( double )nb_cells
        << endl;
    cout << "percentage of voxels discarded because of value: "
        << ( double )nb_wrong_val/( double )nb_cells << std::endl
        << "percentage of voxels discarded because of strength: "
        << ( double )nb_wrong_eval/( double )nb_cells << std::endl
        << "percentage of voxels discarded because of both: "
        << ( double )nb_wrong_both/( double )nb_cells << std::endl;

    cout << "number of segments = " << all_edges.size()
        << ", number of connected components: " << components.size()
        << std::endl;



    /*
    std::cout << "filtering connecting components" << std::endl;
    std::vector< double > strengths;
    std::vector< double > lengths;
    std::vector< unsigned int > ids;
    for ( unsigned int i=0 ; i<components.size() ; i++ )
    {
        // skip trivial stuff
        if ( components[i].size() < 3 ) continue;

        double mean = 0;
        for ( unsigned int j=0 ; j<components[i].size() ; j++ )
        {
            double eval = the_wrapper->eigenvalue( all_face_points[ components[i][j] ], 1 );
            mean += eval;
        }

        if ( components[i].size() )
        {
            mean /= ( double )components[i].size();
            strengths.push_back( fabs( mean ) );
            lengths.push_back( ( double )components[i].size() );
            ids.push_back( i );
        }
        std::cout << "connected component #" << i << ": length=" << components[i].size()
            << ", average strength=" << fabs( mean ) << std::endl;
    }
    display_stats( lengths, "ridge line lengths: " );
    display_stats( strengths, "ridge line strengths: " );

    std::cout << "ridge strength threshold? " << std::flush;
    std::cin >> answer;
    std::vector< std::vector< unsigned int > > filtered;
    for ( unsigned int i=0 ; i<ids.size() ; i++ )
    {
        if ( components[ids[i]].size()>2 && strengths[i]>answer )
            filtered.push_back( components[ids[i]] );
    }
    components.swap( filtered );
*/

}
