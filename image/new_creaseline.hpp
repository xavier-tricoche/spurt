#ifndef __NEW_CREASE_LINE_HPP__
#define __NEW_CREASE_LINE_HPP__

#include <spurt/image/probe.hpp>
#include <nvis/math/fixed_vector.hpp>
#include <teem/nrrd.h>
#include <set>
#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <sstream>
#include <string>

namespace spurt
{
    template< typename T >
    void display_stats( const std::vector< T >& data, const std::string& msg );
    
    namespace crease
    {
        using namespace std;
        using namespace nvis;
        using namespace gage_interface;

        typedef vec3 Vector;
        typedef vec3 Point;
        typedef pair< Point, Point > Edge;
        typedef pair< unsigned int, unsigned int > EdgeId;

        extern vector< Edge > problematic_edges;
        extern vector< Point > all_face_points;
        extern vector< EdgeId > all_edges;
        extern vector< vector< unsigned int > > components;
        extern vector< double > crease_strength;
        extern vector< double > grad_dot_evec;

        extern double threshold;
        extern double eps;
        extern unsigned int subdiv;
        extern unsigned int upsample;

        // 0: scalar field
        // 1: FA of tensor field
        // 2: mode of tensor field
        extern int crease_kind;

        struct FaceId
        {
            FaceId( unsigned int i0, unsigned int i1, 
                unsigned int i2, unsigned int i3 );
            FaceId( const FaceId& fid );

            int operator<( const FaceId& fid ) const;

            vector< unsigned int > is;
        };

        ostream& operator<< ( ostream& os, const FaceId& fid );
        
        class MeasureWrapper
        {
        public:
            MeasureWrapper( const Nrrd* nrrd, int aniso=0 );

            double confidence( const Point& p ) const;
            double value( const Point& p ) const;
            Vector eigenvector( const Point& p, unsigned int idx ) const;
            Vector gradient( const Point& p ) const;
            double eigenvalue( const Point& p, unsigned int idx ) const;

        private:
            scalar_wrapper *_scal_wrap;
            tensor_wrapper *_tens_wrap;
            mutable std::vector< Vector > _evecs;
            mutable std::vector< double > _evals;
            mutable Vector _grad;
            int _measure;
        };

        Vector basis_vector( const Vector& e );
        
        Vector artificial_evec_med( double u, const Vector& ev0, 
            const Vector& ev1, const Edge& edge, unsigned int idx );
            
        Vector artificial_evec_min( double u, 
            const Vector& ev0, const Vector& ev1, 
            const Vector& ev0_ref, const Vector& ev1_ref, 
            const Edge& edge, unsigned int idx );

        double zero_crossing( const Vector& ev0, const Vector& ev1,
            const Vector& g0, const Vector& g1 );
            
        double zero_crossing_eref( const Edge& edge, unsigned int idx );
        
        double zero_crossing_emed( const Edge& edge ); // may need additional parameters
        
        double zero_crossing_emin( const Edge& edge ); // may need additional parameters 

        // track orientation of well-defined eigenvector along an edge
        // returned boolean value indicates whether tracking was successful. 
        bool track_eigenvector_orientation( Vector& out, const Vector& in,
            const Edge& edge, unsigned int idx );

        // return eigenvector with same orientation as input eigenvector
        void next_evec( Vector& out, const Vector& in, 
            const Edge& edge, unsigned int idx );

        bool orient_face( vector< Vector >& evecs, const vector< Point >& face,
            unsigned int idx );

        bool extract_point_with_valid_emedium( Point& out, const vector< Vector >& evecs,
            const vector< Point >& face, unsigned int idx );

        bool extract_point_without_valid_emedium( Point& out, const vector< Point >& face,
            unsigned int idx );

        void connect_segments();
        
        void extract_lines( const Nrrd* nrrd, bool ridge );
    };
};
};


inline
    spurt::crease::FaceId::FaceId( unsigned int i0, unsigned int i1, 
    unsigned int i2, unsigned int i3 )
    : is(4)
{
    is[0] = i0;
    is[1] = i1;
    is[2] = i2;
    is[3] = i3;
    sort( is.begin(), is.end() );
}

inline
    spurt::crease::FaceId::FaceId( const FaceId& fid )
    : is(4)
{
    for ( unsigned int i=0 ; i<4 ; i++ )
    {
        is[i] = fid.is[i];
    }
}

inline
    int spurt::crease::FaceId::operator<( const FaceId& fid ) const
{
    return ( is[0] < fid.is[0] ||
        ( is[0] == fid.is[0] && 
        ( is[1] < fid.is[1] ||
        ( is[1] == fid.is[1] && 
        ( is[2] < fid.is[2] ||
        ( is[2] == fid.is[2] &&
        is[3] < fid.is[3] ) ) ) ) ) );
}

inline
    std::ostream& spurt::crease::operator<<( ostream& os, const spurt::crease::FaceId& fid )
{
    os << "[ " << fid.is[0] << ", " << fid.is[1] << ", " 
        << fid.is[2] << ", " << fid.is[3] << "]";
    return os;
}

#endif
