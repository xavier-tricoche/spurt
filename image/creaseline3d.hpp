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

        typedef pair< vec3, vec3 > Edge;
        typedef pair< unsigned int, unsigned int > EdgeId;
        typedef std::pair< nvis::vec3, double > value;
        
        extern vector< Edge > problematic_edges;
        extern vector< vec3 > all_face_points;
        extern vector< EdgeId > all_edges;
        extern vector< vector< unsigned int > > components;
        extern vector< double > crease_strength;
        extern vector< double > grad_dot_evec;
        extern vector< double > measure_value;
        
        // for debugging purposes
        extern vector< vector< vector< value > > > nonorientable_faces;
        extern vector< vector< nvis::vec3 > > nonorientable_face_points;

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

            double confidence( const vec3& p ) const;
            double value( const vec3& p ) const;
            vec3 eigenvector( const vec3& p, unsigned int idx ) const;
            vec3 gradient( const vec3& p ) const;
            double eigenvalue( const vec3& p, unsigned int idx ) const;
            vec6 hessian( const vec3& p ) const;
            
            int what_measure() const;

        private:
            scalar_wrapper *_scal_wrap;
            tensor_wrapper *_tens_wrap;
            mutable std::vector< vec3 > _evecs;
            mutable std::vector< double > _evals;
            mutable vec3 _grad;
            int _measure;
        };

        vec3 basis_vector( const vec3& e );
        
        vec3 artificial_evec_med( double u, const vec3& ev0, 
            const vec3& ev1, const Edge& edge, unsigned int idx );
            
        vec3 artificial_evec_min( double u, 
            const vec3& ev0, const vec3& ev1, 
            const vec3& ev0_ref, const vec3& ev1_ref, 
            const Edge& edge, unsigned int idx );

        double zero_crossing( const vec3& ev0, const vec3& ev1,
            const vec3& g0, const vec3& g1 );
            
        double zero_crossing_evec_grad( const Edge& edge, unsigned int idx );
        
        double zero_crossing_emed( const Edge& edge ); // may need additional parameters
        
        double zero_crossing_emin( const Edge& edge ); // may need additional parameters 
        
        

        // track orientation of well-defined eigenvector along an edge
        // returned boolean value indicates whether tracking was successful. 
        bool track_eigenvector_orientation( vec3& out, const vec3& in,
            const Edge& edge, unsigned int idx );

        // return eigenvector with same orientation as input eigenvector
        void next_evec( vec3& out, const vec3& in, 
            const Edge& edge, unsigned int idx );

        bool orient_face( vector< vec3 >& evecs, const vector< vec3 >& face,
            unsigned int idx );

        bool extract_point_with_valid_emedium( vec3& out, const vector< vec3 >& evecs,
            const vector< vec3 >& face, unsigned int idx );

        bool extract_point_without_valid_emedium( vec3& out, const vector< vec3 >& face,
            unsigned int idx );

        void connect_segments();
        
        void extract_lines( const Nrrd* nrrd, bool ridge, unsigned int aniso=0 );
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
