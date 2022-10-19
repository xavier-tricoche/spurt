#ifndef __NEW_CREASE_LINE_HPP__
#define __NEW_CREASE_LINE_HPP__

#include <xavier/image/probe.hpp>
#include <nvis/math/fixed_vector.hpp>
#include <teem/nrrd.h>
#include <set>
#include <map>
#include <vector>
#include <list>
#include <iostream>
#include <sstream>

namespace xavier
{
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
        extern vector< list< unsigned int > > components;

        extern double threshold;
        extern double eps;
        extern unsigned int subdiv;

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

        struct Grid
        {
            Grid( unsigned int size[3], double min[3], double max[3] );

            unsigned int id( unsigned int i, unsigned int j, unsigned int k );
            vec3 operator()( unsigned int i, unsigned int j, unsigned int k );

            double _min[3], _max[3];
            unsigned int _size[3];
            double _d[3];
        };

        class MeasureWrapper
        {
        public:
            MeasureWrapper( const scalar_wrapper& scal );
            MeasureWrapper( const tensor_wrapper& tens, int aniso );
            
            Vector eigenvector( const Point& p, unsigned int idx ) const;
            Vector gradient( const Point& p ) const;
            double eigenvalue( const Point& p, unsigned int idx ) const;
            
        private:
            scalar_wrapper _scal_wrap;
            tensor_wrapper _tens_wrap;
            std::vector< Vector > _evecs;
            std::vector< double > _evals;
            Vector _grad;
            double _evecs_array[9];
            double _evals_array[3];
            double _grad_array[3];    
            int _aniso;
        };

        Vector eigenvector( const Point& p, const scalar_wrapper& gH, 
            unsigned int idx );

        Vector gradient( const Point& p, const scalar_wrapper& gH );

        Vector basis_vector( const Vector& e );

        double zero_crossing( const Vector& ev0, const Vector& ev1,
            const Vector& g0, const Vector& g1 );

        // track orientation of well-defined eigenvector along an edge
        // returned boolean value indicates whether tracking was successful. 
        bool track_eigenvector_orientation( Vector& out, const Vector& in,
            const Edge& edge, const scalar_wrapper& gH, unsigned int idx );

        // return eigenvector with same orientation as input eigenvector
        void next_evec( Vector& out, const Vector& in, 
            const Edge& edge, const scalar_wrapper& gH, unsigned int idx );

        bool orient_face( vector< Vector >& evecs, const vector< Point >& face,
            const scalar_wrapper& gH, unsigned int idx );

        bool extract_point_with_valid_emedium( Point& out, const vector< Vector >& evecs,
            const vector< Point >& face,
            const scalar_wrapper& gH, unsigned int idx );

        bool extract_point_without_valid_emedium( Point& out, const vector< Point >& face,
            const scalar_wrapper& gH, unsigned int idx );

        void connect_segments();

        void extract_lines( const Nrrd* nrrd, bool ridge );

    };
};
};


inline
    xavier::crease::FaceId::FaceId( unsigned int i0, unsigned int i1, 
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
    xavier::crease::FaceId::FaceId( const FaceId& fid )
    : is(4)
{
    for ( unsigned int i=0 ; i<4 ; i++ )
    {
        is[i] = fid.is[i];
    }
}

inline
    int xavier::crease::FaceId::operator<( const FaceId& fid ) const
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
    std::ostream& xavier::crease::operator<<( ostream& os, const xavier::crease::FaceId& fid )
{
    os << "[ " << fid.is[0] << ", " << fid.is[1] << ", " 
        << fid.is[2] << ", " << fid.is[3] << "]";
    return os;
}     

inline
    xavier::crease::Grid::Grid( unsigned int size[3], 
    double min[3], double max[3] )
{
    for ( unsigned int i=0 ; i<3 ; i++ )
    {
        if ( false )
        { 
            // physical coordinates
            _min[i] = min[i];
            _max[i] = max[i];
            _size[i] = size[i];
            _d[i] = ( _max[i]-_min[i] )/( double )( _size[i]-1 );
        }
        else
        {
            // raster coordinates
            _min[i] = 0;
            _max[i] = size[i]-1;
            _size[i] = size[i]; 
            _d[i] = 1;
        }
    }
}

inline
    unsigned int xavier::crease::Grid::id( unsigned int i, unsigned int j, unsigned int k )
{
    return i + _size[0]*( j+_size[1]*k );
}

inline
    nvis::vec3 xavier::crease::Grid::operator()( unsigned int i, unsigned int j, unsigned int k )
{
    return vec3( _min[0]+i*_d[0],
        _min[1]+j*_d[1],
        _min[2]+k*_d[2] );
}

#endif
