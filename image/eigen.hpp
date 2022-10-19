#ifndef __EIGEN_HPP__
#define __EIGEN_HPP__

#include <math/fixed_vector.hpp>
#include <complex>

namespace xavier
{
  void eigensystem( nvis::vec2& evec_min, nvis::vec2& evec_max,
            double& lmin, double& lmax,
            const nvis::vec3& Hessian );

  void eigensystem( nvis::vec2& evec_min, nvis::vec2& evec_max,
            std::complex< double >& lmin, 
            std::complex< double >& lmax,
            const nvis::vec4& Jacobian );

  nvis::vec2 singularsystem( const nvis::vec4& M );
};


inline
void xavier::eigensystem( nvis::vec2& evec_min, nvis::vec2& evec_max,
              double& lmin, double& lmax,
              const nvis::vec3& Hessian )
{
    static nvis::vec4 _matrix;
    
    const double& hxx = Hessian[0];
    const double& hxy = Hessian[1];
    const double& hyy = Hessian[2];
    
    double trace = hxx+hyy;
    double det = hxx*hyy - hxy*hxy;
    
    double delta = trace*trace - 4*det;
    lmin = 0.5*( trace - sqrt( delta ) );
    lmax = 0.5*( trace + sqrt( delta ) );
    
    _matrix[1] = _matrix[2] = hxy;
    _matrix[0] = hxx-lmin;
    _matrix[3] = hyy-lmin;
    evec_min = singularsystem( _matrix );
    
    _matrix[0] = hxx-lmax;
    _matrix[3] = hyy-lmax;
    evec_max = singularsystem( _matrix );
}


inline
void xavier::eigensystem( nvis::vec2& evec_min, nvis::vec2& evec_max,
              std::complex< double >& lmin, 
              std::complex< double >& lmax,
              const nvis::vec4& Jacobian )
{
    static double tr, det, delta;
    static nvis::vec4 _matrix;
    
    tr = Jacobian[0]+Jacobian[3];
    det = Jacobian[0]*Jacobian[3]-Jacobian[1]*Jacobian[2];
    delta = tr*tr-4*det;
    _matrix[1] = Jacobian[1];
    _matrix[2] = Jacobian[2];
    
    if ( delta >= 0 ) {
        lmin = std::complex< double >( 0.5*( tr-sqrt(delta) ), 0 );
        _matrix[0] = Jacobian[0]-lmin.real();
        _matrix[3] = Jacobian[3]-lmin.real();
        evec_min = singularsystem( _matrix );
        
        lmax = std::complex< double >( 0.5*( tr+sqrt(delta) ), 0 );
        _matrix[0] = Jacobian[0]-lmax.real();
        _matrix[3] = Jacobian[3]-lmax.real();
        evec_max = singularsystem( _matrix );
    } else {
        lmin = std::complex< double >( 0.5*tr, -0.5*sqrt(-delta) );
        lmax = std::complex< double >( 0.5*tr, 0.5*sqrt(-delta) );
    }
}


inline
nvis::vec2 xavier::singularsystem( const nvis::vec4& M )
{
    static nvis::vec2 v;
    
    double a = M[0]*M[0] + M[1]*M[1];
    double b = M[2]*M[2] + M[3]*M[3];
    
    if ( a>b ) {
        v[0] = -M[1];
        v[1] = M[0];
    } else {
        v[0] = -M[3];
        v[1] = M[2];
    }
    
    double norm = nvis::norm( v );
    if ( norm>0 ) {
        v /= norm;
    }
    
    return v;
}


#endif
