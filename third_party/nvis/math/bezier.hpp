#ifndef __bezier_hpp
#define __bezier_hpp

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <math/array.hpp>

namespace nvis
{
    // --- evaluation -----------------------------------------

    template<typename T, size_type N> 
    inline T eval_bezier( const fixed_vector<double,N>& s,
			  array_ref<T,N> pts )
    {
	const size_type Np = pts.size();
	
	T tmp[Np];
	
	for( size_type i=0; i<Np; ++i )
	    tmp[i] = eval_bezier( subv<0,N-1>(s), pts[i] );
	
	return eval_bezier( s[N-1], array_ref<T,1>( tmp, Np ) );
    };
    
    template<typename T>
    inline T eval_bezier( const fixed_vector<double,1>& s, 
			  array_ref<T,1> pts )
    {
	const size_type Np = pts.size();
	
	if( Np == 1 )
	    return pts[0];
	
	T tmp[Np];
	
	for( size_type i=Np-1; i>=1; --i )
	    tmp[i] = pts[i-1] * (1-s[0]) + pts[i]*s[0]; 
	
	for( size_type l=2; l<Np; ++l )
	    for( size_type i=Np-1; i>=l; --i )
		tmp[i] = tmp[i-1] * (1-s[0]) + tmp[i] * s[0];
	
	return tmp[Np-1];
    }
    
    // --- derivatives ----------------------------------------
    
    template<typename T, size_type N>
    inline T derive_bezier( const fixed_vector<double,N>& s,
			    array_ref<T,N> pts )
    {
	const size_type Np = pts.size();
	
	T tmp[Np];
	
	for( size_type i=0; i<Np; ++i )
	    tmp[i] = derive_bezier( subv<0,N-1>(s), pts[i] );
	
	return eval_bezier( s[N-1], array_ref<T,1>( tmp, Np ) );
    }
    
    template<typename T>
    inline T derive_bezier( const fixed_vector<double,1>& s, 
			    array_ref<T,1> pts )
    {
	const size_type Np = pts.size();
	
	if( Np == 1 )
	    return T( 0 );
	
	T tmp[Np-1];
	
	for( unsigned int i=0; i<Np-1; ++i )
	    tmp[i] = (pts[i+1] - pts[i]); 
	
	return (Np-1)*eval_bezier( s[0], array_ref<T,1>( tmp, Np-1 ) );
    }
    
    template<typename T, size_type N>
    inline T derive_bezier( const size_type& dim,
			    fixed_vector<double,N> s,
			    array_ref<T,N> pts )
    {
	fixed_vector<double,N> ss(s);
	
	std::swap( ss[0], ss[dim] );
	return derive_bezier( ss, pts.permute( N-1-dim, N-1 ) );
	
    }

    // --- splitting ----------------------------------------
    
    template<typename T>
    inline void split_bezier( const double& s,
			      array_ref<T,1> pts, 
			      array_ref<T,1> left, 
			      array_ref<T,1> right )
    {
	const size_type Np = pts.size();
	
	for( size_type i=0; i<Np; ++i )
	    left[i] = pts[i]; 
	
	right[Np-1] = pts[Np-1];
	
	for( size_type l=1; l<Np; ++l )
	{
	    for( size_type i=Np-1; i>=l; --i )
		left[i] = left[i-1] * (1-s) + left[i] * s;
	    
	    right[Np-l-1] = left[Np-1];
	}
    }
    
    
    template<typename T, size_type N>
    inline void split_bezier( const double& s,
			      array_ref<T,N> pts,
			      array_ref<T,N> left,
			      array_ref<T,N> right )
    {
	for( size_type i=0; i<pts.size(); ++i )
	    split_bezier( s, pts[i], left[i], right[i] );
    }
    
    template<typename T, size_type N>
    void split_bezier( size_type dim,
		       const double& s,
		       array_ref<T,N> pts,
		       array_ref<T,N> left,
		       array_ref<T,N> right )
    {
	dim = N-1 - dim;
	
	split_bezier( s, 
		      pts.permute( dim, N-1 ), 
		      left.permute( dim, N-1 ), 
		      right.permute( dim, N-1 ) );
    }
    
    // --- inversion ----------------------------------------------------

    template<typename T, size_type N>
    bool bezier_invert_subdiv( const T& v,
			       array_ref<T,N> opts,
			       array_ref<T,N> pts, 
			       fixed_vector<double,N>& p0,
			       fixed_vector<double,N>& p1,
			       const double& tolerance,
			       size_type dim = 0 )
    {
	array<T,N> split1( pts.shape() );
	array<T,N> split2( pts.shape() );
	
	bounding_box<T> bbox1, bbox2;
	
	split_bezier( dim, 0.5, pts, split1, split2 );

	while( true )
	{
	    bbox1.set( split1.begin(), split1.end() );
	    bbox2.set( split2.begin(), split2.end() );
	    
	    bool in1 = bbox1.inside(v);
	    bool in2 = bbox2.inside(v);
	    
	    if( in1 && !in2 )
	    {
		p1[dim] = 0.5*( p0[dim]+p1[dim] );
	    }
	    else if( in2 && !in1 )
	    {
		p0[dim] = 0.5*(p0[dim]+p1[dim]);
		split1.swap( split2 );
	    }
	    else if( in1 && in2 )
	    {
		fixed_vector<double,N> p0save = p0, p1save = p1;
		
		p1[dim] = 0.5*(p0[dim]+p1[dim]);
		
		if( bezier_invert_subdiv( v, opts, split1, p0, p1, tolerance, (dim+1)%N ) )
		    return true;
		
		p0 = p0save;
		p1 = p1save;
		
		p0[dim] = 0.5*(p0[dim]+p1[dim]);
		
		return bezier_invert_subdiv( v, opts, split2, p0, p1, tolerance, (dim+1)%N );
	    }
	    else
		return false;
	    
	    T f = eval_bezier(0.5*p0+0.5*p1, opts);
	    
	    if( norm(f-v) < tolerance )
		return true;

	    dim = (dim+1) % N;
	    split_bezier( dim, 0.5, split1, split1, split2 );
	    
	}
	
	return true;
    }


    namespace { 

	// TODO: clean this up, roll out a fixed_matrix class and implement invert()
	// until then, this must remain here

	template<typename T> 
	fixed_vector<T,2> inverse_mult( const fixed_vector<T,2> dd[2],
					const fixed_vector<T,2>& v )
	{
	    double norm = dd[0][0] * dd[1][1] - dd[0][1]*dd[1][0];

	    fixed_vector<T,2> res;
	    
	    res[0] = (dd[1][1]*v[0] - dd[1][0]*v[1]);
	    res[1] = (dd[0][0]*v[1] - dd[0][1]*v[0]);

	    return res/norm;
	}

	template<typename T> 
	fixed_vector<T,3> inverse_mult( const fixed_vector<T,3>& ds, 
					const fixed_vector<T,3>& dt, 
					const fixed_vector<T,3>& du, 
					const fixed_vector<T,3>& v )
	{
	    fixed_vector<T,3> s, res;
	    
	    s[0] = (dt[1]*du[2] - dt[2]*du[1]);
	    s[1] = (dt[0]*du[2] - dt[2]*du[0]);
	    s[2] = (dt[0]*du[1] - dt[1]*du[0]);
	    
	    double det = ds[0]*s[0] + ds[1]*s[1] + ds[2]*s[2];
	    
	    res[0] = (s[0]*v[0] + s[1]*v[1] + s[2]*v[2]) / det;
	    
	    s[0] = (ds[1]*du[2] - ds[2]*du[1]);
	    s[1] = (ds[0]*du[2] - ds[2]*du[0]);
	    s[2] = (ds[0]*du[1] - ds[1]*du[0]);
	    
	    res[1] = (s[0]*v[0] + s[1]*v[1] + s[2]*v[2]) / det;
	    
	    s[0] = (ds[1]*dt[2] - ds[2]*dt[1]);
	    s[1] = (ds[0]*dt[2] - ds[2]*dt[0]);
	    s[2] = (ds[0]*dt[1] - ds[1]*dt[0]);
	    
	    res[2] = (s[0]*v[0] + s[1]*v[1] + s[2]*v[2]) / det;
	    
	    return res;
	}
	
    } // namespace

    template<typename T, size_type N> 
    bool bezier_invert_newton( const T& v, 
			       array_ref<T,N> pts, 
			       fixed_vector<double,N>& p,
			       const double& tolerance )
    {
	T f2, f = eval_bezier( p, pts );
	fixed_vector<double,N> p2;
	
	size_type niter;
	const size_type maxiter = 20;
	
	T deriv[N];

	for( niter=0; niter<maxiter; ++niter )
	{
	    for( size_type i=0; i<N; ++i )
		deriv[i] = derive_bezier( i, p, pts );

	    fixed_vector<double,N> h = inverse_mult( deriv, v-f );
	    
// 	    double lambda = 1.0;
	    
// 	    do
// 	    {
// 		p2 = p+lambda*h;
// 		f2 = eval_bezier( p2, pts );
		
// 		lambda *= 0.5;

// 		std::cout << "lambda iter: " << (p2<=fixed_vector<double,N>(0.)) << '\t' << (p2>=fixed_vector<double,N>(1.)) << '\n';
// 	    }
// 	    while( ( norm(v-f2) > norm(v-f) ) 
// 		   && lambda > 0.01 );
	    
// 	    std::cout << "newton: " << p << '\t' << h << '\n';

	    p = p + h;

	    p = max( p, fixed_vector<double,N>( 0. ) );
	    p = min( p, fixed_vector<double,N>( 1. ) );

	    f = eval_bezier( p, pts );
	    
	    if( norm(v-f) < tolerance )
		break;
	}
	
	if( niter == maxiter )
	    return false;
	
	assert( all(p>=T(0.)) && all(p<=T(1.)) );

	return true;
    }

    template<typename T, size_type N>
    bool bezier_invert( const T& v, array_ref<T,N> pts,
			fixed_vector<double,N>& param,
			const double& tolerance = 1e-9 )
    {
 	const double init_tol = 0.125;
	
	fixed_vector<double,N> p0( 0.0 ), p1( 1.0 );
	
	bool found = bezier_invert_subdiv( v, pts, pts, p0, p1, init_tol );
	
	if( found )
	{
// 	    std::cout << "found by subdiv\n";

	    param = 0.5*p0 + 0.5*p1;

	    found = bezier_invert_newton( v, pts, param, tolerance );
	    
	    if( !found )
	    {
//  		std::cout << "need deep subdiv\n";

		p0 = fixed_vector<double,2>( 0.0 );
		p1 = fixed_vector<double,2>( 1.0 );
		
		found = bezier_invert_subdiv( v, pts, pts, p0, p1, tolerance );
		param = 0.5*p0 + 0.5*p1;

//  		std::cout << "deep subdiv found: " << found << '\n';
	    }
//  	    else
//  		std::cout << "found by newton\n";


	}
	return found;
    }


} // namespace nvis

#endif //__bezier_hpp


