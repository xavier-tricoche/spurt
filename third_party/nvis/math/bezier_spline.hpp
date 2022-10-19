#ifndef __bezier_spline_hpp
#define __bezier_spline_hpp

#include <vector>
#include <iosfwd>
#include <cmath>
#include <stdexcept>
#include <iostream>

namespace {
#if 1
template<unsigned int D, typename Iter>
typename Iter::value_type decasteljau( double t, Iter pi )
{
    if( D == 0 )
	return *pi;

    typename Iter::value_type tmp[D+1];

    std::copy( pi, pi+D+1, tmp );

    for( size_t l=1; l<D+1; ++l )
	for( size_t i=D; i>=l; --i )
	    tmp[i] = tmp[i-1] * (1-t) + tmp[i] * t;

    return tmp[D];
}

template<unsigned int D, typename Iter>
typename Iter::value_type decasteljau_derivative( double t, Iter pi )
{
    if( D == 0 )
	    return *pi;

    std::vector<typename Iter::value_type> tmp( D );

    for( size_t i=0; i<D; ++i )
        tmp[i] = pi[i+1] - pi[i];

    return decasteljau<D-1>( t, tmp.begin() );
}

#else
template<unsigned int D, typename Iter>
typename Iter::value_type decasteljau( double t, Iter pi )
{
    if( D == 0 )
	    return *pi;

    typename Iter::value_type tmp[D];

    for( size_t i=0; i<D; ++i )
        tmp[i] = *pi * (1-t) + *(++pi) * t;

    for( size_t l=D-1; l>0; --l )
        for( size_t i=0; i<l; ++i )
            tmp[i] = tmp[i]*(1-t) + tmp[i+1]*t;

    return tmp[0];
}

template<unsigned int D, typename Iter>
typename Iter::value_type decasteljau_derivative( double t, Iter pi )
{
    if( D == 0 )
	return *pi;

    typename Iter::value_type tmp[D];

    for( size_t i=0; i<D; ++i )
        tmp[i] = *pi * (1-t) + *(++pi) * t;

    for( size_t l=D-2; l>0; --l )
        for( size_t i=0; i<l; ++i )
            tmp[i] = tmp[i]*(1-t) + tmp[i+1]*t;

    return tmp[1]-tmp[0];
}
#endif
}

// --------------------------------------------------------------------------

namespace nvis {

template<typename T, unsigned int D>
class bezier_spline
{
    typedef std::vector<double> knot_vector;
    typedef std::vector<T>      ctrl_vector;

public:

    bezier_spline()
    {
    }

    bezier_spline( knot_vector& knots,
		   ctrl_vector& ctrlp )
    {
        _knots.swap( knots );
        _ctrlp.swap( ctrlp );

        consistency_check();
    }

    template<typename KIter, typename CPIter>
    bezier_spline( KIter  kbegin,  KIter  kend,
                   CPIter cpbegin, CPIter cpend ) :
        _knots( kbegin, kend ), _ctrlp( cpbegin, cpend )
    {
        consistency_check();
    }

    bezier_spline( std::istream& in )
    {
        read( in );

        if( in )
            consistency_check();
    }

    // ---

    const double& t_min() const
    {
    	return _knots.front();
    }

    const double& t_max() const
    {
    	return _knots.back();
    }

    // ---

    const knot_vector& knots() const
    {
    	return _knots;
    }

    const ctrl_vector& ctrlp() const
    {
    	return _ctrlp;
    }

    // ---

    void reverse()
    {
    	std::reverse( _knots.begin(), _knots.end() );
    	std::reverse( _ctrlp.begin(), _ctrlp.end() );
    }

    // ---

    T operator()( double t ) const
    {
    	if(_knots.size() == 1)
            return _ctrlp[0];

        size_t s = interval( t );

    	t -= _knots[s];
    	t /= _knots[s+1] - _knots[s];

        return decasteljau<D>( t, _ctrlp.begin()+D*s );
    }

    T derivative( double t ) const
    {

        if(_knots.size() == 1)
        {
            //std::cout << "SIZE --> X " << _knots.size() << " " << t << std::endl;
            return (_ctrlp[0]);
        }

        size_t s = interval( t );

        t -= _knots[s];
        t /= _knots[s+1] - _knots[s];

        return decasteljau_derivative<D>( t, _ctrlp.begin()+D*s ) / (_knots[s+1]-_knots[s]);
    }

    // ---

    void write( std::ostream& out ) const
    {
    	unsigned int tmp;

    	tmp = D;
    	out.write( (const char*)&tmp, sizeof(tmp) );

    	tmp = _knots.size();
    	out.write( (const char*)&tmp, sizeof(tmp) );
    	out.write( (const char*)&_knots.front(), tmp * sizeof(double) );

    	tmp = 4*tmp - 3;
    	out.write( (const char*)&_ctrlp.front(), tmp * sizeof(T) );
    }

    void read( std::istream& in )
    {
    	unsigned int tmp;

    	in.read( (char*)&tmp, sizeof(tmp) );

    	if( !in )
    	    return;

    	if( D != tmp )
    	    throw std::runtime_error( "degree mismatch while reading bezier_spline" );

    	in.read( (char*)&tmp, sizeof(tmp) );

    	if( !in )
    	    return;

    	_knots.resize( tmp );
    	in.read( (char*)&_knots.front(), tmp * sizeof(double) );

    	if( !in )
    	    return;

    	tmp = 4*tmp - 3;

    	_ctrlp.resize( tmp );
    	in.read( (char*)&_ctrlp.front(), tmp * sizeof(T) );
    }

    void consistency_check() const
    {
    	if( _knots.size()-1 != (_ctrlp.size()-1) / D )
    	    throw std::logic_error( "knots/control points inconsistent" );

    	if( _knots.size() < 2 )
    	    throw std::logic_error( "less than two knots" );
    }

    size_t interval( double t ) const
    {
    	knot_vector::const_iterator ki =
    	    std::lower_bound( _knots.begin(), _knots.end(), t );

    	if( ki == _knots.begin() )
    	    return 0;
    	else if( ki == _knots.end() )
    	    return _knots.size()-2;
    	else
    	    return ki - _knots.begin() - 1;
    }

    std::vector<double>  _knots;
    std::vector<T>       _ctrlp;
};

// --------------------------------------------------------------------------

template<typename KIter, typename CPIter>
bezier_spline<typename KIter::value_type,1>
make_linear_spline( KIter  kbegin, KIter kend, CPIter cpbegin, CPIter cpend )
{
	return bezier_spline<typename KIter::value_type,1>( kbegin, kend, cpbegin, cpend );
}

// --------------------------------------------------------------------------

template<typename T>
bezier_spline<T,1> make_linear_spline( const std::vector<double>& knots,
				                       const std::vector<T>& ctrlp )
{
	return bezier_spline<T,1>( knots.begin(), knots.end(), ctrlp.begin(), ctrlp.end() );
}

// --------------------------------------------------------------------------

template<typename T>
bezier_spline<T,3> make_catmull_rom_spline( const std::vector<double>& knots,
                                            const std::vector<T>& ctrlp )
{
	assert( knots.size() == ctrlp.size() );
	const int N = knots.size()-1;

	std::vector<T> m( N+1 );

	m[0] = (ctrlp[1]-ctrlp[0])/(knots[1]-knots[0]);
	m[N] = (ctrlp[N]-ctrlp[N-1])/(knots[N]-knots[N-1]);

	for( int i=1; i<N; ++i )
		m[i] = ( (ctrlp[i+1]-ctrlp[i])/(knots[i+1]-knots[i]) +
				 (ctrlp[i]-ctrlp[i-1])/(knots[i]-knots[i-1]) )/2;

	std::vector<T> crctrl( 3*N+1 );

	for( int i=0; i<N; i++ )
	{
		crctrl[3*i+0] = ctrlp[i];
		crctrl[3*i+1] = ctrlp[i] + m[i]/3*(knots[i+1]-knots[i]);
		crctrl[3*i+2] = ctrlp[i+1] - m[i+1]/3*(knots[i+1]-knots[i]);
	}

	crctrl.back() = ctrlp.back();

	return bezier_spline<T,3>( knots.begin(), knots.end(), crctrl.begin(), crctrl.end() );
}

// --------------------------------------------------------------------------

template<typename T>
bezier_spline<T,3> make_catmull_rom_spline_arclen( const std::vector<T>& ctrlp )
{
	const int N = ctrlp.size()-1;

	std::vector<T> m( N+1 );
    std::vector<double> knots( ctrlp.size() );
    knots[0] = 0.0;

    for( int i=1; i<=N; ++i )
        knots[i] = knots[i-1] + norm( ctrlp[i] - ctrlp[i-1] );

    for( int i=0; i<=N; ++i )
        knots[i] /= knots[N];

	m[0] = (ctrlp[1]-ctrlp[0])/(knots[1]-knots[0]);
	m[N] = (ctrlp[N]-ctrlp[N-1])/(knots[N]-knots[N-1]);

	for( int i=1; i<N; ++i )
		m[i] = ( (ctrlp[i+1]-ctrlp[i])/(knots[i+1]-knots[i]) +
				 (ctrlp[i]-ctrlp[i-1])/(knots[i]-knots[i-1]) ) / 2;

	std::vector<T> crctrl( 3*N+1 );

	for( int i=0; i<N; i++ )
	{
		crctrl[3*i+0] = ctrlp[i];
		crctrl[3*i+1] = ctrlp[i] + m[i]/3*(knots[i+1]-knots[i]);
		crctrl[3*i+2] = ctrlp[i+1] - m[i+1]/3*(knots[i+1]-knots[i]);
	}

	crctrl.back() = ctrlp.back();

	return bezier_spline<T,3>( knots.begin(), knots.end(), crctrl.begin(), crctrl.end() );
}

// --------------------------------------------------------------------------

} // namespace nvis

template<typename T, unsigned int D>
std::ostream& operator<<( std::ostream& out, const nvis::bezier_spline<T,D>& bs )
{
    bs.write( out );
    return out;
}

#endif // __bezier_spline_hpp
