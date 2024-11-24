#ifndef __array_hpp
#define __array_hpp

#include <nvis-math/fixed_vector.hpp>
#include <boost/iterator/iterator_facade.hpp>

namespace nvis
{
    template<typename T, size_type N> class array_ref;
    template<typename T, size_type N> class array;

    // ---------------------------------------------------------
    
    template<typename T, size_type N> class array_access_nd
    {
    public: 
	
	typedef T                            element;
	
	typedef fixed_vector<size_type,N>    size_tuple;
	typedef array<T,N-1>              value_type;
	typedef array_ref<T,N-1>          reference;
    
    protected: 

	reference access_outer( T* data, 
				const size_type&  index,
				const size_tuple& size, 
				const size_tuple& strides ) const
	{
	    return reference( data + index*strides[0],
			      subv<1,N-1>( size ), 
			      subv<1,N-1>( strides ) );
	}
	
	reference access_inner( T* data,
				const size_type& index,
				const size_tuple& size,
				const size_tuple& strides ) const
	{
	    return reference( data + index*strides[N-1],
			      subv<0,N-1>( size ), 
			      subv<0,N-1>( strides ) );
	}
    };

    template<typename T> class array_access_nd<T,1>
    {
    public: 
	
	typedef T                            element;
	typedef T                            value_type;
	typedef fixed_vector<size_type,1>    size_tuple;
	typedef T&                           reference;
	
	reference access_outer( T* data, 
				const size_type&  index,
				const size_tuple& size, 
				const size_tuple& strides ) const
	{
	    return *( data + index*strides[0] );
	}
	
	reference access_inner( T* data,
				const size_type& index,
				const size_tuple& size,
				const size_tuple& strides ) const
	{
	    return *( data + index*strides[0] );
	}
    };

    // ---------------------------------------------------------

    template<typename T, size_type N> class array_ref : 
	public array_access_nd<T,N>
    {
    public:
	
	typedef array_access_nd<T,N>             base_type;
	
	typedef typename base_type::element      element;
	typedef typename base_type::size_tuple   size_tuple;
	typedef typename base_type::value_type   value_type;
	typedef typename base_type::reference    reference;
	typedef T*                               pointer;
	
	// -----------------------------------------------------
	
	class iterator: 
	    public boost::iterator_facade< iterator, element,
					   boost::single_pass_traversal_tag >
	{
	public:
	    
	    iterator( const array_ref& array )
	    {
		    // array* runs fortran order, which is suboptimal for 
		// iterator traversal, so go C order by reversing shape and strides
		_strides = nvis::reverse(array._strides);
		_shape   = nvis::reverse(array._shape);
		_data    = array._base;
		
		for( int i=0; i<N; ++i )
		{
		    _stack[i] = _data;
		    _last[i]  = _data + _shape[i] * _strides[i];
		}
		
		_stride = _strides[0];
	    }
	    
	    iterator()
	    {
		_data = 0;
	    }
	    
	private:
	    
	    friend class boost::iterator_core_access;
	    
	    bool equal( const iterator& other ) const
	    {
		return _data == other._data;
	    }
	    
	    typename iterator::reference dereference() const
	    {
		return *_data;
	    }
	    
	    void increment()
	    {
		_data += _stride;
		
		if( _data != _last[0] )
		    return;
		
		int j;
		for( j=1; j<N; ++j )
		{
		    _data =  _stack[j];
		    _data += _strides[j];
		    
		    if( _data != _last[j] )
			break;
		}
		
		if( j==N )
		{
		    _data = 0;
		    return;
		}
		
		_stack[j] = _data;
		
		for( --j; j>=0; --j )
		{
		    _stack[j] = _data;
		    _last[j]  = _data + _shape[j]*_strides[j]; 
		}
	    }
	    
	    size_tuple _strides, _shape;
	    size_type  _stride;
	    T          *_data;
	    T          *_stack[N], *_last[N];
	};
    
	// ------------------------------------------------

	array_ref() : _base(0), _shape(0), _strides(0)
	{
	}

	array_ref( pointer base, 
		   const size_tuple& shape ) : 
	    _base(base), _shape(shape)
	{
	    compute_strides();
	}

	array_ref( const array_ref& r ) : 
	    _base(r._base), _shape(r._shape), _strides(r._strides)
	 {
	 }

	array_ref( pointer base,
		      const size_tuple& shape,
		      const size_tuple& strides ) : 
	    _base(base), _shape(shape), _strides(strides)
	{
	}

	reference operator[]( const size_type& i )
	{
	    return base_type::access_outer( _base, i, _shape, _strides );
	}
	
	element& operator()( const size_tuple& i )
	{
	    return _base[ inner(i,_strides) ];
	}
	
	const size_tuple& shape() const
	{
	    return _shape;
	}
	
	const size_type num_elements() const
	{
	    return prod(_shape);
	}
	
	array_ref permute( const size_type dim1,
			   const size_type dim2 ) const
	{
	    size_tuple strides = _strides, shape = _shape;

	    std::swap( shape[dim1],   shape[dim2]   );
	    std::swap( strides[dim1], strides[dim2] );

	    return array_ref( _base, shape, strides );
	}
	
	array_ref sub( const size_tuple& base, 
		       const size_tuple& shape, 
		       const size_tuple& strides = size_tuple( 1 ) )
	{
	    return array_ref( _base + inner( base, _strides ),
			      shape, strides * _strides );
	}

	array_ref permute( const size_tuple& dim ) const
	{
	    size_tuple strides, shape;

	    for( size_type i=0; i<N; ++i )
	    {
		strides[i] = _strides[dim[i]];
		shape[i]   = _shape[dim[i]]; 
	    }

	    return array_ref( _base, shape, strides );
	}

	array_ref shift() const
	{

	    return array_ref( _base, 
				 nvis::shift(_shape), 
				 nvis::shift(_strides) );
	}


	array_ref reverse() const
	{
	    return array_ref( _base, 
				 nvis::reverse(_shape), 
				 nvis::reverse(_strides) );
	}

	void swap( array_ref& other )
	{
	    std::swap( _strides, other._strides );
	    std::swap( _shape,   other._shape   );
	    std::swap( _base,    other._base    );
	}

	const size_type& size() const
	{
	    return _shape[0];
	}
	
	static size_type dimension()
        {
	    return N;
	}
	
	iterator begin() const
	{
	    return iterator( *this );
	}
	
	iterator end() const
	{
	    return iterator();
	}
	
	pointer base()
	{
	    return _base;
	}
	
	const size_tuple& strides() const
	{
	    return _strides;
	}

    protected:
	
	void compute_strides()
	{
	    _strides[N-1] = 1;
	    
	    for( size_type s=N-1; s>0; --s )
		_strides[s-1] = _strides[s]*_shape[s];
	}

    public:

	pointer     _base;
	size_tuple  _shape;
	size_tuple  _strides;
    };

    // ---------------------------------------------------------------
    
    template<typename T, size_type N> class array : 
	public array_ref<T,N>
    {
    public:

	typedef array_ref<T,N>                   base_type;
	
	typedef typename base_type::element      element;
	typedef typename base_type::size_tuple   size_tuple;
	typedef typename base_type::value_type   value_type;
	typedef typename base_type::reference    reference;
	typedef T*                               pointer;
	typedef pointer                          iterator;

	array( const base_type& other ) : 
	    base_type( new T[other.num_elements()], other.shape() )
        {
	    std::copy( other.begin(), other.end(), base_type::begin() );
	}

	array( const size_tuple& shape ) : 
	    base_type( new T[prod(shape)], shape )
        {
	}
	
	~array()
        {
	    delete[] base_type::_base;
	}

	iterator begin()
	{
	    return base_type::_base;
	}

	iterator end()
	{
	    return base_type::_base + base_type::num_elements();
	}
    };

    // ---------------------------------------------------------------

    template<typename T, size_type N> 
    std::ostream& operator<<( std::ostream& out, array_ref<T,N> a )
    {
	if( N == 1 )
	    for( int i=0; i<a.size(); ++i )
		out << a[i] << '\t';
	else
	    for( int i=0; i<a.size(); ++i )
		out << a[i] << '\n';
	
	return out;
    }

} // namespace nvis

#endif // __array_hpp
