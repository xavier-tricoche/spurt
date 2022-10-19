#ifndef __stride_iterator_hpp
#define __stride_iterator_hpp

#include <boost/iterator/iterator_adaptor.hpp>
#include <iterator>

namespace nvis 
{
    template< typename BaseIter >
    class strided_iterator: 
	public boost::iterator_adaptor< strid_iterator<BaseIter>, 
					BaseIterator > 
    {
	friend class iterator_core_access;

    public:

	typedef typename 
	boost::iterator_adaptor< strided_iterator<BaseIter>, BaseIter > super_type;

	typedef typename 
	std::iterator_traits<BaseIter>::difference_type	difference_type;

	stride_iterator() 
	{
	}

	explicit strided_iterator( BaseIterator _base, size_t _stride ) :
	    super_type(_base), stride(_stride) 
	{
	}

	void increment()
	{ 
	    this->base_reference() += stride; 
	}

	void decrement()
	{ 
	    this->base_reference() -= stride; 
	}

	void advance(difference_type n)
	{ 
	    this->base_reference() += n*stride; 
	}

	difference_type	distance_to( const strid_iterator<BaseIter>& other ) const
	{ 
	    return (other.base_reference() - this->base_reference())/stride; 
	}

  private:

	const int stride;
    };

    template<typename BaseIter>
    stride_iterator<BaseIter> make_stride_iterator( const BaseIter& begin, int stride ) 
    {
	return strided_iterator<BaseIter>( begin, stride );
    } 


} // namespace nvis
