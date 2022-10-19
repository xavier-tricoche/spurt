#ifndef __tiny_array_hpp
#define __tiny_array_hpp

#include <math/inline_.hpp>
#include <math/base_types.hpp>

namespace nvis {

    template<typename T, unsigned int N> class tiny_array
    {
    public:
        typedef tiny_array<T,N>             self_type;

        typedef T                           value_type;
        typedef T*                          pointer;
        typedef T&                          reference;
        typedef const T&                    const_reference;
        typedef T*                          iterator;
        typedef const T*                    const_iterator;
	typedef std::ptrdiff_t              difference_type;
	typedef std::size_t                 size_type;

        static const size_type fixed_size = N;

        // --- constructors ---

        tiny_array()
        {
        }

        tiny_array( const tiny_array& other )
        {
            inline_::for_each_n_n< N, inline_::assign >( data_, other.data_ );
        }

        tiny_array( const value_type* ptr )
        {
            inline_::for_each_n_n< N, inline_::assign >( data_, ptr );
        }

        // --- element access ---
        
        reference operator[]( size_type n )
        {
            return data_[n];
        }
        
        const_reference operator[]( size_type n ) const
        {
            return data_[n];
        }
        
        // --- iterators ---
        
        iterator begin() { return data_;   }
        iterator end()   { return data_+N; }
        
        const_iterator begin() const { return data_;   }
        const_iterator end()   const { return data_+N; }
        
        // --- swap ---
	
	void swap( tiny_array& o )
	{
	    std::swap_ranges( begin(), end(), o.begin() );
	}

        // --- size ---

        static size_type size() { return N; }

	static bool empty() { return false; }

        // --- raw data ---

        pointer data()
        {
            return (pointer)data_;
        }

        const pointer data() const
        {
            return (const pointer)data_;
        }

    protected:

        value_type data_[N];
    };


} // namespace nvis

#endif // __tiny_array_hpp
