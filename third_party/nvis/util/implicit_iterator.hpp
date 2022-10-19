#ifndef __implicit_iterator_hpp
#define __implicit_iterator_hpp

#include <boost/iterator/iterator_facade.hpp>
#include <math/base_types.hpp>

namespace nvis {

    template<typename T> struct square_bracket
    {
        typedef typename T::value_type result_type;

        result_type operator()( const T& t, size_type i )
        {
            return t[i];
        }
    };

    template< typename T, typename Op = square_bracket<T> >
    class implicit_iterator: 
        public boost::iterator_facade< implicit_iterator<T,Op>, 
                                       typename Op::result_type,
                                       boost::bidirectional_traversal_tag,
                                       typename Op::result_type >
    {
    public: 
        
        typedef typename Op::result_type     value_type;
        typedef Op                           op_type;
        typedef T                            container_type;
        
        explicit implicit_iterator( const container_type& c, size_type i ) : 
            container(&c), index(i), cached( op_type()( c, i ) )
        {
        }

        implicit_iterator() : container(0), index(0)
        {
        }
    
    private:

        void cache()
        {
            cached = op_type()( *container, index );
        }

        // --- iterator_facade implementation ---
        
        friend class boost::iterator_core_access;
        
        const value_type& dereference() const
        {
            return cached;
        }
        
        void increment() { index++; cache(); }
        void decrement() { index--; cache(); }
        void advance( size_type n ) { index+=n; cache(); }
        
        bool equal( const implicit_iterator& other ) const
        {
            return other.index == index;
        }
        
        std::ptrdiff_t distance_to( const implicit_iterator& other ) const
        {
            return other.index - index;
        }

        // --- data members ---
        
        const container_type*   container;
        size_type               index;
        value_type              cached;
    };

} // namespace nvis

#endif // __implicit_iterator_hpp
