#ifndef __vector_array_hpp
#define __vector_array_hpp

#include "fixed_vector.hpp"
#include "var_vector.hpp"
#include "array.hpp"

#include <util/implicit_iterator.hpp>

namespace nvis {

    template<typename T> class vector_array: public array<T>
    {
    public:

        typedef var_vector<T>                             value_type;
        typedef implicit_iterator< vector_array<T> >      iterator;

        // --- constructor ---

        vector_array( size_type dim, size_type size = 0 ) : 
            array<T>(size*dim), dim_(dim)
        {
        }

        vector_array( size_type dim, array<T>& r ) : 
            array<T>( r ), dim_(dim)
        {
            assert( r.size() % dim_ == 0 );
        }

        // --- resize ---

        void resize( size_type size, bool preserve = true )
        {
            array<T>::resize( size*dim_ );
        }

        // --- element access ---

        value_type operator[]( size_type index ) const
        {
            return value_type( dim_, array<T>::data() + dim_*index );
        }

        // --- iterator ---

        iterator begin() const { return iterator( this, 0 ); }
        iterator end()   const { return iterator( this, size() );  }

        // --- inspectors ---

        size_type size() const
        {
            return array<T>::size() / dim_;
        }

        size_type dim() const
        {
            return dim_;
        }

    private:

        const size_type dim_;
    };

    // --------------------------------------------------------------------

    template<typename T, std::size_t N> class fixed_vector_array : 
        public array<T>
    {
    public:
        
        typedef fixed_vector<T,N>                         value_type;
        typedef fixed_vector<T,N>*                        pointer;
        typedef const fixed_vector<T,N>*                  const_pointer;
        typedef fixed_vector<T,N>&                        reference;
        typedef const fixed_vector<T,N>&                  const_reference;
        typedef fixed_vector<T,N>*                        iterator;
        typedef const fixed_vector<T,N>*                  const_iterator;
        typedef const fixed_vector<T,N>&                  value_arg;
        
        // --- constructors ---
        
        fixed_vector_array( size_type s ) : array<T>( N*s )
        {
        }
        
        fixed_vector_array( size_type s, value_arg v ) :
            array<T>( N*s )
        {
            std::fill( begin(), end(), v );
        }
        
        fixed_vector_array( array<T>& a ) : array<T>( a )
        {
            assert( a.size() % N == 0 );
        }

        // --- resize ---
        
        void resize( size_type s, bool preserve = true )
        {
            array<T>::resize( N*s, preserve );
        }
    
        // --- element access (fixed_vector) ---
        
        reference operator[]( size_type n )
        {
            return data()[n];
        }
        
        const_reference operator[]( size_type n ) const
        {
            return data()[n];
        }
        
        // --- iterators ---
    
        iterator begin() { return data(); }
        iterator end()   { return data()+size(); }

        const_iterator begin() const { return data(); }
        const_iterator end()   const { return data()+size(); }
        
        // --- size ---
        
        size_type size() const { return array<T>::size()/N; }

    private:

        pointer data()
        {
            return (pointer)array<T>::data();
        }
        
        const_pointer data() const
        {
            return (const_pointer)array<T>::data();
        }
    };
    
    // ---------------------------------------------------------------

    template<typename T> 
    std::ostream& operator<<( std::ostream& out, const vector_array<T>& a )
    {
        typename vector_array<T>::iterator i;

        out << a.size() << " =   [";
        
        for( i=a.begin(); i!=a.end(); ++i )
            out << "\t" << *i << ( i==a.end()-1 ? " ]" : "," ) << '\n';

        return out;
    }

    // ---------------------------------------------------------------

    template<typename T, std::size_t N> 
    std::ostream& operator<<( std::ostream& out, const fixed_vector_array<T,N>& a )
    {
        typename fixed_vector_array<T,N>::const_iterator i;

        out << a.size() << " =   [";
        
        for( i=a.begin(); i!=a.end(); ++i )
            out << "\t" << *i << ( i==a.end()-1 ? " ]" : "," ) << '\n';

        return out;
    }

    // ---------------------------------------------------------------

    typedef vector_array<double>         vec_array;
    typedef vector_array<int>           ivec_array;
    typedef vector_array<unsigned int>  uvec_array;

    // ---------------------------------------------------------------

    typedef fixed_vector_array<double,2>         vec2_array;
    typedef fixed_vector_array<int,2>           ivec2_array;
    typedef fixed_vector_array<unsigned int,2>  uvec2_array;

    typedef fixed_vector_array<double,3>         vec3_array;
    typedef fixed_vector_array<int,3>           ivec3_array;
    typedef fixed_vector_array<unsigned int,3>  uvec3_array;

    typedef fixed_vector_array<double,4>         vec4_array;
    typedef fixed_vector_array<int,4>           ivec4_array;
    typedef fixed_vector_array<unsigned int,4>  uvec4_array;


} // namespace nvis;


#endif // __vector_array_hpp
