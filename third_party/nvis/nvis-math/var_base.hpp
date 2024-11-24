#ifndef __var_base_hpp
#define __var_base_hpp

#include <algorithm>
#include <iosfwd>

#include <boost/call_traits.hpp>
#include <math/base_types.hpp>

namespace nvis {
    
    template<typename T> class var_base
    {
    public:

        typedef typename boost::call_traits<T>::param_type value_arg;

        typedef T                           value_type;
        typedef T*                          pointer;
        typedef const T*                    const_pointer;
        typedef T&                          reference;
        typedef const T&                    const_reference;
        typedef T*                          iterator;
        typedef const T*                    const_iterator;
        
        // --- constructors ---
        
        var_base() : size_(0), data_(0)
        {
        }
        
        var_base( size_type s ) : 
            size_(s), data_( new value_type[s] )
        {
        }

        var_base( size_type s, value_arg v ) :
            size_(s), data_( new value_type[s] )
        {
            std::fill( begin(), end(), v );
        }    

        var_base( size_type s, const_pointer data ) :
            size_(s), data_( new value_type[s] )
        {
            std::copy( data, data+s, begin() );
        }

        var_base( const var_base& v ) : 
            size_(v.size()), data_( new value_type[v.size()] )
        {
            std::copy( v.begin(), v.end(), begin() );
        }

        ~var_base()
        {
            delete[] data_;
        }
            
        // --- assignment operator ---

        var_base& operator=( const var_base& v )
        {
            if( size_ != v.size() )
                resize( v.size(), false );

            std::copy( v.begin(), v.end(), begin() );
        }

        // --- resize ---

        void resize( size_type size, bool preserve = true, 
                     value_arg v = value_type() )
        {
            if (size != size_) 
            {
                pointer data = new value_type[size];

                if( preserve )  
                {
                    std::copy( data_, data_+std::min(size, size_), data );
                    std::fill( data + std::min( size, size_ ), data+size, v );
                }

                delete [] data_;

                size_ = size;
                data_ = data;
            }
        }

        // --- element access ---
        
        reference operator[]( unsigned int n )
        {
            return data_[n];
        }
        
        const_reference operator[]( unsigned int n ) const
        {
            return data_[n];
        }
        
        // --- iterators ---
        
        iterator begin() { return data_;   }
        iterator end()   { return data_+size_; }
        
        const_iterator begin() const { return data_;   }
        const_iterator end()   const { return data_+size_; }
        
        // --- size ---

        size_type size() const { return size_; }

    protected:
        
        pointer      data_;
        size_type    size_;
    };

} // namespace nvis


#endif // __var_base_hpp
