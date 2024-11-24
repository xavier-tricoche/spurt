#ifndef __var_vector_hpp
#define __var_vector_hpp

#include <algorithm>
#include <iosfwd>
#include <cassert>
#include <cmath>

#include <boost/call_traits.hpp>
#include <boost/operators.hpp>
#include <boost/preprocessor.hpp>

#include <math/var_base.hpp>
//#include <math/fixed_vector.hpp>

#define VAR_VECTOR_MAX_CONSTRUCTOR_ARGS 6

#define VAR_VECTOR_CTOR_AUX(z,n,_) base::data_[n] = v##n;

#define VAR_VECTOR_CTOR(z,n,_)                                   \
        var_vector( BOOST_PP_ENUM_PARAMS( n, value_arg v ) ) :   \
            base( n )                                            \
        {                                                        \
            BOOST_PP_REPEAT( n, VAR_VECTOR_CTOR_AUX, _ )         \
        }                                                        

namespace nvis {
    
    template<typename T> 
    class var_vector:
        public var_base<T>,
        boost::additive< var_vector<T>,
        boost::multiplicative< var_vector<T>, T > >
    {
    public:

        typedef var_base<T>                                 base;
        typedef typename boost::call_traits<T>::param_type  value_arg;

        typedef T                           value_type;
        typedef T*                          pointer;
        typedef const T*                    const_pointer;
        typedef T&                          reference;
        typedef const T&                    const_reference;
        typedef T*                          iterator;
        typedef const T*                    const_iterator;
        
        // --- constructors ---
        
        template<unsigned int N> var_vector( const fixed_vector<T,N>& v ) : 
            base( N, v.begin() )
        {
        }

        var_vector()
        {
        }
        
        var_vector( size_type s ) : base(s)
        {
        }

        var_vector( size_type s, value_arg v ) : base( s, v )
        {
        }    

        var_vector( size_type s, const_pointer data ) : base( s, data )
        {
        }

           
        BOOST_PP_REPEAT_FROM_TO( 2, VAR_VECTOR_MAX_CONSTRUCTOR_ARGS, 
                                 VAR_VECTOR_CTOR, _ ) 
            
        // handled by base class:
        // copy constr. + assignment op., destructor,
        // element access, iterators

        // --- arithmetic operators ---
        
        var_vector& operator+=( const var_vector& rhs )
        {
            assert( rhs.size() == base::size_ );

            for( size_type i=0; i<base::size_; ++i )
                base::data_[i] += rhs[i];
            
            return *this;
        }

        var_vector& operator-=( const var_vector& rhs )
        {
            assert( rhs.size() == base::size_ );

            for( size_type i=0; i<base::size_; ++i )
                base::data_[i] -= rhs[i];
            
            return *this;
        }

        var_vector& operator*=( value_arg rhs )
        {
            for( size_type i=0; i<base::size_; ++i )
                base::data_[i] *= rhs;
            
            return *this;
        }

        var_vector& operator/=( value_arg rhs )
        {
            for( size_type i=0; i<base::size_; ++i )
                base::data_[i] /= rhs;
            
            return *this;
        }
    };

    // -----------------------------------------------------------------

    template<typename T>
    typename var_vector<T>::value_type 
    inner( const var_vector<T>& v1, 
           const var_vector<T>& v2 )
    {
        typename var_vector<T>::value_type ip = 0;

        assert( v1.size() == v2.size() );

        for( size_type i=0; i<v1.size(); ++i )
            ip += v1[i] * v2[i];
        
        return ip;
    }

    template<typename T>
    var_vector<T> cross( const var_vector<T>& v1, 
                         const var_vector<T>& v2 )
    {
        assert( v1.size() == 3 && v2.size() == 3 );

        var_vector<T> cp( 3 );
        
        cp[0] = v1[1]*v2[2] - v2[1]*v1[2];
        cp[1] = v1[2]*v2[0] - v2[2]*v1[0];
        cp[2] = v1[0]*v2[1] - v2[0]*v1[1];
        
        return cp;
    }

    template<typename T>
    typename var_vector<T>::value_type
    norm( const var_vector<T>& v )
    {
        return sqrt( inner(v,v) );
    }

    // -----------------------------------------------------------------

    template<typename T>
    std::ostream& operator<<( std::ostream& out, const var_vector<T>& v )
    {
        typename var_vector<T>::const_iterator i;

        out << "[ ";

        if( v.size() )
        {
            for( i=v.begin(); i!=v.end()-1; ++i )
                out << *i << ", ";

            out << *i;
        }

        out << " ]";

        return out;
    }

    // -----------------------------------------------------------------

    typedef var_vector<double>         vec;
    typedef var_vector<int>           ivec;
    typedef var_vector<unsigned int>  uvec;

} // namespace nvis


#undef VAR_VECTOR_CTOR
#undef VAR_VECTOR_CTOR_AUX

#endif // __var_vector_hpp
