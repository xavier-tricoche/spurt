#ifndef __bounding_box_hpp
#define __bounding_box_hpp

#include <math/fixed_vector.hpp>
#include <limits>
#include <algorithm>

// nvis::bounding_box
// Copyright Christoph Garth 2005
namespace nvis {

    template<typename T>
    struct _less {
        bool operator()(const T& a, const T& b) const {
            return a<b;
        }
    };

    template<typename T, size_t N>
    struct _less<spurt::fixed_vector<T,N>> {
        typedef spurt::fixed_vector<T,N> value_type;
        bool operator()(const value_type& a, const value_type& b) const {
            return spurt::all(a < b);
        }
    };

    template<typename T> class bounding_box
    {
        typedef T                                value_type;
        typedef typename T::value_type           inner_type;
        typedef std::numeric_limits<inner_type>  limits_type;
        typedef _less<T>                         less_type;

    public:

        // --- constructors ---

        bounding_box( const value_type& min, const value_type& max ) :
        min_(min), max_(max)
        {
        }

        bounding_box()
        {
            reset();
        }

        template<typename Iter>
        bounding_box( Iter begin, Iter end )
        {
            set( begin, end );
        }

        // ---

        void reset()
        {
            min_ = value_type(  limits_type::max() );
            max_ = value_type( -limits_type::max() );
        }

        // --- main functionality ---

        void add( const value_type& p )
        {
            min_ = std::min( min_, p, less_type() );
            max_ = std::max( max_, p, less_type() );
        }

        template<typename Iter>
        void add( Iter begin, Iter end )
        {
            while( begin != end )
                add( *(begin++) );
        }

        template<typename Iter>
        void set( Iter begin, Iter end )
        {
            reset();
            add( begin, end );
        }

        bool inside( const value_type& p ) const
        {
            return all( p >= min_ ) && all( p <= max_ );
        }

        bool empty() const {
            return any( min_ >= max_ );
        }

        // --- inspection ---

        value_type& min()
        {
            return min_;
        }

        const value_type& min() const
        {
            return min_;
        }

        value_type& max()
        {
            return max_;
        }

        const value_type& max() const
        {
            return max_;
        }

        value_type center() const
        {
            return (min_+max_)/2;
        }

        value_type size() const
        {
            return max_ - min_;
        }

        typename value_type::value_type volume() const
        {
            return prod( max_ - min_ );
        }

        void scale( const double& factor )
        {
            value_type c = center();

            min_ = c + factor * (min_ - c);
            max_ = c + factor * (max_ - c);
        }

        // --- normalizer ---

        class normalizer
        {
            const value_type min_, ext_;

        public:
            normalizer( const bounding_box& box ) :
            min_( box.min() ),
            ext_( box.size() )
            {
            }

            value_type operator()( const value_type& v )
            {
                return (v-min_)/ext_;
            }
        };

        normalizer normalize() const
        {
            return normalizer( *this );
        }

        // --- random point ---

        value_type random() const
        {
            value_type p = min_;
            value_type s = max_ - min_;

            for( int i=0; i<p.size(); ++i )
                p[i] += drand48() * s[i];

            return p;
        }

        // --- predicates ---

        struct inside_predicate
        {
            const bounding_box& box;

            inside_predicate( const bounding_box& _box ) :
            box(_box)
            {
            }

            bool operator()( const value_type& p )
            {
                return box.inside( p );
            }
        };

        struct outside_predicate
        {
            const bounding_box& box;

            outside_predicate( const bounding_box& _box ) :
            box(_box)
            {
            }

            bool operator()( const value_type& p )
            {
                return !box.inside( p );
            }
        };

    private:

        value_type min_, max_;
    };

    // --- ostream operator ---

    template<typename T> std::ostream& operator<<( std::ostream& out,
    const bounding_box<T>& bb )
    {
        return out << "[\n\tmin =  " << bb.min()
            << ",\n\tmax =  " << bb.max()
                << ",\n\tsize = " << bb.size()
                    << ",\n\tvol =  " << bb.volume()
                        << "\n]\n";
    }

    // --- predicate generators ---

    template<typename Box>
    typename Box::inside_predicate inside( const Box& box )
    {
        return typename Box::inside_predicate( box );
    }

    template<typename Box>
    typename Box::outside_predicate outside( const Box& box )
    {
        return typename Box::outside_predicate( box );
    }

} // namespace nvis

namespace spurt {
    // --- commonly used types ---

    typedef nvis::bounding_box<spurt::vec2> bbox2;
    typedef nvis::bounding_box<spurt::vec3> bbox3;
}

#endif // __bounding_box_hpp
