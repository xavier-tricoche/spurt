#ifndef __inline_hpp
#define __inline_hpp

namespace nvis { namespace inline_ {

    template<int M, int N> struct forward_loop_type
    {
        template<typename A, typename T1, typename T2>
        void iterate( const A& actor, T1& t1, T2& t2 )
        {
            actor( N-M, t1, t2 );
            forward_loop_type<M-1,N>().iterate( actor, t1, t2 );
        }

        template<typename A, typename T1, typename T2, typename T3>
        void iterate( const A& actor, T1& t1, T2& t2, T3& t3 )
        {
            actor( N-M, t1, t2, t3 );
            forward_loop_type<M-1,N>().iterate( actor, t1, t2, t3 );
        }

        template<typename A, typename T1, typename T2, typename T3, typename T4>
        void iterate( const A& actor, T1& t1, T2& t2, T3& t3, T4& t4 )
        {
            actor( N-M, t1, t2, t3, t4 );
            forward_loop_type<M-1,N>().iterate( actor, t1, t2, t3, t4 );
        }
    };

    // ----------------------------------------------
    
    template<int N> struct forward_loop_type<0,N>
    {
        template<typename A, typename T1, typename T2>
        void iterate( const A&, T1&, T2& )
        {
        }

        template<typename A, typename T1, typename T2, typename T3>
        void iterate( const A&, T1&, T2&, T3& )
        {
        }

        template<typename A, typename T1, typename T2, typename T3, typename T4>
        void iterate( const A&, T1&, T2&, T3&, T4& )
        {
        }
    };
    
    // ----------------------------------------------
    
    template<int N,typename A,typename T1, typename T2> 
    inline void loop( const A& actor, T1& t1, T2& t2 )
    {
        forward_loop_type<N,N>().iterate( actor, t1, t2 );
    }

    template<int N,typename A,typename T1, typename T2, typename T3 > 
    inline void loop( const A& actor, T1& t1, T2& t2, T3& t3 )
    {
        forward_loop_type<N,N>().iterate( actor, t1, t2, t3 );
    }

    template<int N,typename A,typename T1, typename T2, typename T3, typename T4 > 
    inline void loop( const A& actor, T1& t1, T2& t2, T3& t3, T4& t4 )
    {
        forward_loop_type<N,N>().iterate( actor, t1, t2, t3, t4 );
    }

    // ----------------------------------------------

    template<typename Op> struct apply_n_1
    {
        template<typename T1, typename T2> 
        void operator()( int n, T1& t, const T2& p ) const
        {
            Op()( t[n], p );
        }

    };

    template<typename Op> struct apply_1_n
    {
        template<typename T1, typename T2> 
        void operator()( int n, T1& t, const T2& p ) const
        {
            Op()( t, p[n] );
        }

    };

    template<typename Op> struct apply_n_n
    {
        template<typename T1, typename T2> 
        void operator()( int n, T1& t1, const T2& t2 ) const
        {
            Op()( t1[n], t2[n] );
        }
    };
    

    template<typename Op> struct apply_1_n_n
    {
        template<typename T1, typename T2, typename T3> 
        void operator()( int n, T1& t1, const T2& t2, const T3& t3 ) const
        {
            Op()( t1, t2[n], t3[n] );
        }
    };

    template<typename Op> struct apply_n_n_n
    {
        template<typename T1, typename T2, typename T3> 
        void operator()( int n, T1& t1, const T2& t2, const T3& t3 ) const
        {
            Op()( t1[n], t2[n], t3[n] );
        }
    };

    // ----------------------------------------------

    template<int N, typename Op, typename T1, typename T2> 
    inline void for_each_n_1( T1& t1, const T2& t2 )
    {
        forward_loop_type<N,N>().iterate( apply_n_1<Op>(), t1, t2 );
    }

    template<int N, typename Op, typename T1, typename T2> 
    inline void for_each_1_n( T1& t1, const T2& t2 )
    {
        forward_loop_type<N,N>().iterate( apply_1_n<Op>(), t1, t2 );
    }

    template<int N, typename Op, typename T1, typename T2> 
    inline void for_each_n_n( T1& t1, const T2& t2 )
    {
        forward_loop_type<N,N>().iterate( apply_n_n<Op>(), t1, t2 );
    }

    template<int N, typename Op, typename T1, typename T2, typename T3> 
    inline void for_each_1_n_n( T1& t1, const T2& t2, const T3& t3 )
    {
        forward_loop_type<N,N>().iterate( apply_1_n_n<Op>(), t1, t2, t3 );
    }

    template<int N, typename Op, typename T1, typename T2, typename T3> 
    inline void for_each_n_n_n( T1& t1, const T2& t2, const T3& t3 )
    {
        forward_loop_type<N,N>().iterate( apply_n_n_n<Op>(), t1, t2, t3 );
    }

    // ----------------------------------------------

    struct assign
    {
        template<typename T1, typename T2>
        void operator()( T1& to, const T2& val ) const
        {
            to = val;
        }
    };

    // comparison ops

    struct min
    {
	template<typename T1, typename T2>
	void operator()( T1& to, const T2& val )
	{
	    if( val < to )
		to = val;
	}
    };

    struct max
    {
	template<typename T1, typename T2>
	void operator()( T1& to, const T2& val )
	{
	    if( val > to )
		to = val;
	}
    };

    // inplace ops
    
    struct inplace_mul
    {
        template<typename T1, typename T2>
        void operator()( T1& to, const T2& val ) const
        {
            to *= val;
        }
    };

    struct inplace_div
    {
        template<typename T1, typename T2>
        void operator()( T1& to, const T2& val ) const
        {
            to /= val;
        }
    };


    struct inplace_add
    {
        template<typename T1, typename T2>
        void operator()( T1& to, const T2& val ) const
        {
            to += val;
        }
    };

    struct inplace_sub
    {
        template<typename T1, typename T2>
        void operator()( T1& to, const T2& val ) const
        {
            to -= val;
        }
    };

    struct inplace_mul_add
    {
        template<typename T1, typename T2, typename T3>
        void operator()( T1& to, const T2& val1, const T3& val2 ) const
        {
            to += val1*val2;
        }
    };

} } // namespace nvis, namespace inline_

// ----------------------------------------------

#endif // __inline_hpp
