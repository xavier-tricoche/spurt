#pragma once

#include <algorithm>
#include <array>
#include <complex>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <random>
#include <tuple>
#include <type_traits>

#include <math/better_type.hpp>
#include <math/basic_math.hpp>
#include <math/data_block.hpp>
#include <Eigen/Eigen>

namespace spurt 
{
    // small_vector_interface: 
    // minimal fixed size array interface
    template <typename Storage_>
    class small_vector_interface
    {
    public:
        typedef typename Storage_::value_type scalar_type;
        typedef typename Storage_::value_type value_type;
        typedef Storage_ storage_type;
        typedef typename storage_type::iterator iterator;
        typedef typename storage_type::const_iterator const_iterator;
        
        typedef small_vector_interface<storage_type> self_type;
        static constexpr size_t _size_ = storage_type::size;
        
        typedef typename storage_type::eigen_vector_type eigen_type;
        typedef typename storage_type::eigen_vector_map_type eigen_map_type;
        typedef typename storage_type::const_eigen_vector_map_type const_eigen_map_type;
        
        typedef internal::random_filler<value_type> random_filler_type;
        
        template<typename OtherStorage, typename = typename std::enable_if<OtherStorage::size == _size_>::type >
        using matching_vector = small_vector_interface<OtherStorage>;
        
        template<typename OtherStorage, typename = typename std::enable_if<OtherStorage::size == _size_>::type >
        using matching_array = std::array<typename OtherStorage::value_type, _size_> ;
        
        template<typename T>
        using datablock = DataBlock<T, _size_>;
        
        typedef matching_vector<datablock<bool>> bool_vector;
        
        template<typename OtherStorage>
        using right_S_type = typename better_type<value_type, typename OtherStorage::value_type>::type;
        
        template<typename T>
        using right_T_type = typename better_type<value_type, T>::type; 
        
        template<typename T>
        using datablock_view = DataBlockView<T, _size_, 1>;
        template<typename T>
        using const_datablock_view = ConstDataBlockView<T, _size_, 1>;
        
        template<typename OtherStorage>
        using right_S_storage = DataBlock<right_S_type<OtherStorage>, _size_>;
        template<typename T>
        using right_T_storage = DataBlock<right_T_type<T>, _size_>;
        
        template<typename OtherStorage, typename = typename std::enable_if<OtherStorage::size == _size_>::type>
        using right_S_output = matching_vector< right_S_storage<OtherStorage> >;
        template<typename T, typename = typename internal::must_be_scalar<T>::type>
        using right_T_output = matching_vector< right_T_storage<T> >;
        
        //---------------------------------
        //
        // Constructors and Initializers
        //
        //---------------------------------
        
        iterator begin() { return m_storage.begin(); }
        iterator end() { return m_storage.end(); }
        
        const_iterator begin() const { return m_storage.begin(); }
        const_iterator end() const { return m_storage.end(); }

        small_vector_interface(scalar_type val=0) : m_storage() {
            std::fill(begin(), end(), val);
        }
        
        small_vector_interface(storage_type storage) : m_storage(storage) {}
        
        small_vector_interface(const self_type& other) {
            std::copy(other.begin(), other.end(), begin());
        }

        template<typename OtherStorage>
        small_vector_interface(const matching_vector<OtherStorage>& other) 
            : m_storage() 
        {
            std::copy(other.begin(), other.end(), begin());
        }

        template<typename OtherStorage>
        small_vector_interface(const matching_array<OtherStorage>& other)
            : m_storage() 
        {
            std::copy(other.begin(), other.end(), begin());
        }
        
        template<typename T>
        small_vector_interface(const T* ptr) : m_storage() {
            std::copy(ptr, ptr+_size_, begin());
        }

        template<typename T, typename = typename internal::must_be_scalar<T>::type >
        small_vector_interface(std::initializer_list<T> vals)
            : m_storage() 
        {
            std::copy(vals.begin(), vals.end(), begin());
        }

        template<typename T1, typename T2, typename = typename internal::must_all_be_scalar<T1, T2>::type>
        small_vector_interface(T1 v0, T2 v1)
            : m_storage() 
        {
            std::fill(begin(), end(), 0);
            static_assert(_size_ >= 2, 
                          "Invalid initializer for array of size < 2");
            m_storage[0] = v0;
            m_storage[1] = v1;
        }

        template <typename T1, typename T2, typename T3,
                  typename = typename internal::must_all_be_scalar<T1, T2, T3>::type>
        small_vector_interface(T1 v0, T2 v1, T3 v2)
            : m_storage() 
        {
            std::fill(begin(), end(), 0);
            static_assert(_size_ >= 3,
                          "Invalid initializer for array of size < 3");
            m_storage[0] = v0;
            m_storage[1] = v1;
            m_storage[2] = v2;
        }

        template <typename T1, typename T2, typename T3, typename T4,
                typename = typename internal::must_all_be_scalar<T1, T2, T3, T4>::type>
        small_vector_interface(T1 v0, T2 v1, T3 v2, T4 v3)
            : m_storage() 
        {
            std::fill(begin(), end(), 0);
            static_assert(_size_ >= 4,
                          "Invalid initializer for array of size < 4");
            m_storage[0] = v0;
            m_storage[1] = v1;
            m_storage[2] = v2;
            m_storage[3] = v3;
        }
        
        const storage_type& storage() const { return m_storage; }
        
        const value_type* data() const { 
            return static_cast<const value_type*>(m_storage.data);
        }
        
        size_t size() const { return _size_; }

        template <typename T, typename = 
        typename std::enable_if<std::is_scalar<T>::value>::type>
        self_type& operator=(const T& val) 
        {
            std::fill(begin(), end(), static_cast<scalar_type>(val));
            return *this;
        }
        
        self_type& operator=(const self_type& other) {
            std::copy(other.begin(), other.end(), begin());
            return *this;
        }
        
        template<typename OtherStorage, 
                 typename = typename std::enable_if<!std::is_same<storage_type, OtherStorage>::value>::type>
        self_type& operator=(const matching_vector<OtherStorage>& other)
        {
            std::copy(other.begin(), other.end(), begin());
            return *this;
        }

        // miscellaneous
        static self_type random(value_type a=random_filler_type::default_min, 
                                value_type b=random_filler_type::default_max) {
            random_filler_type filler(a, b);
            self_type r;
            filler.fill(r.begin(), r.end());
            return r;
        }
        
        scalar_type norm() const {
            return as_const_eigen().norm();
        }

        //---------------------------------
        //
        //       Accessors
        //
        //---------------------------------
        value_type& operator[](size_t i) { return m_storage[i]; }
        const value_type& operator[](size_t i) const { return m_storage[i]; }
        
        value_type& operator()(size_t i) { return m_storage[i]; }
        const value_type& operator()(size_t i) const { return m_storage[i]; }

        //---------------------------------
        //
        //       Comparison
        //
        //---------------------------------
        template <typename OtherStorage>
        bool_vector operator<(const matching_vector<OtherStorage> &other) const
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { r[i] = ((*this)[i] < other[i]); }
            return r;
        }
        template <typename OtherStorage>
        bool_vector operator<=(const matching_vector<OtherStorage> &other) const
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { r[i] = ((*this)[i] <= other[i]); }
            return r;
        }

        template<typename T, typename = typename internal::must_be_scalar<T>::type >
        bool_vector operator<(T val) const 
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { r[i] = ((*this)[i] < val); }
            return r;
        }

        template<typename T, typename = typename internal::must_be_scalar<T>::type >
        bool_vector operator<=(T val) const 
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { r[i] = ((*this)[i] <= val); }
            return r;
        }

        template <typename OtherStorage>
        bool_vector operator>(const matching_vector<OtherStorage> &other) const
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { r[i] = ((*this)[i] > other[i]); }
            return r;
        }

        template <typename OtherStorage>
        bool_vector operator>=(const matching_vector<OtherStorage> &other) const
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { r[i] = ((*this)[i] >= other[i]); }
            return r;
        }

        template <typename T, typename = typename internal::must_be_scalar<T>::type >
        bool_vector operator>(T val) const
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { r[i] = (m_storage[i] > val); }
            return r;
        }

        template <typename T, typename = typename internal::must_be_scalar<T>::type >
        bool_vector operator>=(T val) const
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { r[i] = (m_storage[i] >= val); }
            return r;
        }

        template <typename OtherStorage>
        bool_vector operator==(const matching_vector<OtherStorage>& other) const 
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { 
                r[i] = (m_storage[i] == static_cast<scalar_type>(other[i]));
            }
            return r;
        }

        template <typename T, typename = typename internal::must_be_scalar<T>::type>
        bool_vector operator==(T val) const
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { 
                r[i] = (m_storage[i] == static_cast<scalar_type>(val));
            }
            return r;
        }

        template <typename OtherStorage>
        bool_vector operator!=(const matching_vector<OtherStorage>& other) const 
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { 
                r[i] = (m_storage[i] != static_cast<scalar_type>(other[i]));
            }
            return r;
        }

        template <typename T, typename = typename internal::must_be_scalar<T>::type>
        bool_vector operator!=(T val) const
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { 
                r[i] = (m_storage[i] != static_cast<scalar_type>(val));
            }
            return r;
        }

        //---------------------------------
        //
        //     Multiplication operators
        //
        //---------------------------------        
        template<typename OtherStorage>
        self_type& operator*=(const matching_vector<OtherStorage>& other)
        {
            for (size_t i=0; i<_size_; ++i) 
            {
                m_storage[i] *= other[i];
            }
            return (*this);
        }
        
        template <typename T, typename = typename internal::must_be_scalar<T>::type >
        self_type &operator*=(T val)
        {
            scalar_type rhs = static_cast<scalar_type>(val);
            for (size_t i=0; i<_size_; ++i)
            {
                m_storage[i] *= rhs;
            }
            return (*this);
        }
        
        template <typename OtherStorage>
        right_S_output<OtherStorage>
        operator*(const matching_vector<OtherStorage> &other) const
        {
            typedef right_S_type<OtherStorage> out_type; 
            right_S_output<OtherStorage> r;
            for (size_t i=0; i<_size_; ++i) {
                r[i] = static_cast<out_type>(m_storage[i]) * static_cast<out_type>(other[i]);
            }
            return r;
        }
        
        template <typename T, typename = typename internal::must_be_scalar<T>::type>
        right_T_output<T>
        operator*(T val) const
        {
            typedef right_T_type<T> out_type; 
            right_T_output<T> r;
            for (size_t i=0; i<_size_; ++i) {
                r[i] = static_cast<out_type>(m_storage[i]) * static_cast<out_type>(val);
            }
            return r;
        }

        //---------------------------------
        //
        //       Division operators
        //
        //---------------------------------
        template <typename OtherStorage>
        self_type &operator/=(const matching_vector<OtherStorage> &other)
        {
            for (size_t i=0; i<_size_; ++i) 
            {
                m_storage[i] /= other[i];
            }
            return (*this);
        }
        
        template <typename T, typename = typename internal::must_be_scalar<T>::type >
        self_type &operator/=(const T &val)
        {
            for (size_t i=0; i<_size_; ++i) 
            {
                m_storage[i] /= val;
            }
            return (*this);
        }
       
        template <typename OtherStorage>
        right_S_output<OtherStorage>
        operator/(const matching_vector<OtherStorage> &other) const
        {
            typedef right_S_type<OtherStorage> out_type;
            right_S_output<OtherStorage> r;
            for (int i = 0; i < _size_; ++i)
            {
                out_type v0 = static_cast<out_type>(m_storage[i]);
                out_type v1 = static_cast<out_type>(other[i]);
                r[i] = v0 / v1;
            }
            return r;
        }
        
        template <typename T, typename = typename internal::must_be_scalar<T>::type>
        right_T_output<T> operator/(const T &val) const
        {
            typedef right_T_type<T> out_type;
            right_T_output<T> r;
            out_type rhs = static_cast<out_type>(val);
            for (int i = 0; i < _size_; ++i)
            {
                out_type v0 = static_cast<out_type>(m_storage[i]);
                r[i] = v0 / rhs;
            }
            return r;
        }

        //---------------------------------
        //
        //       Addition operators
        //
        //---------------------------------
        template <typename OtherStorage>
        self_type &operator+=(const matching_vector<OtherStorage> &other)
        {
            for (int i=0; i<_size_; ++i) 
            {
                m_storage[i] += static_cast<scalar_type>(other[i]);
            }
            return (*this);
        }
        template <typename T, typename = typename internal::must_be_scalar<T>::type >
        self_type &operator+=(const T &val)
        {
            scalar_type rhs = static_cast<scalar_type>(val);
            for (int i = 0; i < _size_; ++i)
            {
                m_storage[i] += rhs;
            }
            return (*this);
        }
        template <typename OtherStorage>
        right_S_output<OtherStorage>
        operator+(const matching_vector<OtherStorage> &other) const
        {
            typedef right_S_type<OtherStorage> out_type;
            right_S_output<OtherStorage> r;
            for (int i=0; i<_size_; ++i) 
            {
                out_type v0 = static_cast<out_type>(m_storage[i]);
                out_type v1 = static_cast<out_type>(other[i]);
                r[i] = v0+v1;
            }
            return r;
        }
        template <typename T, typename = typename internal::must_be_scalar<T>::type>
        right_T_output<T>
        operator+(const T &val) const
        {
            typedef right_T_type<T> out_type;
            right_T_output<T> r;
            out_type rhs = static_cast<out_type>(val);
            for (int i=0; i<_size_; ++i) 
            {
                r[i] = static_cast<out_type>(m_storage[i]) + rhs;
            }
            return r;
        }
        const self_type& operator+() const
        {
            return *this;
        }

        //---------------------------------
        //
        //     Subtraction operators
        //
        //---------------------------------
        template <typename OtherStorage>
        self_type &operator-=(const matching_vector<OtherStorage> &other)
        {
            for (int i = 0; i < _size_; ++i)
            {
                m_storage[i] -= static_cast<scalar_type>(other[i]);
            }
            return (*this);
        }
        template <typename T, typename = typename internal::must_be_scalar<T>::type>
        self_type &operator-=(const T &val)
        {
            scalar_type rhs = static_cast<scalar_type>(val);
            for (int i = 0; i < _size_; ++i)
            {
                m_storage[i] -= rhs;
            }
            return (*this);
        }
        template <typename OtherStorage>
        right_S_output<OtherStorage>
        operator-(const matching_vector<OtherStorage> &other) const
        {
            typedef right_S_type<OtherStorage> out_type;
            right_S_output<OtherStorage> r;
            for (int i = 0; i < _size_; ++i)
            {
                out_type v0 = static_cast<out_type>(m_storage[i]);
                out_type v1 = static_cast<out_type>(other[i]);
                r[i] = v0 - v1;
            }
            return r;
        }
        template <typename T, typename = typename internal::must_be_scalar<T>::type>
        right_T_output<T> operator-(const T &val) const
        {
            typedef right_T_type<T> out_type;
            right_T_output<T> r;
            out_type rhs = static_cast<out_type>(val);
            for (int i=0; i<_size_; ++i) 
            {
                r[i] = static_cast<out_type>(m_storage[i]) - rhs;
            }
            return r;
        }
        self_type operator-() const {
            self_type r;
            r.as_eigen() = -as_const_eigen();
            return r;
        }
    
        const_eigen_map_type as_const_eigen() const
        {
            return const_eigen_map_type(m_storage.data);
        }

        eigen_map_type as_eigen() 
        {
            return eigen_map_type(m_storage.data);
        }
        
        storage_type m_storage;
    };
    

    template<typename T, size_t N>
    using small_vector = small_vector_interface< DataBlock<T, N> >;

    // standard array-type / linear algebra-type vector manipulation
    template <typename Storage1, typename Storage2, 
              typename = typename internal::must_be_3d<Storage1, Storage2>::type>
    small_vector_interface<typename internal::better_storage<Storage1, Storage2>::type > 
    cross(const small_vector_interface<Storage1> &a, const small_vector_interface<Storage2> &b)
    {
        typedef typename internal::better_storage<Storage1, Storage2>::type out_storage;
        typedef typename out_storage::value_type out_type;
        small_vector_interface<out_storage> x;
        x[0] = static_cast<out_type>(a[1])*static_cast<out_type>(b[2]) - 
               static_cast<out_type>(a[2])*static_cast<out_type>(b[1]);
        x[1] = static_cast<out_type>(a[2])*static_cast<out_type>(b[0]) - 
               static_cast<out_type>(a[0])*static_cast<out_type>(b[2]);
        x[2] = static_cast<out_type>(a[0])*static_cast<out_type>(b[1]) - 
               static_cast<out_type>(a[1])*static_cast<out_type>(b[0]);
        return x;
    }

    template <typename Storage1, typename Storage2, 
              typename = typename internal::same_size<Storage1, Storage2>::type>
    typename internal::better_storage<Storage1, Storage2>::value_type
    inner(const small_vector_interface<Storage1> &a, const small_vector_interface<Storage2> &b)
    {
        typedef typename internal::better_storage<Storage1, Storage2>::value_type type;
        type sum = 0;
        for (size_t i=0; i<Storage1::size; ++i) 
        {
            sum += static_cast<type>(a[i])*static_cast<type>(b[i]);
        }
        return sum;
    }

    template <typename Storage>
    typename Storage::value_type
    norm_square(const small_vector_interface<Storage> &a)
    {
        return inner(a, a);
    }

    template <typename Storage>
    typename Storage::value_type
    norm(const small_vector_interface<Storage> &a)
    {
        return sqrt(norm_square(a));
    }
    
    template <typename Storage>
    typename Storage::value_type
    l1_norm(const small_vector_interface<Storage>& a) 
    {
        typename Storage::value_type r = 0;
        std::for_each(a.begin(), a.end(), [&](auto v) { 
            r += std::abs(v);
        });
        return r;
    }
    
    template <typename Storage>
    typename Storage::value_type
    linf_norm(const small_vector_interface<Storage>& a)
    {
        typename Storage::value_type r=0;
        std::for_each(a.begin(), a.end(), [&](auto v) {
            r = std::max(r, std::abs(v));
        });
        return r;
    }

    template <typename T, typename Storage>
    small_vector_interface<typename internal::better_storage<Storage, T>::type>
    operator+(const T &val, const small_vector_interface<Storage> &a)
    {
        return a + val;
    }

    template <typename T, typename Storage>
    small_vector_interface<typename internal::better_storage<Storage, T>::type > 
    operator-(const T &val, const small_vector_interface<Storage> &a)
    {
        return a - val;
    }

    template <typename T, typename Storage>
    small_vector_interface<typename internal::better_storage<Storage, T>::type> 
    operator*(const T &val, const small_vector_interface<Storage> &a)
    {
        return a * val;
    }

    // Reductions
    template <typename Storage>
    bool all(const small_vector_interface<Storage> &a)
    {
        return std::all_of(a.begin(), a.end(), internal::trivial_test());
    }

    template<typename Storage, typename Predicate>
    bool all(const small_vector_interface<Storage>& a, const Predicate& pred=Predicate())
    {
        return std::all_of(a.begin(), a.end(), pred);
    }

    template<typename Storage, 
             typename = typename std::enable_if<std::is_same<typename Storage::value_type, bool>::value>::type>
    bool any(const small_vector_interface<Storage>& a) 
    {
        return std::any_of(a.begin(), a.end(), internal::trivial_test());
    }

    template <typename Storage, typename Predicate>
    bool any(const small_vector_interface<Storage> &a, const Predicate &pred = Predicate())
    {
        return std::any_of(a.begin(), a.end(), pred);
    }

    template <typename Storage, typename = typename std::enable_if<std::is_same<typename Storage::value_type, bool>::value>::type>
    bool none(const small_vector_interface<Storage> &a)
    {
        return std::none_of(a.begin(), a.end(), internal::trivial_test());
    }

    template <typename Storage, typename Predicate>
    bool none(const small_vector_interface<Storage> &a, const Predicate &pred = Predicate())
    {
        return std::none_of(a.begin(), a.end(), pred);
    }

    template<typename Storage>
    typename Storage::value_type min(const small_vector_interface<Storage>& a) 
    {
        return *std::min_element(a.begin(), a.end());
    }

    template <typename Storage>
    typename Storage::value_type max(const small_vector_interface<Storage> &a)
    {
        return *std::max_element(a.begin(), a.end());
    }
    
    template<typename Storage>
    typename Storage::value_type sum(const small_vector_interface<Storage>& v)
    {
        typedef typename Storage::value_type value_type;
        value_type r = 0;
        for (auto it=v.begin(); it!=v.end() ; ++it) 
        {
            r += *it;
        }
        return r;
    }

    template<typename Storage, typename ReduceOperator, 
             typename ReturnType=typename ReduceOperator::return_type>
    ReturnType reduce(const small_vector_interface<Storage>& a, 
                      const ReduceOperator& reduce=ReduceOperator(),
                      const ReturnType& init=ReturnType(0))
    {
        ReturnType r = init;
        std::for_each(a.begin(), a.end(), [&](auto v) { reduce(r, v); });
        return r;
    }
    
    template<typename Storage>
    typename Storage::value_type product(const small_vector_interface<Storage>& v)
    {
        typedef typename Storage::value_type value_type;
        return reduce(v, [&](value_type& v0, value_type v) -> value_type { 
            v0 *= v; return v0; 
        }, value_type(1));
    }
    
    template<typename Storage>
    typename Storage::value_type mean(const small_vector_interface<Storage>& v)
    {
        typedef typename Storage::value_type value_type;
        return reduce(v, [&](value_type& v0, value_type v) -> value_type { 
            v0 += v; return v0; 
        }, value_type(0)) / static_cast<value_type>(v.size());
    }

    // broadcasters
    template<typename Storage, typename UnaryOperator>
    small_vector_interface<Storage>& apply(small_vector_interface<Storage>& a, const UnaryOperator& op=UnaryOperator())
    {
        for (typename Storage::value_type& val : a)
        {
            val = op(val);
        }
        return a;
    }

    template<typename Storage, typename UnaryOperator>
    small_vector_interface<typename internal::output_storage<Storage>::type> 
    map(const small_vector_interface<Storage>& a, const UnaryOperator& op=UnaryOperator())
    {
        small_vector_interface<Storage> r = a;
        for (typename Storage::value_type& val : r)
        {
            val = op(val);
        }
        return r;
    }

    // unary mappings
    template<typename Storage>
    small_vector_interface<typename internal::output_storage<Storage>::type> 
    square(const small_vector_interface<Storage>& a)
    {
        return map(a, [&](auto v) { return v*v; });
    }
    
    template<typename Storage>
    small_vector_interface<typename internal::output_storage<Storage>::type>
    floor(const small_vector_interface<Storage>& a)
    {
        return map(a, [&](auto v) { return std::floor(v); });
    }
    
    template<typename Storage>
    small_vector_interface<typename internal::output_storage<Storage>::type>
    ceil(const small_vector_interface<Storage>& a)
    {
        return map(a, [&](auto v) { return std::ceil(v); });
    }

    template<typename Storage>
    small_vector_interface<typename internal::output_storage<Storage>::type> 
    abs(const small_vector_interface<Storage>& a) 
    {
        return map(a, [&](auto v) { return std::abs(v); });
    }

    template<typename Storage>
    small_vector_interface<typename internal::output_storage<Storage>::type> 
    round(const small_vector_interface<Storage>& a) 
    {
        return map(a, [&](auto v) { return std::round(v); });
    }
    
    template<typename Storage>
    small_vector<bool, Storage::size> isinf(const small_vector_interface<Storage>& a)
    {
        return map(a, [&](auto v) { return std::isinf(v); });
    }
    
    template<typename Storage>
    small_vector<bool, Storage::size> isnan(const small_vector_interface<Storage>& a)
    {
        return map(a, [&](auto v) { return std::isnan(v); });
    }
    
    template<typename Storage>
    small_vector<bool, Storage::size> isinvalid(const small_vector_interface<Storage>& a)
    {
        return map(a, [&](auto v) { return std::isnan(v) || std::isinf(v); });
    }
    
    template<typename Storage, 
             typename = typename std::enable_if<internal::is_complex<typename Storage::value_type>::value>::type>    
    small_vector_interface<typename internal::real_storage<Storage>::type> 
    real(const small_vector_interface<Storage>& v) 
    {
        small_vector_interface<typename internal::real_storage<Storage>::type> r;
        for (size_t i=0; i<Storage::size; ++i) r[i] = v[i].real();
        return r;
    }
    
    template<typename Storage,
             typename = typename std::enable_if<internal::is_complex<typename Storage::value_type>::value>::type>
    small_vector_interface<typename internal::real_storage<Storage>::type> 
    imaginary(const small_vector_interface<Storage>& v) 
    {
        small_vector_interface<typename internal::real_storage<Storage>::type> r;
        for (size_t i=0; i<Storage::size; ++i) r[i] = v[i].imag();
        return r;
    }

    template<typename Storage>
    std::ostream& operator<<(std::ostream& os, const small_vector_interface<Storage>& a)
    {
        os << "[";
        for (size_t i=0; i<Storage::size-1; ++i) {
            os << a[i] << ", ";
        }
        os << a[Storage::size-1] << "]";
        return os;
    }

    // convenience typedefs
    typedef small_vector<double, 1> vec1;
    typedef small_vector<double, 2> vec2;
    typedef small_vector<double, 3> vec3;
    typedef small_vector<double, 4> vec4;
    typedef small_vector<double, 5> vec5;
    typedef small_vector<double, 6> vec6;

    typedef small_vector<double, 1> dvec1;
    typedef small_vector<double, 2> dvec2;
    typedef small_vector<double, 3> dvec3;
    typedef small_vector<double, 4> dvec4;
    typedef small_vector<double, 5> dvec5;
    typedef small_vector<double, 6> dvec6;

    typedef small_vector<float, 1> fvec1;
    typedef small_vector<float, 2> fvec2;
    typedef small_vector<float, 3> fvec3;
    typedef small_vector<float, 4> fvec4;
    typedef small_vector<float, 5> fvec5;
    typedef small_vector<float, 6> fvec6;
    
    typedef small_vector< std::complex<double>, 1> cvec1;
    typedef small_vector< std::complex<double>, 2> cvec2;
    typedef small_vector< std::complex<double>, 3> cvec3;
    typedef small_vector< std::complex<double>, 4> cvec4;
    typedef small_vector< std::complex<double>, 5> cvec5;
    typedef small_vector< std::complex<double>, 5> cvec6;
    
    typedef small_vector<std::complex<double>, 1> dcvec1;
    typedef small_vector<std::complex<double>, 2> dcvec2;
    typedef small_vector<std::complex<double>, 3> dcvec3;
    typedef small_vector<std::complex<double>, 4> dcvec4;
    typedef small_vector<std::complex<double>, 5> dcvec5;
    typedef small_vector<std::complex<double>, 5> dcvec6;
    
    typedef small_vector<std::complex<float>, 1> fcvec1;
    typedef small_vector<std::complex<float>, 2> fcvec2;
    typedef small_vector<std::complex<float>, 3> fcvec3;
    typedef small_vector<std::complex<float>, 4> fcvec4;
    typedef small_vector<std::complex<float>, 5> fcvec5;
    typedef small_vector<std::complex<float>, 5> fcvec6;

    typedef small_vector<int, 1> ivec1;
    typedef small_vector<int, 2> ivec2;
    typedef small_vector<int, 3> ivec3;
    typedef small_vector<int, 4> ivec4;
    typedef small_vector<int, 5> ivec5;
    typedef small_vector<int, 6> ivec6;

    typedef small_vector<long int, 1> lvec1;
    typedef small_vector<long int, 2> lvec2;
    typedef small_vector<long int, 3> lvec3;
    typedef small_vector<long int, 4> lvec4;
    typedef small_vector<long int, 5> lvec5;
    typedef small_vector<long int, 6> lvec6;

    typedef small_vector<unsigned int, 1> uvec1;
    typedef small_vector<unsigned int, 2> uvec2;
    typedef small_vector<unsigned int, 3> uvec3;
    typedef small_vector<unsigned int, 4> uvec4;
    typedef small_vector<unsigned int, 5> uvec5;
    typedef small_vector<unsigned int, 6> uvec6;

    typedef small_vector<size_t, 1> svec1;
    typedef small_vector<size_t, 2> svec2;
    typedef small_vector<size_t, 3> svec3;
    typedef small_vector<size_t, 4> svec4;
    typedef small_vector<size_t, 5> svec5;
    typedef small_vector<size_t, 6> svec6;

    struct lexicographical_order
    {
        template <typename Storage1, typename Storage2, 
                  typename = typename std::enable_if<internal::same_size<Storage1, Storage2>::value>::type>
        bool operator()(const small_vector_interface<Storage1> &v1,
                        const small_vector_interface<Storage2> &v2) const
        {
            for (auto i=0; i<Storage1::size; ++i) 
            {
                if (v1[i] < v2[i])
                    return true;
                else if (v1[i] > v2[i])
                    return false;
            }
            return false;
        }

        template <typename Storage1, typename Storage2, 
                  typename = typename std::enable_if<Storage1::size == Storage2::size>::type >
        static bool equal(const small_vector_interface<Storage1> &v1, 
                          const small_vector_interface<Storage2> &v2)
        {
            for (auto i = 0; i < Storage1::size; ++i)
            {
                if (v1[i] != v2[i])
                    return false;
            }
            return true;
        }

        template <typename Storage>
        static bool is_zero(const small_vector_interface<Storage> &v)
        {
            for (auto i = 0; i < Storage::size; ++i)
            {
                if (v[i] != 0)
                    return false;
            }
            return true;
        }
    };

    struct eps_lexicographical_order
    {
        eps_lexicographical_order(double eps) : m_eps(eps) {}

        template <typename Storage1, typename Storage2,
                  typename = typename std::enable_if<Storage1::size == Storage2::size>::type>
        bool operator()(const small_vector_interface<Storage1> &v1,
                        const small_vector_interface<Storage2> &v2) const
        {
            for (auto i = 0; i < Storage1::size; ++i)
            {
                if (v1[i] < v2[i] - m_eps)
                    return true;
                else if (v1[i] > v2[i] + m_eps)
                    return false;
            }
            return false;
        }

        template <typename Storage1, typename Storage2,
                  typename = typename std::enable_if<Storage1::size == Storage2::size>::type>
        bool equal(const small_vector_interface<Storage1> &v0, 
                   const small_vector_interface<Storage2> &v1) const
        {
            for (size_t i = 0; i < Storage1::size; ++i)
            {
                if (fabs(v0[i] - v1[i]) > m_eps)
                    return false;
            }
            return true;
        }

        template <typename Storage>
        bool equal_zero(const small_vector_interface<Storage> &v) const
        {
            for (int i = 0; i < Storage::size; ++i)
            {
                if (fabs(v[i]) > m_eps)
                    return false;
            }
            return true;
        }

        double m_eps;
    };

} // namespace spurt