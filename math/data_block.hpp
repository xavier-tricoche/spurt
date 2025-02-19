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
#include <Eigen/Eigen>

namespace spurt 
{
    
    namespace internal
    {
        struct trivial_test 
        {
            bool operator()(bool a) { return a; }
        };

        template<typename T, typename Enable=void>
        struct must_be_scalar {};

        template <typename T >
        struct must_be_scalar<T, typename std::enable_if<std::is_scalar<T>::value>::type> {
            typedef void type;
        };
    
        template<typename T1, typename T2=void, typename T3=void, typename T4=void, typename T5=void>
        struct must_all_be_scalar {};
    
        template<typename T1>
        struct must_all_be_scalar<T1, typename std::enable_if<std::is_scalar<T1>::value>::type> 
        {
            typedef void type;
        };
    
        template<typename T1, typename T2>
        struct must_all_be_scalar<T1, T2, typename std::enable_if<std::is_scalar<T1>::value &&
                                                                  std::is_scalar<T2>::value>::type> 
        {
            typedef void type;
        };
    
        template<typename T1, typename T2, typename T3>
        struct must_all_be_scalar<T1, T2, T3, typename std::enable_if<std::is_scalar<T1>::value &&
                                                                      std::is_scalar<T2>::value &&
                                                                      std::is_scalar<T3>::value>::type > 
        {
            typedef void type;
        };
    
        template<typename T1, typename T2, typename T3, typename T4>
        struct must_all_be_scalar<T1, T2, T3, T4, typename std::enable_if<std::is_scalar<T1>::value &&
                                                                          std::is_scalar<T2>::value &&
                                                                          std::is_scalar<T3>::value &&
                                                                          std::is_scalar<T4>::value>::type > 
        {
            typedef void type;
        };
    
        template<typename Value_, typename Pointer_ = Value_*, typename Reference_ = Value_&>
        class BasicIterator
        {
        public:
            typedef BasicIterator<Value_, Pointer_, Reference_> self_type;
            typedef Value_ value_type;
            typedef Pointer_ pointer;
            typedef Reference_ reference;
            typedef long difference_type;
            typedef std::random_access_iterator_tag iterator_category; 
        
            explicit BasicIterator(pointer init=nullptr) : m_ptr(init) {}
            self_type& operator++() { ++m_ptr; return *this; }
            self_type operator++(int) { self_type r=*this; ++m_ptr; return r; }
            self_type& operator+=(difference_type n) { m_ptr += n; return *this; }
            self_type operator+(difference_type n) const { self_type r=*this; r+=n; return r; }
            self_type& operator--() { --m_ptr; return *this; }
            self_type operator--(int) { self_type r=*this; --m_ptr; return r; }
            self_type& operator-=(difference_type n) { m_ptr -= n; return *this; }
            self_type operator-(difference_type n) const { self_type r = *this; r-=n; return r; }
            difference_type operator-(self_type other) const { return m_ptr - other.m_ptr; }
            bool operator==(self_type other) const { return m_ptr == other.m_ptr; }
            bool operator!=(self_type other) const { return m_ptr != other.m_ptr; }
            bool operator<(self_type other) const { return m_ptr < other.m_ptr; }
            bool operator<=(self_type other) const { return m_ptr <= other.m_ptr; }
            bool operator>(self_type other) const { return m_ptr > other.m_ptr; }
            bool operator>=(self_type other) const { return m_ptr >= other.m_ptr; }
            reference operator*() const { return *m_ptr; }
        
            pointer m_ptr; 
        };
    
        template< typename Value_, typename Pointer_, typename Reference_, size_t Stride_=1>
        class StridedIterator : public BasicIterator<Value_, Pointer_, Reference_>
        {
        public:
            typedef BasicIterator<Value_, Pointer_, Reference_> base_type;
            static constexpr size_t stride = Stride_;
            using typename base_type::value_type;
            using typename base_type::pointer;
            using typename base_type::reference;
            using typename base_type::difference_type;
            using typename base_type::iterator_category;
            typedef StridedIterator<Value_, Pointer_, Reference_, Stride_> self_type;
    
            explicit StridedIterator(pointer init=nullptr)
                : base_type(init) {}
        
            self_type& operator++() { base_type::m_ptr += stride; return *this; }
            self_type operator++(int) { self_type r(base_type::m_ptr); ++(*this); return r; }
            self_type& operator+=(difference_type n) { base_type::m_ptr += n*stride; return *this; }
            self_type operator+(difference_type n) const { self_type r=*this; r+=n; return r; }
            self_type& operator--() { base_type::m_ptr -= stride; return *this; }
            self_type operator--(int) { self_type r(base_type::m_ptr); --(*this); return r; }
            self_type& operator-=(difference_type n) { base_type::m_ptr -= n*stride; return *this; }
            self_type operator-(difference_type n) const { self_type r = *this; r-=n; return r; }
            difference_type operator-(self_type other) const { return (base_type::m_ptr - static_cast<base_type>(other).m_ptr)/stride; }
            bool operator==(self_type other) const { return base_type::m_ptr == static_cast<base_type>(other).m_ptr; }
            bool operator!=(self_type other) const { return !(*this == other); }
            bool operator<(self_type other) const { return base_type::m_ptr < static_cast<base_type>(other).m_ptr; }
            bool operator<=(self_type other) const { return base_type::m_ptr <= static_cast<base_type>(other).m_ptr; }
            bool operator>(self_type other) const { return base_type::m_ptr > static_cast<base_type>(other).m_ptr; }
            bool operator>=(self_type other) const { return base_type::m_ptr >= static_cast<base_type>(other).m_ptr; }
            reference operator*() const { return static_cast<reference>(*base_type::m_ptr); }
        };
    
        template<typename Value_, typename Pointer_, typename Reference_, size_t M, size_t N>
        class OuterIterator : public BasicIterator<Value_, Pointer_, Reference_>
        {
        public:
            typedef BasicIterator<Value_, Pointer_, Reference_> base_type;
            static constexpr long nrows = M;
            static constexpr long ncols = N;
            using typename base_type::value_type;
            using typename base_type::pointer;
            using typename base_type::reference;
            using typename base_type::difference_type;
            using typename base_type::iterator_category;
            typedef OuterIterator<Value_, Pointer_, Reference_, M, N> self_type;
        private:
            difference_type m_i, m_j;
            difference_type fwd(difference_type n) {
                difference_type q = n/ncols;
                difference_type r = n%ncols;
                difference_type i = m_i + q;
                difference_type j = m_j + r;
                if (j >= ncols) {
                    j = j%ncols;
                    i++;
                }
                difference_type delta = (j-m_j)*nrows + (i-m_i);
                m_i = i;
                m_j = j;
                return delta;
            }
            difference_type bwd(difference_type n) {
                difference_type q = n/ncols;
                difference_type r = n%ncols;
                difference_type i = m_i - q;
                difference_type j = m_j - r;
                if (j < 0) {
                    j = j % ncols;
                    i--;
                }
                difference_type delta = (j-m_j)*nrows + (i-m_i);
                m_i = i;
                m_j = j;
                return delta;
            }
        public:
            explicit OuterIterator(pointer init=nullptr, difference_type i=0, difference_type j=0) 
                : base_type(init), m_i(i), m_j(j) {}
        
            self_type operator++() {
                base_type::m_ptr += fwd(1);
                return *this;
            }
            self_type operator++(int) {
                self_type r(base_type::m_ptr, m_i, m_j);
                ++(*this);
                return r;
            }
            self_type operator--() {
                base_type::m_ptr += bwd(1);
            }
            self_type operator--(int) {
                self_type r = *this; --r; return r;
            }
            self_type& operator+=(difference_type n) { 
                if (n>0) 
                    base_type::m_ptr += fwd(n); 
                else 
                    base_type::m_ptr += bwd(-n); 
                return *this;
            }
            self_type operator+(difference_type n) const { self_type r=*this; r+=n; return r; }
            self_type& operator-=(difference_type n) {
                if (n>0)
                    base_type::m_ptr += bwd(n);
                else
                    base_type::m_ptr += fwd(-n);
                return *this;
            }
            self_type operator-(long n) const { self_type r = *this; r-=n; return r; }
            difference_type operator-(self_type other) const { return (m_j-other.m_j)*nrows + (m_i-other.m_i); }
            bool operator==(self_type other) const { return base_type::operator==(static_cast<base_type>(other)); }
            bool operator!=(self_type other) const { return !(*this == other); }
            bool operator<(self_type other) const { return base_type::m_ptr < static_cast<base_type>(other).m_ptr; }
            bool operator<=(self_type other) const { return base_type::m_ptr <= static_cast<base_type>(other).m_ptr; }
            bool operator>(self_type other) const { return base_type::m_ptr > static_cast<base_type>(other).m_ptr; }
            bool operator>=(self_type other) const { return base_type::m_ptr >= static_cast<base_type>(other).m_ptr; }
            reference operator*() const { return static_cast<reference>(*base_type::m_ptr); }
        };
    } // namespace internal
    
    // Basic data storage model for 1D and 2D contiguous data (aka glorified C-array)
    template<typename T, size_t Size_>
    class DataBlock
    {
    public:
        typedef T value_type;
        static constexpr size_t size = Size_;
        typedef internal::BasicIterator<T> iterator;
        typedef internal::BasicIterator<T, const T*, const T&> const_iterator;
        
        // interface to matching Eigen library objects
        typedef Eigen::Vector<T, Size_> eigen_vector_type;
        typedef Eigen::Map<eigen_vector_type> eigen_vector_map_type;
        typedef Eigen::Map<const eigen_vector_type> const_eigen_vector_map_type;
        template<size_t M, size_t N, typename = typename std::enable_if<M*N==size>::type>
        using eigen_matrix_type = Eigen::Matrix<T, M, N>;
        template<size_t M, size_t N, typename = typename std::enable_if<M*N==size>::type>
        using eigen_matrix_map_type = Eigen::Map<eigen_matrix_type<M, N>>;
        template<size_t M, size_t N, typename = typename std::enable_if<M*N==size>::type>
        using const_eigen_matrix_map_type = Eigen::Map<const eigen_matrix_type<M,N>>;
        typedef DataBlock<T, Size_> self_type;
        
        DataBlock() {}
        
        DataBlock(const self_type& other) {
            std::copy(other.begin(), other.end(), begin());
        }

        template<typename OtherStorage, 
                 typename = typename std::enable_if<size == OtherStorage::size>::type>
        DataBlock(const OtherStorage& other) {
            std::copy(other.begin(), other.end(), begin());
        }
        
        iterator begin() { return iterator(data); }
        iterator end() { return iterator(data + size); }
        const_iterator begin() const { return const_iterator(data); }
        const_iterator end() const { return const_iterator(data + size); }
        
        value_type& operator[](size_t i) { assert(i<size); return data[i]; }
        const value_type& operator[](size_t i) const { assert(i<size); return data[i]; }
        
        eigen_vector_map_type as_eigen() { return eigen_vector_map_type(data); }
        const_eigen_vector_map_type as_const_eigen() const { return const_eigen_vector_map_type(data); }
        template<size_t M, size_t N, typename = typename std::enable_if<M*N==size>::type>
        eigen_matrix_map_type<M, N> as_eigen_matrix() { return eigen_matrix_map_type<M,N>(data); }
        template<size_t M, size_t N, typename = typename std::enable_if<M*N==size>::type>
        const_eigen_matrix_map_type<M, M> as_const_eigen_matrix() { return const_eigen_matrix_map_type<M,N>(data); }
        
        value_type data[size];
    };
    
    // 1D strided view of a C-array
    template<typename T, size_t Size_, size_t Stride_>
    class DataBlockView
    {
    public:
        typedef T value_type;
        static constexpr size_t size = Size_;
        static constexpr size_t stride = Stride_;
        typedef internal::StridedIterator<T, T*, T&, stride> iterator;
        typedef internal::StridedIterator<T, const T*, const T&, stride> const_iterator;
        
        // interface to matching Eigen library objects
        typedef Eigen::Vector<T, Size_> eigen_vector_type;
        typedef Eigen::Map<eigen_vector_type, Eigen::Unaligned, Eigen::InnerStride<Stride_> > eigen_vector_map_type;
        typedef Eigen::Map<const eigen_vector_type, Eigen::Unaligned, Eigen::InnerStride<Stride_> > const_eigen_vector_map_type;
        // these should not be needed...
        template<size_t M, size_t N, typename = typename std::enable_if<M*N == size>::type>
        using eigen_matrix_type = Eigen::Matrix<T, M, N>;
        template<size_t M, size_t N, typename = typename std::enable_if<M*N == size>::type>
        using eigen_matrix_map_type = Eigen::Map<Eigen::Matrix<T, M, N>>;
        
        typedef DataBlockView<T,Size_,Stride_> self_type;
    
        iterator begin() { 
            return iterator(data); 
        }
        iterator end() { return iterator(data + stride*size); }
        const_iterator begin() const { 
            return const_iterator(data); 
        }
        const_iterator end() const { return const_iterator(data + stride*size); }
            
        DataBlockView(value_type* _data=nullptr) : data(_data) {}
        DataBlockView(const self_type& other) : data(other.data) {}

        template<typename OtherStorage, 
                 typename = typename std::enable_if<size == OtherStorage::size &&
                 !std::is_const<OtherStorage>::value>::type>
        DataBlockView(OtherStorage& other) : data(other.data) {}
        
        template<typename OtherStorage, 
                 typename = typename std::enable_if<size == OtherStorage::size>::type>
        self_type& operator=(const OtherStorage& other) {
            if (data == nullptr) {
                throw std::runtime_error("Illegal lvalue reference to nonallocated memory");
            }
            std::copy(other.begin(), other.end(), begin());
            return *this;
        }
        
        value_type& operator[](size_t i) { assert(i<size); return *(data + i*stride); }
        const value_type& operator[](size_t i) const { assert(i<size); return *(data + i*stride); }
        
        eigen_vector_map_type as_eigen() { return eigen_vector_map_type(data); }
        const_eigen_vector_map_type as_const_eigen() const { return const_eigen_vector_map_type(data); }
        
        value_type* data;
    };
    
    // 1D strided view of a const C-array
    template<typename T, size_t Size_, size_t Stride_>
    class ConstDataBlockView
    {
    public:
        typedef T value_type;
        static constexpr size_t size = Size_;
        static constexpr size_t stride = Stride_;
        typedef internal::StridedIterator<T, const T*, const T&, stride> iterator;
        typedef internal::StridedIterator<T, const T*, const T&, stride> const_iterator;
        
        // interface to matching Eigen library objects
        typedef const Eigen::Vector<T, Size_> eigen_vector_type;
        typedef Eigen::Map<const Eigen::Vector<T, Size_>, Eigen::Unaligned, Eigen::InnerStride<Stride_> > eigen_vector_map_type;
        typedef const Eigen::Vector<T, Size_> const_eigen_vector_type;
        typedef Eigen::Map<const Eigen::Vector<T, Size_>, Eigen::Unaligned, Eigen::InnerStride<Stride_> > const_eigen_vector_map_type;
        typedef ConstDataBlockView<T, Size_, Stride_> self_type;
        typedef DataBlockView<T, Size_, Stride_> unconst_self_type;
        // these should not be needed...
        template<size_t M, size_t N, typename = typename std::enable_if<M*N == size>::type>
        using eigen_matrix_type = Eigen::Matrix<T, M, N>;
        template<size_t M, size_t N, typename = typename std::enable_if<M*N == size>::type>
        using eigen_matrix_map_type = Eigen::Map<Eigen::Matrix<T, M, N>>;
        
        const_iterator begin() const { return const_iterator(data); }
        const_iterator end() const { return const_iterator(data + stride*size); }
            
        ConstDataBlockView(const value_type* _data=nullptr) : data(_data) {}
        ConstDataBlockView(const self_type& other) : data(other.data) {}
        ConstDataBlockView(unconst_self_type other): data(const_cast<const value_type*>(other.data)) {}
        
        const value_type& operator[](size_t i) const { return *(data + i*stride); }
        const_eigen_vector_map_type as_const_eigen() const { return const_eigen_vector_map_type(data); }
        
        const value_type* data;
    };
    
    namespace internal
    {
    
        template<typename Storage1, typename Storage2, typename Enable = void> struct  better_storage {};
    
        template<typename Storage1, typename Storage2>
        struct  better_storage<Storage1, Storage2, typename std::enable_if<Storage1::size==Storage2::size>::type >
        {
            typedef typename spurt::better_type<typename Storage1::value_type, typename Storage2::value_type>::type value_type;
            typedef DataBlock<value_type, Storage1::size> type;
        };
    
        template<typename Storage, typename T>
        struct  better_storage<Storage, T, typename must_be_scalar<T>::type>
        {
            typedef typename spurt::better_type<typename Storage::value_type, T>::type value_type;
            typedef DataBlock<value_type, Storage::size> type;
        };
    
        template<typename Storage1, typename Storage2, typename Enable = void>
        struct same_size {};
    
        template<typename Storage1, typename Storage2>
        struct same_size<Storage1, Storage2, 
                         typename std::enable_if<Storage1::size == Storage2::size>::type >
        {
            typedef void type;
            static constexpr bool value = true;
        };
    
        template<typename Storage1, typename Storage2, typename Enable = void>
        struct must_be_3d {};
    
        template<typename Storage1, typename Storage2>
        struct must_be_3d<Storage1, Storage2, 
                          typename std::enable_if<Storage1::size == 3 && Storage2::size == 3>::type>
        {
            typedef void type;
            static constexpr bool value = true;
        };
    
        template<typename T>
        struct is_complex : public std::false_type {};
    
        template<typename T>
        struct is_complex<std::complex<T>> : public std::true_type {};
    
        template<typename T>
        struct is_complex<const std::complex<T>> : public std::true_type {};
    
        template<typename Storage>
        struct output_storage {};
    
        template<typename T, size_t N>
        struct output_storage<DataBlock<T,N>> {
            typedef DataBlock<T, N> type;
        };
    
        template<typename T, size_t N, size_t Stride>
        struct output_storage<DataBlockView<T,N,Stride>> {
            typedef DataBlock<T, N> type;
        };
    
        template<typename T, size_t N, size_t Stride>
        struct output_storage<ConstDataBlockView<T,N,Stride>> {
            typedef DataBlock<T, N> type;
        };
    
        template<typename Storage, typename Enable = void>
        struct real_storage {};
    
        template<typename T, size_t N>
        struct real_storage<DataBlock<T, N>, 
                            typename std::enable_if<!is_complex<T>::value>::type>
        {
            typedef T value_type;
            typedef DataBlock<T, N> type;
        };
    
        template<typename T, size_t N>
        struct real_storage<DataBlock<T, N>, 
                            typename std::enable_if<is_complex<T>::value>::type>
        {
            typedef typename T::value_type value_type;
            typedef DataBlock<value_type, N> type;
        };
    
        template<typename T, size_t N, size_t Stride>
        struct real_storage<DataBlockView<T, N, Stride>,
                            typename std::enable_if<is_complex<T>::value>::type>
        {
            typedef typename T::value_type value_type;
            typedef DataBlock<value_type, N> type;
        };
    
        template<typename T, size_t N, size_t Stride>
        struct real_storage<ConstDataBlockView<T, N, Stride>,
                            typename std::enable_if<is_complex<T>::value>::type>
        {
            typedef typename T::value_type value_type;
            typedef DataBlock<value_type, N> type;
        };
    
        template<typename T, typename Enable = void>
        struct random_filler {};
    
        template<typename T>
        struct random_filler<T, typename std::enable_if<std::is_integral<T>::value>::type>
        {
            typedef T value_type;
            typedef std::uniform_int_distribution<value_type> dist_type;
            static constexpr value_type default_min = std::numeric_limits<value_type>::min();
            static constexpr value_type default_max = std::numeric_limits<value_type>::max();
            random_filler(value_type min_=default_min, value_type max_=default_max) {
                m_dist = dist_type(min_, max_);
                std::random_device rd;
                m_gen = std::mt19937(rd());
            }
        
            value_type operator()() {
                return m_dist(m_gen);
            }
        
            template<typename Iterator>
            void fill(Iterator begin, Iterator end) {
                for (Iterator it=begin; it!=end; ++it) {
                    *it = m_dist(m_gen);
                }
            }
        
            dist_type m_dist;
            std::mt19937 m_gen;
        };
    
        template<typename T>
        struct random_filler<T, typename std::enable_if<std::is_floating_point<T>::value>::type>
        {
            typedef T value_type;
            typedef std::uniform_real_distribution<value_type> dist_type;
            static constexpr value_type default_min = value_type(0);
            static constexpr value_type default_max = value_type(1);
            random_filler(value_type min_=default_min, value_type max_=default_max) {
                m_dist = dist_type(min_, max_);
                std::random_device rd;
                m_gen = std::mt19937(rd());
            }
        
            value_type operator()() {
                return m_dist(m_gen);
            }
        
            template<typename Iterator>
            void fill(Iterator begin, Iterator end) {
                for (Iterator it=begin; it!=end; ++it) {
                    *it = m_dist(m_gen);
                }
            }
        
            value_type m_min, m_max;
            dist_type m_dist;
            std::mt19937 m_gen;
        };
    
        template<typename T>
        struct random_filler<std::complex<T>, typename std::enable_if<std::is_integral<T>::value>::type>
        {
            typedef std::complex<T> value_type;
            typedef T scalar_type;
            typedef random_filler<T> base_type;
            static constexpr value_type default_min = value_type(base_type::default_min, base_type::default_min);
            static constexpr value_type default_max = value_type(base_type::default_max, base_type::default_max);
            random_filler(value_type min_ = default_min, value_type max_ = default_max) 
                 : m_real_filler(min_.real(), max_.real()), m_imag_filler(min_.imag(), max_.imag()) {}
        
            value_type operator()() {
                return value_type(m_real_filler(), m_imag_filler());
            }
        
            template<typename Iterator>
            void fill(Iterator begin, Iterator end) {
                for (Iterator it=begin; it!=end; ++it) {
                    *it = (*this)();
                }
            }
        
            base_type m_real_filler, m_imag_filler;
        };
    
        template<typename T>
        struct random_filler<std::complex<T>, typename std::enable_if<std::is_floating_point<T>::value>::type>
        {
            typedef std::complex<T> value_type;
            typedef T scalar_type;
            typedef random_filler<T> base_type;
            static constexpr value_type default_min = value_type(base_type::default_min, base_type::default_min);
            static constexpr value_type default_max = value_type(base_type::default_max, base_type::default_max);
            random_filler(value_type min_ = default_min, value_type max_ = default_max) 
                 : m_real_filler(min_.real(), max_.real()), m_imag_filler(min_.imag(), max_.imag()) {}
        
            value_type operator()() {
                return value_type(m_real_filler(), m_imag_filler());
            }
        
            template<typename Iterator>
            void fill(Iterator begin, Iterator end) {
                for (Iterator it=begin; it!=end; ++it) {
                    *it = (*this)();
                }
            }
        
            base_type m_real_filler, m_imag_filler;
        };
    
    } // namespace internal

} // namespace spurt