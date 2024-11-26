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

namespace 
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
    
    // Basic data storage model for 1D and 2D contiguous data (aka glorified C-array)
    template<typename T, size_t Size_>
    class DataBlock
    {
    public:
        typedef T value_type;
        static constexpr size_t size = Size_;
        typedef BasicIterator<T> iterator;
        typedef BasicIterator<T, const T*, const T&> const_iterator;
        
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
        
        iterator begin() { return iterator(data); }
        iterator end() { return iterator(data + size); }
        const_iterator begin() const { return const_iterator(data); }
        const_iterator end() const { return const_iterator(data + size); }
        
        value_type& operator[](size_t i) { assert(i<size); return data[i]; }
        const value_type& operator[](size_t i) const { assert(i<size); return data[i]; }
        
        eigen_vector_map_type as_eigen() { return eigen_vector_map_type(data); }
        const_eigen_vector_map_type as_const_eigen() const { return const_eigen_vector_map_type(data); }
        template<size_t M, size_t N, typename = typename std::enable_if<M*N==size>::type>
        eigen_matrix_map_type<M, N> as_eigen_matrix() { return eigen_matrix_map_type(data); }
        template<size_t M, size_t N, typename = typename std::enable_if<M*N==size>::type>
        const_eigen_matrix_map_type<M, M> as_const_eigen_matrix() { return const_eigen_matrix_map_type(data); }
        
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
        typedef StridedIterator<T, T*, T&, stride> iterator;
        typedef StridedIterator<T, const T*, const T&, stride> const_iterator;
        
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
                 typename = typename std::enable_if<size == OtherStorage::size>::type>
        self_type& operator=(const OtherStorage& other) {
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
        typedef StridedIterator<T, const T*, const T&, stride> iterator;
        typedef StridedIterator<T, const T*, const T&, stride> const_iterator;
        
        // interface to matching Eigen library objects
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
    
}

namespace spurt {

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
        
        typedef random_filler<value_type> random_filler_type;
        
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
        template<typename T, typename = typename must_be_scalar<T>::type>
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

        template<typename T, typename = typename must_be_scalar<T>::type >
        small_vector_interface(std::initializer_list<T> vals)
            : m_storage() 
        {
            std::copy(vals.begin(), vals.end(), begin());
        }

        template<typename T1, typename T2, typename = typename must_all_be_scalar<T1, T2>::type>
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
                  typename = typename must_all_be_scalar<T1, T2, T3>::type>
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
                typename = typename must_all_be_scalar<T1, T2, T3, T4>::type>
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

        template<typename T, typename = typename must_be_scalar<T>::type >
        bool_vector operator<(T val) const 
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { r[i] = ((*this)[i] < val); }
            return r;
        }

        template<typename T, typename = typename must_be_scalar<T>::type >
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

        template <typename T, typename = typename must_be_scalar<T>::type >
        bool_vector operator>(T val) const
        {
            bool_vector r;
            for (size_t i=0; i<_size_; ++i) { r[i] = (m_storage[i] > val); }
            return r;
        }

        template <typename T, typename = typename must_be_scalar<T>::type >
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

        template <typename T, typename = typename must_be_scalar<T>::type>
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

        template <typename T, typename = typename must_be_scalar<T>::type>
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
        
        template <typename T, typename = typename must_be_scalar<T>::type >
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
        
        template <typename T, typename = typename must_be_scalar<T>::type>
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
        
        template <typename T, typename = typename must_be_scalar<T>::type >
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
        
        template <typename T, typename = typename must_be_scalar<T>::type>
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
        template <typename T, typename = typename must_be_scalar<T>::type >
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
        template <typename T, typename = typename must_be_scalar<T>::type>
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
        template <typename T, typename = typename must_be_scalar<T>::type>
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
        template <typename T, typename = typename must_be_scalar<T>::type>
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

    template <typename Storage1, typename Storage2, 
              typename = typename must_be_3d<Storage1, Storage2>::type>
    small_vector_interface<typename  better_storage<Storage1, Storage2>::type > 
    cross(const small_vector_interface<Storage1> &a, const small_vector_interface<Storage2> &b)
    {
        typedef typename better_storage<Storage1, Storage2>::type out_storage;
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
              typename = typename same_size<Storage1, Storage2>::type>
    typename better_storage<Storage1, Storage2>::value_type
    inner(const small_vector_interface<Storage1> &a, const small_vector_interface<Storage2> &b)
    {
        typedef typename better_storage<Storage1, Storage2>::value_type type;
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

    template <typename T, typename Storage>
    small_vector_interface<typename better_storage<Storage, T>::type>
    operator+(const T &val, const small_vector_interface<Storage> &a)
    {
        return a + val;
    }

    template <typename T, typename Storage>
    small_vector_interface<typename better_storage<Storage, T>::type > 
    operator-(const T &val, const small_vector_interface<Storage> &a)
    {
        return a - val;
    }

    template <typename T, typename Storage>
    small_vector_interface<typename better_storage<Storage, T>::type> 
    operator*(const T &val, const small_vector_interface<Storage> &a)
    {
        return a * val;
    }

    // Reductions
    template <typename Storage>
    bool all(const small_vector_interface<Storage> &a)
    {
        return std::all_of(a.begin(), a.end(), trivial_test());
    }

    template<typename Storage, typename Predicate>
    bool all(const small_vector_interface<Storage>& a, const Predicate& pred=Predicate())
    {
        return std::all_of(a.begin(), a.end(), pred);
    }

    template<typename Storage, typename = typename std::enable_if<std::is_same<typename Storage::value_type, bool>::value>::type>
    bool any(const small_vector_interface<Storage>& a) 
    {
        return std::any_of(a.begin(), a.end(), trivial_test());
    }

    template <typename Storage, typename Predicate>
    bool any(const small_vector_interface<Storage> &a, const Predicate &pred = Predicate())
    {
        return std::any_of(a.begin(), a.end(), pred);
    }

    template <typename Storage, typename = typename std::enable_if<std::is_same<typename Storage::value_type, bool>::value>::type>
    bool none(const small_vector_interface<Storage> &a)
    {
        return std::none_of(a.begin(), a.end(), trivial_test());
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
        return reduce(v, [&](value_type& v0, value_type v) -> value_type { v0 *= v; return v0; }, value_type(1));
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
    small_vector_interface<typename output_storage<Storage>::type> 
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
    small_vector_interface<typename output_storage<Storage>::type> 
    square(const small_vector_interface<Storage>& a)
    {
        return map(a, [&](auto v) { return v*v; });
    }
    
    template<typename Storage>
    small_vector_interface<typename output_storage<Storage>::type>
    floor(const small_vector_interface<Storage>& a)
    {
        return map(a, [&](auto v) { return std::floor(v); });
    }
    
    template<typename Storage>
    small_vector_interface<typename output_storage<Storage>::type>
    ceil(const small_vector_interface<Storage>& a)
    {
        return map(a, [&](auto v) { return std::ceil(v); });
    }

    template<typename Storage>
    small_vector_interface<typename output_storage<Storage>::type> 
    abs(const small_vector_interface<Storage>& a) 
    {
        return map(a, [&](auto v) { return std::abs(v); });
    }

    template<typename Storage>
    small_vector_interface<typename output_storage<Storage>::type> 
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
             typename = typename std::enable_if<is_complex<typename Storage::value_type>::value>::type>    
    small_vector_interface<typename real_storage<Storage>::type> 
    real(const small_vector_interface<Storage>& v) 
    {
        small_vector_interface<typename real_storage<Storage>::type> r;
        for (size_t i=0; i<Storage::size; ++i) r[i] = v[i].real();
        return r;
    }
    
    template<typename Storage,
             typename = typename std::enable_if<is_complex<typename Storage::value_type>::value>::type>
    small_vector_interface<typename real_storage<Storage>::type> 
    imaginary(const small_vector_interface<Storage>& v) 
    {
        small_vector_interface<typename real_storage<Storage>::type> r;
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

    // small_matrix:
    // a matrix interface to a small_vector_interface using Eigen for
    // linear algebra computation. Storage is *column-major* but
    // iterators are available for both column-wise and row-wise traversal
    template<typename T, size_t M, size_t N>
    class small_matrix : public small_vector_interface<DataBlock<T, M*N> >
    {
    public:
        typedef T scalar_type;
        typedef T value_type;
        typedef DataBlock<T, M*N> storage_type;
        typedef small_vector_interface< storage_type > base_type;
        typedef random_filler<value_type> random_filler_type;
        typedef DataBlockView< T, M, 1> column_view_type;
        typedef ConstDataBlockView< T, M, 1> const_column_view_type;
        typedef DataBlockView< T, N, M > row_view_type;
        typedef ConstDataBlockView< T, N, M > const_row_view_type;
        typedef small_vector_interface< column_view_type > column_type;
        typedef small_vector_interface< const_column_view_type > const_column_type;
        typedef small_vector_interface< row_view_type > row_type;
        typedef small_vector_interface< const_row_view_type > const_row_type;
        static constexpr size_t _size_ = M*N;
        static constexpr size_t nrows = M;
        static constexpr size_t ncols = N;
        typedef small_matrix<scalar_type, nrows, ncols> self_type;
        typedef small_matrix<scalar_type, ncols, nrows> self_transpose_type;
        typedef Eigen::Matrix<scalar_type, nrows, ncols> self_matrix_type;
        typedef Eigen::Map<self_matrix_type> self_map_type;
        typedef Eigen::Map<const self_matrix_type> const_self_map_type;
        
        typedef typename base_type::iterator columnwise_iterator;
        typedef typename base_type::const_iterator const_columnwise_iterator;
        
        typedef OuterIterator<scalar_type, scalar_type*, scalar_type&, M, N> rowwise_iterator;
        typedef OuterIterator<scalar_type, const scalar_type*, const scalar_type&, M, N> const_rowwise_iterator;
        
        template<typename Iterator = columnwise_iterator>
        Iterator begin() { return Iterator(base_type::m_storage.data); }
        
        template<typename Iterator = columnwise_iterator>
        Iterator end() { return Iterator(base_type::m_storage.data + _size_); } 
        
        template<typename ConstIterator = const_columnwise_iterator>
        ConstIterator begin() const { return ConstIterator(base_type::m_storage.data); }
        
        template<typename ConstIterator = const_columnwise_iterator>
        ConstIterator end() const { return ConstIterator(base_type::m_storage.data + _size_); }
        
        template<typename T1>
        using right_type = typename better_type<scalar_type, T1>::type;
        
        template<typename T1>
        using right_duplicate = small_matrix<right_type<T1>, nrows, ncols>;
        
        template<typename T1, size_t P>
        using right_output = small_matrix<right_type<T1>, nrows, P>;
        
        template<typename OtherStorage>
        using right_output_vector = small_vector_interface<DataBlock<right_type<typename OtherStorage::value_type>, ncols>>;
        
        template <typename T1, size_t P>
        using rhs_matrix = small_matrix<T1, ncols, P>;
        
        template <typename T1, size_t P>
        using rhs_eigen_matrix = Eigen::Matrix<T1, ncols, P>;
        
        template<typename T1, size_t P>
        using rhs_map_type = Eigen::Map<rhs_eigen_matrix<T1,P>>;
        
        template<typename T1>
        using similar_matrix = small_matrix<T1, nrows, ncols>;
        
        template<typename Storage, 
                 typename = typename std::enable_if<Storage::size == _size_>::type>
        using similar_vector = small_vector_interface<Storage>;
        
        template<typename Storage, 
                 typename = typename std::enable_if<Storage::size == N>::type>
        using rhs_vector = small_vector_interface<Storage>;
                 
        template<typename T1>
        using similar_array = std::array<T1, nrows*ncols>;
        
        small_matrix(scalar_type val=0, bool linalg_mode=true) 
        : base_type(val), m_linalg_mode(linalg_mode) {}

        template<typename OtherStorage>
        small_matrix(const similar_vector<OtherStorage>& other, bool linalg_mode=true)
        : base_type(other), m_linalg_mode(linalg_mode) {}

        template<typename T1=scalar_type>
        small_matrix(const similar_array<T1>& other, bool linalg_mode=true) 
            : base_type(other), m_linalg_mode(linalg_mode) {}
        
        template<typename T1=scalar_type>
        small_matrix(std::initializer_list<T> vals) : base_type(vals), m_linalg_mode(true) {}

        template<typename T1>
        small_matrix(const similar_matrix<T1>& other) 
            : base_type(static_cast<const similar_vector<storage_type>&>(other)), m_linalg_mode(other.m_linalg_mode) {}

        static self_type identity() {
            self_type r;
            r.as_eigen() = self_matrix_type::Identity();
            return r;
        }
        
        static self_type random(value_type min=random_filler_type::default_min,
                                value_type max=random_filler_type::default_max) {
            random_filler_type filler(min, max);
            self_type r;
            filler.fill(r.template begin<columnwise_iterator>(), r.template end<columnwise_iterator>());
            return r;
        }
        
        self_transpose_type transpose() const {
            self_transpose_type r;
            r.as_eigen() = as_const_eigen().transpose();
            return r;
        }

        self_type& linalg(bool active=true) 
        {
            m_linalg_mode = active;
            return *this;
        }

        bool is_linalg() const 
        {
            return m_linalg_mode;
        }

        const scalar_type& operator()(int i, int j) const 
        {
            return static_cast<const base_type&>(*this)[j*nrows+i];
        }
        scalar_type& operator()(int i, int j) 
        {
            return static_cast<base_type&>(*this)[j*nrows+i];
        }

        column_type column(int j)
        {
            return column_type(column_view_type(&((*this)(0,j))));
        }
        
        const_column_type column(int j) const 
        {
            return const_column_type(const_column_view_type(&((*this)(0,j))));
        }

        row_type row(int i)
        {
            return row_type(row_view_type(&((*this)(i,0))));
        }
        
        const_row_type row(int i) const
        {
            return const_row_type(const_row_view_type(&((*this)(i,0))));
        }

        template <typename T1, size_t P>
        right_output<T1, P> operator*(const rhs_matrix<T1, P> &rhs)
        {
            typedef right_type<T1> out_type;
            typedef right_output<T1, P> out_matrix;
            typedef Eigen::Map<Eigen::Matrix<out_type, nrows, P>> out_map_type;
            if (!m_linalg_mode)
                return self_type(base_type::operator*(rhs));
            else
            {
                out_matrix r;
                r.as_eigen() = as_const_eigen().template cast<out_type>() * 
                               rhs.as_const_eigen().template cast<out_type>();
                return r;
            }
        }

        template <typename Storage>
        right_output_vector<Storage>
        operator*(const rhs_vector<Storage> &rhs) const
        {
            typedef right_type<typename Storage::value_type> out_type;
            right_output_vector<Storage> r;
            r.as_eigen() = as_const_eigen().template cast<out_type>() * 
                           rhs.as_const_eigen().template cast<out_type>();
            return r;
        }

        self_map_type as_eigen()
        {
            return self_map_type(&((*this)(0,0)));
        }

        const_self_map_type as_const_eigen() const {
            return const_self_map_type(&((*this)(0,0)));
        }
        
        bool m_linalg_mode;
    };

    template<typename T, size_t M, size_t N>
    std::ostream& operator<<(std::ostream& os, const small_matrix<T, M, N>& m) 
    {
        os << '\n' << m.as_const_eigen() << '\n';
        return os;
    }

    template<typename T, size_t M>
    using small_square_matrix = small_matrix<T, M, M>;
    
    typedef small_square_matrix<double, 2> mat2;
    typedef small_square_matrix<double, 3> mat3;
    typedef small_square_matrix<double, 4> mat4;
    typedef small_square_matrix<double, 5> mat5;
    typedef small_square_matrix<double, 6> mat6;
    
    typedef small_square_matrix<double, 2> dmat2;
    typedef small_square_matrix<double, 3> dmat3;
    typedef small_square_matrix<double, 4> dmat4;
    typedef small_square_matrix<double, 5> dmat5;
    typedef small_square_matrix<double, 6> dmat6;
    
    typedef small_square_matrix<float, 2> fmat2;
    typedef small_square_matrix<float, 3> fmat3;
    typedef small_square_matrix<float, 4> fmat4;
    typedef small_square_matrix<float, 5> fmat5;
    typedef small_square_matrix<float, 6> fmat6;
    
    typedef small_square_matrix<std::complex<double>, 2> cmat2;
    typedef small_square_matrix<std::complex<double>, 3> cmat3;
    typedef small_square_matrix<std::complex<double>, 4> cmat4;
    typedef small_square_matrix<std::complex<double>, 5> cmat5;
    typedef small_square_matrix<std::complex<double>, 6> cmat6;
    
    typedef small_square_matrix<std::complex<double>, 2> dcmat2;
    typedef small_square_matrix<std::complex<double>, 3> dcmat3;
    typedef small_square_matrix<std::complex<double>, 4> dcmat4;
    typedef small_square_matrix<std::complex<double>, 5> dcmat5;
    typedef small_square_matrix<std::complex<double>, 6> dcmat6;
    
    typedef small_square_matrix<std::complex<float>, 2> fcmat2;
    typedef small_square_matrix<std::complex<float>, 3> fcmat3;
    typedef small_square_matrix<std::complex<float>, 4> fcmat4;
    typedef small_square_matrix<std::complex<float>, 5> fcmat5;
    typedef small_square_matrix<std::complex<float>, 6> fcmat6;
    
    typedef small_square_matrix<int, 2> imat2;
    typedef small_square_matrix<int, 3> imat3;
    typedef small_square_matrix<int, 4> imat4;
    typedef small_square_matrix<int, 5> imat5;
    typedef small_square_matrix<int, 6> imat6;
    
    typedef small_square_matrix<long int, 2> lmat2;
    typedef small_square_matrix<long int, 3> lmat3;
    typedef small_square_matrix<long int, 4> lmat4;
    typedef small_square_matrix<long int, 5> lmat5;
    typedef small_square_matrix<long int, 6> lmat6;
    
    typedef small_square_matrix<size_t, 2> smat2;
    typedef small_square_matrix<size_t, 3> smat3;
    typedef small_square_matrix<size_t, 4> smat4;
    typedef small_square_matrix<size_t, 5> smat5;
    typedef small_square_matrix<size_t, 6> smat6;

    template<typename T, size_t M, size_t N>
    small_matrix<T, N, M> transpose(const small_matrix<T, M, N>& m)
    {
        return m.transpose();
    }
    
    template<typename T, size_t N>
    small_matrix<T, N, N> make_symmetric(const small_matrix<T, N, N>& m, bool lower=true)
    {
        small_matrix<T, N, N> r;
        for (size_t i=0; i<N; ++i) {
            r(i,i) = m(i,i);
            for (size_t j=i+1; j<N; ++j) {
                if (lower) {
                    r(i,j) = r(j,i) = m(j,i);
                }
                else {
                    r(i,j) = r(j,i) = m(i,j);
                }
            }
        }
        return r;
    }
    
    template<typename T=double>
    small_square_matrix<T, 3> from_dti(const double* t) 
    {
        small_square_matrix<T, 3> r;
        r(0,0) = t[1];
        r(1,0) = r(0,1) = t[2];
        r(2,0) = r(0,2) = t[3];
        r(1,1) = t[4];
        r(2,1) = r(1,2) = t[5];
        r(2,2) = t[6];
        return r;
    }

    template<typename T, size_t M>
    T trace(const small_square_matrix<T, M>& m)
    {
        return m.as_const_eigen().trace();
    }

    template<typename T, size_t M>
    T determinant(const small_square_matrix<T, M>& m) 
    {
        return m.as_const_eigen().determinant();
    }

    template<typename T, size_t M>
    small_square_matrix<T, M> inverse(const small_square_matrix<T, M>& m) {
        small_square_matrix<T, M> r;
        r.as_eigen() = m.as_const_eigen().inverse();
        return r;
    }
    
    template<typename Storage1, typename Storage2, typename = typename same_size<Storage1, Storage2>::type>
    small_square_matrix<typename better_type<typename Storage1::value_type, typename Storage2::value_type>::type, Storage1::size> 
    outer(const small_vector_interface<Storage1>& v0, const small_vector_interface<Storage2>& v1) {
        typedef typename better_type<typename Storage1::value_type, typename Storage2::value_type>::type out_type;
        small_square_matrix<out_type, Storage1::size> r;
        r.as_eigen() = v0.as_const_eigen().template cast<out_type>() *
                       v1.as_const_eigen().transpose().template cast<out_type>();
        return r;
    }

    // Sort Eigenvalues
    // "Decreasing" sorting order:
    // 1. real eigenvalues in decreasing order
    // 2. complex conjugate eigenvalues in decreasing order of real parts
    // 3. complex conjugate eigenvalues in decreasing order of imaginary parts
    // Example: [ 1, 5, 2+3j, 2-3j, 4-2j, 7, 4+2j ] ->
    //          [ 7, 5, 1, 4+2j, 4-2j, 2+3j, 2-3j ]
    template<typename T1, typename T2, size_t M>
    void sort_complex_eigenvalues(
            small_vector<std::complex<T1>, M>& eigenvalues, 
            small_matrix<std::complex<T2>, M, M>& eigenvectors)
    {
        typedef std::complex<T1> val_t;
        typedef small_vector<val_t, M> vec_t;
        typedef small_matrix<std::complex<T2>, M, M> mat_t;

        std::vector<int> indices(M);
        struct Order {
            bool operator()(const val_t& v0, const val_t& v1) 
            {
                if (v0.imag() == 0 && v1.imag() != 0) return true;
                else if (v0.imag() != 0 && v1.imag() == 0) return false;
                else if (v0.real() > v1.real()) return true;
                else if (v0.real() < v1.real()) return false;
                else return (v0.imag() >= v1.imag());
            }
        };
        spurt::argsort(indices, eigenvalues, Order());
        mat_t new_evecs(eigenvectors);
        vec_t new_evals(eigenvalues);
        for (int i=0; i<M; ++i)
        {
            eigenvalues[i] = new_evals[indices[i]];
            eigenvectors.as_eigen().col(i) = new_evecs.as_const_eigen().col(indices[i]);
        }
    }

    template <typename T1, typename T2, size_t M>
    void sort_real_eigenvalues(
        small_vector<T1, M> &eigenvalues,
        small_matrix<T2, M, M> &eigenvectors)
    {
        typedef T1 val_t;
        typedef small_vector<val_t, M> vec_t;
        typedef small_matrix<T2, M, M> mat_t;

        std::vector<int> indices(M);
        struct Order
        {
            bool operator()(const val_t &v0, const val_t &v1)
            {
                return v0 >= v1;
            }
        };
        spurt::argsort(indices, eigenvalues, Order());
        mat_t new_evecs(eigenvectors);
        vec_t new_evals(eigenvalues);
        for (int i = 0; i < M; ++i)
        {
            eigenvalues[i] = new_evals[indices[i]];
            eigenvectors.as_eigen().col(i) = new_evecs.as_const_eigen().col(indices[i]);
        }
    }

    template<typename T, typename Storage, typename T2, size_t M,
             typename = typename std::enable_if< Storage::size == M &&
                                                 is_complex<typename Storage::value_type>::value &&
                                                 std::is_scalar<T2>::value >::type >
    void eigensystem(small_vector_interface<Storage>& eigenvalues, 
                     small_square_matrix<std::complex<T2>, M>& eigenvectors, 
                     const small_square_matrix<T, M>& m) 
    {
        typedef Eigen::Matrix<T, M, M> eigen_type;
        typedef typename Storage::value_type evals_type;
        
        Eigen::EigenSolver<eigen_type> solver(m.as_const_eigen());
        eigenvalues.as_eigen() = solver.eigenvalues().template cast<evals_type>();
        eigenvectors.as_eigen() = solver.eigenvectors().template cast<std::complex<T2>>();
        sort_complex_eigenvalues(eigenvalues, eigenvectors);
    }

    template<typename T, typename Storage, size_t M,
             typename = typename std::enable_if< Storage::size == M &&
                                                 is_complex<typename Storage::value_type>::value>::type >
    void complex_eigenvalues(small_vector_interface<Storage>& evals,
                             const small_square_matrix<T, M>& m)
    {
        typedef Eigen::Matrix<T, M, M> eigen_type;
        typedef typename Storage::value_type scalar;
        Eigen::EigenSolver<eigen_type> solver(m.as_const_eigen());
        evals.as_eigen() = solver.eigenvalues().template cast<scalar>();
        std::sort(evals.begin(), evals.end(), [&](scalar a, scalar b) {
            if (a.imag() == 0 && b.imag() != 0) return true;
            else if (a.imag() != 0 && b.imag() == 0) return false;
            else if (a.real() > b.real()) return true;
            else if (a.real() < b.real()) return false;
            else return (a.imag() >= b.imag());
        });
    }

    template<typename T, typename Storage, typename T2, size_t M,
             typename = typename std::enable_if<Storage::size == M &&
                                                !is_complex<typename Storage::value_type>::value>::type>
    void sym_eigensystem(small_vector_interface<Storage>& eigenvalues,
                         small_square_matrix<T2, M>& eigenvectors,
                         const small_square_matrix<T, M>& m)
    {
        typedef Eigen::Matrix<T, M, M> eigen_type;
        Eigen::SelfAdjointEigenSolver<eigen_type> solver(m.as_const_eigen());
        eigenvalues.as_eigen() = solver.eigenvalues().template cast<typename Storage::value_type>();
        eigenvectors.as_eigen() = solver.eigenvectors().template cast<T2>();
        sort_real_eigenvalues(eigenvalues, eigenvectors);
    }

    template <typename Storage, typename T, size_t M, 
              typename = typename std::enable_if<Storage::size == M &&
                                                 !is_complex<typename Storage::value_type>::value>::type>
    void real_eigenvalues(small_vector_interface<Storage> &evals,
                          const small_square_matrix<T, M> &m)
    {
        typedef Eigen::Matrix<T, M, M> eigen_type;
        typedef typename Storage::value_type scalar;
        Eigen::SelfAdjointEigenSolver<eigen_type> solver(m.as_const_eigen());
        evals.as_eigen() = solver.eigenvalues().template cast<scalar>();
        std::sort(evals.begin(), evals.end(), 
                  [&](scalar a, scalar b) {
                    return a >= b;
                  });
    }

    struct lexicographical_order
    {
        template <typename Storage1, typename Storage2, 
                  typename = typename std::enable_if<same_size<Storage1, Storage2>::value>::type>
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