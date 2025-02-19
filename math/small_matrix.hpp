#pragma once

#include <math/small_vector.hpp>

namespace spurt {
    
    namespace internal {
        template<typename Storage, size_t M, size_t N, typename Enable=void>
        struct fits {};
    
        template<typename Storage, size_t M, size_t N>
        struct fits<Storage, M, N, typename std::enable_if<Storage::size == M*N>::type>
        {
            static constexpr bool value = true;
            typedef void type;
        };
        
        template<typename Matrix>
        struct is_square {
            static constexpr bool value = Matrix::nrows == Matrix::ncols;
        };
        
        template<typename Matrix, typename Enable=void>
        struct must_be_square {};
        
        template<typename Matrix>
        struct must_be_square<Matrix, typename std::enable_if<is_square<Matrix>::value>::type>
        {
            typedef void type;
        };
        
    } // internal
    
    // small_matrix:
    // a matrix interface to a DataBlock or DataBlock view using Eigen for
    // linear algebra computation. Storage is *column-major* but
    // iterators are available for both column-wise and row-wise traversal
    template< typename Storage_, size_t M, size_t N, 
              typename = typename internal::fits<Storage_, M, N>::type >
    class small_matrix_interface : public small_vector_interface<Storage_>
    {
    public:
        static constexpr size_t _size_ = M*N;
        static constexpr size_t nrows = M;
        static constexpr size_t ncols = N;
        static constexpr bool is_square = (M == N);
        
        typedef typename Storage_::value_type                     scalar_type;
        typedef typename Storage_::value_type                     value_type;
        typedef Storage_                                          storage_type;
        typedef small_matrix_interface<storage_type, M, N>        self_type;
        typedef small_vector_interface<storage_type>              base_type;
        
        // Iterators
        typedef typename base_type::iterator                      columnwise_iterator;
        typedef typename base_type::const_iterator                const_columnwise_iterator;
        typedef internal::OuterIterator<scalar_type, 
            scalar_type*, scalar_type&, M, N>                     rowwise_iterator;
        typedef internal::OuterIterator<scalar_type, 
            const scalar_type*, const scalar_type&, M, N>         const_rowwise_iterator;
        typedef columnwise_iterator                               iterator;
        typedef const_columnwise_iterator                         const_iterator;
        
        // Views
        typedef internal::random_filler<value_type>               random_filler_type;
        typedef DataBlockView< value_type, M, 1>                  column_view_type;
        typedef ConstDataBlockView< value_type, M, 1>             const_column_view_type;
        typedef DataBlockView< value_type, N, M >                 row_view_type;
        typedef ConstDataBlockView< value_type, N, M >            const_row_view_type;
        
        typedef small_vector_interface< column_view_type >        column_type;
        typedef small_vector_interface< const_column_view_type>   const_column_type;
        typedef small_vector_interface< row_view_type >           row_type;
        typedef small_vector_interface< const_row_view_type>      const_row_type;
        
        typedef small_matrix_interface<storage_type, N, M>        self_transpose_type;
        typedef Eigen::Matrix<scalar_type, M, N>                  self_matrix_type;
        typedef Eigen::Map<self_matrix_type>                      self_map_type;
        typedef Eigen::Map<const self_matrix_type>                const_self_map_type;
        
        // Shorthands for compatible types
        template<typename T>
        using datablock = DataBlock<T, _size_>;
        
        template<typename T>
        using like_vector = small_vector_interface< datablock<T> >;
        
        template<typename T>
        using like_matrix = small_matrix_interface<datablock<T>, nrows, ncols>;
        
        template<typename T>
        using like_array = std::array<T, _size_> ;
        
        template<typename OtherStorage = storage_type, 
                 typename = typename internal::fits<OtherStorage, _size_, 1>::type>
        using like_vector_S = small_vector_interface<OtherStorage>;
        
        template<typename OtherStorage, 
                 typename = typename internal::fits<OtherStorage, _size_, 1>::type>
        using like_matrix_S = small_matrix_interface<OtherStorage, nrows, ncols>;
        
        template<typename OtherStorage, 
                 typename = typename internal::fits<OtherStorage, _size_, 1>::type >
        using like_array_S = std::array<typename OtherStorage::value_type, _size_> ;

        template<typename OtherStorage, 
                 size_t P, typename = typename std::enable_if<OtherStorage::size == ncols*P>::type>
        using rhs_matrix_S = small_matrix_interface<OtherStorage, ncols, P>;
        
        // Storage for output types
        typedef like_vector<datablock<bool>> bool_vector;
        typedef like_matrix<datablock<bool>> bool_matrix;
        
        // Iterator ranges
        template<typename Iterator = iterator>
        Iterator begin() { return Iterator(base_type::m_storage.data); }
        
        template<typename Iterator = iterator>
        Iterator end() { return Iterator(base_type::m_storage.data + _size_); } 
        
        template<typename ConstIterator = const_iterator>
        ConstIterator begin() const { return ConstIterator(base_type::m_storage.data); }
        
        template<typename ConstIterator = const_iterator>
        ConstIterator end() const { return ConstIterator(base_type::m_storage.data + _size_); }
        
        columnwise_iterator cbegin() { return columnwise_iterator(base_type::m_storage.data); }
        columnwise_iterator cend() { return columnwise_iterator(base_type::m_storage.data + _size_); }
        
        rowwise_iterator rbegin() { return rowwise_iterator(base_type::m_storage.data); }
        rowwise_iterator rend() { return rowwise_iterator(base_type::m_storage.data + _size_); }
        
        template<typename T>
        using auto_type = typename better_type<scalar_type, T>::type;
        
        template<typename OtherStorage>
        using auto_type_S = auto_type<typename OtherStorage::value_type>;
        
        template<typename T>
        using auto_copy = small_matrix_interface<DataBlock<auto_type<T>, _size_>, nrows, ncols>;
        
        template<typename OtherStorage>
        using auto_copy_S = auto_copy<typename OtherStorage::value_type>;
        
        // template<typename OtherStorage, size_t P>
        // using auto_output_S = auto_output<typename OtherStorage::value_type, P>;

        template<typename OtherStorage, size_t P,
                 typename = typename std::enable_if<OtherStorage::size == P*ncols>::type>
        using auto_output_S = small_matrix_interface<DataBlock<auto_type<typename OtherStorage::value_type>, nrows*P>, nrows, P>;
        
        template<typename OtherStorage>
        using auto_output_vector = 
            small_vector_interface<DataBlock<auto_type_S<OtherStorage>, M>>;
        
        template< typename OtherStorage, size_t P, 
                  typename = typename std::enable_if<OtherStorage::size == ncols*P>::type>
        using rhs_matrix = small_matrix_interface<OtherStorage, ncols, P>;
        
        template <typename T, size_t P>
        using rhs_eigen_matrix = Eigen::Matrix<T, ncols, P>;
        
        template<typename T, size_t P>
        using rhs_map_type = Eigen::Map<rhs_eigen_matrix<T,P>>;
        
        template<typename Storage, 
                 typename = typename std::enable_if<Storage::size == N>::type>
        using rhs_vector = small_vector_interface<Storage>;
        
        // Constructors
        small_matrix_interface(scalar_type val=0) 
        : base_type(val) {}
        
        small_matrix_interface(storage_type storage) : base_type(storage) {}

        template<typename OtherStorage>
        small_matrix_interface(const like_vector_S<OtherStorage>& other)
        : base_type(other) {}

        template<typename T=scalar_type>
        small_matrix_interface(const like_array<T>& other) 
            : base_type(other) {}
        
        template<typename T=scalar_type>
        small_matrix_interface(std::initializer_list<T> vals) : base_type(vals) {}

        template<typename OtherStorage>
        small_matrix_interface(const like_matrix_S<OtherStorage>& other) 
            : base_type(static_cast<const like_vector_S<OtherStorage>&>(other)) {}

        template<typename OtherStorage>
        self_type& operator=(const like_matrix_S<OtherStorage>& other) {
            if (this->m_storage.data == nullptr)
                throw std::runtime_error("Illegal lvalue assignment to unallocated memory");
            std::copy(other.begin(), other.end(), begin());
            return *this;
        }
              
        like_vector<value_type>& as_vector() {
            return static_cast< like_vector<value_type>& >(*this);
        }
        
        const like_vector<value_type>& as_const_vector() const 
        {
            return static_cast<const like_vector<value_type>& >(*this);
        }

        // Initializers
        static self_type identity() {
            self_type r;
            r.as_eigen() = self_matrix_type::Identity();
            return r;
        }
        
        static self_type
        random(value_type min=random_filler_type::default_min,
               value_type max=random_filler_type::default_max) {
            random_filler_type filler(min, max);
            self_type r;
            filler.fill(r.template begin<columnwise_iterator>(), r.template end<columnwise_iterator>());
            return r;
        }
        
        // Basic manipulations
        self_transpose_type transpose() const {
            self_transpose_type r;
            r.as_eigen() = as_const_eigen().transpose();
            return r;
        }

        // Accessors
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

        // Matrix multiplication
        // When square, it can either be a LA multiplication or a 
        // coefficient-wise. LA multiplication assumed
        /*
        template <typename OtherStorage, size_t P>
        auto_output_S<OtherStorage, P> operator*(const rhs_matrix<OtherStorage> &rhs) const
        {
            typedef auto_type_S<OtherStorage> out_type;
            typedef auto_output_S<out_type> out_mat_type;
            out_mat_type r;
            // matrix-matrix multiplication with same-sized matrix
            r.as_eigen() = as_const_eigen().template cast<out_type>() * 
                           rhs.as_const_eigen().template cast<out_type>();
            return r;
        }
        
        // When not square, it must be a coefficient-wise multiplication
        template <typename OtherStorage, 
                  typename = typename std::enable_if<!is_square>::type>
        auto_copy_S<OtherStorage> operator*(const like_matrix_S<OtherStorage> &rhs) const
        {
            typedef auto_type_S<OtherStorage> out_type;
            typedef like_matrix<out_type> out_mat_type;
            out_mat_type r;
            // matrix-matrix multiplication with same-sized matrix
            r.as_eigen() = 
                (as_const_eigen().array().template cast<out_type>() * 
                 rhs.as_const_eigen().array().template cast<out_type>()).matrix();
            return r;
        }*/

        template <typename OtherStorage, size_t P>
        auto_output_S<OtherStorage, P> operator*(const rhs_matrix<OtherStorage, P> &rhs) const
        {
            typedef auto_output_S<OtherStorage, P> out_matrix_type;
            typedef typename out_matrix_type::value_type out_type;
            out_matrix_type r;

            r.as_eigen() = as_const_eigen().template cast<out_type>() * 
                           rhs.as_const_eigen().template cast<out_type>();
            return r;
        }


        // Multiplication by a vector must be LA matrix vector product
        template <typename OtherStorage>
        auto_output_vector<OtherStorage>
        operator*(const rhs_vector<OtherStorage> &rhs) const
        {
            typedef auto_type_S<OtherStorage> out_type;
            auto_output_vector<OtherStorage> r;
            r.as_eigen() = as_const_eigen().template cast<out_type>() * 
                           rhs.as_const_eigen().template cast<out_type>();
            return r;
        }
        
        
        // Multiplication by a scalar is a broadcast operation but it must be
        // replicated here to ensure that a matrix (and not its 
        // small_vector_interface cast) is returned.
        template <typename T>
        auto_copy<T>
        operator*(T val) const
        {
            return base_type::operator*(val);
        }
        
        // Matrix Division
        template <typename OtherStorage>
        auto_output_S<OtherStorage, ncols> 
        operator/(const small_matrix_interface<OtherStorage, ncols, ncols> &rhs) const
        {
            typedef auto_type_S<OtherStorage> out_type;
            typedef like_matrix<out_type> out_mat_type;
            out_mat_type r;
            // matrix-matrix multiplication with same-sized matrix
            r.as_eigen() = as_const_eigen().template cast<out_type>() * 
                               rhs.as_const_eigen().template cast<out_type>().inverse();
            return r;
        }

        // Row-wise division by a vector
        template <typename OtherStorage, 
                  typename = typename std::enable_if<OtherStorage::size == nrows>::type>
        auto_output_S<OtherStorage, ncols>
        operator/(const small_vector_interface<OtherStorage> &rhs) const
        {
            auto_output_S<OtherStorage, ncols> r;
            for (int i=0; i<nrows; ++i) 
            {
                r.row(i) = row(i)/rhs[i];
            }
            return r;
        }
        
        // Division by a scalar is a broadcast operation but it must be
        // replicated here to ensure that a matrix (and not its 
        // small_vector_interface cast) is returned.
        template <typename T>
        auto_copy<T>
        operator/(T val) const
        {
            return base_type::operator/(val);
        }
        
        // Matrix Addition
        template <typename OtherStorage>
        auto_output_S<OtherStorage, ncols> 
        operator+(const small_matrix_interface<OtherStorage, nrows, ncols> &rhs) const
        {
            typedef auto_type_S<OtherStorage> out_type;
            typedef auto_output_S<OtherStorage, ncols> out_mat_type;
            out_mat_type r;
            // matrix-matrix multiplication with same-sized matrix
            r.as_eigen() = as_const_eigen().template cast<out_type>() * 
                               rhs.as_const_eigen().template cast<out_type>().inverse();
            return r;
        }

        // Row-wise addition by a vector
        template <typename OtherStorage, 
                  typename = typename std::enable_if<OtherStorage::size == nrows>::type>
        auto_output_S<OtherStorage, ncols>
        operator+(const small_vector_interface<OtherStorage> &rhs) const
        {
            auto_output_S<OtherStorage, ncols> r;
            for (int i=0; i<nrows; ++i) 
            {
                r.row(i) = row(i) + rhs[i];
            }
            return r;
        }
        
        // Addition by a scalar is a broadcast operation but it must be
        // replicated here to ensure that a matrix (and not its 
        // small_vector_interface cast) is returned.
        template <typename T>
        auto_copy<T>
        operator+(T val) const
        {
            return base_type::operator+(val);
        }
        
        // Matrix subtraction
        template <typename OtherStorage>
        auto_output_S<OtherStorage, ncols> 
        operator-(const small_matrix_interface<OtherStorage, nrows, ncols> &rhs) const
        {
            return base_type::operator-(rhs);
        }

        // Row-wise addition by a vector
        template <typename OtherStorage, 
                  typename = typename std::enable_if<OtherStorage::size == nrows>::type>
        auto_output_S<OtherStorage, ncols>
        operator-(const small_vector_interface<OtherStorage> &rhs) const
        {
            auto_output_S<OtherStorage, ncols> r;
            for (int i=0; i<nrows; ++i) 
            {
                r.row(i) = row(i) - rhs[i];
            }
            return r;
        }
        
        // Subtraction by a scalar is a broadcast operation but it must be
        // replicated here to ensure that a matrix (and not its 
        // small_vector_interface cast) is returned.
        template <typename T>
        auto_copy<T>
        operator-(T val) const
        {
            return base_type::operator-(val);
        }
            
        self_map_type as_eigen()
        {
            return self_map_type(&((*this)(0,0)));
        }

        const_self_map_type as_const_eigen() const {
            return const_self_map_type(&((*this)(0,0)));
        }
    };

    // short-hand for memory-owning small_matrix_interface, squared by default
    template<typename T, size_t M, size_t N=M>
    using small_matrix = small_matrix_interface<DataBlock<T, M*N>, M, N>;
    
    typedef small_matrix<double, 2> mat2;
    typedef small_matrix<double, 3> mat3;
    typedef small_matrix<double, 4> mat4;
    typedef small_matrix<double, 5> mat5;
    typedef small_matrix<double, 6> mat6;
    
    typedef small_matrix<double, 2> dmat2;
    typedef small_matrix<double, 3> dmat3;
    typedef small_matrix<double, 4> dmat4;
    typedef small_matrix<double, 5> dmat5;
    typedef small_matrix<double, 6> dmat6;
    
    typedef small_matrix<float, 2> fmat2;
    typedef small_matrix<float, 3> fmat3;
    typedef small_matrix<float, 4> fmat4;
    typedef small_matrix<float, 5> fmat5;
    typedef small_matrix<float, 6> fmat6;
    
    typedef small_matrix<std::complex<double>, 2> cmat2;
    typedef small_matrix<std::complex<double>, 3> cmat3;
    typedef small_matrix<std::complex<double>, 4> cmat4;
    typedef small_matrix<std::complex<double>, 5> cmat5;
    typedef small_matrix<std::complex<double>, 6> cmat6;
    
    typedef small_matrix<std::complex<double>, 2> dcmat2;
    typedef small_matrix<std::complex<double>, 3> dcmat3;
    typedef small_matrix<std::complex<double>, 4> dcmat4;
    typedef small_matrix<std::complex<double>, 5> dcmat5;
    typedef small_matrix<std::complex<double>, 6> dcmat6;
    
    typedef small_matrix<std::complex<float>, 2> fcmat2;
    typedef small_matrix<std::complex<float>, 3> fcmat3;
    typedef small_matrix<std::complex<float>, 4> fcmat4;
    typedef small_matrix<std::complex<float>, 5> fcmat5;
    typedef small_matrix<std::complex<float>, 6> fcmat6;
    
    typedef small_matrix<int, 2> imat2;
    typedef small_matrix<int, 3> imat3;
    typedef small_matrix<int, 4> imat4;
    typedef small_matrix<int, 5> imat5;
    typedef small_matrix<int, 6> imat6;
    
    typedef small_matrix<long int, 2> lmat2;
    typedef small_matrix<long int, 3> lmat3;
    typedef small_matrix<long int, 4> lmat4;
    typedef small_matrix<long int, 5> lmat5;
    typedef small_matrix<long int, 6> lmat6;
    
    typedef small_matrix<size_t, 2> smat2;
    typedef small_matrix<size_t, 3> smat3;
    typedef small_matrix<size_t, 4> smat4;
    typedef small_matrix<size_t, 5> smat5;
    typedef small_matrix<size_t, 6> smat6;

    template<typename Storage, size_t M, size_t N, 
             typename = typename std::enable_if<Storage::size == M*N>::type >
    std::ostream& operator<<(std::ostream& os, const small_matrix_interface<Storage, M, N>& m) 
    {
        os << '\n' << m.as_const_eigen() << '\n';
        return os;
    }
    
    template<typename Storage, size_t M, size_t N, typename = typename std::enable_if<Storage::size == M*N>::type>
    small_matrix<typename Storage::value_type, N, M> transpose(const small_matrix_interface<Storage, M, N>& m)
    {
        return m.transpose();
    }
    
    template<typename Storage, size_t N, typename = typename std::enable_if<Storage::size == N*N>::type>
    small_matrix<typename Storage::value_type, N, N> make_symmetric(const small_matrix_interface<Storage, N, N>& m, bool lower=true)
    {
        small_matrix<typename Storage::value_type, N, N> r;
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
    small_matrix<T, 3, 3> from_dti(const double* t) 
    {
        small_matrix<T, 3, 3> r;
        r(0,0) = t[1];
        r(1,0) = r(0,1) = t[2];
        r(2,0) = r(0,2) = t[3];
        r(1,1) = t[4];
        r(2,1) = r(1,2) = t[5];
        r(2,2) = t[6];
        return r;
    }

    template<typename Storage, size_t M, size_t N>
    typename Storage::value_type trace(const small_matrix_interface<Storage, M, N>& m)
    {
        return m.as_const_eigen().trace();
    }

    template<typename T, size_t M>
    T determinant(const small_matrix<T, M>& m) 
    {
        return m.as_const_eigen().determinant();
    }

    template<typename Storage, size_t M>
    small_matrix<typename Storage::value_type, M> inverse(const small_matrix_interface<Storage, M, M>& m) {
        small_matrix<typename Storage::value_type, M> r;
        r.as_eigen() = m.as_const_eigen().inverse();
        return r;
    }

    template<typename Storage, size_t M, size_t N,
             typename = typename std::enable_if<(M>N)>::type >
    small_matrix<typename Storage::value_type, N, M>
    moore_penrose_pseudoinverse(const small_matrix_interface<Storage, M, N>& m) {
        return inverse(transpose(m)*m)*transpose(m);
    }

    template<typename Storage, size_t M, size_t N,
             typename std::enable_if<(M<N)>::type >
    small_matrix<typename Storage::value_type, N, M>
    moore_penrose_pseudoinverse(const small_matrix_interface<Storage, M, N>& m) {
        return transpose(m)*inverse(m*transpose(m));
    }

    template<typename Storage, size_t M>
    small_matrix<typename Storage::value_type, M, M>
    moore_penrose_pseudoinverse(const small_matrix_interface<Storage, M, M>& m) {
        return inverse(m);
    }
    
    template<typename Storage1, typename Storage2>
    small_matrix<typename better_type<typename Storage1::value_type, 
                                      typename Storage2::value_type>::type, 
                 Storage1::size, Storage2::size> 
    outer(const small_vector_interface<Storage1>& v0, const small_vector_interface<Storage2>& v1) {
        typedef typename better_type<typename Storage1::value_type,
                                     typename Storage2::value_type>::type out_type;
        small_matrix<out_type, Storage1::size, Storage2::size> r;
        r.as_eigen() = v0.as_const_eigen().template cast<out_type>() *
                       v1.as_const_eigen().transpose().template cast<out_type>();
        return r;
    }
    
    template<typename Storage1, typename Storage2, size_t M, size_t N, 
             typename = typename std::enable_if<Storage1::size == M*N && Storage2::size == M*N >::type>
    typename better_type<typename Storage1::value_type, typename Storage2::value_type>::type
    contraction(const small_matrix_interface<Storage1, M, N>& m1, 
                const small_matrix_interface<Storage2, M, N>& m2) {
        typedef typename Storage1::value_type value_type1;
        typedef typename Storage2::value_type value_type2;
        typedef typename better_type<value_type1, value_type2>::type out_type;
        out_type r = 0;
        auto it1 = m1.begin();
        for (auto it2=m2.begin(); it2!=m2.end(); ++it1, ++it2) {
            r += (*it1) * (*it2);
        }
        return r;
    }
            
    
    template<typename T, typename OtherStorage, size_t M, size_t N,
             typename = typename std::enable_if<std::is_scalar<T>::value>::type>
    small_matrix<typename better_type<T, typename OtherStorage::value_type>::type, M, N> operator*(T val, const small_matrix_interface<OtherStorage, M, N>& m)
    {
        return m * val;
    }
    
    template<typename T, typename OtherStorage, size_t M, size_t N,
             typename = typename std::enable_if<std::is_scalar<T>::value>::type>
    small_matrix<typename better_type<T, typename OtherStorage::value_type>::type, M, N> operator+(T val, const small_matrix_interface<OtherStorage, M, N>& m)
    {
        return m + val;
    }
    
    template<typename T, typename OtherStorage, size_t M, size_t N,
             typename = typename std::enable_if<std::is_scalar<T>::value>::type>
    small_matrix<typename better_type<T, typename OtherStorage::value_type>::type, M, N> operator-(T val, const small_matrix_interface<OtherStorage, M, N>& m)
    {
        return m - val;
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
                                                 internal::is_complex<typename Storage::value_type>::value &&
                                                 std::is_scalar<T2>::value >::type >
    void eigensystem(small_vector_interface<Storage>& eigenvalues, 
                     small_matrix<std::complex<T2>, M>& eigenvectors, 
                     const small_matrix<T, M>& m) 
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
                                                 internal::is_complex<typename Storage::value_type>::value>::type >
    void complex_eigenvalues(small_vector_interface<Storage>& evals,
                             const small_matrix<T, M>& m)
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
                                                !internal::is_complex<typename Storage::value_type>::value>::type>
    void sym_eigensystem(small_vector_interface<Storage>& eigenvalues,
                         small_matrix<T2, M>& eigenvectors,
                         const small_matrix<T, M>& m)
    {
        typedef Eigen::Matrix<T, M, M> eigen_type;
        Eigen::SelfAdjointEigenSolver<eigen_type> solver(m.as_const_eigen());
        eigenvalues.as_eigen() = solver.eigenvalues().template cast<typename Storage::value_type>();
        eigenvectors.as_eigen() = solver.eigenvectors().template cast<T2>();
        sort_real_eigenvalues(eigenvalues, eigenvectors);
    }

    template <typename Storage, typename T, size_t M, 
              typename = typename std::enable_if<Storage::size == M &&
                                                 !internal::is_complex<typename Storage::value_type>::value>::type>
    void real_eigenvalues(small_vector_interface<Storage> &evals,
                          const small_matrix<T, M> &m)
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

    template<typename T, size_t M>
    inline small_matrix<T, M, M> deviatoric(const small_matrix<T, M, M>& tensor) {
        double lambda = trace(tensor) / static_cast<double>(tensor.nrows);
        return tensor - lambda*small_matrix<T, M, M>::identity();
    }
    
    template < typename T> 
    inline T FA(const small_matrix<T, 3>& tensor)
    {
        small_matrix<T, 3, 3> D = deviatoric(tensor);
        return sqrt(1.5) * norm(D) / norm(tensor);
    }
    
    template < typename T >
    inline small_vector<T, 3> westin_aniso(const small_matrix<T, 3>& tensor)
    {
        typedef small_matrix<T, 3> mat_t;
        typedef small_vector<T, 3> vec_t;
        mat_t D = deviatoric(tensor);
        
        vec_t lambdas;
        mat_t evecs;
        sym_eigensystem(lambdas, evecs, D);
        
        double cl = lambdas[0] - lambdas[1];
        double cp = 2.*(lambdas[1] - lambdas[2]);
        double cs = 3.*lambdas[2];
        vec_t ans(cl, cp, cs);
        return ans / trace(D);
    }
    
} // namespace spurt
