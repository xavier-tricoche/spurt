#pragma once

#include <math/basic_math.hpp>
#include <math/types.hpp>
#include <algorithm>

namespace spurt {
    
    namespace internal {
        template<typename Storage1, typename Storage2, typename Storage3>
        struct better_type3 {
            typedef typename better_type<typename Storage1::type, 
                                         typename better_type<typename Storage2::type,
                                                              typename Storage3::type>::type>::type type;
        };
    }
    
    // small tensor:
    // a tensor interface to a DataBlock or DataBlock view using Eigen for
    // linear algebra computation. Storage is *column-major* but
    // iterators are available for column-wise, row-wise, and layer-wise traversals
    template< typename Storage_, size_t M, size_t N=M, size_t P=N,
              typename = typename internal::fits<Storage_, M, N*P>::type >
    class small_tensor_interface : public small_vector_interface<Storage_>
    {
    public:
        static constexpr size_t _size_ = M*N*P;
        static constexpr size_t nrows = M;
        static constexpr size_t ncols = N;
        static constexpr size_t nlays = P;
        static constexpr bool is_cube = ((M == N) && (N == P));
        
        typedef typename Storage_::value_type                     scalar_type;
        typedef typename Storage_::value_type                     value_type;
        typedef Storage_                                          storage_type;
        typedef small_tensor_interface<storage_type, M, N, P>     self_type;
        typedef small_vector_interface<storage_type>              base_type;
        
        // Iterators
        typedef typename base_type::iterator                      iterator;
        typedef typename base_type::const_iterator                const_iterator;
        typedef typename base_type::iterator                      columnwise_iterator;
        typedef typename base_type::const_iterator                const_columnwise_iterator;
        typedef internal::OuterIterator<scalar_type, 
             scalar_type*, scalar_type&, M, N>                    rowwise_iterator;
        typedef internal::OuterIterator<scalar_type,
             const scalar_type*, const scalar_type&, M, N>        const_rowwise_iterator;
        typedef internal::OuterIterator<scalar_type, 
             scalar_type*, scalar_type&, M*N, P>                  layerwise_iterator;
        typedef internal::OuterIterator<scalar_type,
             const scalar_type*, const scalar_type&, M*N, P>      const_layerwise_iterator;
        
        // Views
        typedef internal::random_filler<value_type>               random_filler_type;
        typedef DataBlockView< value_type, M, 1>                  column_view_type;
        typedef ConstDataBlockView< value_type, M, 1>             const_column_view_type;
        typedef DataBlockView< value_type, N, M >                 row_view_type;
        typedef ConstDataBlockView< value_type, N, M >            const_row_view_type;
        typedef DataBlockView< value_type, M*N, 1>                layer_view_type;
        typedef ConstDataBlockView< value_type, M*N, 1>           const_layer_view_type;
        
        typedef small_vector_interface< column_view_type >        column_type;
        typedef small_vector_interface< const_column_view_type >  const_column_type;
        typedef small_vector_interface< row_view_type >           row_type;
        typedef small_vector_interface< const_row_view_type >     const_row_type;
        typedef small_matrix_interface< layer_view_type, M, N >         layer_type;
        typedef small_matrix_interface< const_layer_view_type, M, N >   const_layer_type;
        
        // Shorthands for compatible types
        template<typename T>
        using datablock = DataBlock<T, _size_>;
        
        template<typename T>
        using like_vector = small_vector_interface< datablock<T> >;
        
        template<typename T>
        using like_matrix = small_matrix_interface<datablock<T>, nrows, ncols*nlays>;
        
        template<typename T>
        using like_array = std::array<T, _size_> ;
        
        template<typename T>
        using rhs_matrix = small_matrix_interface<datablock<T>, ncols, nlays>;
        
        template<typename OtherStorage = storage_type, 
                 typename = typename internal::fits<OtherStorage, _size_, 1>::type>
        using like_vector_S = small_vector_interface<OtherStorage>;
        
        template<typename OtherStorage, 
                 typename = typename internal::fits<OtherStorage, _size_, 1>::type>
        using like_matrix_S = small_matrix_interface<OtherStorage, nrows, ncols*nlays>;
        
        template<typename OtherStorage = storage_type, 
                 typename = typename internal::fits<OtherStorage, _size_, 1>::type >
        using like_array_S = std::array<typename OtherStorage::value_type, _size_> ;
                 
        template<typename OtherStorage,
                 typename = typename std::enable_if<OtherStorage::size == nrows*nlays>::type >
        using rhs_matrix_S = small_matrix_interface<OtherStorage, ncols, nlays>;
                 
                 
        // Iterator ranges
        template<typename Iterator = iterator>
        Iterator begin() { return Iterator(base_type::m_storage.data); }

        template<typename Iterator = iterator>
        Iterator end() { return Iterator(base_type::m_storage.data + _size_); } 

        columnwise_iterator cbegin() { return columnwise_iterator(base_type::m_storage.data); }
        columnwise_iterator cend() { return columnwise_iterator(base_type::m_storage.data + _size_); }

        const_columnwise_iterator cbegin() const { return const_columnwise_iterator(base_type::m_storage.data); }
        const_columnwise_iterator cend() const { return const_columnwise_iterator(base_type::m_storage.data + _size_); }

        rowwise_iterator rbegin() { return rowwise_iterator(base_type::m_storage.data); }
        rowwise_iterator rend() { return rowwise_iterator(base_type::m_storage.data + _size_); }

        const_rowwise_iterator rbegin() const { return const_rowwise_iterator(base_type::m_storage.data); }
        const_rowwise_iterator rend() const { return const_rowwise_iterator(base_type::m_storage.data + _size_); }

        layerwise_iterator lbegin() { return layerwise_iterator(base_type::m_storage.data); }
        layerwise_iterator lend() { return layerwise_iterator(base_type::m_storage.data + _size_); }

        const_layerwise_iterator lbegin() const { return const_layerwise_iterator(base_type::m_storage.data); }
        const_layerwise_iterator lend() const { return const_layerwise_iterator(base_type::m_storage.data + _size_); }
        
        // Constructors
        small_tensor_interface(scalar_type val=0) 
        : base_type(val) {}

        template<typename OtherStorage>
        small_tensor_interface(const like_vector_S<OtherStorage>& other)
        : base_type(other) {}

        template<typename T=scalar_type>
        small_tensor_interface(const like_array<T>& other) 
            : base_type(other) {}
        
        template<typename T=scalar_type>
        small_tensor_interface(std::initializer_list<T> vals) : base_type(vals) {}

        template<typename OtherStorage>
        small_tensor_interface(const like_matrix_S<OtherStorage>& other) 
            : base_type(static_cast<const like_vector_S<OtherStorage>&>(other)) {}
              
        like_vector<value_type>& as_vector() {
            return static_cast< like_vector<value_type>& >(*this);
        }
        
        const like_vector<value_type>& as_const_vector() const 
        {
            return static_cast<const like_vector<value_type>& >(*this);
        }
        
        static self_type
        random(value_type min=random_filler_type::default_min,
               value_type max=random_filler_type::default_max) {
            random_filler_type filler(min, max);
            self_type r;
            filler.fill(r.begin(), r.end());
            return r;
        }

        // Accessors
        const scalar_type& operator()(int row, int col, int lay) const 
        {
            return static_cast<const base_type&>(*this)[row + nrows*(col + nlays*lay)];
        }
        scalar_type& operator()(int row, int col, int lay) 
        {
            return static_cast<base_type&>(*this)[row + nrows*(col + nlays*lay)];
        }

        column_type column(int col, int lay=0)
        {
            return column_type(column_view_type(&((*this)(0, col, lay))));
        }
        
        const_column_type column(int col, int lay=0) const 
        {
            return const_column_type(const_column_view_type(&((*this)(0, col, lay))));
        }

        row_type row(int row, int lay=0)
        {
            return row_type(row_view_type(&((*this)(row, 0, lay))));
        }
        
        const_row_type row(int row, int lay=0) const
        {
            return const_row_type(const_row_view_type(&((*this)(row, 0, lay))));
        }
        
        layer_type layer(int lay)
        {
            return layer_type(layer_view_type(&((*this)(0, 0, lay))));
        }
        
        const_layer_type layer(int lay) const
        {
            return const_layer_type(const_layer_view_type(&((*this)(0, 0, lay))));
        }
        
        layer_type operator[](int k) {
            return layer(k);
        }
        
        const_layer_type operator[](int k) const {
            return layer(k);
        }
        
    };
    
    template<typename T, size_t M, size_t N, size_t P>
    using small_tensor = small_tensor_interface<DataBlock<T, M*N*P>, M, N, P>;
    
    template<typename Storage1, typename Storage2, size_t M, size_t N, size_t P,
             typename = typename std::enable_if<Storage1::size == M*N*P && Storage2::size == N*P >::type>
    small_vector<typename better_type<typename Storage1::value_type, typename Storage2::value_type>::type, M>
    contraction(const small_tensor_interface<Storage1, M, N, P>& t, 
                const small_matrix_interface<Storage2, N, P>& m) {
        typedef typename Storage1::value_type value_type1;
        typedef typename Storage2::value_type value_type2;
        typedef typename better_type<value_type1, value_type2>::type out_type;
        small_vector<out_type, M> r = 0;
        for (int i=0; i<M; ++i) {
            for (int j=0; j<N; ++j) {
                for (int k=0; k<P; ++k) {
                    r[i] += t(i,j,k) * m(j,k);
                }
            }
        }
        return r;
    }
    
    template<typename Storage1, typename Storage2, typename Storage3>
    small_tensor_interface<
        DataBlock<typename internal::better_type3<Storage1, Storage2, Storage3>::type, 
                  Storage1::size*Storage2::size*Storage3::size>,
        Storage1::size, Storage2::size, Storage3::size >
    outer(const small_vector_interface<Storage1>& v1,
          const small_vector_interface<Storage2>& v2,
          const small_vector_interface<Storage3>& v3) {
          const size_t& M = Storage1::size;
          const size_t& N = Storage2::size;
          const size_t& P = Storage3::size;
          
          typedef typename internal::better_type3<Storage1, Storage2, Storage3>::type out_type;
          small_tensor<out_type, M, N, P> r = 0;
          for (int k=0; k<P; ++k) {
              for (int j=0; j<N; ++j) {
                  for (int i=0; i<M; ++i) {
                      r(i,j,k) = v1[i]*v2[j]*v3[k];
                  }
              }
          }
          return r;
      }

} // namespace spurt












