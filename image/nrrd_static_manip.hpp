#include <functional>

#include <boost/integer/static_min_max.hpp>

#include <image/nrrd_wrapper.hpp>
#include <Eigen/Core>
#include <Eigen/SVD>

#include <util/timer.hpp>

extern double toverhead;

namespace spurt {
    size_t nrrd_size(const Nrrd* A, bool skip_first_col=false) {
        size_t s=1;
        for (int i=(skip_first_col ? 1 : 0); i<A->dim; ++i) {
            s*=A->axis[i].size;
        }
        return s;
    }
    
    bool matching_sizes(const Nrrd* A, const Nrrd* B, bool skip_first=true) {
        if (A->dim!=B->dim) return false;
        for (int d=(skip_first ? 1 : 0); d<A->dim; ++d) {
            if (A->axis[d].size!=B->axis[d].size) return false;
        }
        return true;
    }
    
    template<typename T, int Nrows, int Ncols, bool SameRawType=false>
    struct nrrd_matrix_wrapper {
        typedef T value_t;
        typedef std::array<value_t, Nrows*Ncols>     vec_t;
        typedef Eigen::Matrix<value_t, Nrows, Ncols> matrix_t;
        typedef Eigen::Matrix<value_t, Ncols, Nrows> transpose_t;
        typedef Eigen::Matrix<value_t, Nrows, 1> col_t;
        typedef Eigen::Matrix<value_t, 1, Ncols> row_t;
        static constexpr int nrows=Nrows;
        static constexpr int ncols=Ncols;
        
        nrrd_matrix_wrapper(const Nrrd* nin) {
            if (nin->axis[0].size!=Nrows*Ncols) {
                throw std::runtime_error("matrix dimensions mismatch");
            }
            spurt::to_vector<value_t>(_val_array, nin);
            _mat_array=reinterpret_cast<matrix_t*>(&_val_array[0]);
        }
        
        const matrix_t& operator[](unsigned long int n) const {
            return _mat_array[n];
        }
        
        const transpose_t& transpose(unsigned long int n) const {
            return _mat_array[n];
        }
        
        const col_t& column(unsigned long int n, int c) const {
            return _mat_array[n].col(c);
        }
        
        const row_t& row(unsigned long int n, int r) const {
            return _mat_array[n].row(r);
        }
        
        std::vector<T> _val_array;
        matrix_t* _mat_array;
    };
    
    namespace detail {
        
        template<typename T>
        void copyAxisInfo(Nrrd* nout, const Nrrd* nin, int info) {
            T attr[NRRD_DIM_MAX];
            nrrdAxisInfoGet_nva(nin, info, attr);
            if (nin->dim==nout->dim)
                nrrdAxisInfoSet_nva(nout, info, attr);
            else
                nrrdAxisInfoSet_nva(nout, info, &attr[1]);
        }
        
        template<typename T, int Nrows, int Ncols>
        T* __matrix_transpose(const Nrrd* nin) {
            typedef T value_t;           
            typedef nrrd_matrix_wrapper<T, Nrows, Ncols> wrapper_t;
            typedef typename wrapper_t::matrix_t rval_matrix_t;
            typedef Eigen::Matrix<T, Ncols, Nrows> lval_matrix_t;

            // nvis::timer tovh;
            //
            assert(Nrows*Ncols==nin->axis[0].size);
            size_t nmats=nrrd_size(nin)/(Nrows*Ncols);
            value_t* data=(value_t*)calloc(Nrows*Ncols*nmats, sizeof(value_t));
            lval_matrix_t* res=reinterpret_cast<lval_matrix_t*>(data);
            wrapper_t matrices(nin);
            //
            // toverhead += tovh.elapsed();
            
            for (size_t n=0; n<nmats; ++n) {
                res[n]=matrices[n].transpose();
            }
            return data;
        }
        
        template<typename T, int Nrows, int Ncols=Nrows, int P=Ncols>
        T* __matrix_product(const Nrrd* nin1, const Nrrd* nin2) {
            typedef T value_t;
            typedef nrrd_matrix_wrapper<T, Nrows, Ncols> leftw_t;
            typedef typename leftw_t::matrix_t left_mat_t;
            typedef nrrd_matrix_wrapper<T, Ncols, P> rightw_t;
            typedef typename rightw_t::matrix_t right_mat_t;
            typedef nrrd_matrix_wrapper<T, Nrows, P> lvalw_t;
            typedef typename lvalw_t::matrix_t lval_mat_t;
            
            // nvis::timer tovh;
            //
            assert(Nrows*Ncols==nin1->axis[0].size);
            assert(Ncols*P==nin2->axis[0].size);
            assert(matching_sizes(nin1, nin2));
            size_t nmats=nrrd_size(nin1)/(Nrows*Ncols);
            value_t* data=(value_t*)calloc(Nrows*P*nmats, sizeof(value_t));
            lval_mat_t* res=reinterpret_cast<lval_mat_t*>(data);
            leftw_t left_matrices(nin1);
            rightw_t right_matrices(nin2);
            //
            // toverhead += tovh.elapsed();
            
            for (size_t n=0; n<nmats; ++n) {
                res[n]=left_matrices[n]*right_matrices[n];
            }
            return data;
        }
        
        template<typename T, int Nrows, int Ncols=Nrows, int P=Ncols>
        T* __matrix_transpose_product(const Nrrd* nin1, const Nrrd* nin2) {
            typedef T value_t;
            typedef nrrd_matrix_wrapper<T, Ncols, Nrows> leftw_t;
            typedef typename leftw_t::transpose_t trans_mat_t;
            typedef nrrd_matrix_wrapper<T, Ncols, P> rightw_t;
            typedef typename rightw_t::matrix_t right_mat_t;
            typedef nrrd_matrix_wrapper<T, Nrows, P> lvalw_t;
            typedef typename lvalw_t::matrix_t lval_mat_t;
            
            // nvis::timer tovh;
            //
            assert(Nrows*Ncols==nin1->axis[0].size);
            assert(Ncols*P==nin2->axis[0].size);
            assert(matching_sizes(nin1, nin2));
            
            size_t nmats=nrrd_size(nin1)/(Nrows*Ncols);
            value_t* data=(value_t*)calloc(Nrows*P*nmats, sizeof(value_t));
            lval_mat_t* res=reinterpret_cast<lval_mat_t*>(data);
            leftw_t left_matrices(nin1);
            rightw_t right_matrices(nin2);
            //
            // toverhead += tovh.elapsed();
            
            for (size_t n=0; n<nmats; ++n) {
                res[n]=left_matrices[n].transpose()*right_matrices[n];
            }
            return data;
        }
        
        template<typename T, int N>
        T* __matrix_inverse(const Nrrd* nin) {
            typedef T value_t;            
            typedef nrrd_matrix_wrapper<T, N, N> wrapper_t;
            typedef typename wrapper_t::matrix_t matrix_t;

            // nvis::timer tovh;
            //
            assert(N*N==nin->axis[0].size);
            size_t nmats=nrrd_size(nin)/(N*N);
            value_t* data=(value_t*)calloc(N*N*nmats, sizeof(value_t));
            matrix_t* res=reinterpret_cast<matrix_t*>(data);
            wrapper_t matrices(nin);
            //
            // toverhead += tovh.elapsed();
            
            for (size_t n=0; n<nmats; ++n) {
                res[n]=matrices[n].inverse();
            }
            return data;
        }
        
        template<typename T, int Nrows, int Ncols>
        void __matrix_SVD(const Nrrd* nin, 
                          T*& sv, /* singular values */
                          T*& leftsv, /* left singular vectors */
                          T*& rightsv) /* right singular vectors */ {
                              
            using boost::static_signed_min;
            constexpr int Ndiag=static_signed_min<Nrows,Ncols>::value; // number of singular values
            
            typedef T value_t;
            typedef nrrd_matrix_wrapper<T, Nrows, Ncols> wrapper_t;
            typedef typename wrapper_t::matrix_t matrix_t;
            typedef Eigen::Matrix<T, Nrows, Nrows> matrixU_t;
            typedef Eigen::Matrix<T, Ncols, Ncols> matrixV_t;
            typedef Eigen::Matrix<T, Ndiag, 1> sigma_t;
            
            // nvis::timer tovh;
            //
            assert(Nrows*Ncols==nin->axis[0].size);
            size_t nmats=nrrd_size(nin)/(Nrows*Ncols);
            sv=(value_t*)calloc(P*nmats, sizeof(value_t));
            leftsv=(value_t*)calloc(Nrows*Nrows*nmats, sizeof(value_t));
            rightsv=(value_t*)calloc(Ncols*Ncols*nmats, sizeof(value_t));
            wrapper_t matrices(nin);
            matrixU_t* left_singular_vecs=
                reinterpret_cast<matrixU_t*>(leftsv);
            matrixV_t* right_singular_vecs=
                reinterpret_cast<matrixV_t*>(rightsv);
            sigma_t* singular_vals=
                reinterpret_cast<sigma_t*>(sv);
            //
            // toverhead += tovh.elapsed();
            
            for (size_t n=0; n<nmats; ++n) {
                Eigen::JacobiSVD<matrix_t> svd(matrices[n], Eigen::ComputeFullU | Eigen::ComputeFullV); 
                singular_vals[n]=svd.singularValues();
                left_singular_vecs[n]=svd.matrixU();
                right_singular_vecs[n]=svd.matrixV();
            }
        }
        
    } // detail
    
    template<typename T, int Size>
    Nrrd* wrap_and_copy_header(const Nrrd* src, T* data) {
        // nvis::timer tovh;
        //
        Nrrd *nout = nrrdNew();
        int dim=(Size<=1 ? src->dim-1 : src->dim);
        size_t sizes[src->dim];
        if (Size>1) {
            sizes[0]=Size;
            for (int i=1; i<src->dim; ++i) sizes[i]=src->axis[i].size;
        }
        else {
            for (int i=0; i<src->dim-1; ++i) sizes[i]=src->axis[i+1].size;
        }
        if (nrrdWrap_nva(nout, data, nrrd_type_id<T>::id, dim, sizes)) {
            throw std::runtime_error(error_msg("copy_header"));
        }
        detail::copyAxisInfo<double>(nout, src, nrrdAxisInfoSpacing);
        detail::copyAxisInfo<double>(nout, src, nrrdAxisInfoMin);
        detail::copyAxisInfo<double>(nout, src, nrrdAxisInfoMax);
        detail::copyAxisInfo<int>(nout, src, nrrdAxisInfoCenter);
        detail::copyAxisInfo<int>(nout, src, nrrdAxisInfoKind);
        detail::copyAxisInfo<char*>(nout, src, nrrdAxisInfoLabel);
        detail::copyAxisInfo<double>(nout, src, nrrdAxisInfoThickness);
        detail::copyAxisInfo<char*>(nout, src, nrrdAxisInfoUnits);
        double d_array[NRRD_DIM_MAX];
        nrrdSpaceOriginGet(src, d_array);
        if (Size>1) nrrdSpaceOriginSet(nout, d_array);
        else nrrdSpaceOriginSet(nout, &d_array[1]);
        //
        // toverhead += tovh.elapsed();
        
        return nout;
    }
    
    template<typename T, int Nrows, int Ncols1=Nrows, int Ncols2=Ncols1>
    struct nrrd_matrix_product {
        typedef T value_t;
        Nrrd* operator()(const Nrrd* A, const Nrrd* B) {
            value_t* data=detail::__matrix_product<value_t, Nrows, Ncols1, Ncols2>(A, B);
            return wrap_and_copy_header<value_t, Nrows*Ncols2>(A, data);
        }
    };
    
    template<typename T, int Nrows, int Ncols=Nrows>
    struct nrrd_matrix_transpose {
        typedef T value_t;
        Nrrd* operator()(const Nrrd* A) {
            value_t* data=detail::__matrix_transpose<value_t, Nrows, Ncols>(A);
            return wrap_and_copy_header<value_t, Nrows*Ncols>(A, data);
        }
    };
    
    template<typename T, int Nrows, int Ncols1=Nrows, int Ncols2=Ncols1>
    struct nrrd_matrix_transpose_product {
        typedef T value_t;
        Nrrd* operator()(const Nrrd* A, const Nrrd* B) {
            value_t* data=detail::__matrix_transpose_product<value_t, Ncols1, Nrows, Ncols2>(A, B);
            return wrap_and_copy_header<value_t, Nrows*Ncols2>(A, data);
        }
    };
    
    template<typename T, int N>
    struct nrrd_matrix_invert {
        typedef T value_t;
        Nrrd* operator()(const Nrrd* A) {
            value_t* data=detail::__matrix_inverse<value_t, N>(A);
            return wrap_and_copy_header<value_t, N*N>(A, data);
        }
    };
    
    template<typename T, int Nrows, int Ncols=Nrows>
    struct nrrd_matrix_svd {
        typedef T value_t;
        void operator()(const Nrrd* A, Nrrd*& singvals, Nrrd*& leftvecs, Nrrd*& rightvecs) {
            value_t *svals_data, *leftvec_data, *rightvec_data;
            using boost::static_signed_min;
            constexpr int P=static_signed_min<Nrows,Ncols>::value;
            detail::__matrix_SVD<T, Nrows, Ncols>(A, svals_data, leftvec_data, rightvec_data);
            singvals=wrap_and_copy_header<value_t,P>(A, svals_data);
            leftvecs=wrap_and_copy_header<value_t,Nrows*Nrows>(A, leftvec_data);
            rightvecs=wrap_and_copy_header<value_t,Ncols*Ncols>(A, rightvec_data);
        }
    };

} // spurt
