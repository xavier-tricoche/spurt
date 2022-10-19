#include <array>
#include <iomanip>
#include <map>

#include <image/nrrd_manip.hpp>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>

#include <util/timer.hpp>

#include <misc/progress.hpp>


namespace xavier { namespace nrrd_manip {

template<typename T>
using mat_t=Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using col_t=Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T>
using row_t=Eigen::Matrix<T, 1, Eigen::Dynamic>;
template<typename T>
using eigensolver_t=Eigen::EigenSolver< mat_t<T> >;
template<typename T>
using sym_eigensolver_t=Eigen::SelfAdjointEigenSolver< mat_t<T> >;

template<typename T>
Nrrd* transpose(const Nrrd* nin, int Nrows, int Ncols, bool show_progress=false);

template<typename T>
Nrrd* product(const Nrrd* nin1, const Nrrd* nin2,
              int Nrows, int Ncols1, int Ncols2, bool show_progress=false);

template<typename T>
Nrrd* sum(const Nrrd* nin1, const Nrrd* nin2, int Nrows, int Ncols, bool show_progress=false);

template<typename T>
Nrrd* subtraction(const Nrrd* nin1, const Nrrd* nin2, int Nrows, int Ncols, bool show_progress=false);

template<typename T>
Nrrd* trans_product(const Nrrd* nin1, const Nrrd* nin2,
                    int Nrows, int Ncols1, int Ncols2, bool show_progress=false);

template<typename T>
Nrrd* cauchy_green(const Nrrd* nin, int Nrows, bool show_progress=false);

template<typename T>
Nrrd* inverse(const Nrrd* nin, int Nrows, bool show_progress=false);

template<typename T>
void SVD(Nrrd*& sings, Nrrd*& left, Nrrd*& right,
         const Nrrd* nin, int Nrows, int Ncols, bool show_progress=false);

template<typename T>
void Eigendecomposition(Nrrd*& evals, Nrrd*& evecs,
                        const Nrrd* nin, int Nrows, bool issym=false,
                        bool show_progress=false);

template<typename T>
Nrrd* FTLE(const Nrrd* nin, int Nrows, bool show_progress=false);

namespace detail {

template<typename T>
void copyAxisInfo(Nrrd* nout, const Nrrd* nin, int info);

template<typename T>
Nrrd* wrap_and_copy_header(const Nrrd* src, T* data, int valdim);

template<typename T>
Nrrd* mat2nrrd(const std::vector< mat_t<T> >& mats,
               const std::array<size_t, 3>& res,
               int Nrows, int Ncols, bool show_progress=false);

template<typename T>
Nrrd* vec2nrrd(const std::vector< col_t<T> >& vecs,
               const std::array<size_t, 3>& res,
               int Nrows, bool show_progress=false);

template<typename T>
void nrrd2mat(const Nrrd* nin, std::vector< mat_t<T> >& mats,
              std::array<size_t, 3>& res,
              int Nrows, int Ncols, bool show_progress=false);

template<typename T, typename T1=T>
void do_transpose(T* out_ptr, T1* in_ptr,
                  size_t N, int Nrows, int Ncols,
                  bool show_progress=false);

template<typename T, typename T1=T>
void do_SVD(T* svals_ptr, T* left_ptr, T* right_ptr, T1* in_ptr,
            size_t N, int Nrows, int Ncols,
            bool show_progress=false);

template<typename TypeOut_, typename TypeIn_, typename Solver_>
void do_Eigendecomposition(TypeOut_* evals_ptr, TypeOut_* evecs_ptr, TypeIn_* in_ptr,
                           size_t N, int Nrows, bool issym=false,
                           bool show_progress=false);

template<typename T, typename T1>
void do_FTLE(T* out_ptr, T1* in_ptr, size_t N, int Nrows, bool show_progress=false);

template<typename T, typename T1=T>
void do_inverse(T* out_ptr, T1* in_ptr, size_t N,
                int Nrows, bool show_progress=false);

template<typename T, typename T1=T, typename T2=T1>
void do_sum(T* out_ptr, T1* a_ptr, T2* b_ptr,
            size_t N, int Nrows, int Ncols, bool show_progress=false);

template<typename T, typename T1=T, typename T2=T1>
void do_subtraction(T* out_ptr, T1* a_ptr, T2* b_ptr,
                    size_t N, int Nrows, int Ncols, bool show_progress=false);

template<typename T, typename T1=T, typename T2=T1>
void do_product(T* out_ptr, T1* a_ptr, T2* b_ptr, size_t N,
                int Nrows, int Ncols1, int Ncols2, bool show_progress=false);

template<typename T, typename T1=T, typename T2=T1>
void do_trans_product(T* out_ptr, T1* a_ptr, T2* b_ptr,
                      size_t N,
                      int Nrows, int Ncols1, int Ncols2, bool show_progress=false);

template<typename T, typename T1=T>
void do_cauchy_green(T* out_ptr, T1* in_ptr, size_t N, int Nrows, bool show_progress=false);

template<typename Type_>
void sort_eigenvalues(Type_* eigenvalues, Type_* eigenvectors, size_t N, size_t Nrows);
} // detail

} // nrrd_manip

} // xavier


namespace xavier { namespace nrrd_manip {
template<typename T>
Nrrd* transpose(const Nrrd* nin, int Nrows, int Ncols, bool show_progress) {

    size_t N=xavier::nrrd_utils::nrrd_size(nin, true);
    T* out_array=(T*)malloc(N*Ncols*Nrows*sizeof(T));

    if (nin->type==nrrdTypeFloat)
        detail::do_transpose<T, float>(out_array, (float*)nin->data, N,
            Nrows, Ncols, show_progress);
    else if (nin->type==nrrdTypeDouble)
        detail::do_transpose<T, double>(out_array, (double*)nin->data, N,
            Nrows, Ncols, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin->type));
    }

    return detail::wrap_and_copy_header<T>(nin, out_array, Nrows*Ncols);
    // size_t dims[4];
    // dims[0]=Nrows*Ncols;
    // for (int i=1; i<nin->dim; ++i) dims[i]=nin->axis[i].size;
    // Nrrd *nout = nrrdNew();
    // if (nrrdWrap_nva(nout, out_array, nrrd_type_id<T>::id,
    //                  nin->dim, &dims[0])) {
    //     std::cout << "error in nrrdWrap" << '\n';
    //     throw std::runtime_error("nrrd wrap error in transpose");
    // }
    // return nout;
}

template<typename T>
Nrrd* product(const Nrrd* nin1, const Nrrd* nin2,
              int Nrows, int Ncols1, int Ncols2, bool show_progress) {
    size_t N=xavier::nrrd_utils::nrrd_size(nin1, true);
    T* out_array=(T*)malloc(N*Ncols2*Nrows*sizeof(T));

    if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeFloat)
        detail::do_product<T, float, float>(out_array, (float*)nin1->data,
            (float*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeFloat)
        detail::do_product<T, double, float>(out_array, (double*)nin1->data,
            (float*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeDouble)
        detail::do_product<T, float, double>(out_array, (float*)nin1->data,
            (double*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeDouble)
        detail::do_product<T, double, double>(out_array, (double*)nin1->data,
            (double*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin1->type)+" or "+std::to_string(nin2->type));
    }

    // size_t dims[4];
    // dims[0]=Nrows*Ncols2;
    // for (int i=1; i<nin1->dim; ++i) dims[i]=nin1->axis[i].size;
    // Nrrd *nout = nrrdNew();
    // if (nrrdWrap_nva(nout, out_array, nrrd_type_id<T>::id,
    //                  nin1->dim, &dims[0])) {
    //     std::cout << "error in nrrdWrap" << '\n';
    //     throw std::runtime_error("nrrd wrap error in product_no_convert");
    // }
    // return nout;

    return detail::wrap_and_copy_header<T>(nin1, out_array, Nrows*Ncols2);
}

template<typename T>
Nrrd* sum(const Nrrd* nin1, const Nrrd* nin2, int Nrows, int Ncols, bool show_progress) {
    size_t N=xavier::nrrd_utils::nrrd_size(nin1, true);
    T* out_array=(T*)malloc(N*Nrows*Ncols*sizeof(T));

    if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeFloat)
        detail::do_sum<T, float, float>(out_array, (float*)nin1->data, (float*)nin2->data, N, Nrows, Ncols, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeFloat)
        detail::do_sum<T, double, float>(out_array, (double*)nin1->data, (float*)nin2->data, N, Nrows, Ncols, show_progress);
    else if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeDouble)
        detail::do_sum<T, float, double>(out_array, (float*)nin1->data, (double*)nin2->data, N, Nrows, Ncols, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeDouble)
        detail::do_sum<T, double, double>(out_array, (double*)nin1->data, (double*)nin2->data, N, Nrows, Ncols, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin1->type)+" or "+std::to_string(nin2->type));
    }
    return detail::wrap_and_copy_header<T>(nin1, out_array, Nrows*Ncols);
}

template<typename T>
Nrrd* subtraction(const Nrrd* nin1, const Nrrd* nin2, int Nrows, int Ncols, bool show_progress) {
    size_t N=xavier::nrrd_utils::nrrd_size(nin1, true);
    T* out_array=(T*)malloc(N*Nrows*Ncols*sizeof(T));

    if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeFloat)
        detail::do_subtraction<T, float, float>(out_array, (float*)nin1->data, (float*)nin2->data, N, Nrows, Ncols, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeFloat)
        detail::do_subtraction<T, double, float>(out_array, (double*)nin1->data, (float*)nin2->data, N, Nrows, Ncols, show_progress);
    else if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeDouble)
        detail::do_subtraction<T, float, double>(out_array, (float*)nin1->data, (double*)nin2->data, N, Nrows, Ncols, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeDouble)
        detail::do_subtraction<T, double, double>(out_array, (double*)nin1->data, (double*)nin2->data, N, Nrows, Ncols, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin1->type)+" or "+std::to_string(nin2->type));
    }
    return detail::wrap_and_copy_header<T>(nin1, out_array, Nrows*Ncols);
}

template<typename T>
Nrrd* trans_product(const Nrrd* nin1, const Nrrd* nin2,
                    int Nrows, int Ncols1, int Ncols2, bool show_progress) {
    size_t N=xavier::nrrd_utils::nrrd_size(nin1, true);
    T* out_array=(T*)malloc(N*Ncols1*Ncols2*sizeof(T));

    if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeFloat)
        detail::do_trans_product<T, float, float>(out_array, (float*)nin1->data,
            (float*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeFloat)
        detail::do_trans_product<T, double, float>(out_array, (double*)nin1->data,
            (float*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeDouble)
        detail::do_trans_product<T, float, double>(out_array, (float*)nin1->data,
            (double*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeDouble)
        detail::do_trans_product<T, double, double>(out_array, (double*)nin1->data,
            (double*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin1->type)+" or "+std::to_string(nin2->type));
    }

    return detail::wrap_and_copy_header<T>(nin1, out_array, Ncols1*Ncols2);
}

template<typename T>
Nrrd* cauchy_green(const Nrrd* nin, int Nrows, bool show_progress) {
    size_t N=xavier::nrrd_utils::nrrd_size(nin, true);
    T* out_array=(T*)malloc(N*Nrows*Nrows*sizeof(T));

    if (nin->type==nrrdTypeFloat)
        detail::do_cauchy_green<T, float>(out_array, (float*)nin->data, N, Nrows, show_progress);
    else if (nin->type==nrrdTypeDouble)
        detail::do_cauchy_green<T, double>(out_array, (double*)nin->data, N, Nrows, show_progress);
    else {
        throw std::runtime_error("Invalid data type: " + std::to_string(nin->type));
    }
    return detail::wrap_and_copy_header<T>(nin, out_array, Nrows*Nrows);
}

template<typename T>
Nrrd* inverse(const Nrrd* nin, int Nrows, bool show_progress) {
    size_t N=xavier::nrrd_utils::nrrd_size(nin, true);
    T* out_array=(T*)malloc(N*Nrows*Nrows*sizeof(T));

    if (nin->type==nrrdTypeFloat)
        detail::do_inverse<T, float>(out_array, (float*)nin->data, N, Nrows, show_progress);
    else if (nin->type==nrrdTypeDouble)
        detail::do_inverse<T, double>(out_array, (double*)nin->data, N, Nrows, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin->type));
    }

    return detail::wrap_and_copy_header<T>(nin, out_array, Nrows*Nrows);
}

template<typename T>
void SVD(Nrrd*& sings, Nrrd*& left, Nrrd*& right,
         const Nrrd* nin, int Nrows, int Ncols, bool show_progress) {
    size_t N=xavier::nrrd_utils::nrrd_size(nin, true);
    int Ndiag=std::min(Nrows, Ncols);

    T* svals_array=(T*)malloc(N*Ndiag*sizeof(T));
    T* left_array=(T*)malloc(N*Nrows*Ndiag*sizeof(T));
    T* right_array=(T*)malloc(N*Ncols*Ndiag*sizeof(T));

    if (nin->type==nrrdTypeFloat)
        detail::do_SVD<T, float>(svals_array, left_array, right_array,
                         (float*)nin->data, N, Nrows, Ncols, show_progress);
    else if (nin->type==nrrdTypeDouble)
        detail::do_SVD<T, double>(svals_array, left_array, right_array,
                          (double*)nin->data, N, Nrows, Ncols, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin->type));
    }

    size_t dims[4];
    for (int i=1; i<nin->dim; ++i) dims[i]=nin->axis[i].size;

    sings=detail::wrap_and_copy_header<T>(nin, svals_array, Ndiag);

    left=detail::wrap_and_copy_header<T>(nin, left_array, Nrows*Ndiag);

    right=detail::wrap_and_copy_header<T>(nin, right_array, Ncols*Ndiag);
}

template<typename Type_>
void general_Eigendecomposition(Nrrd*& evals, Nrrd*& evecs, const Nrrd* nin,
                               int Nrows, bool show_progress) {
    size_t N=xavier::nrrd_utils::nrrd_size(nin, true);
    typedef std::complex<Type_> value_t;
    typedef Type_ base_t;

    value_t* evals_ptr=(value_t*)malloc(N*Nrows*sizeof(value_t));
    value_t* evecs_ptr=(value_t*)malloc(N*Nrows*Nrows*sizeof(value_t));

    if (nin->type==nrrdTypeFloat)
        detail::do_Eigendecomposition<value_t, float, eigensolver_t<float> >(
                        evals_ptr, evecs_ptr, (float*)nin->data,
                        N, Nrows, false, show_progress);
    else if (nin->type==nrrdTypeDouble)
        detail::do_Eigendecomposition<value_t, double, eigensolver_t<double> >(
                        evals_ptr, evecs_ptr, (double*)nin->data,
                        N, Nrows, false, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin->type));
    }
    size_t dims[4];
    for (int i=1; i<nin->dim; ++i) dims[i]=nin->axis[i].size;

    evals=detail::wrap_and_copy_header<Type_>(nin, reinterpret_cast<Type_*>(evals_ptr), 2*Nrows);

    evecs=detail::wrap_and_copy_header<Type_>(nin, reinterpret_cast<Type_*>(evecs_ptr), Nrows*Nrows*2);
}

template<typename Type_>
void symmetric_Eigendecomposition(Nrrd*& evals, Nrrd*& evecs, const Nrrd* nin,
                                  int Nrows, bool show_progress) {
    size_t N=xavier::nrrd_utils::nrrd_size(nin, true);
    typedef Type_ value_t;

    value_t* evals_ptr=(value_t*)malloc(N*Nrows*sizeof(value_t));
    value_t* evecs_ptr=(value_t*)malloc(N*Nrows*Nrows*sizeof(value_t));

    if (nin->type==nrrdTypeFloat)
        detail::do_Eigendecomposition<value_t, float, sym_eigensolver_t<float> >(
                        evals_ptr, evecs_ptr, (float*)nin->data,
                        N, Nrows, true, show_progress);
    else if (nin->type==nrrdTypeDouble)
        detail::do_Eigendecomposition<value_t, double, sym_eigensolver_t<double> >(
                        evals_ptr, evecs_ptr, (double*)nin->data,
                        N, Nrows, true, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin->type));
    }

    evals=detail::wrap_and_copy_header<Type_>(nin, evals_ptr, Nrows);

    evecs=detail::wrap_and_copy_header<Type_>(nin, evecs_ptr, Nrows*Nrows);
}

template<typename T>
Nrrd* FTLE(const Nrrd* nin, int Nrows, bool show_progress) {
    size_t N=xavier::nrrd_utils::nrrd_size(nin, true);

    T* out_array=(T*)malloc(N*sizeof(T));
    if (nin->type==nrrdTypeFloat) {
        detail::do_FTLE(out_array, (float*)nin->data, N, Nrows, show_progress);
    }
    else if (nin->type==nrrdTypeDouble) {
        detail::do_FTLE(out_array, (double*)nin->data, N, Nrows, show_progress);
    }

    return detail::wrap_and_copy_header<T>(nin, out_array, 1);
}

template<typename T>
void Eigendecomposition(Nrrd*& evals, Nrrd*& evecs,
                        const Nrrd* nin, int Nrows, bool issym,
                        bool show_progress) {

    if (issym)
        symmetric_Eigendecomposition<T>(evals, evecs, nin, Nrows, show_progress);
    else
        general_Eigendecomposition<T>(evals, evecs, nin, Nrows, show_progress);
}

} // nrrd_manip
} // namespace xavier

namespace xavier { namespace nrrd_manip { namespace detail {
    template<typename T>
    Nrrd* mat2nrrd(const std::vector< mat_t<T> >& mats,
                   const std::array<size_t, 3>& res,
                   int Nrows, int Ncols, bool show_progress) {

        size_t N=res[0]*res[1]*res[2];
        int ncoef=Ncols*Nrows;
        assert(N==mats.size());
        T* raster=(T*)calloc(ncoef*N, sizeof(T));
        xavier::ProgressDisplay progress(show_progress);
        progress.start(N, "Matrices to Raster");
        for (size_t i=0; i<N; ++i) {
            for (int j=0; j<ncoef; ++j) {
                raster[i*ncoef+j]=mats[i].data()[j];
            }
            progress.update(i);
        }
        progress.end();
        std::vector<size_t> dims(4);
        dims[0]=ncoef;
        dims[1]=res[0];
        dims[2]=res[1];
        dims[3]=res[2];
        Nrrd *nout = nrrdNew();
        if (nrrdWrap_nva(nout, raster, xavier::nrrd_utils::nrrd_value_traits_from_type<T>::index,
                         4, &dims[0])) {
            throw std::runtime_error("nrrd wrap error in mat2nrrd");
        }
        return nout;
    }

    template<typename T>
    Nrrd* vec2nrrd(const std::vector< col_t<T> >& vecs,
                   const std::array<size_t, 3>& res,
                   int Nrows, bool show_progress) {

        size_t N=res[0]*res[1]*res[2];
        int ncoef=Nrows;
        assert(N==vecs.size());
        T* raster=(T*)calloc(ncoef*N, sizeof(T));
        xavier::ProgressDisplay progress(show_progress);
        progress.start(N, "Vectors to Raster");
        for (size_t i=0; i<N; ++i) {
            for (int j=0; j<ncoef; ++j) {
                raster[i*ncoef+j]=vecs[i].data()[j];
            }
            progress.update(i);
        }
        progress.end();
        std::vector<size_t> dims(4);
        dims[0]=ncoef;
        dims[1]=res[0];
        dims[2]=res[1];
        dims[3]=res[2];
        Nrrd *nout = nrrdNew();
        if (nrrdWrap_nva(nout, raster, xavier::nrrd_utils::nrrd_value_traits_from_type<T>::index,
                         4, &dims[0])) {
            throw std::runtime_error("nrrd wrap error in vec2nrrd");
        }
        return nout;
    }

    template<typename T>
    void nrrd2mat(const Nrrd* nin, std::vector< mat_t<T> >& mats,
                  std::array<size_t, 3>& res,
                  int Nrows, int Ncols, bool show_progress) {

        assert(Nrows*Ncols==nin->axis[0].size && nin->dim==4);

        for (int i=0; i<3; ++i) {
            res[i]=nin->axis[i+1].size;
        }

        size_t N=res[0]*res[1]*res[2];
        int ncoef=Ncols*Nrows;
        mat_t<T> refval(Nrows, Ncols);
        mats.resize(N, refval);

        xavier::nrrd_utils::nrrd_data_wrapper<T> raster(nin);
        xavier::ProgressDisplay progress(show_progress);
        progress.start(N, "Raster to Matrices");
        for (size_t i=0; i<N; ++i) {
            mats[i].resize(Nrows, Ncols); // does nothing if already set
            for (int j=0; j<ncoef; ++j) {
                mats[i].data()[j]=raster[i*ncoef+j];
            }
            progress.update(i);
        }
        progress.end();
    }

    template<typename T, typename T1>
    void do_transpose(T* out_ptr, T1* in_ptr,
                      size_t N, int Nrows, int Ncols,
                      bool show_progress) {
        typedef Eigen::Map< mat_t<T> > out_map_t;
        typedef Eigen::Map< mat_t<T1> > in_map_t;

        assert(N>0);
        xavier::ProgressDisplay progress(show_progress);
        progress.start(N, "Transpose no alloc");
        int nterms=Nrows*Ncols;
        out_map_t out_matrix(out_ptr, Ncols, Nrows);
        in_map_t in_matrix(in_ptr, Nrows, Ncols);
        for (size_t i=0; i<N; ++i, out_ptr+=nterms, in_ptr+=nterms) {
            new (&out_matrix) out_map_t(out_ptr, Ncols, Nrows);
            new (&in_matrix) in_map_t(in_ptr, Nrows, Ncols);
            out_matrix=in_matrix.transpose().template cast<T>();
            progress.update(i);
        }
        progress.end();
    }

    template<typename T, typename T1, typename T2>
    void do_sum(T* out_ptr, T1* a_ptr, T2* b_ptr,
                size_t N, int Nrows, int Ncols, bool show_progress) {
        typedef Eigen::Map< mat_t< T > > out_map_t;
        typedef Eigen::Map< mat_t< T1 > > a_map_t;
        typedef Eigen::Map< mat_t< T2 > > b_map_t;

        int ncoef=Nrows*Ncols;
        out_map_t sum(out_ptr, Nrows, Ncols);
        a_map_t A(a_ptr, Nrows, Ncols);
        b_map_t B(b_ptr, Nrows, Ncols);

        xavier::ProgressDisplay progress(show_progress);
        progress.start(N, "Sum");
        for (size_t i=0; i<N ; ++i, out_ptr+=ncoef, a_ptr+=ncoef, b_ptr+=ncoef) {
            new (&sum) out_map_t(out_ptr, Nrows, Ncols);
            new (&A) a_map_t(a_ptr, Nrows, Ncols);
            new (&B) b_map_t(b_ptr, Nrows, Ncols);
            sum=A.template cast<T>() + B.template cast<T>();
            progress.update(i);
        }
        progress.end();
    }

    template<typename T, typename T1, typename T2>
    void do_subtraction(T* out_ptr, T1* a_ptr, T2* b_ptr,
                        size_t N, int Nrows, int Ncols, bool show_progress) {
        typedef Eigen::Map< mat_t< T > > out_map_t;
        typedef Eigen::Map< mat_t< T1 > > a_map_t;
        typedef Eigen::Map< mat_t< T2 > > b_map_t;

        int ncoef=Nrows*Ncols;
        out_map_t subtraction(out_ptr, Nrows, Ncols);
        a_map_t A(a_ptr, Nrows, Ncols);
        b_map_t B(b_ptr, Nrows, Ncols);

        xavier::ProgressDisplay progress(show_progress);
        progress.start(N, "Sum");
        for (size_t i=0; i<N ; ++i, out_ptr+=ncoef, a_ptr+=ncoef, b_ptr+=ncoef) {
            new (&subtraction) out_map_t(out_ptr, Nrows, Ncols);
            new (&A) a_map_t(a_ptr, Nrows, Ncols);
            new (&B) b_map_t(b_ptr, Nrows, Ncols);
            subtraction=A.template cast<T>() - B.template cast<T>();
            progress.update(i);
        }
        progress.end();
    }

    template<typename T, typename T1, typename T2>
    void do_product(T* out_ptr, T1* a_ptr, T2* b_ptr, size_t N,
                    int Nrows, int Ncols1, int Ncols2, bool show_progress) {

        typedef Eigen::Map< mat_t< T > > out_map_t;
        typedef Eigen::Map< mat_t< T1 > > a_map_t;
        typedef Eigen::Map< mat_t< T2 > > b_map_t;

        int ncoef_out=Nrows*Ncols2;
        int ncoef_a=Nrows*Ncols1;
        int ncoef_b=Ncols1*Ncols2;

        out_map_t res(out_ptr, Nrows, Ncols2);
        a_map_t A(a_ptr, Nrows, Ncols1);
        b_map_t B(b_ptr, Ncols1, Ncols2);

        xavier::ProgressDisplay progress(show_progress);
        progress.start(N, "Product");
        for (size_t i=0; i<N; ++i, out_ptr+=ncoef_out, a_ptr+=ncoef_a, b_ptr+=ncoef_b) {
            new (&res) out_map_t(out_ptr, Nrows, Ncols2);
            new (&A) a_map_t(a_ptr, Nrows, Ncols1);
            new (&B) b_map_t(b_ptr, Ncols1, Ncols2);

            res=A.template cast<T>()*B.template cast<T>();
            progress.update(i);
        }
        progress.end();
    }

    template<typename T, typename T1, typename T2>
    void do_trans_product(T* out_ptr, T1* a_ptr, T2* b_ptr,
                          size_t N, int Nrows, int Ncols1, int Ncols2,
                          bool show_progress) {
        typedef Eigen::Map< mat_t< T > > out_map_t;
        typedef Eigen::Map< mat_t< T1 > > a_map_t;
        typedef Eigen::Map< mat_t< T2 > > b_map_t;

        int ncoef_out=Ncols1*Ncols2;
        int ncoef_a=Nrows*Ncols1;
        int ncoef_b=Ncols1*Ncols2;

        out_map_t tprod(out_ptr, Ncols1, Ncols2);
        a_map_t A(a_ptr, Nrows, Ncols1);
        b_map_t B(b_ptr, Ncols1, Ncols2);

        xavier::ProgressDisplay progress(show_progress);
        progress.start(N, "Transpose+Product");
        for (size_t i=0; i<N; ++i, out_ptr+=ncoef_out, a_ptr+=ncoef_a, b_ptr+=ncoef_b) {
            new (&tprod) out_map_t(out_ptr, Ncols1, Ncols2);
            new (&A) a_map_t(a_ptr, Nrows, Ncols1);
            new (&B) b_map_t(b_ptr, Ncols1, Ncols2);

            tprod=A.transpose().template cast<T>()*B.template cast<T>();
            progress.update(i);
        }
        progress.end();
    }

    template<typename T, typename T1>
    void do_cauchy_green(T* out_ptr, T1* in_ptr, size_t N, int Nrows,
                         bool show_progress) {

        typedef Eigen::Map< mat_t<T> > out_map_t;
        typedef Eigen::Map< mat_t<T1> > in_map_t;

        int ncoef=Nrows*Nrows;
        out_map_t cg(out_ptr, Nrows, Nrows);
        in_map_t J(in_ptr, Nrows, Nrows);

        xavier::ProgressDisplay progress(show_progress);
        progress.start(N, "Cauchy-Green");
        for (size_t i=0; i<N; ++i, out_ptr+=ncoef, in_ptr+=ncoef) {
            new (&cg) out_map_t(out_ptr, Nrows, Nrows);
            new (&J) in_map_t(in_ptr, Nrows, Nrows);

            cg=(J.transpose()*J).template cast<T>();
            progress.update(i);
        }
        progress.end();
    }

    template<typename T, typename T1>
    void do_SVD(T* svals_ptr, T* left_ptr, T* right_ptr, T1* in_ptr,
                size_t N, int Nrows, int Ncols,
                bool show_progress) {
        typedef Eigen::Map< mat_t< T > > out_matrix_map_t;
        typedef Eigen::Map< mat_t< T1 > > in_matrix_map_t;
        typedef Eigen::Map< col_t< T> > out_vector_map_t;

        int Ndiag=std::min(Nrows, Ncols);
        typedef Eigen::JacobiSVD< mat_t<T1> > svd_t;
        xavier::ProgressDisplay progress(show_progress);

        out_matrix_map_t left(left_ptr, Nrows, Ndiag);
        out_matrix_map_t right(right_ptr, Ncols, Ndiag);
        out_vector_map_t svals(svals_ptr, Ndiag, 1);
        in_matrix_map_t matrix(in_ptr, Nrows, Ncols);

        int ncoef_left=Nrows*Ndiag;
        int ncoef_right=Ncols*Ndiag;
        int ncoef_svals=Ndiag;
        int ncoef_matrix=Nrows*Ncols;

        progress.start(N, "Singular Value Decomposition");
        for (size_t i=0; i<N; ++i, svals_ptr+=ncoef_svals,
                left_ptr+=ncoef_left, right_ptr+=ncoef_right, in_ptr+=ncoef_matrix) {
            new (&left) out_matrix_map_t(left_ptr, Nrows, Ndiag);
            new (&right) out_matrix_map_t(right_ptr, Ncols, Ndiag);
            new (&svals) out_vector_map_t(svals_ptr, Ndiag, 1);
            new (&matrix) in_matrix_map_t(in_ptr, Nrows, Ncols);

            svd_t svd(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
            svals=svd.singularValues().template cast<T>();
            left=svd.matrixU().template cast<T>();
            right=svd.matrixV().template cast<T>();
            progress.update(i);
        }
        progress.end();
    }

    template<typename Type_>
    struct LessThan : std::less<Type_> {};

    template<typename Type_>
    struct LessThan< std::complex<Type_> > {
        typedef std::complex<Type_> value_t;
        bool operator()(const value_t& a, const value_t& b) const {
            if (a.real() < b.real()) return true;
            else if (a.real() > b.real()) return false;
            else return a.imag() < b.imag();
        }
    };

    template< typename Type_, typename Compare_>
    void sort_eigensolution(Type_* eigenvalues, Type_* eigenvectors,
                            size_t N, size_t Nrows, bool show_progress);

    template<typename TypeOut_, typename TypeIn_, typename EigenSolver_>
    void do_Eigendecomposition(TypeOut_* evals_ptr,
                               TypeOut_* evecs_ptr, TypeIn_* in_ptr,
                               size_t N, int Nrows, bool issym,
                               bool show_progress) {
        typedef Eigen::Map< mat_t< TypeOut_ > > out_matrix_map_t;
        typedef Eigen::Map< mat_t< TypeIn_ > > in_matrix_map_t;
        typedef EigenSolver_ solver_t;
        xavier::ProgressDisplay progress(show_progress);

        out_matrix_map_t evals(evals_ptr, Nrows, 1);
        out_matrix_map_t evecs(evecs_ptr, Nrows, Nrows);
        in_matrix_map_t matrix(in_ptr, Nrows, Nrows);

        int ncoef_evals=Nrows;
        int ncoef_evecs=Nrows*Nrows;
        int ncoef_matrix=Nrows*Nrows;

        TypeOut_* _cur_evals=evals_ptr;
        TypeOut_* _cur_evecs=evecs_ptr;
        TypeIn_* _cur_in=in_ptr;

        progress.start(N, "Eigendecomposition");
        for (size_t i=0; i<N; ++i, _cur_evals+=ncoef_evals,
                _cur_evecs+=ncoef_evecs, _cur_in+=ncoef_matrix) {
            new (&evecs) out_matrix_map_t(_cur_evecs, Nrows, Nrows);
            new (&evals) out_matrix_map_t(_cur_evals, Nrows, 1);
            new (&matrix) in_matrix_map_t(_cur_in, Nrows, Nrows);
            solver_t eigen(matrix);
            evals=eigen.eigenvalues().template cast<TypeOut_>();
            evecs=eigen.eigenvectors().template cast<TypeOut_>();
            progress.update(i);
        }
        progress.end();

        sort_eigensolution<TypeOut_, LessThan<TypeOut_> >(evals_ptr, evecs_ptr, N, Nrows, show_progress);
    }

    template<typename T, typename T1>
    void do_FTLE(T* out_ptr, T1* in_ptr, size_t N, int Nrows, bool show_progress) {
        typedef Eigen::Map< mat_t< T1 > > in_map_t;
        typedef mat_t< T > mat_t;

        typedef sym_eigensolver_t<T> solver_t;

        int ncoef=Nrows*Nrows;

        T* ftle=out_ptr;

        mat_t cg;

        xavier::ProgressDisplay progress(show_progress);
        progress.start(N, "FTLE computation");
        for (size_t i=0; i<N ; ++i, ++ftle, in_ptr+=ncoef) {
            in_map_t A(in_ptr, Nrows, Nrows);
            cg=(A.transpose()*A).template cast<T>();
            solver_t eigen(cg);
            *ftle=eigen.eigenvalues().maxCoeff();

            progress.update(i);
        }
        progress.end();
    }

    template<typename T, typename T1>
    void do_inverse(T* out_ptr, T1* in_ptr, size_t N,
                    int Nrows, bool show_progress) {
        typedef Eigen::Map< mat_t<T> > out_map_t;
        typedef Eigen::Map< mat_t<T1> > in_map_t;

        xavier::ProgressDisplay progress(show_progress);

        int ncoef=Nrows*Nrows;
        out_map_t inv(out_ptr, Nrows, Nrows);
        in_map_t mat(in_ptr, Nrows, Nrows);
        progress.start(N, "Inverse");
        for (size_t i=0; i<N; ++i, out_ptr+=ncoef, in_ptr+=ncoef) {
            new (&inv) out_map_t(out_ptr, Nrows, Nrows);
            new (&mat) in_map_t(in_ptr, Nrows, Nrows);
            inv=mat.lu().inverse().template cast<T>();
            progress.update(i);
        }
        progress.end();
    }

    template<typename Type_>
    inline void swap_arrays(Type_* array0, Type_* array1, Type_* tmp, size_t stride) {
        // tmp=array0
        std::memcpy(tmp, array0, stride*sizeof(Type_));
        // array0=array1
        std::memcpy(array0, array1, stride*sizeof(Type_));
        // array1=tmp
        std::memcpy(array1, tmp, stride*sizeof(Type_));
    }

    template<typename Type_>
    void reorder(Type_* array, const std::vector<size_t>& order, const size_t& stride)  {
        size_t nleft=order.size(); // # elements left to process
        Type_* tmp=(Type_*)calloc(stride, sizeof(Type_));
        Type_* swap=(Type_*)calloc(stride, sizeof(Type_));
        const size_t size_in_bytes=stride*sizeof(Type_);
        // std::cout << "stride=" << stride << ", size_in_bytes=" << size_in_bytes << '\n';
        for (size_t s=0, d; nleft>0; ++s) { // iterate over all elements
            // std::cout << "s=" << s << '\n';
            for (d=order[s]; d>s; d=order[d]) ; // identify circular swap chain
            if (d==s) { // chain is closed
                --nleft;
                // exchange all chain elements circularly
                // tmp=v[s];
                // std::cout << "tmp=array[" << s << "]=" << array[s*stride] << '\n';
                std::memcpy(tmp, array+s*stride, size_in_bytes);
                while (d=order[d], d!=s) {
                    // std::swap(v[d], tmp)
                    // std::cout << "array[" << d << "]=" << array[d*stride] << " <-> tmp=" << *tmp << "\n";
                    swap_arrays<Type_>(tmp, array+d*stride, swap, stride);--nleft;
                }
                // v[s]=tmp
                // std::cout << "array[" << s << "]=" << array[s*stride] << "=tmp=" << *tmp << '\n';
                std::memcpy(array+s*stride, tmp, size_in_bytes);
            }
        }
    }

    template<typename Type_, typename Compare_>
    bool check_sorted(Type_* array) {
        Compare_ cmp;
        if (cmp(array[0], array[1]) && cmp(array[1], array[2])) return true;

        std::cout << "ERROR: " << array[0] << ">=" << array[1] << " or "
            << array[1] << ">=" << array[2] << '\n';
        return false;
    }

    template< typename Type_, typename Compare_>
    void sort_eigensolution(Type_* eigenvalues, Type_* eigenvectors,
                            size_t N, size_t Nrows,
                            bool show_progress) {
        // eigenvalues is [Nrows x N] array
        // eigenvectors is [Nrows x Nrows x N] array
        std::vector<size_t> order(Nrows);
        Type_* evals_ptr=eigenvalues;
        Type_* evecs_ptr=eigenvectors;
        size_t stride=Nrows*Nrows;
        std::map<Type_, int, Compare_> order_map;
        xavier::ProgressDisplay progress(show_progress);

        progress.start(N, "Sort eigenvalues and eigenvectors");
        for (size_t n=0; n<N; ++n, evals_ptr+=Nrows, evecs_ptr+=stride) {
            // determine order of Nrows eigenvalues
            order_map.clear();
            for (int i=0; i<Nrows; ++i) {
                order_map[evals_ptr[i]]=i;
            }
            // reorder eigenvalues and eigenvectors accordingly
            std::vector<size_t> order(Nrows);
            int i=0;
            for (auto it=order_map.begin(); it!=order_map.end(); ++it, ++i) {
                order[it->second]=i;
            }
            /*std::vector<Type_> __evals;
            for (int i=0; i<Nrows; ++i) {
                __evals.push_back(evals_ptr[i]);
            }*/
            reorder<Type_>(evals_ptr, order, 1);
            reorder<Type_>(evecs_ptr, order, Nrows);

            /*if (!check_sorted<Type_, Compare_>(evals_ptr)) {
                std::cout << "ERROR: map ordering was:\n";
                for (auto it=order_map.begin(); it!=order_map.end(); ++it) {
                    std::cout << it->first << " -> " << it->second << '\n';
                }
                std::cout << "corresponding order array:\n";
                std::copy(order.begin(), order.end(), std::ostream_iterator<size_t>(std::cout, ", "));
                std::cout << "\n";
                std::vector<Type_> __evals_after;
                for (int i=0; i<Nrows; ++i) {
                    __evals_after.push_back(evals_ptr[i]);
                }
                std::cout << "before sorting, eigenvalue array is\n";
                std::copy(__evals.begin(), __evals.end(), std::ostream_iterator<Type_>(std::cout, ", "));
                std::cout << "\nafter sorting, eigenvalue array is\n";
                std::copy(__evals_after.begin(), __evals_after.end(), std::ostream_iterator<Type_>(std::cout, ", "));
                std::cout << "\n";


                throw std::runtime_error("invalid ordering in sort_eigensolutions");
            }*/
            progress.update(n);
        }
        progress.end();
    }

    template<typename T>
    Nrrd* wrap_and_copy_header(const Nrrd* src, T* data, int valdim) {
        Nrrd *nout = nrrdNew();
        int dim=(valdim<=1 ? src->dim-1 : src->dim);
        size_t sizes[src->dim];
        if (valdim>1) {
            sizes[0]=valdim;
            for (int i=1; i<src->dim; ++i) sizes[i]=src->axis[i].size;
        }
        else {
            for (int i=0; i<src->dim-1; ++i) sizes[i]=src->axis[i+1].size;
        }
        if (nrrdWrap_nva(nout, data, xavier::nrrd_utils::nrrd_value_traits_from_type<T>::index, dim, sizes)) {
            throw std::runtime_error(xavier::nrrd_utils::error_msg("wrap_and_copy_header"));
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
        if (valdim>1) nrrdSpaceOriginSet(nout, d_array);
        else nrrdSpaceOriginSet(nout, &d_array[1]);

        return nout;
    }

    template<typename T>
    void copyAxisInfo(Nrrd* nout, const Nrrd* nin, int info) {
        T attr[NRRD_DIM_MAX];
        nrrdAxisInfoGet_nva(nin, info, attr);
        if (nin->dim==nout->dim)
            nrrdAxisInfoSet_nva(nout, info, attr);
        else
            nrrdAxisInfoSet_nva(nout, info, &attr[1]);
    }

} // detail
} // nrrd_manip
} // xavier
