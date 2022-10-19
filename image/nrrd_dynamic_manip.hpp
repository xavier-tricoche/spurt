#include <array>
#include <iomanip>

#include <image/nrrd_manip.hpp>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <Eigen/Dense>

#include <util/timer.hpp>

#include <misc/progress.hpp>


namespace spurt { namespace nrrd_manip {
    
template<typename T>
using mat_t=Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
template<typename T>
using col_t=Eigen::Matrix<T, Eigen::Dynamic, 1>;
template<typename T>
using row_t=Eigen::Matrix<T, 1, Eigen::Dynamic>;

template<typename T>
Nrrd* mat2nrrd(const std::vector< mat_t<T> >& mats,
               const std::array<size_t, 3>& res, 
               int Nrows, int Ncols, bool show_progress=false) {

    size_t N=res[0]*res[1]*res[2];
    int ncoef=Ncols*Nrows;
    assert(N==mats.size());
    T* raster=(T*)calloc(ncoef*N, sizeof(T));
    spurt::ProgressDisplay progress(show_progress);
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
    if (nrrdWrap_nva(nout, raster, nrrd_type_id<T>::id, 
                     4, &dims[0])) {
        throw std::runtime_error("nrrd wrap error in mat2nrrd");
    }
    return nout;
}

template<typename T>
Nrrd* vec2nrrd(const std::vector< col_t<T> >& vecs,
               const std::array<size_t, 3>& res, 
               int Nrows, bool show_progress=false) {

    size_t N=res[0]*res[1]*res[2];
    int ncoef=Nrows;
    assert(N==vecs.size());
    T* raster=(T*)calloc(ncoef*N, sizeof(T));
    spurt::ProgressDisplay progress(show_progress);
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
    if (nrrdWrap_nva(nout, raster, nrrd_type_id<T>::id, 
                     4, &dims[0])) {
        throw std::runtime_error("nrrd wrap error in vec2nrrd");
    }
    return nout;
}

template<typename T>
void nrrd2mat(const Nrrd* nin, std::vector< mat_t<T> >& mats,
              std::array<size_t, 3>& res, 
              int Nrows, int Ncols, bool show_progress=false) {

    assert(Nrows*Ncols==nin->axis[0].size && nin->dim==4);
    
    for (int i=0; i<3; ++i) {
        res[i]=nin->axis[i+1].size;
    }

    size_t N=res[0]*res[1]*res[2];
    int ncoef=Ncols*Nrows;
    mat_t<T> refval(Nrows, Ncols);
    mats.resize(N, refval);
    
    nrrd_data_wrapper<T> raster(nin);
    spurt::ProgressDisplay progress(show_progress);
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

template<typename T>
void compute_transpose(std::vector< mat_t<T> >& res,
                       const std::vector< mat_t<T> >& As, bool show_progress=false) {
    assert(!As.empty());
    mat_t<T> refval(As[0].cols(), As[0].rows());
    res.resize(As.size(), refval);
    spurt::ProgressDisplay progress(show_progress);
    progress.start(As.size(), "Transpose");
    for (size_t i=0; i<As.size(); ++i) {
        res[i]=As[i].transpose();
        progress.update(i);
    }
    progress.end();
}

template<typename T, typename T1=T>
void compute_transpose_no_copy(T* out_ptr, T1* in_ptr,
                               size_t N, int Nrows, int Ncols, 
                               bool show_progress=false) {
    assert(N>0);
    spurt::ProgressDisplay progress(show_progress);
    progress.start(N, "Transpose no alloc");
    int nterms=Nrows*Ncols;
    T* out_arr=out_ptr;
    T1* in_arr=in_ptr;
    std::vector<int> offsets(nterms);
    for (int i=0; i<nterms; ++i) {
        int m=i/Ncols;
        int n=i%Ncols;
        offsets[i]=m+Nrows*n; 
    }
    std::copy(offsets.begin(), offsets.end(), std::ostream_iterator<int>(std::cout, ", "));
    std::cout << "\n";
    for (size_t i=0; i<N; ++i, out_arr+=nterms, in_arr+=nterms) {
        for (int j=0; j<nterms; ++j) {
            out_arr[j]=in_arr[offsets[j]];
        }
        progress.update(i);
    }
    progress.end();
}

template<typename T, typename T1=T>
void compute_transpose_map(T* out_ptr, T1* in_ptr,
                           size_t N, int Nrows, int Ncols, 
                           bool show_progress=false) {
    typedef Eigen::Map< mat_t<T> > out_map_t;
    typedef Eigen::Map< mat_t<T1> > in_map_t;                           
                               
    assert(N>0);
    spurt::ProgressDisplay progress(show_progress);
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

template<typename T>
void compute_SVD(std::vector< col_t<T> >& values,
                 std::vector< mat_t<T> >& leftvs,
                 std::vector< mat_t<T> >& rightvs,
                 const std::vector< mat_t<T> >& As, bool show_progress=false) {
    assert(!As.empty());
    int Nrows=As[0].rows();
    int Ncols=As[0].cols();
    int Ndiag=std::min(Nrows, Ncols);
    col_t<T> colref(Ndiag);
    mat_t<T> leftref(Nrows, Nrows);
    mat_t<T> rightref(Ncols, Ncols);
    values.resize(As.size(), colref);
    leftvs.resize(As.size(), leftref);
    rightvs.resize(As.size(), rightref);
    typedef Eigen::JacobiSVD< mat_t<T> > svd_t;
    spurt::ProgressDisplay progress(show_progress);
    progress.start(As.size(), "Singular Value Decomposition");
    for (size_t i=0; i<As.size(); ++i) {
        svd_t svd(As[i], Eigen::ComputeFullU | Eigen::ComputeFullV);
        values[i]=svd.singularValues();
        leftvs[i]=svd.matrixU();
        rightvs[i]=svd.matrixV();
        progress.update(i);
    }
    progress.end();
}

template<typename T, typename T1=T>
void compute_SVD_map(T* svals_ptr, T* left_ptr, T* right_ptr, T1* in_ptr,
                     size_t N, int Nrows, int Ncols, 
                     bool show_progress=false) {
    typedef Eigen::Map< mat_t< T > > out_matrix_map_t;
    typedef Eigen::Map< mat_t< T1 > > in_matrix_map_t;     
    typedef Eigen::Map< col_t< T> > out_vector_map_t;            
                         
    int Ndiag=std::min(Nrows, Ncols);
    typedef Eigen::JacobiSVD< mat_t<T1> > svd_t;
    spurt::ProgressDisplay progress(show_progress);
    
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

template<typename T>
void compute_inverse(std::vector< mat_t<T> >& res,
                     const std::vector<mat_t<T> >& As, bool show_progress=false) {
    assert(!As.empty() && As[0].rows()==As[0].cols());
    res.resize(As.size(), As[0]);
    spurt::ProgressDisplay progress(show_progress);
    progress.start(As.size(), "Inverse");
    
    for (size_t i=0; i<As.size(); ++i) {
        res[i]=As[i].lu().inverse();
        progress.update(i);
    } 
    progress.end();
}

template<typename T, typename T1=T>
void compute_inverse_map(T* out_ptr, T1* in_ptr, size_t N, int Nrows, bool show_progress=false) {
    typedef Eigen::Map< mat_t<T> > out_map_t;
    typedef Eigen::Map< mat_t<T1> > in_map_t;

    spurt::ProgressDisplay progress(show_progress);
    
    int ncoef=Nrows*Nrows;
    out_map_t inv(out_ptr, Nrows, Nrows);
    in_map_t mat(in_ptr, Nrows, Nrows);
    progress.start(N, "Sum");
    for (size_t i=0; i<N; ++i, out_ptr+=ncoef, in_ptr+=ncoef) {
        new (&inv) out_map_t(out_ptr, Nrows, Nrows);
        new (&mat) in_map_t(in_ptr, Nrows, Nrows);
        inv=mat.lu().inverse().template cast<T>();
        
        progress.update(i);
    }
    progress.end();
}

template<typename T>
void compute_sum(std::vector< mat_t<T> >& res,
                 const std::vector< mat_t<T> >& As,
                 const std::vector< mat_t<T> >& Bs, bool show_progress=false) {
    assert(!As.empty() && As.size()==Bs.size());
    res.resize(As.size(), As[0]);
    spurt::ProgressDisplay progress(show_progress);
    progress.start(As.size(), "Sum");
    for (size_t i=0; i<As.size(); ++i) {
        res[i]=As[i]+Bs[i];
        progress.update(i);
    } 
    progress.end();
}

template<typename T, typename T1=T, typename T2=T1>
void compute_sum_map(T* out_ptr, T1* a_ptr, T2* b_ptr, 
                     size_t N, int Nrows, int Ncols, bool show_progress=false) {
    typedef Eigen::Map< mat_t< T > > out_map_t;
    typedef Eigen::Map< mat_t< T1 > > a_map_t;
    typedef Eigen::Map< mat_t< T2 > > b_map_t;
    
    int ncoef=Nrows*Ncols;
    out_map_t sum(out_ptr, Nrows, Ncols);
    a_map_t A(a_ptr, Nrows, Ncols);
    b_map_t B(b_ptr, Nrows, Ncols);
    
    spurt::ProgressDisplay progress(show_progress);
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

template<typename T>
void compute_product(std::vector< mat_t<T> >& res,
                 const std::vector< mat_t<T> >& As,
                 const std::vector< mat_t<T> >& Bs, bool show_progress=false) {
    assert(As.size()==Bs.size());
    res.resize(As.size(), mat_t<T>(As[0].rows(), Bs[0].cols()));
    spurt::ProgressDisplay progress(show_progress);
    progress.start(As.size(), "Product");
    for (size_t i=0; i<As.size(); ++i) {
        res[i]=As[i]*Bs[i];
        progress.update(i);
    }
    progress.end();
}

template<typename T, typename T1=T, typename T2=T1>
void compute_product_no_conversion(T* out_ptr, T1* a_ptr, T2* b_ptr, size_t N,
                                   int Nrows, int Ncols1, int Ncols2, bool show_progress=false) {
    // assert(As.size()==Bs.size());
    // res.resize(As.size(), mat_t<T>(As[0].rows(), Bs[0].cols()));
    spurt::ProgressDisplay progress(show_progress);
    progress.start(N, "Product");
    int ncoef_out=Nrows*Ncols2;
    int ncoef_a=Nrows*Ncols1;
    int ncoef_b=Ncols1*Ncols2;
    for (size_t i=0; i<N; ++i, out_ptr+=ncoef_out, a_ptr+=ncoef_a, b_ptr+=ncoef_b) {
        for (int r=0; r<Nrows ; ++r) {
            for (int c=0; c<Ncols2; ++c) {
                out_ptr[c*Nrows+r]=0;
                for (int l=0; l<Ncols1; ++l) {
                    out_ptr[c*Nrows+r]+=a_ptr[r+l*Nrows]*b_ptr[l+c*Ncols2];
                }
            }
        }
        progress.update(i);
    }
    progress.end();
}

template<typename T, typename T1=T, typename T2=T1>
void compute_product_map(T* out_ptr, T1* a_ptr, T2* b_ptr, size_t N,
                         int Nrows, int Ncols1, int Ncols2, bool show_progress=false) {

    typedef Eigen::Map< mat_t< T > > out_map_t;
    typedef Eigen::Map< mat_t< T1 > > a_map_t;
    typedef Eigen::Map< mat_t< T2 > > b_map_t; 
                             
    int ncoef_out=Nrows*Ncols2;
    int ncoef_a=Nrows*Ncols1;
    int ncoef_b=Ncols1*Ncols2;
    
    out_map_t res(out_ptr, Nrows, Ncols2);
    a_map_t A(a_ptr, Nrows, Ncols1);
    b_map_t B(b_ptr, Ncols1, Ncols2);
                             
    spurt::ProgressDisplay progress(show_progress);
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

template<typename T>
void compute_tprod(std::vector< mat_t<T> >& res,
                   const std::vector< mat_t<T> >& As,
                   const std::vector< mat_t<T> >& Bs, bool show_progress) {
    assert(As.size()==Bs.size());
    res.resize(As.size(), mat_t<T>(As[0].cols(), Bs[0].cols()));
    spurt::ProgressDisplay progress(show_progress);
    progress.start(As.size(), "Transpose+Product");
    for (size_t i=0; i<As.size(); ++i) {
        res[i]=As[i].transpose()*Bs[i];
        progress.update(i);
    }
    progress.end();
}

template<typename T, typename T1=T, typename T2=T1>
void compute_tprod_map(T* out_ptr, T1* a_ptr, T2* b_ptr, 
                       size_t N, int Nrows, int Ncols1, int Ncols2, bool show_progress) {
    typedef Eigen::Map< mat_t< T > > out_map_t;
    typedef Eigen::Map< mat_t< T1 > > a_map_t;
    typedef Eigen::Map< mat_t< T2 > > b_map_t;
    
    int ncoef_out=Ncols1*Ncols2;
    int ncoef_a=Nrows*Ncols1;
    int ncoef_b=Ncols1*Ncols2;
    
    out_map_t tprod(out_ptr, Ncols1, Ncols2);
    a_map_t A(a_ptr, Nrows, Ncols1);
    b_map_t B(b_ptr, Ncols1, Ncols2);

    spurt::ProgressDisplay progress(show_progress);
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


template<typename T>
Nrrd* transpose(const Nrrd* nin, int Nrows, int Ncols, bool show_progress=false) {

    size_t N=spurt::nrrd_size(nin, true);
    std::vector< mat_t<T> > mats;
    std::array<size_t, 3> dims;
    nrrd2mat<T>(nin, mats, dims, Nrows, Ncols, show_progress);
    
    std::vector< mat_t<T> > res;
    compute_transpose<T>(res, mats, show_progress);
    
    return mat2nrrd<T>(res, dims, Ncols, Nrows, show_progress);
}

template<typename T>
Nrrd* transpose_no_conversion(const Nrrd* nin, int Nrows, int Ncols, bool show_progress=false) {

    size_t N=spurt::nrrd_size(nin, true);
    T* out_array=(T*)malloc(N*Ncols*Nrows*sizeof(T));
    
    if (nin->type==nrrdTypeFloat) 
        compute_transpose_no_copy<T, float>(out_array, (float*)nin->data, N, Nrows, Ncols, show_progress);
    else if (nin->type==nrrdTypeDouble) 
        compute_transpose_no_copy<T, double>(out_array, (double*)nin->data, N, Nrows, Ncols, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin->type));
    }
    
    size_t dims[4];
    dims[0]=Nrows*Ncols;
    for (int i=1; i<nin->dim; ++i) dims[i]=nin->axis[i].size;
    Nrrd *nout = nrrdNew();
    if (nrrdWrap_nva(nout, out_array, nrrd_type_id<T>::id, 
                     nin->dim, &dims[0])) {
        std::cout << "error in nrrdWrap" << '\n';
        throw std::runtime_error("nrrd wrap error in transpose_no_convert");
    }
    return nout;
}

template<typename T>
Nrrd* transpose_map(const Nrrd* nin, int Nrows, int Ncols, bool show_progress=false) {

    size_t N=spurt::nrrd_size(nin, true);
    T* out_array=(T*)malloc(N*Ncols*Nrows*sizeof(T));
    
    if (nin->type==nrrdTypeFloat) 
        compute_transpose_map<T, float>(out_array, (float*)nin->data, N, Nrows, Ncols, show_progress);
    else if (nin->type==nrrdTypeDouble) 
        compute_transpose_map<T, double>(out_array, (double*)nin->data, N, Nrows, Ncols, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin->type));
    }
    
    size_t dims[4];
    dims[0]=Nrows*Ncols;
    for (int i=1; i<nin->dim; ++i) dims[i]=nin->axis[i].size;
    Nrrd *nout = nrrdNew();
    if (nrrdWrap_nva(nout, out_array, nrrd_type_id<T>::id, 
                     nin->dim, &dims[0])) {
        std::cout << "error in nrrdWrap" << '\n';
        throw std::runtime_error("nrrd wrap error in transpose_no_convert");
    }
    return nout;
}

template<typename T>
Nrrd* product(const Nrrd* nin1, const Nrrd* nin2, int Nrows, int Ncols1, int Ncols2, bool show_progress=false) {
    std::vector< mat_t<T> > mats1, mats2;
    std::array<size_t, 3> dims;
    nrrd2mat<T>(nin1, mats1, dims, Nrows, Ncols1, show_progress);
    nrrd2mat<T>(nin2, mats2, dims, Ncols1, Ncols2, show_progress);
    
    std::vector< mat_t<T> > res;
    compute_product<T>(res, mats1, mats2, show_progress);
    
    return mat2nrrd<T>(res, dims, Nrows, Ncols2, show_progress);
}

template<typename T>
Nrrd* product_no_conversion(const Nrrd* nin1, const Nrrd* nin2, int Nrows, int Ncols1, int Ncols2, bool show_progress=false) {
    size_t N=spurt::nrrd_size(nin1, true);
    T* out_array=(T*)malloc(N*Ncols2*Nrows*sizeof(T));
    
    if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeFloat) 
        compute_product_no_conversion<T, float, float>(out_array, (float*)nin1->data, (float*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeFloat) 
        compute_product_no_conversion<T, double, float>(out_array, (double*)nin1->data, (float*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeDouble) 
        compute_product_no_conversion<T, float, double>(out_array, (float*)nin1->data, (double*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeDouble) 
        compute_product_no_conversion<T, double, double>(out_array, (double*)nin1->data, (double*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin1->type)+" or "+std::to_string(nin2->type));
    }
    
    size_t dims[4];
    dims[0]=Nrows*Ncols2;
    for (int i=1; i<nin1->dim; ++i) dims[i]=nin1->axis[i].size;
    Nrrd *nout = nrrdNew();
    if (nrrdWrap_nva(nout, out_array, nrrd_type_id<T>::id, 
                     nin1->dim, &dims[0])) {
        std::cout << "error in nrrdWrap" << '\n';
        throw std::runtime_error("nrrd wrap error in product_no_convert");
    }
    return nout;
}

template<typename T>
Nrrd* product_map(const Nrrd* nin1, const Nrrd* nin2, int Nrows, int Ncols1, int Ncols2, bool show_progress=false) {
    size_t N=spurt::nrrd_size(nin1, true);
    T* out_array=(T*)malloc(N*Ncols2*Nrows*sizeof(T));
    
    if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeFloat) 
        compute_product_map<T, float, float>(out_array, (float*)nin1->data, (float*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeFloat) 
        compute_product_map<T, double, float>(out_array, (double*)nin1->data, (float*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeDouble) 
        compute_product_map<T, float, double>(out_array, (float*)nin1->data, (double*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeDouble) 
        compute_product_map<T, double, double>(out_array, (double*)nin1->data, (double*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin1->type)+" or "+std::to_string(nin2->type));
    }
    
    size_t dims[4];
    dims[0]=Nrows*Ncols2;
    for (int i=1; i<nin1->dim; ++i) dims[i]=nin1->axis[i].size;
    Nrrd *nout = nrrdNew();
    if (nrrdWrap_nva(nout, out_array, nrrd_type_id<T>::id, 
                     nin1->dim, &dims[0])) {
        std::cout << "error in nrrdWrap" << '\n';
        throw std::runtime_error("nrrd wrap error in product_no_convert");
    }
    return nout;
}

template<typename T>
Nrrd* sum_map(const Nrrd* nin1, const Nrrd* nin2, int Nrows, int Ncols, bool show_progress=false) {
    size_t N=spurt::nrrd_size(nin1, true);
    T* out_array=(T*)malloc(N*Nrows*Ncols*sizeof(T));
    
    if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeFloat) 
        compute_sum_map<T, float, float>(out_array, (float*)nin1->data, (float*)nin2->data, N, Nrows, Ncols, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeFloat) 
        compute_sum_map<T, double, float>(out_array, (double*)nin1->data, (float*)nin2->data, N, Nrows, Ncols, show_progress);
    else if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeDouble) 
        compute_sum_map<T, float, double>(out_array, (float*)nin1->data, (double*)nin2->data, N, Nrows, Ncols, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeDouble) 
        compute_sum_map<T, double, double>(out_array, (double*)nin1->data, (double*)nin2->data, N, Nrows, Ncols, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin1->type)+" or "+std::to_string(nin2->type));
    }
    
    size_t dims[4];
    dims[0]=Nrows*Ncols;
    for (int i=1; i<nin1->dim; ++i) dims[i]=nin1->axis[i].size;
    Nrrd *nout = nrrdNew();
    if (nrrdWrap_nva(nout, out_array, nrrd_type_id<T>::id, 
                     nin1->dim, &dims[0])) {
        std::cout << "error in nrrdWrap" << '\n';
        throw std::runtime_error("nrrd wrap error in sum_map");
    }
    return nout;
}


template<typename T>
Nrrd* trans_product(const Nrrd* nin1, const Nrrd* nin2, int Nrows, int Ncols1, int Ncols2, bool show_progress=false) {
    std::vector< mat_t<T> > mats1, mats2;
    std::array<size_t, 3> dims;
    nrrd2mat<T>(nin1, mats1, dims, Nrows, Ncols1, show_progress);
    nrrd2mat<T>(nin2, mats2, dims, Ncols1, Ncols2, show_progress);
    
    std::vector< mat_t<T> > res;
    compute_tprod<T>(res, mats1, mats2, show_progress);
    
    return mat2nrrd<T>(res, dims, Ncols1, Ncols2, show_progress);
}

template<typename T>
Nrrd* trans_product_map(const Nrrd* nin1, const Nrrd* nin2, 
                        int Nrows, int Ncols1, int Ncols2, bool show_progress=false) {
    size_t N=spurt::nrrd_size(nin1, true);
    T* out_array=(T*)malloc(N*Ncols1*Ncols2*sizeof(T));
    
    if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeFloat) 
        compute_tprod_map<T, float, float>(out_array, (float*)nin1->data, (float*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeFloat) 
        compute_tprod_map<T, double, float>(out_array, (double*)nin1->data, (float*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeFloat && nin2->type==nrrdTypeDouble) 
        compute_tprod_map<T, float, double>(out_array, (float*)nin1->data, (double*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else if (nin1->type==nrrdTypeDouble && nin2->type==nrrdTypeDouble) 
        compute_tprod_map<T, double, double>(out_array, (double*)nin1->data, (double*)nin2->data, N, Nrows, Ncols1, Ncols2, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin1->type)+" or "+std::to_string(nin2->type));
    }
    
    size_t dims[4];
    dims[0]=Ncols1*Ncols2;
    for (int i=1; i<nin1->dim; ++i) dims[i]=nin1->axis[i].size;
    Nrrd *nout = nrrdNew();
    if (nrrdWrap_nva(nout, out_array, nrrd_type_id<T>::id, 
                     nin1->dim, &dims[0])) {
        std::cout << "error in nrrdWrap" << '\n';
        throw std::runtime_error("nrrd wrap error in trans_product_map");
    }
    return nout;
}

template<typename T>
Nrrd* invert(const Nrrd* nin, int N, bool show_progress=false) {
    std::vector< mat_t<T> > mats;
    std::array<size_t, 3> dims;
    nrrd2mat<T>(nin, mats, dims, N, N);
    
    std::vector< mat_t<T> > res;
    compute_inverse<T>(res, mats, show_progress);
    
    return mat2nrrd<T>(res, dims, N, N);
}

template<typename T>
Nrrd* invert_map(const Nrrd* nin, int Nrows, bool show_progress=false) {
    size_t N=spurt::nrrd_size(nin, true);
    T* out_array=(T*)malloc(N*Nrows*Nrows*sizeof(T));
    
    if (nin->type==nrrdTypeFloat) 
        compute_inverse_map<T, float>(out_array, (float*)nin->data, N, Nrows, show_progress);
    else if (nin->type==nrrdTypeDouble) 
        compute_inverse_map<T, double>(out_array, (double*)nin->data, N, Nrows, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin->type));
    }
    
    size_t dims[4];
    dims[0]=Nrows*Nrows;
    for (int i=1; i<nin->dim; ++i) dims[i]=nin->axis[i].size;
    Nrrd *nout = nrrdNew();
    if (nrrdWrap_nva(nout, out_array, nrrd_type_id<T>::id, 
                     nin->dim, &dims[0])) {
        std::cout << "error in nrrdWrap" << '\n';
        throw std::runtime_error("nrrd wrap error in invert_map");
    }
    return nout;

}

template<typename T>
void SVD(Nrrd*& sings, Nrrd*& left, Nrrd*& right, 
         const Nrrd* nin, int Nrows, int Ncols, bool show_progress=false) {
             
    int Ndiag=std::min(Nrows, Ncols);
    std::vector< mat_t<T> > mats;
    std::array<size_t, 3> dims;
    nrrd2mat<T>(nin, mats, dims, Nrows, Ncols, show_progress);
        
    std::vector< mat_t<T> > leftmats, rightmats;
    std::vector< col_t<T> > singvecs;
    
    compute_SVD<T>(singvecs, leftmats, rightmats, mats, show_progress);
    
    sings=vec2nrrd<T>(singvecs, dims, Ndiag);
    left=mat2nrrd<T>(leftmats, dims, Nrows, Nrows, show_progress);
    right=mat2nrrd<T>(rightmats, dims, Ncols, Ncols);
}

template<typename T>
void SVD_map(Nrrd*& sings, Nrrd*& left, Nrrd*& right, 
             const Nrrd* nin, int Nrows, int Ncols, bool show_progress=false) {
    size_t N=spurt::nrrd_size(nin, true);
    int Ndiag=std::min(Nrows, Ncols);
    
    T* svals_array=(T*)malloc(N*Ndiag*sizeof(T));
    T* left_array=(T*)malloc(N*Nrows*Ndiag*sizeof(T));
    T* right_array=(T*)malloc(N*Ncols*Ndiag*sizeof(T));
    
    if (nin->type==nrrdTypeFloat) 
        compute_SVD_map<T, float>(svals_array, left_array, right_array,
                                  (float*)nin->data, N, Nrows, Ncols, show_progress);
    else if (nin->type==nrrdTypeDouble) 
        compute_SVD_map<T, double>(svals_array, left_array, right_array,
                                   (double*)nin->data, N, Nrows, Ncols, show_progress);
    else {
        throw std::runtime_error("Invalid data type:"+std::to_string(nin->type));
    }
    
    size_t dims[4];
    for (int i=1; i<nin->dim; ++i) dims[i]=nin->axis[i].size;
    
    dims[0]=Ndiag;
    sings=nrrdNew();
    if (nrrdWrap_nva(sings, svals_array, nrrd_type_id<T>::id, 
                     nin->dim, &dims[0])) {
        std::cout << "error in nrrdWrap" << '\n';
        throw std::runtime_error("nrrd wrap error in SVD_map");
    }
    
    dims[0]=Nrows*Ndiag;
    left=nrrdNew();
    if (nrrdWrap_nva(left, left_array, nrrd_type_id<T>::id, 
                     nin->dim, &dims[0])) {
        std::cout << "error in nrrdWrap" << '\n';
        throw std::runtime_error("nrrd wrap error in SVD_map");
    }
    
    dims[0]=Ncols*Ndiag;
    right=nrrdNew();
    if (nrrdWrap_nva(right, right_array, nrrd_type_id<T>::id, 
                     nin->dim, &dims[0])) {
        std::cout << "error in nrrdWrap" << '\n';
        throw std::runtime_error("nrrd wrap error in SVD_map");
    }
}
} // nrrd_manip 
} // namespace spurt
