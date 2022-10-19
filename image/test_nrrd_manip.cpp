#include <array>
#include <iomanip>

#include <image/nrrd_manip.hpp>
#include <Eigen/Core>
#include <Eigen/SVD>

#include <util/timer.hpp>

double toverhead;

template<int M, int N>
double check_same(const Nrrd* nrrd,
                  const std::vector<Eigen::Matrix<float, M, N> >& vec,
                  const std::string& what="") {
    typedef Eigen::Matrix<float, M, N> mat_t;

    size_t nmats=vec.size();
    assert(M*N==nrrd->axis[0].size ||
           (vec[0].cols()*vec[0].rows()==nrrd->axis[0].size));
    assert(nrrd->dim==4);
    assert(nrrd->axis[1].size*nrrd->axis[2].size*nrrd->axis[3].size==nmats);
    float* raw=(float*)nrrd->data;
    mat_t* array=reinterpret_cast<mat_t*>(raw);
    double err=0;
    for (size_t i=0; i<nmats; ++i) {
        err+=(array[i]-vec[i]).norm()/vec[i].norm();
    }
    std::cout << what << std::setprecision(5) << ": average relative error = "
              << 100.*err/(double)nmats << "%\n";
    return err;
}

template<int Nrows, int Ncols>
Nrrd* mat2nrrd(const std::vector<Eigen::Matrix<float, Nrows, Ncols> >& mats,
               const std::array<size_t, 3>& res) {

    size_t N=res[0]*res[1]*res[2];
    int ncoef=Ncols*Nrows;
    assert(N==mats.size());
    float* raster=(float*)calloc(ncoef*N, sizeof(float));
    for (size_t i=0; i<N; ++i) {
        for (int j=0; j<ncoef; ++j) {
            raster[i*ncoef+j]=mats[i].data()[j];
        }
    }
    std::vector<size_t> dims(4);
    dims[0]=ncoef;
    dims[1]=res[0];
    dims[2]=res[1];
    dims[3]=res[2];
    Nrrd *nout = nrrdNew();
    if (nrrdWrap_nva(nout, raster, nrrdTypeFloat, 4, &dims[0])) {
        throw std::runtime_error("nrrd wrap error in mat2nrrd");
    }
    return nout;
}

template<int Nrows, int Ncols, int Order>
void nrrd2mat(const Nrrd* nin,
              std::vector<Eigen::Matrix<float, Nrows, Ncols, Order> >& mats,
              std::array<size_t, 3>& res) {

    assert(Nrows*Ncols==nin->axis[0].size && nin->dim==4);

    for (int i=0; i<3; ++i) {
        res[i]=nin->axis[i+1].size;
    }

    size_t N=res[0]*res[1]*res[2];
    int ncoef=Ncols*Nrows;
    mats.resize(N);

    float* raster=(float*)nin->data;

    for (size_t i=0; i<N; ++i) {
        for (int j=0; j<ncoef; ++j) {
            mats[i].data()[j]=raster[i*ncoef+j];
        }
    }
}

template<int Nrows, int Ncols, int Order1, int Order2>
void compute_transpose(std::vector<Eigen::Matrix<float, Ncols, Nrows, Order1> >& res,
                       const std::vector<Eigen::Matrix<float, Nrows, Ncols, Order2> >& As) {
    res.resize(As.size());
    for (size_t i=0; i<As.size(); ++i) {
        res[i]=As[i].transpose();
    }
}

template<int Nrows, int Ncols, int Ndiag>
void compute_SVD(std::vector<Eigen::Matrix<float, Ndiag, 1> >& values,
                 std::vector<Eigen::Matrix<float, Nrows, Nrows> >& leftvs,
                 std::vector<Eigen::Matrix<float, Ncols, Ncols> >& rightvs,
                 const std::vector<Eigen::Matrix<float, Nrows, Ncols> >& As) {
     values.resize(As.size());
     leftvs.resize(As.size());
     rightvs.resize(As.size());
     typedef Eigen::Matrix<float, Nrows, Ncols> mat_t;
     typedef Eigen::JacobiSVD<mat_t> svd_t;
     for (size_t i=0; i<As.size(); ++i) {
         svd_t svd(As[i], Eigen::ComputeFullU | Eigen::ComputeFullV);
         values[i]=svd.singularValues();
         leftvs[i]=svd.matrixU();
         rightvs[i]=svd.matrixV();
     }
}

template<int N>
void compute_inverse(std::vector<Eigen::Matrix<float, N, N> >& res,
                     const std::vector<Eigen::Matrix<float, N, N> >& As) {
    res.resize(As.size());
    for (size_t i=0; i<As.size(); ++i) {
        res[i]=As[i].inverse();
    }
}

template<int Nrows, int Ncols>
void compute_sum(std::vector<Eigen::Matrix<float, Ncols, Nrows> >& res,
                 const std::vector<Eigen::Matrix<float, Nrows, Ncols> >& As,
                 const std::vector<Eigen::Matrix<float, Nrows, Ncols> >& Bs) {
    assert(As.size()==Bs.size());
    res.resize(As.size());
    for (size_t i=0; i<As.size(); ++i) {
        res[i]=As[i]+Bs[i];
    }
}

template<int M, int N, int P>
void compute_prod(std::vector<Eigen::Matrix<float, M, P> >& res,
                 const std::vector<Eigen::Matrix<float, M, N> >& As,
                 const std::vector<Eigen::Matrix<float, N, P> >& Bs) {
    assert(As.size()==Bs.size());
    res.resize(As.size());
    for (size_t i=0; i<As.size(); ++i) {
        res[i]=As[i]*Bs[i];
    }
}

template<int M, int N, int P>
void compute_tprod(std::vector<Eigen::Matrix<float, M, P> >& res,
                   const std::vector<Eigen::Matrix<float, N, M> >& As,
                   const std::vector<Eigen::Matrix<float, N, P> >& Bs) {
    assert(As.size()==Bs.size());
    res.resize(As.size());
    for (size_t i=0; i<As.size(); ++i) {
        res[i]=As[i].transpose()*Bs[i];
    }
}

template<int Nrows, int Ncols>
void test_transpose(size_t n=50) {
    typedef Eigen::Matrix<float, Nrows, Ncols> mat_t;
    typedef Eigen::Matrix<float, Ncols, Nrows> trans_t;
    size_t N=n*n*n;

    // generate random matrices of prescribed size
    std::vector<mat_t> matrices(N);
    // std::cout << "test matrices:\n";
    for (size_t i=0; i<N; ++i) {
        matrices[i]=mat_t::Random();
    }

    // convert to Nrrd
    std::array<size_t, 3> dims = {n, n, n};
    Nrrd* nrrd=mat2nrrd(matrices, dims);

    // process each matrix individually
    std::vector<trans_t> transposes;

    nvis::timer timer;
    compute_transpose<Nrows, Ncols>(transposes, matrices);
    double t1=timer.elapsed();

    xavier::nrrd_utils::nrrd_matrix_transpose<float, Nrows, Ncols> trans_op;
    toverhead=0;
    timer.restart();
    Nrrd* res=trans_op(nrrd);
    double t2=timer.elapsed()-toverhead;

    std::cout << "\n\ntest_transpose<" << Nrows << ", " << Ncols << ">:\n";
    double err=check_same<Ncols, Nrows>(res, transposes, "transpose");
    std::cout << "vec solution took " << t1 << " s. | nrrd solution took " << t2 << " s.\n";
}

template<int Nrows, int Ncols1, int Ncols2>
void test_product(size_t n=50) {
    typedef Eigen::Matrix<float, Nrows, Ncols1> left_mat_t;
    typedef Eigen::Matrix<float, Ncols1, Ncols2> right_mat_t;
    typedef Eigen::Matrix<float, Nrows, Ncols2> res_mat_t;
    size_t N=n*n*n;

    // generate random matrices of prescribed size
    std::vector<left_mat_t> left(N);
    std::vector<right_mat_t> right(N);
    for (size_t i=0; i<N; ++i) {
        left[i]=left_mat_t::Random();
        right[i]=right_mat_t::Random();
    }

    // convert to Nrrd
    std::array<size_t, 3> dims={n, n, n};
    Nrrd* nrrd1=mat2nrrd(left, dims);
    Nrrd* nrrd2=mat2nrrd(right, dims);
    // std::cout << "Checking correctness of conversion from vector to nrrd\n";
    // check_same<Nrows, Ncols1>(nrrd1, left);
    // check_same<Ncols1, Ncols2>(nrrd2, right);

    // process each matrix individually
    std::vector<res_mat_t> product;

    nvis::timer timer;
    compute_prod<Nrows, Ncols1, Ncols2>(product, left, right);
    double t1=timer.elapsed();

    xavier::nrrd_utils::nrrd_matrix_product<float, Nrows, Ncols1, Ncols2> prod_op;
    toverhead=0;
    timer.restart();
    Nrrd* res=prod_op(nrrd1, nrrd2);
    double t2=timer.elapsed()-toverhead;
    std::cout << "\n\ntest_product<" << Nrows << ", " << Ncols1 << ", " << Ncols2 << ">:\n";
    double err=check_same<Nrows, Ncols2>(res, product, "product");
    std::cout << "vec solution took " << t1 << " s. | nrrd solution took " << t2 << " s.\n";
}


template<int Nrows, int Ncols>
void test_SVD(size_t n=50) {
    using boost::static_signed_min;
    constexpr int P=static_signed_min<Nrows, Ncols>::value; // nb. sing values

    typedef Eigen::Matrix<float, Nrows, Ncols> mat_t;
    typedef Eigen::Matrix<float, Nrows, Nrows> Lmat_t;
    typedef Eigen::Matrix<float, Ncols, Ncols> Rmat_t;
    typedef Eigen::Matrix<float, P, 1> vec_t;
    size_t N=n*n*n;

    // generate random matrices of prescribed size
    std::vector<mat_t> matrices(N);
    for (size_t i=0; i<N; ++i) {
        matrices[i]=mat_t::Random();
    }

    // convert to Nrrd
    std::array<size_t, 3> dims={n, n, n};
    Nrrd* nrrd=mat2nrrd(matrices, dims);

    // process each matrix individually
    std::vector<vec_t> singvals;
    std::vector<Lmat_t> leftvecs;
    std::vector<Rmat_t> rightvecs;

    nvis::timer timer;
    compute_SVD<Nrows, Ncols, P>(singvals, leftvecs, rightvecs, matrices);
    double t1=timer.elapsed();

    xavier::nrrd_utils::nrrd_matrix_svd<float, Nrows, Ncols> svd_op;
    toverhead=0;
    timer.restart();
    Nrrd *snrrd, *lnrrd, *rnrrd;
    svd_op(nrrd, snrrd, lnrrd, rnrrd);
    double t2=timer.elapsed()-toverhead;

    std::cout << "\n\ntest_SVD<" << Nrows << ", " << Ncols << ">:\n";
    double err=check_same<P, 1>(snrrd, singvals, "singular values");
    err=check_same<Nrows, Nrows>(lnrrd, leftvecs, "left singular vectors");
    err=check_same<Ncols, Ncols>(rnrrd, rightvecs, "right singular vectors");
    std::cout << "vec solution took " << t1 << " s. | nrrd solution took " << t2 << " s.\n";
}

int main(int argc, char* argv[]) {
    int n=20;
    if (argc>1) n=atoi(argv[1]);

    std::cout << "Testing transpose...\n";
    test_transpose<2,2>(n);
    test_transpose<3,2>(n);
    test_transpose<2,3>(n);
    test_transpose<3,3>(n);
    test_transpose<3,4>(n);
    test_transpose<4,4>(n);
    test_transpose<5,4>(n);
    test_transpose<4,7>(n);

    std::cout << "Testing product...\n";
    test_product<2,2,2>(n);
    test_product<3,2,2>(n);
    test_product<3,3,3>(n);
    test_product<2,3,4>(n);
    test_product<4,4,4>(n);
    test_product<5,4,7>(n);
    test_product<2,7,5>(n);

    std::cout << "Testing SVD...\n";
    test_SVD<2,2>(n);
    test_SVD<3,2>(n);
    test_SVD<2,3>(n);
    test_SVD<3,3>(n);
    test_SVD<4,4>(n);
    test_SVD<5,4>(n);
    test_SVD<4,7>(n);

    return 0;
}
