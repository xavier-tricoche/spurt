#ifndef __XAVIER_READ_SPARSE__
#define __XAVIER_READ_SPARSE__

#include <image/nrrd_wrapper.hpp>
#include <Eigen/Sparse>
#include <map>
#include <misc/progress.hpp>
#include <misc/option_parse.hpp>
#include "teem/nrrd.h"

size_t n_iter=0;
double eps = 1.0e-6;
double tol = 1;
std::string name_out;

template<typename Size_, typename Scalar_>
void import_sparse_system(Eigen::SparseMatrix<Scalar_>& matrix,
                          Eigen::Matrix<Scalar_, Eigen::Dynamic, 1>& v,
                          size_t n_rows, Size_* indices, Scalar_* values, Scalar_* rhs) {
    typedef Size_                         size_type;
    typedef Scalar_                       scalar_type;
    typedef Eigen::Triplet<Scalar_>       triplet_type;

    std::vector<triplet_type> triplets;
    size_t val_p=0;
    size_t index_p=0;
    size_t n_columns=0;
    for (size_t row=0 ; row<n_rows ; ++row) {
        v(row) = rhs[row];
        unsigned int n_cols = indices[index_p++];
#if 0
        if (n_cols > 9) std::cout << n_cols << " non-zero terms on row " << row << '\n';
        auto it = number_of_rows_per_number_of_non_zeros.find(n_cols);
        if (it == number_of_rows_per_number_of_non_zeros.end()) {
            number_of_rows_per_number_of_non_zeros[n_cols] = 1;
        }
        else
            it->second = it->second + 1;
#endif
        for (unsigned int j=0; j<n_cols; ++j) {
            size_t col=indices[index_p++];
            triplets.push_back(triplet_type(row, col-1, values[val_p++]));
#if 0
            col_ids.push_back(col);
#endif
        }
        n_columns = std::max<size_t>(n_columns, triplets.back().col()+1);
    }

    matrix.resize(n_rows, n_columns);
    std::cout << "there are " << n_rows << " rows and " << n_columns << " columns\n";

    matrix.setFromTriplets(triplets.begin(), triplets.end());
}

template<typename Scalar_>
std::array<Scalar_ , 9> errors(
    const Eigen::SparseMatrix<Scalar_>& A,
    const Eigen::Matrix<Scalar_, Eigen::Dynamic, 1>& x,
    const Eigen::Matrix<Scalar_, Eigen::Dynamic, 1>& b) {
    typedef Eigen::SparseMatrix<Scalar_>               sparse_matrix_type;
    typedef Eigen::Matrix<Scalar_, Eigen::Dynamic, 1>  vector_type;
    typedef Scalar_                                    scalar_type;

    const size_t nb_output = A.cols();
    const size_t nb_input = A.rows()-A.cols();

    // fit quality
    vector_type fit_res = A.topRows(nb_input)*x-b.topRows(nb_input);
    vector_type laplace_res = A.bottomRows(nb_output)*x-b.bottomRows(nb_output);

    scalar_type fit_err2 = fit_res.norm();
    scalar_type fit_relerr2 = fit_err2 / b.topRows(nb_input).norm();
    scalar_type fit_err1 = fit_res.cwiseAbs().sum();
    scalar_type fit_relerr1 = fit_err1 / b.topRows(nb_input).cwiseAbs().sum();
    scalar_type fit_errinf = fit_res.cwiseAbs().maxCoeff();
    scalar_type fit_relerrinf = fit_errinf / b.topRows(nb_input).cwiseAbs().maxCoeff();

    scalar_type laplace_err2 = laplace_res.norm();
    scalar_type laplace_err1 = laplace_res.cwiseAbs().sum();
    scalar_type laplace_errinf = laplace_res.cwiseAbs().maxCoeff();

    return std::array<scalar_type, 9>({{fit_err2, fit_relerr2, fit_err1, fit_relerr1, fit_errinf, fit_relerrinf, laplace_err2, laplace_err1, laplace_errinf}});
}

template<typename Size_, typename Scalar_>
void run(Nrrd* nrrd_indices, Nrrd* nrrd_values, Nrrd* nrrd_rhs) {
    typedef Eigen::SparseMatrix<Scalar_>               sparse_matrix_type;
    typedef Eigen::Matrix<Scalar_, Eigen::Dynamic, 1>  vector_type;
    typedef Size_                                      size_type;
    typedef Scalar_                                    scalar_type;

    sparse_matrix_type A;
    vector_type b;
    size_type* indices = static_cast<size_type *>(nrrd_indices->data);
    scalar_type* values = static_cast<scalar_type *>(nrrd_values->data);
    scalar_type* rhs = static_cast<scalar_type *>(nrrd_rhs->data);
    size_t nrows = indices[0];
    b.resize(nrows);
    import_sparse_system<size_type, scalar_type>(A, b, nrows, indices+1, values, rhs);

    Eigen::LeastSquaresConjugateGradient<sparse_matrix_type> lscg;

    // Eigen::setNbThreads(6);
    std::cout << "Eigen is using " << Eigen::nbThreads() << " threads\n";

    spurt::ProgressDisplay timer(false), total_timer(false);
    total_timer.start();
    timer.start();
    lscg.compute(A);
    timer.end();
    std::cout << "LSCG initialization took "<< timer.cpu_time() << " ms. (cpu) / " << timer.wall_time() << " ms. (wall)\n";
    vector_type x(A.cols());
    if (n_iter > 0) lscg.setMaxIterations(n_iter);
    if (eps > 0) lscg.setTolerance(eps);

    scalar_type last_error = std::numeric_limits<scalar_type>::max();
    scalar_type error = 0.99*last_error;
    int counter = 0;
    while (error > eps && error < tol*last_error) {
        timer.start();
        if (!counter)
            x = lscg.solve(b);
        else
            x = lscg.solveWithGuess(b, x);
        ++counter;
        timer.end();
        std::array<double, 2> partial_times({{timer.cpu_time(), timer.wall_time()}});
        std::array<double, 2> total_times({{total_timer.instant_cpu_time(), total_timer.instant_wall_time()}});
        size_t iterations = (n_iter == 0 ? lscg.iterations() : n_iter*counter);
        last_error = error;

        scalar_type myerr = (A.transpose()*(A*x-b)).norm() / b.norm();
        std::array<scalar_type, 9> errs = errors(A, x, b);

        error = lscg.error();
        std::cout << iterations << " iteration, est. error: " << error << ", actual error: " << myerr << ", delta: " << -error+last_error << "%, "
        << "fit err: l2=" << errs[0] << " (" << errs[1] << "), l1="
        << errs[2] << " (" << errs[3] << "), linf=" << errs[4] << " (" << errs[5] << "), laplace l2=" << errs[6] << ", l1=" << errs[7] << ", linf=" << errs[8]
        << ", total: cpu: " << total_times[0]/1000. << " s., wall: " << total_times[1]/1000 << " s.\n";
    }


    // vector_type residual = (b-A*x).topRows(A.cols());
    // vector_type laplace = (b-A*x).bottomRows(A.rows()-A.cols());
    // scalar_type l1, l2, linf, zerol1, zerol2, zerolinf;
    // l1 = l2 = linf = zerol1 = zerol2 = zerolinf = 0;
    // for (size_t i=0; i<err.rows(); ++i) {
    //     const scalar_type& bi = b(i);
    //     const scalar_type& erri = err(i);
    //     if (bi != 0) {
    //         scalar_type rel = std::abs(erri/bi);
    //         l1 += rel;
    //         l2 += rel*rel;
    //         linf = std::max(linf, rel);
    //     }
    //     else {
    //         zerol1 += std::abs(erri);
    //         zerol2 += erri*erri;
    //         zerolinf = std::max(zerolinf, std::abs(erri));
    //     }
    // }
    //
    // std::cout << "L2 norm      : " << sqrt(l2) << '\n';
    // std::cout << "L1 norm      : " << l1 << '\n';
    // std::cout << "Linf norm    : " << linf << '\n';
    // std::cout << "zerol2 norm  : " << sqrt(zerol2) << '\n';
    // std::cout << "zerol1 norm  : " << sqrt(zerol1) << '\n';
    // std::cout << "zerolinf norm: " << zerolinf << '\n';

    size_t size = x.rows();
    spurt::nrrd_utils::writeNrrd(&x(0), name_out, spurt::nrrd_utils::nrrd_value_traits_from_type<Scalar_>::index, 1, &size);
}

int main(int argc, const char* argv[]) {

    std::string name_id, name_coef, name_rhs;
    namespace xcl = spurt::command_line;
    int verbose=0;

    xcl::option_traits
        required_group(true, false, "Required parameters"),
        positional_group(true, true, "Positional parameters"),
        optional_group(false, false, "Optional parameters");

    xcl::option_parser parser(argv[0],
        "Solve large sparse least squares problem using Eigen's\npreconditioned LSCG implementation");

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("indices", name_id, "Indices filename (.nrrd)",
                         positional_group);
        parser.add_value("values", name_coef, "Coefficients filename (.nrrd)",
                         positional_group);
        parser.add_value("rhs", name_rhs, "Right hand side filename (.nrrd)",
                         positional_group);
        parser.add_value("output", name_out, "Output filename (.nrrd)",
                         positional_group);
        parser.add_value("niter", n_iter, "Number of iterations per round");
        parser.add_value("epsilon", eps, "Error threshold");
        parser.add_value("tolerance", tol, "Error improvement tolerance");
        parser.add_value("verbose", verbose, verbose,
                         "Verbose level", optional_group);
        parser.parse(argc, argv);
    }
    catch (std::exception& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }


    // name_id = argv[1];
    // name_coef = argv[2];
    // name_rhs = argv[3];

    Nrrd* nin_id = spurt::nrrd_utils::readNrrd(name_id);
    Nrrd* nin_coef = spurt::nrrd_utils::readNrrd(name_coef);
    Nrrd* nin_rhs = spurt::nrrd_utils::readNrrd(name_rhs);

    if (nin_coef->type != nrrdTypeFloat && nin_coef->type != nrrdTypeDouble) {
        std::cerr << "value type " << nin_coef->type << " not supported\n";
        exit(1);
    }

    switch (nin_id->type) {
        case nrrdTypeInt:
            if (nin_coef->type == nrrdTypeFloat) {
                run<int, float>(nin_id, nin_coef, nin_rhs);
            }
            else {
                run<int, double>(nin_id, nin_coef, nin_rhs);
            }
            break;
        case nrrdTypeUInt:
            if (nin_coef->type == nrrdTypeFloat) {
                run<unsigned int, float>(nin_id, nin_coef, nin_rhs);
            }
            else {
                run<unsigned int, double>(nin_id, nin_coef, nin_rhs);
            }
            break;
        case nrrdTypeLLong:
            if (nin_coef->type == nrrdTypeFloat) {
                run<long int, float>(nin_id, nin_coef, nin_rhs);
            }
            else {
                run<long int, double>(nin_id, nin_coef, nin_rhs);
            }
            break;
        case nrrdTypeULLong:
            if (nin_coef->type == nrrdTypeFloat) {
                run<unsigned long int, float>(nin_id, nin_coef, nin_rhs);
            }
            else {
                run<unsigned long int, double>(nin_id, nin_coef, nin_rhs);
            }
            break;
        default:
            std::cerr << "unsupported index type: " << nin_id->type << '\n';
            exit(1);
    };

    nrrdNuke(nin_coef);
    nrrdNuke(nin_id);
    nrrdNuke(nin_rhs);

    return 0;
}


#endif // __XAVIER_READ_SPARSE__
