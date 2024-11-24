#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <boost/lexical_cast.hpp>
#include <math/polynomial.hpp>
#include <misc/progress.hpp>


template<typename T>
struct type_traits {};

template<>
struct type_traits<float> {
    static std::string name() { return "float"; }
};

template<>
struct type_traits<double> {
    static std::string name() { return "double"; }
};

template<typename T, size_t Dim, size_t Order>
void compare(size_t n) {
    typedef T                                                       val_t;
    typedef spurt::small_vector<T, Dim>                             vec_t;
    typedef spurt::polynomial::polynomial_basis<T, Dim, vec_t>      poly_t;
    typedef typename poly_t::monomial_type                          mono_t;
    typedef typename poly_t::derivative_type                        deriv_t;
    typedef spurt::polynomial::alt_polynomial_basis<T, Dim, vec_t, Order> alt_poly_t;
    
    const val_t max_value = 1000.;
    
    spurt::timer _timer;
    double t=0, t_alt=0;
    double err_val=0, err_deriv=0;
    
    std::vector<mono_t> basis_monos;
    std::vector<deriv_t> basis_monos_deriv;
    poly_t::compute_basis(basis_monos, Order);
    poly_t::compute_basis_derivative(basis_monos_deriv, basis_monos);
    
    std::cout << "Polynomial basis of order " << Order 
        << " in " << Dim << " dimensions: [ ";
    for (size_t i=0 ; i<basis_monos.size() ; ++i) {
        std::cout << basis_monos[i] << " ";
    }
    std::cout << "]\n";
    
    std::vector<val_t> values, values_alt;
    std::vector<vec_t> derivatives, derivatives_alt;
    
    // Repeat n times
    for (size_t i=0 ; i<n ; ++i) {
        // Random Dim-D vector
        vec_t x;
        for (size_t j=0 ; j<Dim ; ++j) {
            x[j] = (-0.5 + drand48())*max_value;
        }
        values.clear();
        values_alt.clear();
        derivatives.clear();
        derivatives_alt.clear();
        
        _timer.start();
        values.resize(basis_monos.size());
        derivatives.resize(basis_monos.size());
        for (size_t i=0 ; i<basis_monos.size() ; ++i) {
            values[i] = basis_monos[i](x);
            for (size_t j=0 ; j<Dim ; ++j) {
                derivatives[i][j] = basis_monos_deriv[i][j](x);
            }
        }
        t += _timer.elapsed();
        
        _timer.start();
        alt_poly_t::eval(values_alt, x);
        alt_poly_t::derivative(derivatives_alt, x);
        t_alt += _timer.elapsed();
        
        // basic sanity check, first
        if (values.size() != values_alt.size()) {
            throw std::runtime_error("Wrong number of entries in values: " +
                boost::lexical_cast<std::string>(values.size()) + " != " +
                boost::lexical_cast<std::string>(values_alt.size()));
        }
        
        if (derivatives.size() != derivatives_alt.size()) {
            throw std::runtime_error("Wrong number of entries in derivatives: " +
                boost::lexical_cast<std::string>(derivatives.size()) + " != " +
                boost::lexical_cast<std::string>(derivatives_alt.size()));
        }
        
        // compute delta in value
        double err = 0;
        for (size_t i=0 ; i<values.size() ; ++i) {
            double delta = fabs(values[i] - values_alt[i])/fabs(values_alt[i]);
            // std::cerr << "\n" << delta;
            err += delta;
        }
        if (std::isnan(err)) {
            std::cerr << "\n\nNaN detected in values:\n"
                      << "x=" << x << "\nval(x) = [";
            std::copy(values.begin(), values.end(), 
                      std::ostream_iterator<val_t>(std::cerr, " "));
            std::cerr << "]\n" << "alt_val(x) = [";
            std::copy(values_alt.begin(), values_alt.end(), 
                      std::ostream_iterator<val_t>(std::cerr, " "));
            std::cerr << "]\n";
            exit(-1);
        }
        if (err > 1.0e-8) {
            std::cerr << "\n\nLarge value error detected:\n"
                      << "value 1 = [";
            std::copy(values.begin(), values.end(),
                std::ostream_iterator<val_t>(std::cerr, " "));
            std::cerr << "]\n" << "value 2 = [";
            std::copy(values_alt.begin(), values_alt.end(),
                std::ostream_iterator<val_t>(std::cerr, " "));
            std::cerr << "]\n";
            exit(-1);
        }
        
        err_val += err*err;
        
        // compute delta in derivatives
        err = 0;
        for (size_t i=0 ; i<derivatives.size() ; ++i) {
            if (spurt::norm(derivatives_alt[i])) {
                vec_t delta = derivatives[i] - derivatives_alt[i];
                err += spurt::norm(delta) / spurt::norm(derivatives_alt[i]);
            }
        }
        if (err > 1.0e-8) {
            std::cerr << "\n\nLarge derivative error detected:\n"
                << "deriv = [";
            std::copy(derivatives.begin(), derivatives.end(),
                std::ostream_iterator<vec_t>(std::cerr, " "));
            std::cerr << "]\n" << "alt_deriv = [";
            std::copy(derivatives_alt.begin(), derivatives_alt.end(),
                std::ostream_iterator<vec_t>(std::cerr, " "));
            std::cerr << "]\n";
            exit(-1);
        }
        err_deriv += err*err;
    }
        
    std::cout << "type " << type_traits<T>::name() 
              << " \tdim " << Dim
              << " \torder " << Order
              << " \tProcedural:  " << t << " s. (" 
              << (double)n/t << " Hz)"
              << " \tHardcoded: " << t_alt << " s. (" 
              << (double)n/t_alt << " Hz)"
              << " \tSpeedup: " << t/t_alt
              << " \terror: " << 1/(double)n*sqrt(err_val) << " / "
              << 1/(double)n*sqrt(err_deriv) << '\n';
}

std::string me;
void usage() {
    std::cout 
        << me << ": Compare two methods of evaluation for polynomial bases.\n\n"
        << "USAGE: " << me << " [options]\n\n"
        << "Options:\n"
        << " -h         Print this message\n"
        << " -n <int>   Number of iterations for benchmarking (default: 1e+6)\n"
        << "\n";
    exit(0);
}

bool has_enough(int at, int size, int n=1) {
    return at + n < size;
}

int main(int argc, char* argv[]) {
    srand48(time(0));
    
    me = argv[0];
    size_t n = 1000000;
    
    for (int i=1; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") usage();
        else if (arg == "-n") {
            if (!has_enough(i, argc)) usage();
            n = atoi(argv[++i]);
        }
        else {
            std::cerr << "Unrecognized argument: " << arg << '\n';
            usage();
        }
    }
    
    std::cout << "\n\nTesting float precision...\n";
    std::cout << "\n\t ... in 2 dimensions...\n";
    compare<float, 2, 1>(n);
    compare<float, 2, 2>(n);
    compare<float, 2, 3>(n);
    compare<float, 2, 4>(n);
    compare<float, 2, 5>(n);
    compare<float, 2, 6>(n);
    compare<float, 2, 7>(n);
    compare<float, 2, 8>(n);
    compare<float, 2, 9>(n);
    compare<float, 2, 10>(n);
    compare<float, 2, 11>(n);
    compare<float, 2, 12>(n);
    
    std::cout << "\n\t ... in 3 dimensions...\n";
    compare<float, 3, 1>(n);
    compare<float, 3, 2>(n);
    compare<float, 3, 3>(n);
    compare<float, 3, 4>(n);
    compare<float, 3, 5>(n);
    compare<float, 3, 6>(n);
    compare<float, 3, 7>(n);
    compare<float, 3, 8>(n);
    compare<float, 3, 9>(n);
    compare<float, 3, 10>(n);
    compare<float, 3, 11>(n);
    compare<float, 3, 12>(n);
    
    std::cout << "\n\t ... in 4 dimensions...\n";
    compare<float, 4, 1>(n);
    compare<float, 4, 2>(n);
    compare<float, 4, 3>(n);
    compare<float, 4, 4>(n);
    compare<float, 4, 5>(n);
    compare<float, 4, 6>(n);
    compare<float, 4, 7>(n);
    compare<float, 4, 8>(n);
    compare<float, 4, 9>(n);
    compare<float, 4, 10>(n);
    compare<float, 4, 11>(n);
    compare<float, 4, 12>(n);
    
    std::cout << "\n\nTesting double precision...\n";
    std::cout << "\n\t ... in 2 dimensions...\n";
    compare<double, 2, 1>(n);
    compare<double, 2, 2>(n);
    compare<double, 2, 3>(n);
    compare<double, 2, 4>(n);
    compare<double, 2, 5>(n);
    compare<double, 2, 6>(n);
    compare<double, 2, 7>(n);
    compare<double, 2, 8>(n);
    compare<double, 2, 9>(n);
    compare<double, 2, 10>(n);
    compare<double, 2, 11>(n);
    compare<double, 2, 12>(n);
    
    std::cout << "\n\t ... in 3 dimensions...\n";
    compare<double, 3, 1>(n);
    compare<double, 3, 2>(n);
    compare<double, 3, 3>(n);
    compare<double, 3, 4>(n);
    compare<double, 3, 5>(n);
    compare<double, 3, 6>(n);
    compare<double, 3, 7>(n);
    compare<double, 3, 8>(n);
    compare<double, 3, 9>(n);
    compare<double, 3, 10>(n);
    compare<double, 3, 11>(n);
    compare<double, 3, 12>(n);
    
    std::cout << "\n\t ... in 4 dimensions...\n";
    compare<double, 4, 1>(n);
    compare<double, 4, 2>(n);
    compare<double, 4, 3>(n);
    compare<double, 4, 4>(n);
    compare<double, 4, 5>(n);
    compare<double, 4, 6>(n);
    compare<double, 4, 7>(n);
    compare<double, 4, 8>(n);
    compare<double, 4, 9>(n);
    compare<double, 4, 10>(n);
    compare<double, 4, 11>(n);
    compare<double, 4, 12>(n);
    
    return 0;
}
