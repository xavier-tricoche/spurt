#include <Eigen/Core>
#include <cstdlib>
#include "math_array.hpp"
#include <algorithm>
#include <iostream>

template<typename T>
struct randval {
    T operator()() {
        return static_cast<T>(std::rand())/static_cast<T>(RAND_MAX);
    }
};

int main(int argc, char* argv[]) {

    using namespace spurt;

    vec3 d[4];
    fvec3 f[4];
    ivec3 j[4];

    lexicographical_order lorder;
    eps_lexicographical_order elorder(1.0e-6);

    randval<double> dd;
    randval<float> ff;
    randval<int> ii;

    std::srand(time(nullptr));
    for (int i=0; i<4; ++i) {
        d[i] = vec3(dd(), dd(), dd());
        f[i] = fvec3(ff(), ff(), ff());
        j[i] = ivec3(ff()*100, ff()*100, ff()*100);
    }
    vec2 a(dd(), dd());
    vec4 b(dd(), dd(), dd(), dd());

    vec3 l = d[0] + d[1];
    std::cout << "l=" << l << '\n';
    if (any((d[2]+d[3]) > 0.5)) std::cout << "any TRUE\n";
    else std::cout << "any is FALSE\n";
    fvec3 g = f[0]*f[1];
    fvec3 h = f[2]/f[3];
    std::cout << "the norm of " << f[0] << " * " << f[1] <<  " (=" << g << ") is " << norm(f[0]*f[1]) << '\n';
    std::cout << "the norm of " << f[2] << " / " << f[3] << " (=" << h << ") is " << norm(f[2]/f[3]) << '\n';
    std::cout << "square of " << d[0] << " is " << square(d[0]) << '\n';
    std::cout << "min coefficient in " << f[0] << " is " << min(f[0]) << '\n';
    std::cout << "partial initialization of a 4D vector with 2.3 and 4.5 is " << vec4(2.3, 4.5) << '\n';

    vec5 a5 = 203.45;
    std::cout << "initialization of a 5D vector from constant 203.45 yields " <<  a5 << '\n';


    vec3 k(f[0]);
    std::cout << "Initializing a double vector with a float vector produced " << k << '\n';

    std::cout << "The sum of double vector " << d[0] << " and a float vector " << f[0] << " is " << d[0]+f[0] << '\n';
    std::cout << "The sum of a float vector " << f[1] << " and an int vector " << j[0] << " is " << f[1]+j[0] << '\n';
    std::cout << "The product of a double vector " << d[1] << " and an int vector " << j[1] << " is " << d[1] * j[1] << '\n';

    std::sort(d, d+4, lorder);
    std::cout << "after sorting, the vectors are:\n";
    for (int i=0; i<4; ++i) std::cout << d[i] << '\n';

    d[2][1] *= -1;
    std::cout << "the absolute value of " << d[2] << " is " << abs(d[2]) << '\n';

    std::cout << "the square of " << d[3] << " is " << square(d[3]) << '\n';

    small_square_matrix<double, 2> A(small_vector<double, 4>(dd(), dd(), dd(), dd())), B(small_vector<double, 4>(dd(), dd(), dd(), dd()));

    small_vector<double, 16> av({dd(), dd(), dd(), dd(), dd(), dd(), dd(), dd(), dd(), dd(), dd(), dd(), dd(), dd(), dd(), dd()});

    typedef small_square_matrix<double, 4> sq4;
    static_assert(sq4::base_type::size == 16, "Size of vector is wrong");

    sq4 C(small_vector<double, 16>({dd(), dd(), dd(), dd(), 
                                    dd(), dd(), dd(), dd(), 
                                    dd(), dd(), dd(), dd(), 
                                    dd(), dd(), dd(), dd()}));

    sq4 D(small_vector<float, 16>({ff(), ff(), ff(), ff(),
                                   ff(), ff(), ff(), ff(), 
                                   ff(), ff(), ff(), ff(), 
                                   ff(), ff(), ff(), ff()}));

    std::cout << "A=" << A << ", B=" << B << ", C=" << C << ", D=" << D << '\n';
    auto E = A*B; 
    auto F = C*inverse(D);
    std::cout << "AxB=" << E << '\n';
    std::cout << "C*D^{-1}=" << F << "\n";
    std::cout << "C*D^{-1}*D=" << C*inverse(D)*D << '\n';
    std::cout << "C*D^{-1}*D*C^{-1}=" << C*inverse(D)*D*inverse(C) << '\n';
    
    //", and their sum is " << E+F << '\n'; 
    std::cout << "Computing eigenvalues of C*D^{1}\n";
    small_vector<std::complex<double>, 4> evals; 
    small_matrix<std::complex<float>, 4, 4> evecs;
    eigensystem(evals, evecs, F);
    std::cout << "sorted eigenvalues are:\n" << evals << '\n';
    std::cout << "matching eigenvectors are:" << evecs << '\n';
    
    auto rit = F.begin<sq4::rowwise_iterator>();
    std::cout << "iterating over the matrix F rowwise:\n";
    for (; rit!=F.end<sq4::rowwise_iterator>(); ++rit) {
        std::cout << " " << *rit;
    }
    std::cout << '\n';
    auto cit = F.begin<sq4::columnwise_iterator>();
    std::cout << "iterating over the matrix F columnwise:\n";
    for (; cit!=F.end<sq4::columnwise_iterator>(); ++cit) {
        std::cout << " " << *cit;
    }
    std::cout << '\n';
    
    std::cout << "F is \n" << F << '\n';
    std::cout << "iterating over row 2 of matrix F\n";
    for (auto it=F.row(2).begin(); it!=F.row(2).end(); ++it) std::cout << *it << ", ";
    std::cout << '\n';
    std::cout << "printing row 2 directly:\n" << F.row(2) << '\n';
    std::cout << "iterating over column 1 of matrix F\n";
    for (auto it=F.column(1).begin(); it!=F.column(1).end(); ++it) std::cout << *it << ", ";
    std::cout << '\n';
    std::cout << "printing column 1 directly:\n" << F.column(1) << '\n';
    
    std::cout << "replacing col 3 of F with col 1 of F\n";
    std::copy(F.column(3).begin(), F.column(3).end(), F.column(1).begin());
    std::cout << "F is now:\n" << F << '\n';
    
    std::cout << "Using operator= to replace column 2 of F by column 0 of F\n";
    std::cout << "performing " << F.column(2) << " = " << F.column(0) << '\n';
    
    std::cout << "now\n";
    F.column(2) = F.column(0);
    std::cout << "F is now:\n" << F << '\n';
    
    std::cout << "same thing for row 0 and row 3\n";
    F.row(3) = F.row(0);
    std::cout << "F is now\n" << F << '\n';

    return 0;
}