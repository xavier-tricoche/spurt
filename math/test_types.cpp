#include <Eigen/Core>
#include <cstdlib>
#include "small_vector.hpp"
#include "small_matrix.hpp"
#include <algorithm>
#include <iostream>
// #include <misc/meta_utils.hpp>

#include "small_matrix.hpp"

template<typename T>
std::string whattypeami(T v) { return ""; }

template<>
std::string whattypeami(float f) {
    return "float";
}

template<>
std::string whattypeami(double d) {
    return "double";
}

template<>
std::string whattypeami(int i) {
    return "int";
}

template<>
std::string whattypeami(long i) {
    return "long int";
}

template<>
std::string whattypeami(unsigned long i) {
    return "unsigned long int";
}


template<typename Matrix, typename Iterator>
bool check_rowwise_iterator(const Matrix& m, Iterator begin, Iterator end) {
    typename Matrix::value_type value_t;
    const size_t nrows = m.nrows;
    const size_t ncols = m.ncols;
    Iterator it = begin;
    for (size_t r=0; r<nrows; ++r) {
        for (size_t c=0; c<ncols; ++c, ++it) {
            if (m(r,c) != *it) return false;
        }
    } 
    return true;
}

template<typename Matrix, typename Iterator>
bool check_colwise_iterator(const Matrix& m, Iterator begin, Iterator end) {
    typename Matrix::value_type value_t;
    const size_t nrows = m.nrows;
    const size_t ncols = m.ncols;
    Iterator it = begin;
    for (size_t c=0; c<ncols; ++c) {
        for (size_t r=0; r<nrows; ++r, ++it) {
            if (m(r,c) != *it) return false;
        }
    } 
    return true;
}


int main(int argc, char* argv[]) {

    using namespace spurt;

    vec3 d[4];
    fvec3 f[4];
    ivec3 j[4];

    lexicographical_order lorder;
    eps_lexicographical_order elorder(1.0e-6);

    for (int i=0; i<4; ++i) {
        d[i] = vec3::random(0, 1);
        f[i] = fvec3::random(0, 1);
        j[i] = ivec3::random(0, 100);
    }
    vec2 a = vec2::random(0, 1);
    vec4 b = vec4::random(0, 1);

    vec3 l = d[0] + d[1];
    std::cout << "sum of " << d[0] << " and " << d[1] << " is " << l << '\n';
    if (any((d[2]+d[3]) > 0.5)) std::cout << "any larger than 0.5 is TRUE in " << d[2]+d[3] << "\n";
    else std::cout << "None is larger than 0.5 in " << d[2]+d[3] << "\n";
    fvec3 g = f[0]*f[1];
    fvec3 h = f[2]/f[3];
    std::cout << "the norm of " << f[0] << " * " << f[1] <<  " (=" << g << ") is " << norm(f[0]*f[1]) << '\n';
    std::cout << "norm computed differently is " << norm(g) << '\n';
    std::cout << "the norm of " << f[2] << " / " << f[3] << " (=" << h << ") is " << norm(f[2]/f[3]) << '\n';
    std::cout << "norm computed differently is " << norm(h) << '\n';
    std::cout << "square of " << d[0] << " is " << square(d[0]) << '\n';
    std::cout << "min coefficient in " << f[0] << " is " << min(f[0]) << '\n';
    std::cout << "partial initialization of a 4D vector with 2.3 and 4.5 is " << vec4(2.3, 4.5) << '\n';

    vec5 a5 = 203.45;
    std::cout << "initialization of a 5D vector from constant 203.45 yields " <<  a5 << '\n';


    vec3 k(f[0]);
    std::cout << "Initializing a double vector with a float vector " << f[0] << " produced " << k << '\n';

    std::cout << "The sum of double vector " << d[0] << " and a float vector " << f[0] << " is " << d[0]+f[0] << '\n';
    std::cout << "The sum of a float vector " << f[1] << " and an int vector " << j[0] << " is " << f[1]+j[0] << '\n';
    std::cout << "The product of a double vector " << d[1] << " and an int vector " << j[1] << " is " << d[1] * j[1] << '\n';

    std::sort(d, d+4, lorder);
    std::cout << "after sorting, the vectors are:\n";
    for (int i=0; i<4; ++i) std::cout << d[i] << '\n';

    d[2][1] *= -1;
    std::cout << "the absolute value of " << d[2] << " is " << abs(d[2]) << '\n';

    std::cout << "the square of " << d[3] << " is " << square(d[3]) << '\n';

    /*
    mat2 A = mat2::random(0,1);
    mat2 B = mat2::random(0,1);

    small_vector<double, 16> av = small_vector<double, 16>::random(0,1);
    typedef small_square_matrix<double, 4> sq4;
    static_assert(sq4::base_type::_size_ == 16, "Size of vector is wrong");

    sq4 C(av);

    sq4 D = sq4::random(0, 1);

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
    
    std::cout << "F is \n" << F << '\n';
    
    if (check_rowwise_iterator(F, F.begin<sq4::rowwise_iterator>(), F.end<sq4::rowwise_iterator>()))
    {
        std::cout << "Rowwise iterator check PASSED\n";
    }
    else {
        std::cout << "Rowwise iterator check FAILED\n";
    }
    
    if (check_colwise_iterator(F, F.begin<sq4::colwise_iterator>(), F.end<sq4::colwise_iterator>()))
    {
        std::cout << "Colwise iterator check PASSED\n";
    }
    else {
        std::cout << "Colwise iterator check FAILED\n";
    }
    
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
    
    vec3 real3 = vec3::random(-5. , 5.);
    fvec4 float4 = fvec4::random();
    ivec5 int5 = ivec5::random(-10, 10);
    lvec6 long6 = lvec6::random();
    
    mat4 double44 = mat4::random(-2. , 2.);
    imat5 int55 = imat5::random(-12, 12);
    
    std::cout << "random vec3 in (-5, 5)=" << real3 << '\n';
    std::cout << "random fvec4 in default range=" << float4 << '\n';
    std::cout << "random ivec5 in (-10, 10)=" << int5 << '\n';
    std::cout << "random lvec6 in default range=" << long6 << '\n';
    std::cout << "random mat4 in (-2, 2)=" << double44 << '\n';
    std::cout << "random imat5 in (-12, 12)=" << int55 << '\n';
    
    */
    {
        typedef small_matrix_interface<DataBlock<double, 4>, 2, 2> alt_mat2;
        typedef small_matrix_interface<DataBlock<float, 4>, 2, 2> alt_fmat2;
        typedef small_matrix_interface<DataBlock<double, 16>, 4, 4> alt_mat4;
        typedef small_matrix_interface<DataBlock<float, 16>, 4, 4> alt_fmat4;
        typedef small_matrix_interface<DataBlock<int, 25>, 5, 5> alt_imat5;
        alt_mat2 _A = alt_mat2::random(0, 1);
        alt_mat2 _B = alt_mat2::random(0, 1);

        alt_mat4 _C = alt_mat4::random(0, 1);
        alt_mat4 _D = alt_mat4::random(0, 1);

        std::cout << "A=" << _A << ", B=" << _B << ", C=" << _C << ", D=" << _D << '\n';
        auto _E = _A*_B; 
        auto _F = _C*inverse(_D);
        std::cout << "AxB=" << _E << '\n';
        std::cout << "A+B=" << _A+_B << '\n';
        std::cout << "A-B=" << _A-_B << '\n';
        std::cout << "A/B=" << _A/_B << '\n';
        std::cout << "4. + C=" << 4. + _C << '\n';
        std::cout << "int(4) + C=" << (int(4) + _C) << '\n';
        auto blah = int(4) + _C;
        std::cout << "after storage value is " << blah << '\n';
        std::cout << "C + 4.=" << _C + 4. << '\n';
        std::cout << "C + int(4)=" << _C + int(4) << '\n';
        std::cout << "C * 4.=" << _C*4. << '\n';
        std::cout << "C * int(4)=" << _C*int(4) << '\n';
        std::cout << "4. * C=" << 4.*_C << '\n';
        std::cout << "int(4)*C=" << (int(4)*_C) << '\n';
        auto blih = int(4) * _C;
        std::cout <<"after storage value is " << blih << "\n";
        std::cout << "C / 4.=" << _C/4. << '\n';
        std::cout << "C / int(4)=" << _C/int(4) << '\n';
        std::cout << "b=" << b << '\n';
        std::cout << "C + b=" << _C + b << '\n';
        std::cout << "C - b=" << _C - b << '\n';
        std::cout << "C.*b=" << _C.linalg(false)*b << '\n';
        
        
        std::cout << "C*D^{-1}=" << _F << "\n";
        std::cout << "C*D^{-1}*D=" << _C*inverse(_D)*_D << '\n';
        std::cout << "C*D^{-1}*D*C^{-1}=" << _C*inverse(_D)*_D*inverse(_C) << '\n';
        std::cout << "F is \n" << _F << '\n';
        
        std::cout << "the output of C*D^{-1}*D*C^{-1} has size " << whattypeami((_C*inverse(_D)*_D*inverse(_C))(0,0)) << '\n';
    
        if (check_rowwise_iterator(_F, _F.begin<alt_mat4::rowwise_iterator>(), _F.end<alt_mat4::rowwise_iterator>()))
        {
            std::cout << "Alt rowwise iterator check PASSED\n";
        }
        else {
            std::cout << "Alt rowwise iterator check FAILED\n";
        }
        
        if (check_colwise_iterator(_F, _F.begin<alt_mat4::columnwise_iterator>(), _F.end<alt_mat4::columnwise_iterator>()))
        {
            std::cout << "Alt colwise iterator check PASSED\n";
        }
        else {
            std::cout << "Alt colwise iterator check FAILED\n";
        }
    
        std::cout << "iterating over row 2 of matrix F\n";
        for (auto it=_F.row(2).begin(); it!=_F.row(2).end(); ++it) std::cout << *it << ", ";
        std::cout << '\n';
        std::cout << "printing row 2 directly:\n" << _F.row(2) << '\n';
        std::cout << "iterating over column 1 of matrix F\n";
        for (auto it=_F.column(1).begin(); it!=_F.column(1).end(); ++it) std::cout << *it << ", ";
        std::cout << '\n';
        std::cout << "printing column 1 directly:\n" << _F.column(1) << '\n';
    
        std::cout << "replacing col 3 of F with col 1 of F\n";
        std::copy(_F.column(3).begin(), _F.column(3).end(), _F.column(1).begin());
        std::cout << "F is now:\n" << _F << '\n';
    
        std::cout << "Using operator= to replace column 2 of F by column 0 of F\n";
        std::cout << "performing " << _F.column(2) << " = " << _F.column(0) << '\n';
    
        std::cout << "now\n";
        _F.column(2) = _F.column(0);
        std::cout << "F is now:\n" << _F << '\n';
    
        std::cout << "same thing for row 0 and row 3\n";
        _F.row(3) = _F.row(0);
        std::cout << "F is now\n" << _F << '\n';
    
        alt_mat4 double44 = alt_mat4::random(-2. , 2.);
        alt_imat5 int55 = alt_imat5::random(-12, 12);
    
        std::cout << "random mat4 in (-2, 2)=" << double44 << '\n';
        std::cout << "random imat5 in (-12, 12)=" << int55 << '\n';
    }

    return 0;
}