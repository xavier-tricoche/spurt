#ifndef __EIGENVECTOR_FIELD_HPP__
#define __EIGENVECTOR_FIELD_HPP__

#include <vector>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <teem/ell.h>
#include <math/math.hpp>

template<typename SymTensor7_>
inline void sym7to9(double* nine, const SymTensor7_& seven)
{
    const int ref[9] = { 0, 1, 2, 1, 3, 4, 2, 4, 5 };
    for (int i = 0 ; i < 9 ; ++i) {
        nine[i] = seven[1+ref[i]];
    }
}

template <typename TensorField>
class EigenvectorField {
    int eigensolver(double* eval, double* evec, const nvis::vec3& x) const {
        // if (!bounds().inside(x)) throw std::runtime_error("1");  // out of bounds
        
        double _t[9];
        sym7to9(_t, tfield(x));
        return ell_3m_eigensolve_d(eval, evec, _t, 0); // teem/ell library call
    }
    
public:
    EigenvectorField(const TensorField& field)
        : tfield(field) {}
        
    const nvis::bbox3 bounds() const {
        return tfield.bounds();
    }
    
    nvis::vec3 evec(const nvis::vec3& x, int which) const {
        assert(which >= 0 && which < 3);
        double eval[3], evec[9];
        int ans = eigensolver(eval, evec, x);
        if (!ans || ans == ell_cubic_root_triple) {
            throw std::runtime_error("100"); // singular tensor
        }
        int off = 3*which;
        return nvis::vec3(evec[off+0], evec[off+1], evec[off+2]);
    }
    
    double eval(const nvis::vec3& x, int which) const {
        assert(which >= 0 && which < 3);
        double eval[3], evec[9];
        int ans = eigensolver(eval, evec, x);
        if (!ans || ans == ell_cubic_root_triple) {
            throw std::runtime_error("100"); // singular tensor
        }
        return eval[which];
    }
    
    nvis::vec4 eigen(const nvis::vec3& x, int which) const {
        assert(which >= 0 && which < 3);
        double eval[3], evec[9];
        int ans = eigensolver(eval, evec, x);
        if (!ans || ans == ell_cubic_root_triple) {
            throw std::runtime_error("100"); // singular tensor
        }
        int off = 3*which;
        return nvis::vec4(evec[off+0], evec[off+1], evec[off+2], eval[which]);
    }
    
private:
    const TensorField&  tfield;
};


#endif












