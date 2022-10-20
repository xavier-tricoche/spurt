#ifndef __MAP_APPROXIMATION_HPP__
#define __MAP_APPROXIMATION_HPP__

#include <vector>
#include <iostream>
#include <sstream>
#include <limits>
#include <algorithm>
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/MLS.hpp>
#include <poincare/metric.hpp>

namespace spurt {

template<typename RHS>
class map_approximation {

public:
    map_approximation(const RHS& rhs, const nvis::bbox2& bounds,
                      int order, const nvis::ivec2& res)
        : _order(order), _coef(MLS::dof(2, order), 2) {
        nvis::vec2 step = bounds.size() / nvis::vec2(res[0]-1, res[1]-1);
        _x0 = 0.5*(bounds.min() + bounds.max());
        
        MLS::least_squares<nvis::vec2, nvis::vec2> LS(2, order, 2);
        
        std::vector<nvis::vec2> points(res[0]*res[1]);
        std::vector<nvis::vec2> values(res[0]*res[1]);
        for (int i=0 ; i<res[0] ; ++i) {
            for (int j=0 ; j<res[1] ; ++j) {
                int n = j+res[1]*i;
                points[n] = bounds.min() + step*nvis::vec2(i,j);
                values[n] = rhs(points[n]);
            }
        }
        
        int actual_order = LS(_coef, points, values, _x0);
    }
    
    nvis::vec2 operator()(const nvis::vec2& x) const {
        nvis::vec2 y = x - _x0;
        std::vector<double> b(MLS::dof(2,_order));
        MLS::set_basis(b, y, 2, _order);
        
        // std::cerr << "basis set\n";
        
        nvis::vec2 r(0,0);
        for (int i=0 ; i<b.size() ; ++i) {
            r[0] += b[i]*_coef[0][i];
            r[1] += b[i]*_coef[1][i];
        }
        return r;
    }
    
    nvis::mat2 jacobian(const nvis::vec2& x) const {
        nvis::vec2 y = x - _x0;
        nvis::mat2 J;
        double u = y[0];
        double v = y[1];
        const std::vector<double>& a = _coef[0];
        const std::vector<double>& b = _coef[1];
        if (_order == 3) {
            // Jxx
            J(0,0) = a[1] + 2*a[3]*u + a[4]*v + 3*a[6]*u*u + 2*a[7]*u*v + a[8]*v*v;
            // Jyx
            J(1,0) = b[1] + 2*b[3]*u + b[4]*v + 3*b[6]*u*u + 2*b[7]*u*v + b[8]*v*v;
            // Jxy
            J(0,1) = a[2] + a[4]*u + 2*a[5]*v + a[7]*u*u + 2*a[8]*u*v + 3*a[9]*v*v;
            // Jyy
            J(1,1) = b[2] + b[4]*u + 2*b[5]*v + b[7]*u*u + 2*b[8]*u*v + 3*b[9]*v*v;
        }
        
        return J;
    }
    
private:
    int                                 _order;
    nvis::vec2                          _x0;
    Eigen::MatrixXd                     _coef;
};

}


#endif