#ifndef __DPL_HPP__
#define __DPL_HPP__

#include <limits>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <stdexcept>


class DoublePointLoad {
    
public:
    DoublePointLoad(double length, double width, double depth, double l = 1.)
        : check_inside(true) {
        _bounds.min() = nvis::vec3(-0.5 * length, -0.5 * width, -depth);
        _bounds.max() = nvis::vec3(0.5 * length, 0.5 * width, 0.);
        std::cerr << "bounds = " << _bounds << std::endl;
        const nvis::vec3& origin = _bounds.min();
        center[0] = origin + nvis::vec3(0.25 * length, 0.5 * width, depth);
        center[1] = origin + nvis::vec3(0.75 * length, 0.5 * width, depth);
        twoPi = 2.*M_PI;
        load = l;
        nu = 0.4;
    }
    
    const nvis::bbox3& bounds() const {
        return _bounds;
    }
    
    void set_check_inside(bool check) {
        check_inside = check;
    }
    
    nvis::fixed_vector<double, 7> operator()(const nvis::vec3& x) const {
        if (check_inside && !_bounds.inside(x)) {
            // std::cerr << "DoublePointLoad::tensor(): invalid position = " << x << std::endl;
            throw std::runtime_error("1"); // outside
        }
        
        nvis::fixed_vector<double, 7> t(0.);
        nvis::vec3 s;
        double P, rho, rho2, rho3, rho5;
        double rhoPlusz2, zPlus2rho, txy, txz, tyz;
        t[0] = 1.;
        
        for (unsigned int n = 0 ; n < 2 ; ++n) {
            // points are evaluated in local coordinate system of applied force.
            P = -load;
            nvis::vec3 p = center[n] - x;
            p[0] *= -1.;
            rho = nvis::norm(p);
            if (rho < 1.0e-10) {
                throw std::runtime_error("100"); // singular tensor (for practical purposes)
            }
            
            rho2 = rho * rho;
            rho3 = rho2 * rho;
            rho5 = rho2 * rho3;
            const double& x = p[0];
            const double& y = p[1];
            const double& z = p[2];
            double x2 = x * x;
            double y2 = y * y;
            double z2 = z * z;
            rhoPlusz2 = (rho + z) * (rho + z);
            zPlus2rho = (2.0 * rho + z);
            
            // normal stresses
            s[0] = P / (twoPi * rho2) * (3.0 * z * x2 / rho3 - nu * (z / rho - rho / (rho + z) +
                                         x2 * (zPlus2rho) / (rho * rhoPlusz2)));
            s[1] = P / (twoPi * rho2) * (3.0 * z * y2 / rho3 - nu * (z / rho - rho / (rho + z) +
                                         y2 * (zPlus2rho) / (rho * rhoPlusz2)));
            s[2] = 3.0 * P * z2 * z / (twoPi * rho5);
            
            //shear stresses - negative signs are coordinate transformations
            //that is, equations (in text) are in different coordinate system
            //than volume is in.
            txy = -(P / (twoPi * rho2) * (3.0 * x * y * z / rho3 -
                                          nu * x * y * (zPlus2rho) / (rho * rhoPlusz2)));
            txz = -(3.0 * P * x * z2 / (twoPi * rho5));
            tyz = 3.0 * P * y * z2 / (twoPi * rho5);
            
            t[1] += s[0];  // Component(0,0);
            t[4] += s[1];  // Component(1,1);
            t[6] += s[2];  // Component(2,2);
            t[2] += txy;   // Component(0,1);  real symmetric matrix
            t[3] += txz;   // Component(0,2);
            t[5] += tyz;   // Component(1,2);
        }
        return t;
    }
    
private:
    nvis::vec3 center[2];
    double twoPi, load, poisson_ratio, nu;
    nvis::bbox3 _bounds;
    bool check_inside;
};



#endif


























