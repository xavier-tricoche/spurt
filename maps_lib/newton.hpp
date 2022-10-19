#ifndef __MAPS_LIB_NEWTON_HPP
#define __MAPS_LIB_NEWTON_HPP

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <vector>
#include <assert.h>
#include <iostream>
#include <sstream>
#include <poincare/metric.hpp>

namespace xavier {
template<typename Func>
inline nvis::vec2 cdd(const Func& f, const nvis::vec2& x, double h,  int dim)
{
    nvis::vec2 dx(0);
    dx[dim] = h;
    return 0.5*(f(x+dx)-f(x-dx))/h;
}

template<typename Func>
inline nvis::mat2 richardson(const Func& f, const nvis::vec2& x, double h, double t=2.)
{
    double ht = h/t;
    double t2 = t*t;
    
    nvis::vec2 ddxh = cdd(f, x, h, 0);
    nvis::vec2 ddyh = cdd(f, x, h, 1);
    nvis::vec2 ddxht = cdd(f, x, ht, 0);
    nvis::vec2 ddyht = cdd(f, x, ht, 1);
    
    nvis::vec2 ddx = (t2*ddxht - ddxh)/(t2-1);
    nvis::vec2 ddy = (t2*ddyht - ddyh)/(t2-1);
    
    nvis::mat2 J;
    J(0,0) = ddx[0];
    J(0,1) = ddy[0];
    J(1,0) = ddx[1];
    J(1,1) = ddy[1];
    return J;
}

extern bool record_newton_steps;
extern std::vector<nvis::vec2> newton_steps;

template<typename MAP>
struct rhs_only_wrapper {
    rhs_only_wrapper(const MAP& pmap, const default_metric_type& metric, int period)
        : _map(pmap), _metric(metric), _p(period), _h(0.1) {}
        
    void set_h(double h) {
        _h = h;
    }
    
    nvis::vec2 operator()(const nvis::vec2& x) const {
        nvis::vec2 y = _map.map(x, _p);
        return _metric.displacement(x, y);
    }
    
    nvis::mat2 jacobian(const nvis::vec2& x, double h=0) const {
        if (!h) {
            h = _h;
        }
        // nvis::vec2 dx = nvis::vec2(h, 0);
        // nvis::vec2 dy = nvis::vec2(0, h);
        // nvis::fixed_vector<nvis::vec2, 2> cols;
        // cols[0] = 0.5/h*((*this)(x+dx) - (*this)(x-dx)); // d/dx
        // cols[1] = 0.5/h*((*this)(x+dy) - (*this)(x-dy)); // d/dy
        // return nvis::transpose(nvis::mat2(cols));
        return richardson(*this, x, h, 10.);
    }
    
    int _p;
    double _h;
    const MAP& _map;
    const default_metric_type _metric;
};

template<typename MAP>
struct rhs_wrapper {
    rhs_wrapper(const MAP& pmap, const default_metric_type& metric, int period, double Jeps=0)
        : _map(pmap), _metric(metric), _p(period), _jeps(Jeps),
          _lastx(std::numeric_limits<double>::max(), 0)
    {}
    
    nvis::vec2 operator()(const nvis::vec2& x) const {
        if (nvis::all(x == _lastx)) {
            return _lastv;
        } else {
            try {
                std::pair<nvis::vec2, nvis::mat2> tmp =
                    _map.map_and_jacobian(x, _p, _jeps);
                _lastx = x;
                _lastv = _metric.displacement(x, tmp.first);
                _lastJ = tmp.second;
                _lastJ[0][0] -= 1;
                _lastJ[1][1] -= 1;
            } catch(std::runtime_error& e) {
                _lastx[0] = std::numeric_limits<double>::max();
                throw;
            }
            // std::cerr << "f = " << _lastv << ", J = " << _lastJ << '\n';
            return _lastv;
        }
    }
    
    nvis::mat2 jacobian(const nvis::vec2& x) const {
        if (nvis::all(x == _lastx)) {
            return _lastJ;
        } else {
            try {
                std::pair<nvis::vec2, nvis::mat2> tmp =
                    _map.map_and_jacobian(x, _p, _jeps);
                _lastx = x;
                _lastv = _metric.displacement(x, tmp.first);
                _lastJ = tmp.second;
                _lastJ[0][0] -= 1;
                _lastJ[1][1] -= 1;
            } catch(std::runtime_error& e) {
                _lastx[0] = std::numeric_limits<double>::max();
                throw;
            }
            // std::cerr << "f = " << _lastv << ", J = " << _lastJ << '\n';
            return _lastJ;
        }
    }
    
    nvis::mat2 cd_jacobian(const nvis::vec2& x, const double h=0.05) const {
        nvis::vec2 dx = nvis::vec2(h, 0);
        nvis::vec2 dy = nvis::vec2(0, h);
        nvis::fixed_vector<nvis::vec2, 2> cols;
        cols[0] = 0.5/h*((*this)(x+dx) - (*this)(x-dx)); // d/dx
        cols[1] = 0.5/h*((*this)(x+dy) - (*this)(x-dy)); // d/dy
        return nvis::transpose(nvis::mat2(cols));
    }
    
    int _p;
    double _jeps;
    const MAP& _map;
    const default_metric_type _metric;
    mutable nvis::vec2 _lastx, _lastv;
    mutable nvis::mat2 _lastJ;
};

struct box_constraint {
    box_constraint(const nvis::bbox2& box) : _box(box) {}
    
    nvis::vec2 operator()(const nvis::vec2& from, const nvis::vec2& to) const {
        return to;
        assert(_box.inside(from));
        if (_box.inside(to)) {
            return to;
        }
        std::cerr << "from = " << from << ", to = " << to << " (" << nvis::norm(from-to) << ")" << std::endl;
        std::cerr << "box = " << _box << std::endl;
        nvis::vec2 mid = 0.5*(from + to);
        return (*this)(from, mid);
    }
    
    double size() const {
        return nvis::norm(_box.size());
    }
    
    nvis::bbox2 _box;
};

template<typename RHS>
bool lnsearch(const RHS& rhs, nvis::vec2& x, nvis::vec2& f,
              const nvis::vec2& dd, double maxlength)
{
    double lambda = 1.0;
    const double alpha = 1e-4;
    
    nvis::vec2 xsave = x, fsave = f, d = norm(dd) > maxlength ? dd * maxlength / norm(dd) : dd;
    
    for (unsigned int i = 0; i < 7; ++i) {
        x = xsave + lambda * d;
        f = rhs(x);
        if (norm(f) < (1 - alpha*lambda)*norm(fsave)) {
            return true;
        }
        
        lambda *= 0.5;
    }
    
    return false;
}

template<typename RHS>
bool newton(const RHS& rhs, nvis::vec2& x,
            double eps, size_t maxiter, double maxlength, bool verbose = false)
{
    nvis::vec2 d, f; // d is destination, f is rhs
    nvis::mat2 J;
    
    if (verbose) {
        std::cerr << "entering Newton at " << x << std::endl;
    }
    
    nvis::vec2 best;
    double minnorm = std::numeric_limits<double>::max();
    bool check_improve = false;
    int nb_failed = 0;
    
    unsigned int k;
    try {
        f = rhs(x);
        double dinit = nvis::norm(f);
        for (k = 0; k < maxiter; k++) {
            if (verbose) {
                std::cerr <<"newton: k = " << k << ", norm(f) = " << norm(f) << " (eps=" << eps << ") at " << x << '\n';
            }
            if (record_newton_steps) {
                newton_steps.push_back(x);
            }
            
            double _norm = norm(f);
            if (_norm < eps) {
                check_improve = true;
            }
            if (_norm < minnorm) {
                minnorm = _norm;
                best = x;
                if (check_improve) {
                    nb_failed = 0;
                }
            } else if (check_improve) {
                ++nb_failed;
                if (nb_failed > 5) {
                    break;
                }
            }
            
            // determine local search direction
            J = rhs.jacobian(x);
            d = nvis::solve(J, J * x - f);
            
            // do a relaxation linesearch
            // (updates x and f)
            lnsearch(rhs, x, f, d - x, maxlength);
        }
        
        if (k == maxiter) {
            if (verbose) {
                double dfin = nvis::norm(rhs(x));
                std::cerr << "\t\t initial distance = " << dinit
                          << ", final distance = " << dfin << ". failed.\n";
            }
            // WARNING_MACRO(1, "\t\t initial distance = " << dinit
            // << ", final distance = " << dfin << ". failed.\n");
        }
    } catch (std::runtime_error& e) {
        if (verbose) {
            std::cerr << e.what() << std::endl;
            std::cerr << "\texception caught in Newton. current position is " << x << std::endl;
        }
        // WARNING_MACRO(0, "\texception caught in Newton. current position is " << x << std::endl);
        return false;
    }
    
    if (k == maxiter) {
        x = best;
    }
    return (minnorm < eps);
}


template<typename RHS>
bool newton_in_box(const RHS& rhs, nvis::vec2& x, const nvis::bbox2& box,
                   double eps, size_t maxiter, bool verbose = false)
{
    nvis::vec2 d, f; // d is destination, f is rhs
    nvis::mat2 J;
    box_constraint constraint(box);
    double maxlength = 0.5*constraint.size();
    
    std::ostringstream os;
    
    if (verbose) {
        os << "entering Newton in a box at " << x << std::endl;
//        os << "box = " << box << ", maxiter = " << maxiter << std::endl;
        std::cerr << os.str();
        os.clear();
        os.str("");
    }
    nvis::vec2 best;
    double minnorm = std::numeric_limits<double>::max();
    bool must_improve = false;
    int nb_failed = 0;
    
    unsigned int k;
    try {
        f = rhs(x);
        double dinit = nvis::norm(f);
        minnorm = dinit;
        for (k = 0; k < maxiter; k++) {
            if (verbose) {
                os << "newton: k = " << k << ", norm(f) = " << norm(f)
                   << " (eps=" << eps << ") at " << x << std::endl;
                std::cerr << os.str();
                os.clear();
                os.str("");
            }
            if (record_newton_steps) {
                newton_steps.push_back(x);
            }
            double _norm = norm(f);
            if (_norm < eps) {
                must_improve = true;
            }
            if (_norm < minnorm) {
                minnorm = _norm;
                best = x;
                if (must_improve) {
                    nb_failed = 0;
                }
            } else if (must_improve) {
                ++nb_failed;
                if (nb_failed > 5) {
                    break;
                }
            }
            
            // determine local search direction
            J = rhs.jacobian(x);
            d = nvis::solve(J, J * x - f);
//           if (verbose) {
//               os << "Jacobian(" << x << ") = " << J << std::endl;
//               os << "estimated fixpoint location = " << d << std::endl;
//               os << "norm at that location = " << nvis::norm(rhs(d)) << std::endl;
//               std::cerr << os.str();
//               os.clear();
//               os.str("");
//           }

            // do a relaxation linesearch
            // (updates x and f)
            nvis::vec2 save(x);
            lnsearch(rhs, x, f, d - x, maxlength);
            x = constraint(save, x);
        }
        
        if (k == maxiter) {
            if (verbose) {
                os << "\t\t initial distance = " << dinit
                   << ", final distance = " << minnorm << ". failed.\n";
                std::cerr << os.str();
                os.clear();
                os.str("");
            }
        }
    } catch (std::runtime_error& e) {
        if (verbose) {
            os << e.what() << std::endl;
            os << "\texception caught in Newton. current position is " << x << std::endl;
            std::cerr << os.str();
        }
        return false;
    }
    
    if (k == maxiter) {
        x = best;
    }
    return (minnorm < eps);
}

template<typename RHS>
bool broyden( const RHS& map, nvis::vec2& x, const nvis::mat2& J,
              double eps, unsigned int iter, std::vector<nvis::vec2>& hist )
{
    hist.clear();
    
    try {
        nvis::mat2 A;
        nvis::vec2 xsave, fsave, f, s;
        double d, tmp;
        int i,j,k;
        
        A = J;
        std::swap( A[0][1], A[1][0] );
        
        f = map( x );
        double norm = nvis::norm(f);
        
        // Main loop
        for( k=0; k<iter; k++ ) {
            hist.push_back( x );
            
            std::cout << "broyden: " << x << ' ' << f << ' ' << norm << '\n';
            
            if( norm < eps ) {
                return true;
            }
            
            s = nvis::solve( A, -1*f );
            
            d = nvis::inner( s, s );
            
            fsave = f;
            
            x += s;
            f  = map( x );
            
            // update A
            for( i=0; i<2; i++) {
                tmp = fsave[i]+f[i];
                
                for( j=0; j<2; j++ ) {
                    tmp += A[i][j]*s[j];
                }
                
                for( j=0; j<2; j++ ) {
                    A[i][j] += tmp*s[j]/d;
                }
            }
        }
    } catch(...) {
        std::cout << "broyden: abort\n";
    }
    
    return false;
}

}

#endif




































































