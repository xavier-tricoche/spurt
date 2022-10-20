#ifndef __XAVIER_INTEGRATION_HPP__
#define __XAVIER_INTEGRATION_HPP__

#include <list>

namespace spurt {
template<typename RHS, typename T>
struct euler {

    typedef RHS     rhs_type;
    typedef T       vec_type;
    
    enum state {
        OK = 0,
        LEFT_DOMAIN,
    };
    
    euler(const rhs_type& rhs, double h) : _rhs(rhs), _h(h) {}
    
    state integrate(std::list<vec_type>& sl,
                    double length, bool record = true) const {
                    
        if (!length) {
            return OK;
        }
        
        bool fwd = length > 0;
        vec_type x = (fwd ? sl.back() : sl.front()), f, lastx, lastf = 0;
        lastx = x;
        double h = (fwd ? _h : -_h);
        int n = floor(fabs(length) / _h);
        for (int i = 0 ; i < n ; ++i) {
            if (_rhs(x, f)) {
                lastx = x;
                lastf = f;
                if (record && i) {
                    if (fwd) {
                        sl.push_back(x);
                    } else {
                        sl.push_front(x);
                    }
                }
                x += h * f;
            } else {
                // binary search
                double hmin = 0, hmax = h;
                vec_type _x;
                while (fabs(hmax - hmin) > 0.01*h) {
                    double _h = 0.5 * (hmin * hmax);
                    _x = lastx + _h * lastf;
                    if (!_rhs(_x, f)) {
                        hmax = _h;
                    } else {
                        hmin = _h;
                    }
                }
                if (fwd) {
                    sl.push_back(_x);
                } else {
                    sl.push_front(_x);
                }
                
                return LEFT_DOMAIN;
            }
        }
        double dt = length - n * h;
        if (_rhs(x, f)) {
            if (fwd) {
                sl.push_back(x);
            } else {
                sl.push_front(x);
            }
        } else {
            if (!record) {
                if (fwd) {
                    sl.push_back(x);
                } else {
                    sl.push_front(x);
                }
            }
            return LEFT_DOMAIN;
        }
        return OK;
    }
    
    state flow_map(const vec_type& in, vec_type& out, double length) const {
        std::list<vec_type> res;
        res.push_back(in);
        state s = integrate(res, length, false);
        out = (length > 0 ? res.back() : res.front());
        return s;
    }
    
    state streamline(std::list<vec_type>& io, double length) const {
        return integrate(io, length, true);
    }
    
    
    const rhs_type& _rhs;
    double          _h;
};


}

#endif






