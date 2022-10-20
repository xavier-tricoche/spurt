#ifndef __MAPS_LIB_FIXPOINT_HPP__
#define __MAPS_LIB_FIXPOINT_HPP__

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <iostream>
#include "newton.hpp"
#include <math/angle.hpp>
#include <map>

namespace spurt {

extern bool record_search_steps;
extern std::vector<std::pair<nvis::vec2, nvis::vec2> > search_steps;

struct fixpoint {
    fixpoint() : isolated(true) {}
    fixpoint(const fixpoint& fp)
        : pos(fp.pos), saddle(fp.saddle), K(fp.K), isolated(fp.isolated) {
        evec[0] = fp.evec[0];
        evec[1] = fp.evec[1];
        eval[0] = fp.eval[0];
        eval[1] = fp.eval[1];
    }
    
    nvis::vec2 pos;
    bool saddle;
    nvis::vec2 evec[2];
    double eval[2];
    unsigned int K;
    bool isolated;
};

inline nvis::vec2 nullspace(const nvis::mat2& A)
{
    nvis::vec2 row1(A[0]);
    nvis::vec2 row2(A[1]);
    nvis::vec2 e;
    if (nvis::norm(row1) > nvis::norm(row2)) {
        e = nvis::vec2(-row1[1], row1[0]);
    } else {
        e = nvis::vec2(-row2[1], row2[0]);
    }
    e /= nvis::norm(e);
    return e;
}

inline bool eigen(nvis::vec2 evecs[2], double evals[2], const nvis::mat2& J)
{
    double tr = nvis::trace(J);
    double det = nvis::det(J);
    double delta = tr * tr - 4 * det;
    
    if (delta >= 0) { // real eigenvalues: saddle type
        evals[0] = 0.5 * (tr - sqrt(delta)); // small eigenvalue first
        evals[1] = 0.5 * (tr + sqrt(delta));
        
        nvis::mat2 A(J);
        A[0][0] -= evals[0];
        A[1][1] -= evals[0];
        evecs[0] = nullspace(A);
        
        A = J;
        A[0][0] -= evals[1];
        A[1][1] -= evals[1];
        evecs[1] = nullspace(A);
        
        return true;
    }
    
    return false;
}

std::ostream& operator<<(std::ostream& os, const fixpoint& fp)
{
    os << "fixed point: \tpos=" << fp.pos
       << ", \tperiod=" << fp.K
       << ", \t" << (fp.saddle ? "saddle" : "center");
    if (fp.saddle) {
        os << ", \tev0=" << fp.evec[0] << ",\t ev1=" << fp.evec[1];
    }
    
    return os;
}

void linear_analysis(const nvis::mat2& J, unsigned int period,
                     const nvis::vec2& x, spurt::fixpoint& fp)
{
    fp.pos = x;
    fp.K = period;
    fp.saddle = eigen(fp.evec, fp.eval, J);
}

bool similar(const spurt::fixpoint& fp0, const spurt::fixpoint& fp1)
{
    if ((fp0.saddle && !fp1.saddle) || (fp1.saddle && !fp0.saddle)) {
        return false;
    } else if (!fp0.saddle) {
        return true;
    } else if (fabs(nvis::inner(fp0.evec[0], fp1.evec[0])) < 0.95) {
        return false;
    } else if (fabs(nvis::inner(fp0.evec[1], fp1.evec[1])) < 0.95) {
        return false;
    } else {
        return true;
    }
}

template<typename MAP>
bool linear_analysis(const MAP& map, unsigned int period, const default_metric_type& metric,
                     const nvis::vec2& x, spurt::fixpoint& fp,
                     double hmin, double hmax)
{
    rhs_only_wrapper<MAP> rhs(map, metric, period);
    nvis::mat2 J_coarse = rhs.jacobian(x, hmax);
    nvis::mat2 J_fine = rhs.jacobian(x, hmin);
    spurt::fixpoint fp_coarse, fp_fine;
    linear_analysis(J_coarse, period, x, fp_coarse);
    linear_analysis(J_fine, period, x, fp_fine);
    std::ostringstream os;
    os << "fp_coarse = " << fp_coarse << ", fp_fine = " << fp_fine << std::endl;
    std::cerr << os.str();
    if (similar(fp_coarse, fp_fine)) {
        fp = fp_fine;
        return true;
    } else {
        nvis::mat2 J_mid = rhs.jacobian(x, 0.5*(hmin + hmax));
        spurt::fixpoint fp_mid;
        linear_analysis(J_mid, period, x, fp_mid);
        if (similar(fp_coarse, fp_mid)) {
            fp = fp_mid;
        } else if (similar(fp_mid, fp_fine)) {
            fp = fp_fine;
        } else {
            fp = fp_mid;
        }
        return false;
    }
}

struct map_quad {

    map_quad() {}
    
    nvis::vec2 pos(int i) const {
        switch (i) {
            case 0:
                return _bounds.min();
            case 1:
                return nvis::vec2(_bounds.max()[0], _bounds.min()[1]);
            case 2:
                return _bounds.max();
            case 3:
                return nvis::vec2(_bounds.min()[0], _bounds.max()[1]);
            default: throw std::runtime_error("Invalid index in map_quad::pos()");
        }
    }
    
    nvis::bbox2& bounds() {
        return _bounds;
    }
    
    nvis::vec2& val(int i) {
        return _v[i];
    }
    
    nvis::bbox2     _bounds;
    nvis::vec2      _v[4];
};

void split_box(const nvis::bbox2& box, std::vector<nvis::bbox2>& subs)
{
    nvis::vec2 center = box.center();
    const nvis::vec2& min = box.min();
    const nvis::vec2& max = box.max();
    nvis::bbox2 sub(min, center);
    subs.push_back(sub);
    sub.min() = center;
    sub.max() = max;
    subs.push_back(sub);
    sub.min() = nvis::vec2(center[0], min[1]);
    sub.max() = nvis::vec2(max[0], center[1]);
    subs.push_back(sub);
    sub.min() = nvis::vec2(min[0], center[1]);
    sub.max() = nvis::vec2(center[0], max[1]);
    subs.push_back(sub);
}

template<typename RHS>
void search(const RHS& rhs, map_quad& quad, int depth, bool must_proceed,
            std::vector<nvis::vec2>& found)
{

    // std::cerr << "depth = " << depth << ", bounds = " << quad.bounds() << std::endl;
    
    double theta = 0;
    bool split = false;
    for (int e=0 ; e<4 && !split ; ++e) {
        double dtheta = signed_angle(quad.val(e), quad.val((e+1)%4));
        if (dtheta > 2./3.*M_PI) {
            split = true;
        } else {
            theta += dtheta;
        }
    }
    
    if (!split) {
        long int idx = lrint(0.5*theta/M_PI);
        // if (record_search_steps) {
        // std::cerr << "#0 index = " << idx << "\n";
        // }
        if (idx == 0 && !must_proceed) {
            // std::cerr << "#0 index is zero\n";
            return;
        }
    }
    
    if (depth == 0) {
        // std::cerr << "#1 depth is zero\n";
        
        found.push_back(0.5*(quad.pos(0) + quad.pos(2)));
        return;
    }
    
    // std::cerr << "#1" << std::endl;
    
    nvis::vec2 midv[5];
    try {
        for (int i=0 ; i<4 ; ++i) {
            nvis::vec2 x = 0.5*(quad.pos(i) + quad.pos((i+1)%4));
            // std::cerr << "calling rhs at " << x << std::endl;
            midv[i] = rhs(x);
            if (record_search_steps) {
                search_steps.push_back(std::pair<nvis::vec2, nvis::vec2>(x, midv[i]));
            }
        }
        midv[4] = rhs(0.5*(quad.pos(0) + quad.pos(2)));
    } catch(std::runtime_error& err) {
        // std::cerr << "exception caugh: " << err.what() << std::endl;
        return;
    }
    if (record_search_steps) {
        search_steps.push_back(std::pair<nvis::vec2, nvis::vec2>(0.5*(quad.pos(0) + quad.pos(2)), midv[4]));
    }
    
    // std::cerr << "#2" << std::endl;
    
    map_quad subquads[4];
    subquads[0].bounds().min() = quad.pos(0);
    subquads[0].bounds().max() = 0.5*(quad.pos(0) + quad.pos(2));
    subquads[0].val(0) = quad.val(0);
    subquads[0].val(1) = midv[0];
    subquads[0].val(2) = midv[4];
    subquads[0].val(3) = midv[3];
    
    subquads[1].bounds().min() = 0.5*(quad.pos(0) + quad.pos(1));
    subquads[1].bounds().max() = 0.5*(quad.pos(1) + quad.pos(2));
    subquads[1].val(0) = midv[0];
    subquads[1].val(1) = quad.val(1);
    subquads[1].val(2) = midv[1];
    subquads[1].val(3) = midv[4];
    
    subquads[2].bounds().min() = 0.5*(quad.pos(0) + quad.pos(2));
    subquads[2].bounds().max() = quad.pos(2);
    subquads[2].val(0) = midv[4];
    subquads[2].val(1) = midv[1];
    subquads[2].val(2) = quad.val(2);
    subquads[2].val(3) = midv[2];
    
    subquads[3].bounds().min() = 0.5*(quad.pos(0) + quad.pos(3));
    subquads[3].bounds().max() = 0.5*(quad.pos(2) + quad.pos(3));
    subquads[3].val(0) = midv[3];
    subquads[3].val(1) = midv[4];
    subquads[3].val(2) = midv[2];
    subquads[3].val(3) = quad.val(3);
    
    // std::cerr << "subquad #1\n";
    search(rhs, subquads[0], depth-1, false, found);
    // std::cerr << "subquad #2\n";
    search(rhs, subquads[1], depth-1, false, found);
    // std::cerr << "subquad #3\n";
    search(rhs, subquads[2], depth-1, false, found);
    // std::cerr << "subquad #4\n";
    search(rhs, subquads[3], depth-1, false, found);
}

struct Lt_pos_epsilon {
    Lt_pos_epsilon(double epsilon) : _eps(epsilon), _Lt() {}
    
    bool operator()(const nvis::vec2& x0, const nvis::vec2& x1) const {
        if (nvis::norm(x0-x1) < _eps) {
            return false;
        } else {
            return _Lt(x0, x1);
        }
    }
    
    double _eps;
    nvis::lexicographical_order _Lt;
};

template<typename RHS>
bool find_seed(const RHS& rhs, const default_metric_type& metric, const nvis::bbox2& bounds,
               nvis::vec2& first_guess, int depth)
{

    // std::cerr << "entering find_seed in cell #" << cell_id << std::endl;
    
    std::map<double, nvis::vec2> norm_to_location;
    
    std::vector<nvis::bbox2> cur_res;
    cur_res.push_back(bounds);
    
    try {
        for (int d=0 ; d<depth ; ++d) {
            std::vector<nvis::bbox2> next_res;
            for (int i=0 ; i<cur_res.size() ; ++i) {
                nvis::vec2 x = cur_res[i].center();
                nvis::vec2 v = rhs(x);
                norm_to_location.insert(std::pair<double, nvis::vec2>(nvis::norm(v), x));
                split_box(cur_res[i], next_res);
            }
            std::swap(next_res, cur_res);
        }
    } catch(...) {
        return false;
    }
    
    first_guess = norm_to_location.begin()->second;
    // std::ostringstream os;
    // os << "found best first guess at " << first_guess << " with norm "
    //    << norm_to_location.begin()->first << std::endl;
    // std::cerr << os.str();
    return true;
}

template<typename MAP>
bool meta_newton(const MAP& pmap, const default_metric_type& metric, const nvis::bbox2& bounds,
                 const nvis::vec2& first_guess, int depth,
                 int period, fixpoint& fp, std::vector<nvis::vec2>& iterates,
                 double eps, double Jeps=0, bool verbose=false,
                 size_t maxiter=50, bool prefound=false)
{

    // std::cerr << "meta-newton: cell #" << cell_id << ", guess = " << first_guess << ", period = " << period << '\n';
    
    rhs_only_wrapper<MAP> rhs(pmap, metric, period);
    iterates.resize(period);
    nvis::vec2 x=first_guess, f;
    try {
        f = rhs(x);
    } catch(...) {
        // std::ostringstream os;
        // os << "unable to compute map at seed point" << std::endl;
        // std::cerr << os.str();
        return false;
    }
    
    nvis::mat2 Jinit = rhs.jacobian(x);
    double dinit = nvis::det(Jinit);
    
    if (nvis::norm(f) > eps) {
        if (verbose) {
            std::cerr << "norm at seed point (" << first_guess << ") = " << nvis::norm(f) << " is too large\n";
        }
        bool ok = find_seed(rhs, metric, bounds, x, depth);
        if (!ok) {
            if (verbose) {
                std::cerr << "unable to find seed\n";
            }
            return false;
        } else if (verbose) {
            std::cerr << "initial Jacobian was: " << Jinit << " with determinant " << dinit << '\n';
            std::cerr << "improved seed found\n"
                      << "norm at new seed (" << x << ") is " << nvis::norm(rhs(x)) << '\n';
            for (int k=0 ; k<8 ; ++k) {
                double h = 0.1 / (double)(1 << k);
                nvis::mat2 Jfin = rhs.jacobian(x, h);
                double dfin = nvis::det(Jfin);
                std::cerr << "final Jacobian (h=" << h << ") is " << Jfin << " with determinant " << dfin << '\n';
            }
            
            nvis::mat2 Jrichardson = richardson(rhs, x, 0.1, 10.);
            std::cerr << "Richardson gives us " << Jrichardson << " with determinant "
                      << nvis::det(Jrichardson) << '\n';
                      
        }
        try {
            nvis::vec2 g = rhs(x);
            if (nvis::norm(f) < nvis::norm(g)) {
                x = first_guess;
            }
        } catch(...) {
            if (verbose) {
                std::cerr << "exception caught in metanewton\n";
            }
            return false;
        }
    }
    
    bool found = newton_in_box(rhs, x, bounds, eps, maxiter, verbose);
    // bool found = newton(rhs, x, 1.0e-3, maxiter, 0.5*nvis::norm(cell_bounds.size()), true);
    // std::vector<nvis::vec2> hist;
    // nvis::mat2 J = rhs.jacobian(x);
    // bool found = broyden(rhs, x, J, 1.0e-3, 50, hist);
    if (!found && !prefound) {
        // std::ostringstream os;
        // os << "FAILED" << std::endl;
        // std::cerr << os.str();
        return false;
    }
    if (!bounds.inside(x)) {
        // std::ostringstream os;
        // os << "singularity found in wrong cell (" << idf
        //    << " instead of " << cell_id << ")" << std::endl;
        // std::cerr << os.str();
        
        // something has to be done here
    } else {
        // std::ostringstream os;
        // os << "SUCCESSFUL" << std::endl;
        // std::cerr << os.str();
    }
    
    try {
        if (found) {
            MAP* amap = pmap.clone();
            rhs_wrapper<MAP> jrhs(pmap, metric, period, Jeps);
            iterates[0] = x;
            for (int i=1 ; i<period ; ++i) {
                iterates[i] = metric.modulo(amap->map(iterates[i-1], 1));
            }
            fp.pos = x;
            nvis::mat2 J = jrhs.jacobian(x);
            linear_analysis(J, period, x, fp);
        } else if (prefound) {
            MAP* amap = pmap.clone();
            rhs_wrapper<MAP> jrhs(pmap, metric, period, Jeps);
            fp.pos = x;
            nvis::mat2 J = jrhs.jacobian(x);
            linear_analysis(J, period, x, fp);
            iterates.clear(); // we won't trust iterates from this location
        }
    } catch(...) {
        return false;
    }
    
    return true;
}

template<typename MAP>
bool meta_newton_stdmap(const MAP& pmap, const default_metric_type& metric, const nvis::bbox2& bounds,
                        const nvis::vec2& first_guess, int depth,
                        int period, fixpoint& fp, std::vector<nvis::vec2>& iterates,
                        double eps, bool verbose=false,
                        size_t maxiter=50, bool prefound=false)
{

    if (verbose) {
        std::cerr << "meta-newton: guess = " << first_guess << ", period = " << period << '\n';
    }
    
    rhs_wrapper<MAP> rhs(pmap, metric, period);
    iterates.resize(period);
    nvis::vec2 x=first_guess, f;
    try {
        f = rhs(x);
    } catch(...) {
        if (verbose) {
            std::ostringstream os;
            os << "unable to compute map at seed point" << std::endl;
            std::cerr << os.str();
        }
        return false;
    }
    
    if (nvis::norm(f) > eps) {
        bool ok = find_seed(rhs, metric, bounds, x, depth);
        if (!ok) {
            return false;
        }
        try {
            nvis::vec2 g = rhs(x);
            if (nvis::norm(f) < nvis::norm(g)) {
                x = first_guess;
            }
        } catch(...) {
            return false;
        }
    }
    
    bool found = newton_in_box(rhs, x, bounds, eps, maxiter, verbose);
    // bool found = newton(rhs, x, 1.0e-3, maxiter, 0.5*nvis::norm(cell_bounds.size()), true);
    // std::vector<nvis::vec2> hist;
    // nvis::mat2 J = rhs.jacobian(x);
    // bool found = broyden(rhs, x, J, 1.0e-3, 50, hist);
    if (!found && !prefound) {
        if (verbose) {
            std::ostringstream os;
            os << "FAILED" << std::endl;
            std::cerr << os.str();
        }
        return false;
    }
    if (!bounds.inside(x)) {
        if (verbose) {
            std::ostringstream os;
            os << "left boundaries" << std::endl;
            std::cerr << os.str();
        }
    } else {
        if (verbose) {
            std::ostringstream os;
            os << "SUCCESSFUL" << std::endl;
            std::cerr << os.str();
        }
    }
    
    if (verbose && found) {
        std::cerr << "\n\nfound a critical point of period " << period << " at " << x << '\n'
                  << "processing iterates\n";
    }
    
    try {
        if (found) {
            MAP* amap = pmap.clone();
            rhs_wrapper<MAP> jrhs(pmap, metric, period, 0);
            iterates[0] = x;
            for (int i=1 ; i<period ; ++i) {
                iterates[i] = metric.modulo(amap->map(iterates[i-1], 1));
                if (verbose) {
                    std::cerr << "iterate #" << i << " from " << x << " is " << iterates[i] << std::endl;
                }
            }
            fp.pos = x;
            nvis::mat2 J = jrhs.jacobian(x);
            linear_analysis(J, period, x, fp);
            if (verbose) {
                std::cerr << "found " << fp << std::endl;
            }
        } else if (prefound) {
            MAP* amap = pmap.clone();
            rhs_wrapper<MAP> jrhs(pmap, metric, period, 0);
            fp.pos = x;
            nvis::mat2 J = jrhs.jacobian(x);
            linear_analysis(J, period, x, fp);
            iterates.clear(); // we won't trust iterates from this location
        }
    } catch(...) {
        return false;
    }
    
    return true;
}

template<typename MAP>
bool linear_chain_analysis(const MAP& pmap, const default_metric_type& metric,
                           const nvis::vec2& first_guess, int period,
                           std::vector<fixpoint>& fps, double maxlength,
                           double eps, double Jeps=0, bool verbose=false,
                           size_t maxiter=50)
{
    rhs_wrapper<MAP> rhs(pmap, metric, period, Jeps);
    fps.resize(period);
    
    nvis::vec2 x = first_guess;
    
    bool saddle;
    for (int i=0 ; i<period ; ++i) {
        // find precise location of this fixed point
        bool found = newton(rhs, x, eps, maxiter, maxlength, verbose);
        if (!found) {
            return false;
        }
        
        fixpoint& fp = fps[i];
        nvis::mat2 J = rhs.jacobian(x);
        linear_analysis(J, period, x, fp);
        // check type consistency
        if (!i) {
            saddle = (fp.saddle ? 1 : 0);
        } else if ((saddle && !fp.saddle) || (!saddle && fp.saddle)) {
            std::ostringstream os;
            os << "incompatible linear types for period " << period << " - giving up" << std::endl;
            std::cerr << os.str();
            return false;
        }
        if (i<period-1) {
            x = metric.modulo(pmap.map(x, 1));
        }
    }
    
    return true;
}

}

#endif
