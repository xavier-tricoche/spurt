#ifndef __MAPSTUFF_HPP__
#define __MAPSTUFF_HPP__

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <poincare/map.hpp>
#include <poincare/newton.hpp>
#include <math/math.hpp>
#include <poincare/ls.hpp>


namespace xavier {
struct fixpoint {
    fixpoint() : isolated(true) {}
    nvis::vec2 pos;
    bool saddle;
    nvis::vec2 evec[2];
    unsigned int K;
    bool isolated;
};
}

namespace map_metric {
static nvis::bbox2 bbox;
static bool periodic[2];

inline double width()
{
    return bbox.max()[0] - bbox.min()[0];
}

inline double height()
{
    return bbox.max()[1] - bbox.min()[1];
}

inline double sign(double a)
{
    return (a < 0) ? -1 : 1;
}

inline double mindist(double a, double b)
{
    if (a > 0.5*b) {
        return a - b;
    } else if (a < -0.5*b) {
        return a + b;
    } else {
        return a;
    }
}

inline nvis::vec2 displacement(const nvis::vec2& a, const nvis::vec2& b)
{
    nvis::vec2 dx = (b - a);
    double _x = dx[0];
    double _y = dx[1];
    
    if (periodic[0]) {
        _x = sign(dx[0]) * fmod(fabs(dx[0]), width());
        _x = mindist(_x, width());
    }
    if (periodic[1]) {
        _y = sign(dx[1]) * fmod(fabs(dx[1]), height());
        _y = mindist(_y, height());
    }
    
    return nvis::vec2(_x, _y);
}

inline double distance(const nvis::vec2& a, const nvis::vec2& b)
{
    return nvis::norm(displacement(a, b));
}

inline nvis::vec2 modulo(const nvis::vec2& x)
{
    nvis::vec2 y = x - bbox.min();
    for (int i = 0 ; i < 2 ; ++i)
        if (periodic[i]) {
            double span = bbox.max()[i] - bbox.min()[i];
            y[i] = fmod(y[i], span);
            if (y[i] < 0) {
                y[i] += span;
            }
        }
        
    return y + bbox.min();
}

};

namespace {

inline double cross(const nvis::vec2& x, const nvis::vec2& y)
{
    return x[0]*y[1] - x[1]*y[0];
}

inline double vec_align(const nvis::vec2& v0, const nvis::vec2& v1, const nvis::vec2& r0, const nvis::vec2& r1)
{
    // normalize input vectors
    nvis::vec2 e0 = r0 / nvis::norm(r0);
    nvis::vec2 e1 = r1 / nvis::norm(r1);
    nvis::vec2 f0 = v0 / nvis::norm(v0);
    nvis::vec2 f1 = v1 / nvis::norm(v1);
    
    double a = cross(f0, e0) + cross(f1, e1) - cross(f0, e1) - cross(f1, e0);
    double b = cross(f0, e1) + cross(f1, e0) - 2 * cross(f0, e0);
    double c = cross(f0, e0);
    
    std::complex<double> x[2];
    int nb_roots = xavier::quadratic_equation(a, b, c, x);
    switch (nb_roots) {
        case - 1:
            return -1;
        case 0:
            return -1;
        case 1: {
            return x[0].real();
        }
        case 2: {
            if (x[0].imag()) {
                return -1;
            }
            double t0 = x[0].real();
            double t1 = x[1].real();
            // cheap way to pick one in [0,1] if any is available
            if (fabs(t0 - 0.5) < fabs(t1 - 0.5)) {
                return t0;
            } else {
                return t1;
            }
        }
        default:
            return -1;
    }
}

};

namespace poincare {

template< typename Map >
struct MyWrapper {

    MyWrapper(const Map& map, int p)
        : __map(map), __p(p), ok(true) {}
        
    nvis::vec2 operator()(const nvis::vec2& x) const {
        return map_metric::displacement(x, __map.map(x, __p));
    }
    
    nvis::vec2 step(const nvis::vec2& p0, const nvis::vec2& p1) const {
        return map_metric::displacement(p0, p1);
    }
    
    const Map& __map;
    int __p;
    bool ok;
};

template< typename Map >
unsigned int period(const Map& map, const nvis::vec2& x, unsigned int pmax,
                    double tolerance = 1.0e-6)
{
    std::cout << "period called at " << x << std::endl;
    
    nvis::vec2 y(x);
    unsigned int p;
    double dmin = std::numeric_limits<double>::max();
    unsigned int pmin = 0;
    
    std::vector< nvis::vec2 > hits;
    try {
        map.map(x, hits, pmax);
    } catch (...) {
        return 0;
    }
    for (p = 0 ; p < pmax ; ++p) {
        y = hits[p];
        double d = map_metric::distance(x, y);
        if (d < tolerance) {
            return p + 1;
        } else if (d < dmin) {
            dmin = d;
            pmin = p + 1;
        }
    }
    // unable to determine period
    std::cout << "unable to determine period at " << x
              << " - min dist was " << dmin << " (>" << tolerance << ") for period " << pmin
              << std::endl;
    return 0;
}

inline double avg_dist(const std::vector< nvis::vec2 >& x, unsigned int p)
{
    if (x.size() <= p) {
        return std::numeric_limits<double>::max();
    }
    
    double dist = 0;
    unsigned int i;
    for (i = 0 ; i < x.size() - p ; ++i) {
        dist += map_metric::distance(x[i], x[i+p]);
    }
    return dist / (double)i;
}

template< typename Map >
unsigned int best_period(const Map& map, const nvis::vec2& x, unsigned int pmax,
                         double tolerance = 1.0e-6)
{
    // return period(map, x, pmax, tolerance);
    //
    // std::cout << "period returns " << period(map, x, pmax, tolerance) << std::endl;
    
    nvis::vec2 y(x);
    unsigned int p;
    std::vector< nvis::vec2 > hits;
    
    try {
        map.map(y, hits, 3*pmax);
    } catch (...) {
        return 0;
    }
    
    std::vector< double > dist(pmax, std::numeric_limits<double>::max());
    for (unsigned int p = 0 ; p < pmax ; ++p) {
        dist[p] = avg_dist(hits, p + 1);
    }
    std::vector< unsigned int > sorted;
    xavier::sort_ids(sorted, dist, true);
    
    // std::cout << "best_period returns " << sorted[0] + 1 << std::endl;
    // std::cout << "corresponding avg distance is " << dist[sorted[0]] << std::endl;
    
    return sorted[0] + 1;
}

bool known(const std::vector< xavier::fixpoint >& fps, const nvis::vec2 x,
           unsigned int period, double tol = 1.0e-5)
{
    typedef std::vector< xavier::fixpoint >::const_iterator cst_it_type;
    for (cst_it_type it = fps.begin() ; it != fps.end() ; ++it) {
        double d = map_metric::distance(it->pos, x);
        if (d < tol && !(period % it->K)) {
            return true;
        }
    }
    
    return false;
}

void eigen(double& lmin, nvis::vec2& evmin, const nvis::vec3& H)
{
    double tr = H[0] + H[2];
    double det = H[0] * H[2] - H[1] * H[1];
    /* l^2 - tr*l + det = 0 */
    /* delta = tr*tr - 4*det */
    /* l = 0.5*(tr - sqrt(delta)) */
    double delta = tr * tr - 4.0 * det;
    lmin = 0.5 * (tr - sqrt(delta));
    double a = fabs(H[0] - lmin);
    double b = fabs(H[2] - lmin);
    if (a >= b) {
        evmin[0] = -H[1];
        evmin[1] = H[0] - lmin;
    } else {
        evmin[0] = H[2] - lmin;
        evmin[1] = -H[1];
    }
    evmin /= nvis::norm(evmin);
}

inline nvis::vec2 nullspace(const nvis::vec4& A)
{
    nvis::vec2 row1(A[0], A[1]);
    nvis::vec2 row2(A[2], A[3]);
    nvis::vec2 e;
    if (nvis::norm(row1) > nvis::norm(row2)) {
        e = nvis::vec2(-row1[1], row1[0]);
    } else {
        e = nvis::vec2(-row2[1], row2[0]);
    }
    e /= nvis::norm(e);
    return e;
}

bool eigen(nvis::vec2 evecs[2], const nvis::vec4& J)
{
    double tr = J[0] + J[3];
    double det = J[0] * J[3] - J[1] * J[2];
    
    if (det < 0) {
        double delta = tr * tr - 4 * det;
        double lmin = 0.5 * (tr - sqrt(delta));
        double lmax = 0.5 * (tr + sqrt(delta));
        
        nvis::vec4 A(J);
        A[0] -= lmin;
        A[3] -= lmin;
        evecs[0] = nullspace(A);
        
        A = J;
        A[0] -= lmax;
        A[3] -= lmax;
        evecs[1] = nullspace(A);
        
        return true;
    }
    return false;
}

template< typename Jacobian >
void linear_analysis(const Jacobian& jac, unsigned int period,
                     const nvis::vec2& x, xavier::fixpoint& fp)
{
    nvis::vec4 J = jac(x);
    fp.pos = x;
    fp.K = period;
    fp.saddle = eigen(fp.evec, J);
    if (fp.saddle) {
        std::cout << "found a saddle point at " << fp.pos << "\n";
    } else {
        std::cout << "found a center at " << fp.pos << "\n";
    }
}

template< typename Jacobian >
bool linear_chain_analysis(const Jacobian& jac,
                           const std::vector< nvis::vec2 >& pos,
                           std::vector< xavier::fixpoint >& fps)
{
    unsigned int period = pos.size();
    fps.resize(period);
    
    try {
        linear_analysis(jac, period, pos[0], fps[0]);
        
        unsigned int nb_saddles = 0;
        
        bool saddle = fps[0].saddle;
        if (saddle) {
            ++nb_saddles;
        }
        for (unsigned int i = 1 ; i < pos.size() ; ++i) {
            linear_analysis(jac, period, pos[i], fps[i]);
            if (fps[i].saddle) {
                ++nb_saddles;
            }
        }
        if (fps.size() > 1 && nb_saddles*(fps.size() - nb_saddles) > 0) {
            std::cout << "inconsistent results: "
                      << nb_saddles << " saddles and " << fps.size() - nb_saddles << " centers\n";
            return false;
        }
        
        // analysis successful
        return true;
    } catch (...) {
        std::cout << "caught exception in linear_chain_analysis" << std::endl;
        return false;
    }
}

template< typename Map, typename RHS >
bool compute_iterates(const Map& map, const RHS& rhs, unsigned int period,
                      const nvis::vec2& x,
                      std::vector< nvis::vec2 >& chain,
                      double eps, double hx, double hy, unsigned int niter = 20)
{
    try {
        chain.clear();
        chain.reserve(period);
        nvis::vec2 y = x;
        for (unsigned int i = 0 ; i < period ; ++i) {
            bool found = xavier::Map::meta_newton(rhs, y, eps, 1.0e-9, hx, hy, niter);
            if (!found) {
                std::cout << "unable to converge from iterated position" << std::endl;
                return false;
            }
            nvis::vec2 f = rhs(y);
            std::cout << "after Newton correction, iterate #"
                      << i << " of " << rhs.__p
                      << " is now associated with rhs norm "
                      << nvis::norm(f)
                      << std::endl;
            chain.push_back(y);
            y = map.map(y, 1);
        }
        
        return true;
    } catch (...) {
        std::cout << "exception caught in compute_iterates" << std::endl;
        return false;
    }
}

template< typename Map >
void transport_eigenvector(const Map& map, unsigned int period,
                           std::vector< nvis::vec2 >& out, const nvis::vec2& x,
                           const nvis::vec2& evec, double radius, bool stable)
{
    nvis::vec2 y = x + 0.5 * radius * evec;
    out.clear();
    for (unsigned int i = 1 ; i < period ; ++i) {
        y = map.map(y, (stable ? -1 : 1));
        out.push_back(y);
    }
}


template< typename PMap, typename Map >
bool robust_linear_chain_analysis(const PMap& pmap, const Map& map,
                                  const std::vector< nvis::vec2 >& pos,
                                  std::vector< xavier::fixpoint >& fps,
                                  double dx, double dy)
{
    // static MLS::CLAPACK_helper helper(500, 2, 1, 2, true);
    
    xavier::Map::jacobian_sample_pos.clear();
    xavier::Map::jacobian_sample_vals.clear();
    
    const double eps = 0.01;
    
    unsigned int period = pos.size();
    fps.resize(period);
    
    double radius = std::min(dx, dy);
    
    try {
        // compute index at first position
        std::vector< std::pair< nvis::vec2, nvis::vec2 > > samples;
        int idx = 0;
        unsigned int ref_id = 0;
        for (; ref_id < period ; ++ref_id) {
            idx = xavier::Map::poincare_index(samples, pmap, pos[ref_id], radius, M_PI / 6., 100);
            if (fabs(idx) == 1) {
                break;
            }
        }
        
        if (fabs(idx) != 1) {
            std::cout << "unable to determine valid poincare index along this chain. giving up\n";
            return false;
        }
        
        std::cout << "fix point #" << ref_id << " on chain of period " << period << " has index "
                  << idx << " hence chain is of type " << (idx > 0 ? "center" : "saddle") << std::endl;
                  
        if (idx > 0) {
            std::cout << "found " << period << " center(s)" << std::endl;
            for (unsigned int i = 0 ; i < period ; ++i) {
                fps[i].pos = pos[i];
                fps[i].saddle = false;
                fps[i].K = period;
            }
        } else {
            nvis::vec2 x = pos[ref_id];
            std::vector< nvis::vec2 > ps, fs;
            for (unsigned int i = 0 ; i < samples.size() ; ++i) {
                ps.push_back(samples[i].first);
                fs.push_back(samples[i].second);
            }
            
            // nvis::vec4 J;
            // nvis::vec2 zero;
            // xavier::ls_linear_fit(zero, J, ps, fs, x, helper);
            
            nvis::vec4 J = xavier::ls_jacobian(ps, fs, x);
            nvis::vec2 evecs[2];
            bool saddle = eigen(evecs, J);
            bool passed = false;
            if (saddle) {
                // check that eigenvectors provide reasonable alignment with flow direction
                double alpha = std::min(acos(nvis::inner(evecs[0], evecs[1])),
                                        acos(nvis::inner(evecs[0], -1.*evecs[1])));
                double mincos = cos(0.5 * alpha);
                double cos0, cos1;
                nvis::vec2 x0, x1, f0, f1;
                passed = true;
                double sign[] = { -1, 1};
                
                for (unsigned int i = 0 ; i < 2 ; ++i) {
                    // check stable eigenvector
                    x0 = pos[ref_id] + radius * evecs[i];
                    // x0 = zero + radius * evecs[i];
                    f0 = pmap(x0);
                    f0 /= nvis::norm(f0);
                    cos0 = nvis::inner(evecs[i], sign[i] * f0);
                    
                    x1 = pos[ref_id] - radius * evecs[i];
                    // x1 = zero - radius * evecs[i];
                    f1 = pmap(x1);
                    f1 /= nvis::norm(f1);
                    cos1 = nvis::inner(-1.*evecs[i], sign[i] * f0);
                    if (std::min(cos0, cos1) < mincos) {
                        passed = false;
                        break;
                    }
                }
                
                std::cout << "linear analysis of resulting set of vectors around " << pos[ref_id]
                          << " yields a saddle type: "
                          << (passed ? "check passed" : "check failed") << '\n'
                          << "eigenvectors were: " << evecs[0] << " and " << evecs[1] << '\n';
            }
            
            if (!passed) {
                // looking for position on the circle where map aligns with radial direction
                std::vector< nvis::vec2 > candidates_ev0, candidates_ev1;
                for (unsigned int i = 0 ; i < samples.size() - 1 ; ++i) {
                    double t = vec_align(samples[i].second, samples[i+1].second, samples[i].first - x, samples[i+1].first - x);
                    if (t >= -eps && t <= 1 + eps) {
                        nvis::vec2 ft = (1. - t) * samples[i].second + t * samples[i+1].second;
                        nvis::vec2 r = map_metric::displacement(x, (1. - t) * samples[i].first + t * samples[i+1].first);
                        
                        std::cout << "checking align solution found\n";
                        std::cout << "at " << x + r << " (t=" << t << "), alignment = "
                                  << nvis::inner(r, ft) / nvis::norm(r) / nvis::norm(ft) << "\n";
                                  
                        if (nvis::inner(ft, r) > 0) {
                            candidates_ev1.push_back(x + r);
                        } else {
                            candidates_ev0.push_back(x + r);
                        }
                    }
                }
                
                if (!candidates_ev0.size() || !candidates_ev1.size()) {
                    std::cout << "unable to identify positions aligning with radial directions in both directions\n";
                    return false;
                }
                
                // identify best alignmnent among candidates found for both eigendirections
                double maxdot, mindot;
                int maxdotid = 0, mindotid = 0;
                nvis::vec2 f, r;
                
                f = pmap(candidates_ev0[0]);
                r = candidates_ev0[0] - x;
                f /= nvis::norm(f);
                r /= nvis::norm(r);
                mindot = nvis::inner(f, r);
                
                for (unsigned int i = 1 ; i < candidates_ev0.size() ; ++i) {
                    f = pmap(candidates_ev0[i]);
                    r = candidates_ev0[i] - x;
                    f /= nvis::norm(f);
                    r /= nvis::norm(r);
                    double dot = nvis::inner(f, r);
                    if (dot < mindot) {
                        mindot = dot;
                        mindotid = i;
                    }
                }
                if (mindot > -0.8) {
                    std::cout << "unable to identify minor eigenvector: best alignment was " << mindot
                              << std::endl;
                    return false;
                }
                evecs[0] = candidates_ev0[mindotid] - x;
                evecs[0] /= nvis::norm(evecs[0]);
                
                f = pmap(candidates_ev1[0]);
                r = candidates_ev1[0] - x;
                f /= nvis::norm(f);
                r /= nvis::norm(r);
                maxdot = nvis::inner(f, r);
                for (unsigned int i = 1 ; i < candidates_ev1.size() ; ++i) {
                    f = pmap(candidates_ev1[i]);
                    r = candidates_ev1[i] - x;
                    f /= nvis::norm(f);
                    r /= nvis::norm(r);
                    double dot = nvis::inner(f, r);
                    if (dot > maxdot) {
                        maxdot = dot;
                        maxdotid = i;
                    }
                }
                if (maxdot < 0.8) {
                    std::cout << "unable to identify major eigenvector: best alignment was " << maxdot
                              << std::endl;
                    return false;
                }
                evecs[1] = candidates_ev1[maxdotid] - x;
                evecs[1] /= nvis::norm(evecs[1]);
            }
            
            // eigenvectors have now been found. carry them with the map
            // fps[ref_id].pos = pos[ref_id];
            // x = (passed ? zero : pos[ref_id]);
            fps[ref_id].pos = pos[ref_id];
            fps[ref_id].K = period;
            fps[ref_id].saddle = true;
            fps[ref_id].evec[0] = evecs[0];
            fps[ref_id].evec[1] = evecs[1];
            
            std::vector< nvis::vec2 > evs0, evs1;
            transport_eigenvector(map, period, evs0, x, evecs[0], radius, true);
            transport_eigenvector(map, period, evs1, x, evecs[1], radius, false);
            for (unsigned int i = 0 ; i < evs0.size() ; ++i) {
                unsigned int j = (ref_id + i + 1) % period;
                fps[j].pos = pos[j];
                fps[j].K = period;
                fps[j].saddle = true;
                fps[j].evec[0] = evs0[i] - pos[j];
                fps[j].evec[0] /= nvis::norm(fps[j].evec[0]);
                fps[j].evec[1] = evs1[i] - pos[j];
                fps[j].evec[1] /= nvis::norm(fps[j].evec[1]);
            }
            std::cout << "found " << period << " saddles" << std::endl;
        }
        
        // analysis successful
        return true;
    } catch (...) {
        std::cout << "caught exception in linear_chain_analysis" << std::endl;
        return false;
    }
}



}

#endif










































































































