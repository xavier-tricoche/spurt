#ifndef __INVARIANT_MANIFOLD_HPP__
#define __INVARIANT_MANIFOLD_HPP__

#include <vector>
#include <list>
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <poincare/metric.hpp>
#include "fixpoints.hpp"
#include "map_analysis.hpp"



namespace xavier {
struct fp_chain {
    typedef xavier::fixpoint    fp_type;
    
    fp_chain(const std::vector<fp_type>& _fps, const map_metric& _metric)
        : fixed_points(_fps), metric(_metric) {}
        
    template<typename Iterator>
    fp_chain(Iterator begin, Iterator end, const map_metric& _metric)
        : fixed_points(begin, end), metric(_metric) {}
        
    unsigned int period() const {
        return fixed_points.size();
    }
    
    bool saddle() const {
        return fixed_points.size() && fixed_points[0].saddle;
    }
    
    const fp_type& operator[](unsigned int i) const {
        if (!period()) {
            throw std::runtime_error("empty chain");
        }
        
        unsigned int j = (i % period());
        return fixed_points[j];
    }
    
    unsigned int closest(const nvis::vec2& x) const {
        std::vector<double> dist(period());
        for (int i = 0 ; i < dist.size() ; ++i) {
            dist[i] = metric.distance(x, fixed_points[i].pos);
        }
        return std::distance(dist.begin(), std::min_element(dist.begin(), dist.end()));
    }
    unsigned int closest(const fp_type& fp) const {
        return closest(fp.pos);
    }
    
    double distance(const nvis::vec2& x) const {
        return metric.distance(x, fixed_points[closest(x)].pos);
    }
    double distance(const fp_type& fp) const {
        return distance(fp.pos);
    }
    
    double distance(const fp_chain& c) const {
        if (c.period() != period()) {
            return std::numeric_limits<double>::max();
        }
        
        std::vector<double> dist(period());
        for (int i = 0 ; i < dist.size() ; ++i) {
            dist[i] = distance(c[i]);
        }
        return *std::min_element(dist.begin(), dist.end());
    }
    
    size_t size() const {
        return fixed_points.size();
    }
    
    std::vector<xavier::fixpoint>   fixed_points;
    map_metric                      metric;
};

struct Neighbor {
    Neighbor() : chain_id(-1), p_id(-1) {}
    Neighbor(int cid, int pid) : chain_id(cid), p_id(pid) {}
    Neighbor(const Neighbor& n) : chain_id(n.chain_id), p_id(n.p_id) {}
    int chain_id;
    int p_id;
};

struct Connection {
    Connection() : neigh(), d(std::numeric_limits<double>::max()) {}
    Connection(const Neighbor& n, double dist) : neigh(n), d(dist) {}
    Neighbor neigh;
    double d;
};

struct Lt_Connection {
    bool operator()(const Connection& c0, const Connection& c1) const {
        return (c0.d < c1.d);
    }
};

typedef std::vector<nvis::vec2>     manifold_type;
typedef std::pair<int, int>         fp_index_type;

extern std::vector<manifold_type> broken_manifolds;

struct separatrix {
    fp_index_type   start, end;
    manifold_type   manifold;
    bool            forward;
    double          length;
};

//////////////////////////////////////////////////////////////////////////////
//
// England et al.'s 1D manifold construction technique
// (SIAM J. Appl. Dynamical Systems, 2005)
//
//////////////////////////////////////////////////////////////////////////////

struct manifold_progress {
    manifold_progress() : length(0) {}
    
    double length;
    double tau;
    unsigned int segment;
};

double min_dist(const std::vector<fp_chain>& chains, const xavier::map_metric& metric)
{
    std::vector<double> dist;
    std::vector<nvis::vec2> all_pos;
    
    for (int i=0 ; i<chains.size() ; ++i) {
        for (int j=0 ; j<chains[i].size() ; ++j) {
            all_pos.push_back(chains[i][j].pos);
        }
    }
    
    for (int i=0 ; i<all_pos.size()-1 ; ++i) {
        for (int j=i+1 ; j<all_pos.size() ; ++j) {
            dist.push_back(metric.distance(all_pos[i], all_pos[j]));
        }
    }
    
    return *std::min_element(dist.begin(), dist.end());
}

struct manifold_stop {

    manifold_stop(const std::vector<fp_chain>& all_chains,
                  const xavier::map_metric& metric, double eps)
        : _chains(all_chains), _metric(metric), _eps(eps),
          _left_dist(0.3*min_dist(all_chains, metric)), _left(false) {}
          
    void start(const fp_index_type fpi) {
        _fpi = fpi;
        _start = _chains[fpi.first][fpi.second].pos;
        _end_chain = _end_fp = -1;
        _min_dist_id = fp_index_type(-1, -1);
        _min_dist = std::numeric_limits<double>::max();
        _approaching = false;
        _counter = 0;
        _last = _start;
        _length = 0;
        _left = false;
    }
    
    bool operator()(const nvis::vec2& x) const {
        ++_counter;
        _length += _metric.distance(_last, x);
        // std::cerr << "counter = " << _counter << std::endl;
        if (!_left && _metric.distance(x, _start) < _left_dist) {
            if (_counter > 2000) {
                std::cerr << "unable to leave start. giving up.\n";
                throw std::runtime_error("wrong type suspected");
            }
            return false;
        } else if (!_left) {
            std::cerr << "left start\n";
            _left = true;
        }
        std::map<double, fp_index_type> dist;
        for (int i=0 ; i<_chains.size() ; ++i) {
            for (int j=0 ; j<_chains[i].size() ; ++j) {
                double d = nvis::norm(x-_chains[i][j].pos);
                dist[d] = fp_index_type(i, j);
            }
        }
        // std::cerr << "there are " << _chains.size() << " chains and " << dist.size() << " total fixed points\n";
        
        std::map<double, fp_index_type>::iterator it = dist.begin();
        double d = it->first;
        fp_index_type fp_index = it->second;
        nvis::vec2 pos = _chains[it->second.first][it->second.second].pos;
        if (!_approaching &&
                (fp_index.first != _fpi.first || fp_index.second != _fpi.second)) {
            // nvis::vec2 tangent = _metric.displacement(_last, x);
            //  nvis::vec2 direction = _metric.displacement(x, pos);
            if (d < 0.5*_metric.distance(x, _start)) {
                // std::cerr << "closest to fix point at " << pos << std::endl;
                _approaching = true;
                // double dot = nvis::inner(tangent, direction)/nvis::norm(tangent)/nvis::norm(direction);
                //                 if (dot>0.9) {
                //                 std::cerr << "approaching fixed point "
                //                            << it->second.second << " on chain "
                //                            << it->second.first << " at " << x << std::endl;
                //                  _approaching = true;
                // }
                // else {
                //  std::cerr << "cosine test failed (" << dot << ")\n";
                // }
            }
            return false;
        }
        // double d2 = (++it)->first;
        if (_approaching && d < 0.002*_length) {
            _end_chain = dist.begin()->second.first;
            _end_fp = dist.begin()->second.second;
            // std::cerr << "approaching (d=" << d << ") fixed point " << _end_fp << " on chain " << _end_chain << " while total length is " << _length << std::endl;
            return true;
        } else {
            // std::cerr << "min distance " << d << " is larger than 1\% of " << _length << std::endl;
        }
        return false;
        //         if (d < _eps || _counter > 1000) {
        //  std::cerr << "stopping\n";
        //  _end_chain = dist.begin()->second.first;
        //  _end_fp = dist.begin()->second.second;
        //  return true;
        // }
        // else if (!_approaching) {
        //  if (d < 0.5*nvis::norm(x-_start)) _approaching = true;
        //  _min_dist = d;
        //  _min_dist_id = dist.begin()->second;
        //  return false;
        // }
        //
        // if (_approaching) {
        //  if (d < _min_dist) {
        //      _min_dist = d;
        //      _min_dist_id = dist.begin()->second;
        //      return false;
        //  }
        //  else {
        //      _end_chain = dist.begin()->second.first;
        //      _end_fp = dist.begin()->second.second;
        //      return true;
        //  }
        // }
        //         return false;
    }
    
    fp_index_type connection() const {
        return std::make_pair(_end_chain, _end_fp);
    }
    
    const std::vector<fp_chain>&    _chains;
    xavier::map_metric              _metric;
    double                          _eps;
    double                          _left_dist;
    nvis::vec2                      _start;
    fp_index_type                   _fpi;
    mutable int                     _end_chain;
    mutable int                     _end_fp;
    mutable fp_index_type           _min_dist_id;
    mutable double                  _min_dist;
    mutable bool                    _approaching;
    mutable int                     _counter;
    mutable nvis::vec2              _last;
    mutable double                  _length;
    mutable bool                    _left;
};

struct manifold_stop_period_one {

    manifold_stop_period_one(const xavier::map_metric& metric, double eps)
        : _metric(metric), _eps(eps) {}
        
    void start(const nvis::vec2& x) {
        _start = x;
        _has_left = false;
        _last = _start;
        _approaching = false;
        _min_dist = std::numeric_limits<double>::max();
    }
    
    bool operator()(const nvis::vec2& x) const {
        _length += _metric.distance(x, _last);
        _last = x;
        double dist = _metric.distance(x, _start);
        if (!_has_left && dist > 50*_eps) {
            _has_left = true;
        }
        if (!_has_left) {
            return false;
        }
        if (_has_left && !_approaching && dist < 10*_eps) {
            _approaching = true;
            _nb_backtracked = 0;
        }
        if (_approaching) {
            if (dist<_min_dist) {
                _min_dist = dist;
                _nb_backtracked = 0;
            } else {
                ++_nb_backtracked;
            }
            if (_nb_backtracked > 3) {
                std::cerr << "too many backtrack steps. stop.\n";
                return true;
            }
        }
        if (_approaching && _metric.distance(x, _start) < _eps) {
            // std::cerr << "we are approaching and distance to " << _start
            //           << " is only " << _metric.distance(x, _start) << ". stop.\n";
            return true;
        }
        return false;
    }
    
    xavier::map_metric      _metric;
    double                  _eps;
    nvis::vec2              _start;
    mutable bool            _has_left;
    mutable bool            _approaching;
    mutable double          _length;
    mutable nvis::vec2      _last;
    mutable double          _min_dist;
    mutable int             _nb_backtracked;
};

template<typename Map>
inline nvis::vec2
BVP_Step(const Map& map, int period, const map_metric& metric,
         const nvis::vec2& start, const nvis::vec2& end,    // seeding segment
         const nvis::vec2& prev, const nvis::vec2& cur,     // end segment
         double delta_min, double alpha_max, double delta_alpha_max,    // quality control
         double& tau, map_analysis_param& param)    // where are we along seeding segment
{
    double delta, alpha;
    double t = 1;
    nvis::vec2 x, y, span;
    span = metric.displacement(start, end);
    
    bool verbose = param.verbose;
    verbose = false;
    
    // nvis::vec2 q = (1 - tau) * start + tau * end;
    nvis::vec2 q = metric.modulo(start + tau*span);
    
    if (verbose)
        std::cerr << "BVP called at " << q << " (" << tau << ") on segment "
                  << start << " - " << end << ", prev = " << prev
                  << ", cur = " << cur << '\n';
                  
    while (t > tau) {
        x = metric.modulo(start + t*span);
        if (verbose) {
            std::cerr << "currently trying at " << x << ". ";
        }
        nvis::vec2 _y = map.map(x, period);
        if (verbose) {
            std::cerr << "map(" << x << ") = " << _y;
        }
        y = metric.modulo(_y);
        if (verbose) {
            std::cerr << " (=" << y << ")\n";
            std::cerr << "\tat " << y << " for t=" << t << '\n';
        }
        alpha = metric.angle(prev, cur, y);
        delta = metric.distance(cur, y);
        if (verbose) {
            std::cerr << "delta = " << delta << " (" << delta_min << "), alpha = " << alpha
                      << ", delta*alpha = " << delta* alpha << " (" << delta_alpha_max << "): ";
        }
        if (fabs(delta) < delta_min ||
                (fabs(alpha) < alpha_max && fabs(delta*alpha) < delta_alpha_max)) {
            if (verbose) {
                std::cerr << "PASSED\n";
            }
            break;
        } else {
            if (verbose) {
                std::cerr << "FAILED\n";
            }
            t = 0.5 * (tau + t);
        }
        if (verbose) {
            std::cerr << "delta = " << delta << ", alpha = " << alpha << '\n';
        }
    }
    
    tau = t;
    return y;
}

template<typename MAP, typename STOP>
inline void ManBVP(manifold_type& manifold, manifold_progress& progress,
                   const MAP& map, const STOP& stop, int period, const map_metric& metric,
                   const xavier::fixpoint& saddle, bool fwd, double eps,
                   double max_length, double delta_min, double alpha_max,
                   double delta_alpha_max,  map_analysis_param& param)
{
    unsigned int i0, i1;
    double tau;
    double length;
    
    // std::cerr << "entering ManBVP" << std::endl;
    
    if (progress.length == 0) {
        if (param.verbose)
            // std::cerr << "zero length initially\n";
        {
            manifold.clear();
        }
        manifold.push_back(saddle.pos);
        
        // move along eigenvector until p-step aligns within alpha_max
        // with the eigenvector
        nvis::vec2 evec = (period > 0 ? saddle.evec[1] : saddle.evec[0]);
        nvis::vec2 p = saddle.pos;
        // std::cerr << "starting at " << p << std::endl;
        bool already_aligned = false;
        while (true) {
            p += eps * (fwd ? 1 : -1) * evec;
            nvis::vec2 dir = metric.displacement(p, map.map(p, period));
            dir /= nvis::norm(dir);
            double cosalpha = nvis::inner(dir, (fwd ? 1 : -1)*evec);
            if (!already_aligned && cosalpha > cos(alpha_max)) {
                already_aligned = true;
                // std::cerr << "alignment found (" << dir << " // " << (fwd ? 1 : -1) * evec << " at " << p << std::endl;
            } else if (already_aligned) {
                // std::cerr << "stopping at first step following alignment at " << p << std::endl;
                break;
            }
        }
        manifold.push_back(p);
        length = metric.distance(manifold[0], manifold[1]);
        if (length > max_length) {
            manifold[1] = max_length * (fwd ? 1 : -1)*evec;
            progress.length = max_length;
            return;
        }
        
        // std::cerr << "initial length between " << manifold[0] << " and "
        //               << manifold[1] << " is " << length << '\n';
        
        i0 = 0;
        i1 = 1;
        tau = 0;
    } else {
        length = progress.length;
        i0 = progress.segment;
        i1 = i0 + 1;
        tau = progress.tau;
    }
    // std::cerr << "there are currently " << manifold.size() << " positions in manifold\n";
    
    while (!stop(manifold.back())) {
    
        const nvis::vec2& q0 = manifold[manifold.size()-2];
        const nvis::vec2& q1 = manifold.back();
        const nvis::vec2& p0 = manifold[i0];
        const nvis::vec2& p1 = manifold[i1];
        
        nvis::vec2 next =
            BVP_Step(map, period, metric, p0, p1, q0, q1,
                     delta_min, alpha_max, delta_alpha_max, tau, param);
                     
        // if (metric.distance(manifold.back(), next) > 2) {
        //  std::cerr << "\n\n\n\nHUGE STEP TAKEN\n\n\n\n";
        // }
        
        if (tau == 1) {
            ++i0;
            ++i1;
            tau = 0;
        }
        
        double dl = metric.distance(manifold.back(), next);
        
        length += dl;
        
        if (length > 3*metric.diameter()/(double)fabs(period)) {
            std::cerr << "reached length " << length << " exceeds prescribed upper bound\n";
            break;
        }
        manifold.push_back(next);
        
        // if (param.verbose)
        //     std::cerr << "at " << next << ", length = " << length << '\n';
    }
    
    std::cerr << "stop criterion tested true at " << manifold.back() << std::endl;
    
    progress.length = length;
    progress.segment = i0;
    progress.tau = tau;
}

bool sanity_check(const separatrix& sep, double alpha, const map_metric& metric)
{

    // if (sep.manifold.size() < 10) return false;
    //
    // return true;
    
    double huge = 0.05*metric.diameter();
    nvis::vec2 p0, p1, p2;
    const std::vector<nvis::vec2>& m = sep.manifold;
    if (m.size() < 2) {
        std::cerr << "\n\nnot enough points at " << m[0] << "\n\n\n" << std::endl;
        return false;
    }
    int c = 0;
    p1 = m[c++];
    p2 = m[c++];
    if (metric.distance(p1, p2) > huge) {
        std::cerr << "segment too large (" << metric.distance(p1, p2) << ") at " << p0 << "\n\n\n" << std::endl;
        return false;
    }
    while (c < m.size()) {
        p0 = p1;
        p1 = p2;
        p2 = m[c++];
        if (metric.distance(p1, p2) > huge) {
            std::cerr << "segment too large (" << metric.distance(p1, p2) << ") at " << m[0] << "\n\n\n" << std::endl;
            return false;
        }        //
        // else if (fabs(metric.angle(p0, p1, p2)) > alpha) {
        //     std::cerr << "angle too large (" << metric.angle(p0, p1, p2) << ") at " << m[0] << "\n\n\n" << std::endl;
        //     return false;
        // }
    }
    return true;
}


template< typename MAP >
bool compute_separatrix(std::vector<fp_chain>& all_p_chains,
                        std::vector<separatrix>& separatrices,
                        unsigned int chain_id, const MAP& map,
                        const map_metric& metric, double _eps,
                        map_analysis_param& param)
{
    separatrices.clear();
    param.verbose = true;
    
    if (!all_p_chains.size() || chain_id >= all_p_chains.size()) {
        return true;
    }
    const fp_chain& chain = all_p_chains[chain_id];
    
    if (!chain.saddle()) {
        return true;    // that was easy!
    }
    
    int period = chain.period();
    
    if (param.verbose)
        std::cout << "\n\nprocessing saddle chain of period " << period
                  << " starting at " << chain[0].pos << std::endl;
                  
    double eps = 0.0001 * metric.diameter();
    double init_eps = 0.005*metric.diameter();
    double alpha_max = 0.1;
    double liberal_alpha_max = M_PI/2.;
    double delta_min = (period == 1 ? 1.0e-4 : 1.0e-3) * metric.diameter();
    double delta_length = metric.diameter();
    manifold_type man1, man2, man3, man4;
    
    if (period == 1) {
        manifold_progress prog;
        manifold_stop_period_one stop(metric, _eps);
        
        std::cerr << "\nunstable manifold #0..." << std::endl;
        stop.start(chain[0].pos);
        prog.length = 0;
        ManBVP(man1, prog, map, stop, period, metric, chain[0], true, eps, init_eps, delta_min, alpha_max, 1.0e-5, param);
        separatrices.push_back(separatrix());
        std::copy(man1.begin(), man1.end(), std::back_inserter(separatrices.back().manifold));
        separatrices.back().start = fp_index_type(chain_id, 0);
        separatrices.back().end = fp_index_type(chain_id, 0);
        separatrices.back().forward = true;
        separatrices.back().length = prog.length;
        if (!sanity_check(separatrices.back(), liberal_alpha_max, metric)) {
        
            broken_manifolds.push_back(separatrices.back().manifold);
            
            std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
            separatrices.back().manifold.resize(2);
            nvis::vec2& x0 = separatrices.back().manifold[0];
            nvis::vec2& x1 = separatrices.back().manifold[1];
            x0 = chain[0].pos;
            x1 = x0 + 0.001 * metric.displacement(x0, x1) / metric.distance(x0, x1);
        }
        std::cerr << "length = " << prog.length << ", " << man1.size() << " vertices\n";
        
        std::cerr << "\nunstable manifold #1..." << std::endl;
        stop.start(chain[0].pos);
        prog.length = 0;
        ManBVP(man2, prog, map, stop, period, metric, chain[0], false, eps, init_eps, delta_min, alpha_max, 1.0e-5, param);
        separatrices.push_back(separatrix());
        std::copy(man2.begin(), man2.end(), std::back_inserter(separatrices.back().manifold));
        separatrices.back().start = fp_index_type(chain_id, 0);
        separatrices.back().end = fp_index_type(chain_id, 0);
        separatrices.back().forward = true;
        separatrices.back().length = prog.length;
        if (!sanity_check(separatrices.back(), liberal_alpha_max, metric)) {
            std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
            separatrices.back().manifold.resize(2);
            nvis::vec2& x0 = separatrices.back().manifold[0];
            nvis::vec2& x1 = separatrices.back().manifold[1];
            x0 = chain[0].pos;
            x1 = x0 + 0.001 * metric.displacement(x0, x1) / metric.distance(x0, x1);
        }
        std::cerr << "length = " << prog.length << ", " << man2.size() << " vertices\n";
        
        std::cerr << "\nstable manifold #0..." << std::endl;
        stop.start(chain[0].pos);
        prog.length = 0;
        ManBVP(man3, prog, map, stop, -period, metric, chain[0], true, eps, init_eps, delta_min, alpha_max, 1.0e-5, param);
        separatrices.push_back(separatrix());
        std::copy(man3.begin(), man3.end(), std::back_inserter(separatrices.back().manifold));
        separatrices.back().start = fp_index_type(chain_id, 0);
        separatrices.back().end = fp_index_type(chain_id, 0);
        separatrices.back().forward = false;
        separatrices.back().length = prog.length;
        if (!sanity_check(separatrices.back(), liberal_alpha_max, metric)) {
        
            broken_manifolds.push_back(separatrices.back().manifold);
            
            std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
            separatrices.back().manifold.resize(2);
            nvis::vec2& x0 = separatrices.back().manifold[0];
            nvis::vec2& x1 = separatrices.back().manifold[1];
            x0 = chain[0].pos;
            x1 = x0 + 0.001 * metric.displacement(x0, x1) / metric.distance(x0, x1);
        }
        std::cerr << "length = " << prog.length << ", " << man3.size() << " vertices\n";
        
        std::cerr << "\nstable manifold #1..." << std::endl;
        stop.start(chain[0].pos);
        prog.length = 0;
        ManBVP(man4, prog, map, stop, -period, metric, chain[0], false, eps, init_eps, delta_min, alpha_max, 1.0e-5, param);
        separatrices.push_back(separatrix());
        std::copy(man4.begin(), man4.end(), std::back_inserter(separatrices.back().manifold));
        separatrices.back().start = fp_index_type(chain_id, 0);
        separatrices.back().end = fp_index_type(chain_id, 0);
        separatrices.back().forward = false;
        separatrices.back().length = prog.length;
        if (!sanity_check(separatrices.back(), liberal_alpha_max, metric)) {
        
            broken_manifolds.push_back(separatrices.back().manifold);
            
            std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
            separatrices.back().manifold.resize(2);
            nvis::vec2& x0 = separatrices.back().manifold[0];
            nvis::vec2& x1 = separatrices.back().manifold[1];
            x0 = chain[0].pos;
            x1 = x0 + 0.001 * metric.displacement(x0, x1) / metric.distance(x0, x1);
        }
        std::cerr << "length = " << prog.length << ", " << man4.size() << " vertices\n";
    } else {
        std::cerr << "compute_separatrix called for period " << period << std::endl;
        
        manifold_stop stop(all_p_chains, metric, _eps);
        
        try {
            for (int i=0 ; i<chain.size() ; ++i) {
                manifold_progress prog;
                std::cerr << "\n\nprocessing saddle " << i+1 << " from " << chain.size() << " at " << chain[i].pos << std::endl;
                
                std::cerr << "\nunstable manifold #0..." << std::endl;
                stop.start(fp_index_type(chain_id, i));
                prog.length = 0;
                ManBVP(man1, prog, map, stop, period, metric, chain[i], true, eps, init_eps, delta_min, alpha_max, 1.0e-5, param);
                separatrices.push_back(separatrix());
                std::copy(man1.begin(), man1.end(), std::back_inserter(separatrices.back().manifold));
                separatrices.back().start = fp_index_type(chain_id, i);
                separatrices.back().end = stop.connection();
                separatrices.back().forward = true;
                separatrices.back().length = prog.length;
                if (!sanity_check(separatrices.back(), liberal_alpha_max, metric)) {
                
                    broken_manifolds.push_back(separatrices.back().manifold);
                    
                    std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
                    separatrices.back().manifold.resize(2);
                    nvis::vec2& x0 = separatrices.back().manifold[0];
                    nvis::vec2& x1 = separatrices.back().manifold[1];
                    x0 = chain[i].pos;
                    x1 = x0 + 0.001 * metric.displacement(x0, x1) / metric.distance(x0, x1);
                }
                std::cerr << "length = " << prog.length << ", " << man1.size() << " vertices\n";
                
                std::cerr << "\nunstable manifold #1..." << std::endl;
                stop.start(fp_index_type(chain_id, i));
                prog.length = 0;
                ManBVP(man2, prog, map, stop, period, metric, chain[i], false, eps, init_eps, delta_min, alpha_max, 1.0e-5, param);
                separatrices.push_back(separatrix());
                std::copy(man2.begin(), man2.end(), std::back_inserter(separatrices.back().manifold));
                separatrices.back().start = fp_index_type(chain_id, i);
                separatrices.back().end = stop.connection();
                separatrices.back().forward = true;
                separatrices.back().length = prog.length;
                if (!sanity_check(separatrices.back(), liberal_alpha_max, metric)) {
                
                    broken_manifolds.push_back(separatrices.back().manifold);
                    
                    std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
                    separatrices.back().manifold.resize(2);
                    nvis::vec2& x0 = separatrices.back().manifold[0];
                    nvis::vec2& x1 = separatrices.back().manifold[1];
                    x0 = chain[i].pos;
                    x1 = x0 + 0.001 * metric.displacement(x0, x1) / metric.distance(x0, x1);
                }
                std::cerr << "length = " << prog.length << ", " << man2.size() << " vertices\n";
                
                std::cerr << "\nstable manifold #0..." << std::endl;
                stop.start(fp_index_type(chain_id, i));
                prog.length = 0;
                ManBVP(man3, prog, map, stop, -period, metric, chain[i], true, eps, init_eps, delta_min, alpha_max, 1.0e-5, param);
                separatrices.push_back(separatrix());
                std::copy(man3.begin(), man3.end(), std::back_inserter(separatrices.back().manifold));
                separatrices.back().start = fp_index_type(chain_id, i);
                separatrices.back().end = stop.connection();
                separatrices.back().forward = false;
                separatrices.back().length = prog.length;
                if (!sanity_check(separatrices.back(), liberal_alpha_max, metric)) {
                
                    broken_manifolds.push_back(separatrices.back().manifold);
                    
                    std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
                    separatrices.back().manifold.resize(2);
                    nvis::vec2& x0 = separatrices.back().manifold[0];
                    nvis::vec2& x1 = separatrices.back().manifold[1];
                    x0 = chain[i].pos;
                    x1 = x0 + 0.001 * metric.displacement(x0, x1) / metric.distance(x0, x1);
                }
                std::cerr << "length = " << prog.length << ", " << man3.size() << " vertices\n";
                
                std::cerr << "\nstable manifold #1..." << std::endl;
                stop.start(fp_index_type(chain_id, i));
                prog.length = 0;
                ManBVP(man4, prog, map, stop, -period, metric, chain[i], false, eps, init_eps, delta_min, alpha_max, 1.0e-5, param);
                separatrices.push_back(separatrix());
                std::copy(man4.begin(), man4.end(), std::back_inserter(separatrices.back().manifold));
                separatrices.back().start = fp_index_type(chain_id, i);
                separatrices.back().end = stop.connection();
                separatrices.back().forward = false;
                separatrices.back().length = prog.length;
                if (!sanity_check(separatrices.back(), liberal_alpha_max, metric)) {
                
                    broken_manifolds.push_back(separatrices.back().manifold);
                    
                    std::cerr << "\n\nSEPARATRIX FAILED SANITY CHECK\n";
                    separatrices.back().manifold.resize(2);
                    nvis::vec2& x0 = separatrices.back().manifold[0];
                    nvis::vec2& x1 = separatrices.back().manifold[1];
                    x0 = chain[i].pos;
                    x1 = x0 + 0.001 * metric.displacement(x0, x1) / metric.distance(x0, x1);
                }
                std::cerr << "length = " << prog.length << ", " << man4.size() << " vertices\n";
            }
        } catch(...) {
            return false;
        }
    }
    
    return true;
}

} // xavier

#endif





