#ifndef __SEPARATRIX_HPP__
#define __SEPARATRIX_HPP__

#include <vector>
#include <list>
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <poincare/metric.hpp>
#include <poincare/fixpoints.hpp>

using namespace spurt;

namespace map_analysis {

struct fp_chain {
    typedef spurt::fixpoint    fp_type;
    
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
    
    std::vector<spurt::fixpoint>   fixed_points;
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

typedef std::vector<nvis::vec2>                 manifold_type;
typedef std::pair<unsigned int, unsigned int>   fp_index_type;

struct separatrix {
    fp_index_type   start, end;
    manifold_type   manifold;
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

template<typename Map>
inline nvis::vec2
BVP_Step(const Map& map,
         const nvis::vec2& start, const nvis::vec2& end,    // seeding segment
         const nvis::vec2& prev, const nvis::vec2& cur,     // end segment
         double delta_min, double alpha_max, double delta_alpha_max,    // quality control
         double& tau)   // where are we along seeding segment
{
    const map_metric& metric = map.metric();
    double delta, alpha;
    double t = 1;
    nvis::vec2 x, y;
    
    nvis::vec2 q = (1 - tau) * start + tau * end;
    
    // std::cerr << "BVP called at " << q << " on segment "
    //           << start << " - " << end << ", prev = " << prev
    //           << ", cur = " << cur << '\n';
    
    while (t > tau) {
        x = (1 - t) * start + t * end;
        // std::cerr << "currently trying at " << x << ". ";
        nvis::vec2 _y = map.p_step(x);
        // std::cerr << "map(" << x << ") = " << _y;
        y = metric.modulo(_y);
        // std::cerr << " (=" << y << ")\n";
        // std::cerr << "\tat " << y << " for t=" << t << '\n';
        alpha = metric.angle(prev, cur, y);
        delta = metric.distance(cur, y);
        if (delta < delta_min ||
                (alpha < alpha_max && delta*alpha < delta_alpha_max)) {
            break;
        } else {
            t = 0.5 * (tau + t);
        }
        // std::cerr << "delta = " << delta << ", alpha = " << alpha << '\n';
    }
    
    tau = t;
    return y;
}

template<typename Map>
inline void ManBVP(manifold_type& manifold, manifold_progress& progress,
                   const Map& map,
                   const spurt::fixpoint& saddle, bool fwd, double eps,
                   double max_length,
                   double delta_min, double alpha_max, double delta_alpha_max)
{
    const map_metric& metric = map.metric();
    
    unsigned int i0, i1;
    double tau;
    double length;
    
    if (progress.length == 0) {
        // std::cerr << "zero length initially\n";
        manifold.clear();
        manifold.push_back(saddle.pos);
        
        // move along eigenvector until p-step aligns within alpha_max
        // with the eigenvector
        nvis::vec2 p = saddle.pos;
        while (true) {
            p += eps * (fwd ? 1 : -1) * saddle.evec[1];
            nvis::vec2 dir = map(p);
            dir /= nvis::norm(dir);
            if (nvis::inner(dir, (fwd ? 1 : -1)*saddle.evec[1]) < cos(alpha_max)) {
                break;
            }
        }
        manifold.push_back(p);
        length = metric.distance(manifold[0], manifold[1]);
        
        // std::cerr << "initial length between " << manifold[0] << " and "
        //           << manifold[1] << " is " << length << '\n';
        
        i0 = 0;
        i1 = 1;
        tau = 0;
    } else {
        length = progress.length;
        i0 = progress.segment;
        i1 = i0 + 1;
        tau = progress.tau;
    }
    
    while (length < max_length) {
    
        const nvis::vec2& q0 = manifold[manifold.size()-2];
        const nvis::vec2& q1 = manifold.back();
        const nvis::vec2& p0 = manifold[i0];
        const nvis::vec2& p1 = manifold[i1];
        
        nvis::vec2 next =
            BVP_Step(map, p0, p1, q0, q1,
                     delta_min, alpha_max, delta_alpha_max, tau);
                     
        if (tau == 1) {
            ++i0;
            ++i1;
            tau = 0;
        }
        
        length += metric.distance(manifold.back(), next);
        manifold.push_back(next);
        
        // std::cerr << "at " << next << ", length = " << length << '\n';
    }
    
    progress.length = length;
    progress.segment = i0;
    progress.tau = tau;
}


template< typename Map >
bool compute_separatrix(std::vector<fp_chain>& all_p_chains,
                        std::vector<separatrix>& separatrices,
                        std::vector< bool >& invalid,
                        unsigned int chain_id,
                        const Map& map)
{
    separatrices.clear();
    const map_metric& metric = map.metric();
    
    if (!all_p_chains.size() || chain_id >= all_p_chains.size()) {
        return true;
    }
    const fp_chain& chain = all_p_chains[chain_id];
    
    if (!chain.saddle()) {
        return true;    // that was easy!
    }
    
    std::cout << "\n\nprocessing saddle chain of period " << chain.period()
              << " starting at " << chain[0].pos << std::endl;
              
    unsigned int period = chain.period();
    
    // period 1 is a special case since we have no way to estimate the length
    // of the separatrix
    if (period == 1) {
        // we are going to progress in steps
        // 10% of the bounding box diameter sounds like a reasonable number
        
        double delta_length = 0.1 * metric.diameter();
        double eps = 1.0e-6 * metric.diameter();
        double delta_min = 1.0e-4 * metric.diameter();
        manifold_progress prog1, prog2;
        manifold_type man1, man2;
        
        while (true) {
            ManBVP(man1, prog1, map, chain[0], true, eps, delta_length, delta_min, 0.3, 1.0e-5);
            ManBVP(man2, prog2, map, chain[0], true, eps, delta_length, delta_min, 0.3, 1.0e-5);
        }
    }
}

} // map_analysis

#endif





