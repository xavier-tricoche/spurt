#ifndef __SEPARATRIX_HPP__
#define __SEPARATRIX_HPP__

#include <vector>
#include <list>
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <poincare/newton.hpp>
#include <poincare/map.hpp>
#include <poincare/metric.hpp>

using namespace spurt;

namespace map_analysis {
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

static unsigned int nb_samples_on_separatrix = 100;

// returns distance between x1 and target if the shortest connection between these two
// points is aligned with the direction pointed to by the shortest connection between x0 and x1
// otherwise returns an arbitrarily large number
double directed_distance(const nvis::vec2& x0, const nvis::vec2& x1, const nvis::vec2& target,
                         const map_metric& metric)
{
    using namespace static_data;
    nvis::vec2 dir = metric.displacement(x0, x1);
    nvis::vec2 y = metric.displacement(x1, target);
    if (nvis::inner(dir, y) < 0) {
        return std::numeric_limits<double>::max();
    } else {
        return nvis::norm(y);
    }
}

void make_aperiodic(std::vector< nvis::vec2 >& out, const std::vector< nvis::vec2 >& in,
                    const map_metric& metric)
{
    using namespace static_data;
    out.resize(in.size());
    out[0] = in[0];
    for (unsigned int i = 1 ; i < in.size() ; ++i) {
        out[i] = out[i-1] + metric.displacement(in[i-1], in[i]);
    }
}

void make_periodic(std::vector< nvis::vec2 >& out, const std::vector< nvis::vec2 >& in,
                   const map_metric& metric)
{
    out.resize(in.size());
    out[0] = in[0];
    for (unsigned int i = 1 ; i < in.size() ; ++i) {
        out[i] = metric.modulo(out[i-1] + (in[i] - in[i-1]));
    }
}

template< typename T >
void advance(const T& map, int period,
             const nvis::vec2& x,
             const nvis::vec2& begin, const nvis::vec2& end,
             double close_distance,
             unsigned maxcount,
             std::vector< nvis::vec2 >& steps,
             std::vector< nvis::vec2 >& tail,
             const map_metric& metric)
{
    // std::cout << "advance: period = " << period << std::endl;
    static nvis::bbox2 check_box(metric.bounds().min() - nvis::vec2(1, 1),
                                 metric.bounds().max() + nvis::vec2(1, 1));
                                 
                                 
    bool must_escape = (fabs(period) == 1);
    double last_d = std::numeric_limits<double>::max();
    unsigned int count;
    steps.clear();
    tail.clear();
    steps.push_back(begin);
    steps.push_back(x);
    nvis::vec2 y = x;
    bool converging = false;
    bool approaching = false;
    int first_approach = -1;
    
    for (count = 0 ; count < maxcount ; ++count) {
        y = metric.modulo(map.map(y, period));
        if (must_escape) {
            double d_bwd = metric.distance(y, begin);
            if (d_bwd > 2.0*close_distance) {
                must_escape = false;
            }
        }
        
        // double d_fwd = metric.distance(y, end);
        double d_fwd = directed_distance(steps.back(), y, end, metric);
        
        if (!converging &&  // we have not converged yet
                !must_escape && // we have left the vicinity of the origin
                !approaching &&
                d_fwd < last_d) { // we are moving toward the target)
                
            std::cout << "setting approaching to TRUE at " << y << " (#" << count << ")\n";
            std::cout << "converging = " << (converging ? "TRUE" : "FALSE")
                      << ", must_escape = " << (must_escape ? "TRUE" : "FALSE")
                      << ", approaching = " << (approaching ? "TRUE" : "FALSE")
                      << ", d_fwd=" << d_fwd << ", last_d=" << last_d << std::endl;
                      
            approaching = true;
            first_approach = steps.size();
        }
        
        if (!converging && approaching && d_fwd < 2.0 * close_distance && // we are close to the target
                d_fwd < last_d) {
            // std::cout << "converging at " << y << " with d = " << d_fwd
            //  << ", and close distance = " << close_distance << std::endl;
            converging = true;
            steps.push_back(y); // last position in steps
            if (!check_box.inside(y)) {
                std::cout << "position " << y << " reached after " << count << " steps from " << x
                          << " is out of bounds!" << std::endl;
            }
        }
        
        if (converging && d_fwd > last_d) {
            return;
        }
        
        if (!converging) {
            steps.push_back(y);
            if (!check_box.inside(y)) {
                std::cout << "position " << y << " reached after " << count << " steps from " << x
                          << " is out of bounds!" << std::endl;
            }
        } else {
            tail.push_back(y);
        }
        
        last_d = d_fwd;
    }
    
    std::cout << "unable to reach prescribed vicinity along separatrix in " << maxcount
              << " steps" << std::endl;
              
    // double dmin = metric.distance(steps[first_approach], end);
    double dlast = directed_distance(steps[first_approach-1], steps[first_approach], end,
                                     metric);
    int mini = first_approach;
    for (int i = first_approach + 1 ; i < steps.size() ; ++i) {
        // double d = metric.distance(steps[i], end);
        double d = directed_distance(steps[i-1], steps[i], end, metric);
        if (d < dlast) {
            dlast = d;
            mini = i;
        } else {
            break;
        }
    }
    
    std::vector< nvis::vec2 > fixed_steps;
    int nb_in_tail = std::min((int)mini - first_approach, 5);
    // we will consider 5 points (if available) between the first converging point and the closest point
    // to determine the best alignment between the separatrix and the target eigenvector
    int nb_pts = std::min(mini + 2, (int)steps.size());
    if (nb_pts <= 0 || nb_in_tail <= 0) {
        std::cout << "advance failed between " << begin << " and " << end << std::endl
                  << "steps taken were\n";
        for (unsigned int i = 0 ; i < steps.size() ; ++i) {
            std::cout << steps[i] << " - ";
        }
        std::cout << std::endl;
        tail.clear();
        return;
    }
    // we include 2 points (if available) past the shortest distance for subsequent optimization
    for (int i = 0 ; i < nb_pts ; ++i) {
        if (i < mini - nb_in_tail) {
            fixed_steps.push_back(steps[i]);
        } else {
            tail.push_back(steps[i]);
        }
    }
    fixed_steps.push_back(tail.front());
    std::swap(steps, fixed_steps);
}


template< typename Map >
void refine_separatrix(std::vector< nvis::vec2 >& refined,
                       const spurt::fixpoint& fp0, const spurt::fixpoint& fp1,
                       const std::vector< nvis::vec2 >& sep,
                       const Map& map,
                       const map_metric& metric)
{
    double l = 0;
    for (unsigned int i = 1 ; i < sep.size(); ++i) {
        l += nvis::norm(metric.displacement(sep[i-1], sep[i]));
    }
    
    if (sep.size() < 4 || l < 1.0e-3) {
        refined.clear();
        std::copy(sep.begin(), sep.end(), std::back_inserter(refined));
        // std::cout << "early exit(" << sep.size() << "pts, " << l << " length)" << std::endl;
        return;
    }
    
    // convert input curve to aperiodic representation
    std::vector< nvis::vec2 > aperiodic;
    make_aperiodic(aperiodic, sep, metric);
    
    unsigned int nb_rounds = 5;
    std::vector < nvis::bezier_spline< nvis::vec2, 3 > > bs(2*nb_rounds + 1);
    bs[0] = make_catmull_rom_spline_arclen(aperiodic);
    
    double dt = 1. / (double)nb_samples_on_separatrix;
    
    double ref_dist = 0.1 * l;
    
    for (unsigned int i = 0 ; i < 2*nb_rounds ; ++i) {
        int p;
        nvis::vec2 x, target, origin;
        std::vector< nvis::vec2 > tmp, new_sep, tail;
        double t = (double)(1 + i / 2) / (double)(nb_rounds + 1);
        if (i % 2) {
            x = (1. - t) * sep.back() + t * sep[sep.size()-2];
            p = -fp1.K;
            origin = sep.back();
            target = sep[0];
            new_sep.push_back(sep.back());
        } else {
            x = (1. - t) * sep[0] + t * sep[1];
            p = fp0.K;
            origin = sep[0];
            target = sep.back();
            new_sep.push_back(sep[0]);
        }
        
        advance(map, p, x, origin, target, ref_dist, 100, new_sep, tail, metric);
        new_sep.push_back(target);
        tmp.clear();
        if (i % 2) {
            std::copy(new_sep.rbegin(), new_sep.rend(), std::back_inserter(tmp));
        } else {
            std::copy(new_sep.begin(), new_sep.end(), std::back_inserter(tmp));
        }
        make_aperiodic(aperiodic, tmp, metric);
        
        bs[i+1] = make_catmull_rom_spline_arclen(aperiodic);
    }
    
    std::vector< nvis::vec2 > blend(nb_samples_on_separatrix + 1);
    blend[0] = sep[0];
    blend.back() = sep.back();
    unsigned int count = 1;
    for (double t = dt ; t < 1 ; t += dt, ++count) {
        double wf = exp(-7.0 * t * t);
        double wb = exp(-7.0 * (1 - t) * (1 - t));
        double wc = exp(-7.0 * 0.25);
        double w = 0;
        nvis::vec2 y = wc * bs[0](t);
        w += wc;
        for (unsigned int i = 1 ; i < bs.size() ; ++i) {
            unsigned int j = i - 1;
            if (j % 2) {
                y += wb * bs[i](t);
                w += wb;
            } else {
                y += wf * bs[i](t);
                w += wf;
            }
            
            nvis::vec2 z = bs[i](t);
            if (std::isnan(z[0]) || std::isnan(z[1])) {
                std::cout << "something went south with " << z << std::endl;
                std::cout << "corresponding knots: " << std::endl;
                const std::vector< double >& knots = bs[i].knots();
                for (unsigned int k = 0 ; k < knots.size() ; ++k) {
                    std::cout << "knot " << k << ": " << knots[k] << std::endl;
                }
                std::cout << "corresponding control points:" << std::endl;
                const std::vector< nvis::vec2 >& cps = bs[i].ctrlp();
                for (unsigned int k = 0 ; k < cps.size() ; ++k) {
                    std::cout << "cp " << k << ": " << cps[k] << std::endl;
                }
            }
        }
        
        
        blend[count] = 1. / w * y;
    }
    
    refined.clear();
    make_periodic(refined, blend, metric);
    std::cerr << "separatrix refinement completed\n";
    std::cerr << "refiner version contains " << refined.size() << " points\n";
}

double oriented_distance(unsigned int& minid, unsigned int id,
                         const nvis::vec2& dir, const std::vector< nvis::vec2 >& points,
                         const map_metric& metric)
{
    unsigned int n = points.size();
    minid = id;
    double mind = std::numeric_limits<double>::max();
    for (unsigned int i = 1 ; i < n ; ++i) {
        unsigned int j = (id + i) % n;
        nvis::vec2 step = metric.displacement(points[id], points[j]);
        if (nvis::inner(step, dir) < 0.) {
            continue;
        }
        double d = nvis::norm(step);
        if (d < mind) {
            mind = d;
            minid = j;
        }
    }
    
    std::cout << "closest point from " << points[id] << " in direction " << dir << " is "
              << points[minid] << " at distance " << mind << std::endl;
              
    return mind;
}

void connect_points(std::vector< nvis::vec2 >& curve, const std::vector< nvis::vec2 >& points,
                    const map_metric& metric)
{
    // crappy method with cubic complexity
    unsigned int n = points.size();
    std::vector< std::vector< double > > distances(n);
    for (unsigned i = 0 ; i < n ; ++i) {
        distances[i].resize(n);
        std::fill(distances[i].begin(), distances[i].end(), std::numeric_limits<double>::max());
    }
    for (unsigned int i = 0 ; i < n ; ++i) {
        nvis::vec2 x = points[i];
        for (unsigned int j = i + 1; j < n ; ++j) {
            nvis::vec2 y = points[j];
            double dist = metric.distance(x, y);
            distances[i][j] = distances[j][i] = dist;
        }
    }
    
    std::vector< std::pair< unsigned int, unsigned int > > neighbors(n);
    
    std::vector< unsigned int > sorted;
    for (unsigned int i = 0 ; i < n ; ++i) {
        spurt::sort(distances[i], sorted);
        // points[sorted[0]] is the closest neighbor
        // find closest point in opposite direction
        neighbors[i].first = sorted[0];
        nvis::vec2 dir = metric.displacement(points[i], points[sorted[0]]);
        for (unsigned int j = 1 ; j < n - 1 ; ++j) {
            nvis::vec2 a = metric.displacement(points[i], points[sorted[j]]);
            if (inner(a, dir) < 0) {
                neighbors[i].second = sorted[j];
                break;
            }
        }
    }
    
    curve.clear();
    std::vector< bool > used(n, false);
    unsigned int nused = 0;
    curve.push_back(points[0]);
    used[0] = true;
    ++nused;
    unsigned int cur = 0;
    while (nused < n) {
        unsigned int n0 = neighbors[cur].first;
        unsigned int n1 = neighbors[cur].second;
        unsigned int nnext = (used[n0] ? n1 : n0);
        curve.push_back(points[nnext]);
        used[nnext] = true;
        cur = nnext;
        ++nused;
    }
    
    curve.push_back(curve.front());
}

template< typename PMap, typename Map >
bool trace_rational_manifold(const PMap& pmap, const Map& map, unsigned int period,
                             const nvis::vec2& x, double h, double step_size, double eps,
                             std::vector< nvis::vec2 >& curve,
                             const map_metric& metric)
{
    curve.clear();
    float color[] = { drand48(), drand48(), drand48() };
    
    typedef std::vector< nvis::vec2 >   chain;
    std::vector< chain > chains(period);
    
    std::vector< nvis::vec2 > point_cloud, hits;
    std::vector< std::pair< double, nvis::vec2 > > samples;
    map.map(x, hits, period - 1);
    point_cloud.push_back(x);
    chains[0].push_back(x);
    for (unsigned int i = 0 ; i < hits.size() ; ++i) {
        point_cloud.push_back(hits[i]);
        chains[i+1].push_back(hits[i]);
    }
    
    nvis::vec2 cur, next, minpos, dir, newdir, f;
    double minf = std::numeric_limits<double>::max();
    
    // initialize tangent direction
    for (unsigned int i = 0 ; i < 4 ; ++i) {
        double theta0 = i * M_PI / 4.;
        double theta1 = (i + 1) * M_PI / 4.;
        bool ok = track_map_rotation(pmap, x, theta0, theta1, h, samples, 50, 0.24 * M_PI, next);
        // if (ok) continue;
        next = metric.modulo(next);
        f = pmap(next);
        double nf = nvis::norm(f);
        if (nf < minf) {
            minf = nf;
            minpos = next;
        }
        if (minf < eps) {
            break;
        }
    }
    
    if (minf > eps) {
        std::cout << "min norm around " << x << " = " << minf << " exceeds threshold = " << eps << std::endl;
        std::cout << "no rational surface here\n";
        return false;
    }
    
    next = minpos;
    dir = metric.displacement(x, next);
    dir /= nvis::norm(dir);
    cur = next;
    map.map(cur, hits, period - 1);
    int id = point_cloud.size();
    point_cloud.push_back(cur);
    chains[0].push_back(cur);
    std::copy(hits.begin(), hits.end(), std::back_inserter(point_cloud));
    for (unsigned int i = 0 ; i < period - 1 ; ++i) {
        chains[i+1].push_back(hits[i]);
    }
    
    unsigned int minid;
    double remaining = oriented_distance(minid, id, dir, point_cloud, metric);
    
    unsigned int max_nb_steps = 10.*nvis::norm(metric.bounds().max() -
                                metric.bounds().min()) / (double)period / step_size;
    std::cout << "max_nb_steps = " << max_nb_steps << std::endl;
    
    const double dtheta = M_PI / 20.; // sampling rate of angular domain
    const unsigned int nb_intervals = (unsigned int)floor(M_PI / 2. / dtheta); // number of angular intervals lying in correct half plane
    
    for (unsigned int k = 0 ; k < max_nb_steps && remaining > step_size ; ++k) {
        std::cout << "trace_manifold: iteration " << k << ": " << point_cloud.size()
                  << " points on manifold currently" << std::endl;
                  
        bool found = false;
        double dir_theta = (dir[1] < 0 ? -1 : 1) * acos(dir[0]); // angle of pointing direction
        std::cout << "pointing direction has angle = " << dir_theta*180 / M_PI << std::endl;
        std::cout << "dir vector is " << dir << std::endl;
        
        double init_theta = dir_theta + 0.5 * dtheta;
        std::cout << "dtheta = " << dtheta*180 / M_PI << ", init_theta = " << init_theta*180 / M_PI << "\n";
        
        for (unsigned int i = 0 ; i < nb_intervals ; ++i) {
            unsigned int j = i / 2;
            double theta0, theta1;
            theta0 = init_theta + (i % 2 ? j * dtheta : 2 * M_PI - (j + 1) * dtheta);
            theta1 = theta0 + dtheta;
            std::cout << "considering interval #" << i << " between " << theta0*180 / M_PI
                      << "degrees and " << theta1*180 / M_PI << "degrees\n";
            std::cout << "tracking rotation between vector at " << cur + step_size* nvis::vec2(cos(theta0), sin(theta0))
                      << " = " << pmap(cur + step_size*nvis::vec2(cos(theta0), sin(theta0)))
                      << " and vector at " << cur + step_size* nvis::vec2(cos(theta1), sin(theta1))
                      << " = " << pmap(cur + step_size*nvis::vec2(cos(theta1), sin(theta1))) << "\n";
                      
            bool ok = track_map_rotation(pmap, cur, theta0, theta1, step_size, samples,
                                         20, 0.25 * M_PI, next);
                                         
            if (true) {
                next = metric.modulo(next);
                f = pmap(next);
                double nf = nvis::norm(f);
                if (nf < eps) {
                    found = true;
                    break;
                }
            }
        }
        
        if (!found) {
            std::cout << "unable to proceed from " << cur << "\n";
            return false;
        }
        
        id = point_cloud.size();
        point_cloud.push_back(next);
        chains[0].push_back(next);
        map.map(next, hits, period - 1);
        std::copy(hits.begin(), hits.end(), std::back_inserter(point_cloud));
        for (unsigned int i = 0 ; i < period - 1 ; ++i) {
            chains[i+1].push_back(hits[i]);
        }
        
        dir = metric.displacement(cur, next);
        dir /= nvis::norm(dir);
        cur = next;
        
        remaining = oriented_distance(minid, id, dir, point_cloud, metric);
    }
    
    if (remaining < step_size) {
        // connect_points(curve, point_cloud);
        std::copy(chains[0].begin(), chains[0].end(), std::back_inserter(curve));
        id = minid;
        for (unsigned int i = 0 ; i < period - 1 ; ++i, id += minid) {
            id = id % period;
            std::copy(chains[id].begin(), chains[id].end(), std::back_inserter(curve));
        }
        
        return true;
    }
    
    std::cout << "failed to converge.\n";
    return false;
    
}

template< typename Map >
void refine(std::vector< nvis::vec2 >& finer,
            const spurt::fixpoint& fp0, const spurt::fixpoint& fp1,
            const std::vector< nvis::vec2 >& sep,
            const Map& map,
            const map_metric& metric)
{
    std::cout << "SPLINE_BASED is defined" << std::endl;
    refine_separatrix(finer, fp0, fp1, sep, map, metric);
    return;
}

template< typename Map >
unsigned int find_saddle_neighbors(Neighbor& nfwd, Neighbor& nbwd,
                                   const std::vector< std::vector< spurt::fixpoint > >& all_p_chains,
                                   unsigned int chain_id,
                                   const Map& map,
                                   const map_metric& metric)
{
    unsigned int start_id = 0;
    
    try {
        // chain period. assuming no period doubling
        unsigned int p = all_p_chains[0].size();
        const std::vector< spurt::fixpoint >& saddle_chain = all_p_chains[chain_id];
        nfwd = Neighbor(chain_id, -1);
        nbwd = Neighbor(chain_id, -1);
        if (p == 1) {
            nfwd = Neighbor(chain_id, 0);
            nbwd = Neighbor(chain_id, 0);
            return start_id;
        }
        
        for (; start_id < p ; ++start_id) {
        
            // reset neighbors to invalid values
            nfwd = Neighbor(chain_id, -1);
            nbwd = Neighbor(chain_id, -1);
            
            const spurt::fixpoint& fp = saddle_chain[start_id];
            int d_id;
            
            // determine order:
            // heuristic: saddle #0 is connected to saddle that is approached first at a distance
            // smaller than 10% of the min distance between saddle #0 and any other saddle on any other
            // saddle chain of the same period
            double ref_distance = 0.001; // pixel size
            Connection closest;
            closest.d = std::numeric_limits<double>::max();
            for (unsigned int j = 0 ; j < all_p_chains.size() ; ++j) {
                const std::vector< spurt::fixpoint >& chain = all_p_chains[j];
                if (chain[0].saddle) {
                    for (unsigned int i = 0 ; i < p ; ++i) {
                        if (j == chain_id && i == start_id) {
                            continue;
                        }
                        double d = metric.distance(fp.pos, chain[i].pos);
                        if (d < closest.d) {
                            closest.neigh = Neighbor(j, i);
                            closest.d = d;
                        }
                    }
                }
            }
            ref_distance = 0.1 * closest.d;
            std::cout << "ref distance for connection is " << ref_distance << std::endl;
            
            for (unsigned int dir = 0 ; dir < 2 ; ++dir) {
                std::vector< Connection > connections;
                
                // walk along the manifold
                nvis::vec2 x = fp.pos + 0.0002 * (dir ? 1.0 : -1.0) * fp.evec[1]; // unstable manifold has index 1 per construction
                
                unsigned int c;
                std::vector< nvis::vec2 > steps;
                try {
                    for (c = 0 ; c < 200 ; ++c) {
                        steps.push_back(x);
                        x = map.map(x, p);
                        
                        // compute min distance to points
                        closest.d = std::numeric_limits<double>::max();
                        for (unsigned int j = 0 ; j < all_p_chains.size() ; ++j) {
                            const std::vector< spurt::fixpoint >& chain = all_p_chains[j];
                            if (chain[0].saddle) {
                                for (unsigned int i = 0 ; i < p ; ++i) {
                                    if (j == chain_id && i == start_id) {
                                        continue;
                                    }
                                    double d = metric.distance(x, chain[i].pos);
                                    if (d < closest.d) {
                                        closest.neigh = Neighbor(j, i);
                                        closest.d = d;
                                    }
                                }
                            }
                        }
                        double mindist = closest.d;
                        if (mindist < ref_distance) {
                            if (!dir) {
                                nfwd = closest.neigh;
                            } else {
                                nbwd = closest.neigh;
                            }
                            break;
                        }
                    }
                } catch (...) {
                    std::cout << "caught exception in find_saddle_neighbor while walking along separatrix" << std::endl;
                    std::cout << "the last position we reached was " << steps.back() << std::endl;
                    return 0;
                }
                
                if (c == 200) {
                    std::cout << "failed to find next saddle from " << fp.pos << std::endl;
                }
            }
            
            if (nfwd.p_id >= 0 && nbwd.p_id >= 0) {
                // we were able to connect current saddle to his two neighbors
                break;
            }
        }
    } catch (...) {
        std::cout << "caught exception in find_saddle_neighbor" << std::endl;
        return 0;
    }
    
    return start_id;
}


template< typename Map >
void walk(std::vector< nvis::vec2 >& steps,
          const spurt::fixpoint& fp0, const spurt::fixpoint& fp1,
          const Map& map, double dir,
          const map_metric& metric, bool fwd = true)
{
    try {
        std::cout << "walking between " << fp0.pos << " and " << fp1.pos << std::endl;
        
        nvis::vec2 start_pt = fp0.pos;
        nvis::vec2 target_pt = fp1.pos;
        bool must_escape = false;
        bool converging = false;
        std::vector< nvis::vec2 > remain;
        steps.clear();
        
        int p = fp0.K;
        if (p == 1) {
            must_escape = true; // means that we have to exit saddle's neighborhood
            // before we can consider convergence
            std::cerr << "we must escape\n";
        }
        if (!fwd) {
            p *= -1;    // we will be stepping backward
        }
        
        // walk along the manifold
        steps.push_back(start_pt);
        double ref_distance = 0.1 * metric.distance(start_pt, target_pt);
        if (std::fabs(p) == 1) {
            ref_distance = 0.05 * metric.diameter();
        }
        
        std::cerr << "ref distance = " << ref_distance << '\n';
        
        nvis::vec2 x = start_pt + dir * 0.1 * ref_distance * fp0.evec[(fwd ? 1 : 0)];
        // unstable manifold has index 1 per construction
        steps.push_back(x);
        double last_d = std::numeric_limits<double>::max();
        unsigned int count;
        for (count = 0 ; count < 200 ; ++count) {
            x = metric.modulo(map.map(x, p));
            if (must_escape) {
                double d_bwd = metric.distance(x, fp0.pos);
                if (d_bwd > 2.0*ref_distance) {
                    must_escape = false;
                }
            }
            
            std::cerr << "at " << x << " after " << count + 1 << " steps\n";
            
            double d_fwd = metric.distance(x, target_pt);
            if (!converging && (!must_escape && d_fwd < 2.0 * ref_distance)) {
                std::cout << "we are converging at " << x << std::endl;
                converging = true;
                steps.push_back(x);
            }
            if (converging && d_fwd > last_d) {
                std::cout << "distance stopped to decrease: " << d_fwd << " > " << last_d << std::endl;
                break;
            }
            
            if (!converging) {
                steps.push_back(x);
            } else {
                std::cout << "added a remain position" << std::endl;
                remain.push_back(x);
            }
            
            last_d = d_fwd;
        }
        
        if ((count == 200 && last_d > ref_distance) || !remain.size()) {
            std::cout << "we are unable to pursue (max steps = " << 200 << "): steps so far" << std::endl;
            for (unsigned int i = 0 ; i < steps.size() ; ++i) {
                std::cout << i << ": " << steps[i] << std::endl;
            }
            return;
        }
        
        // close the connection
        nvis::vec2 target_evec = fp1.evec[(fwd ? 0 : 1)];
        nvis::vec2 last_step = metric.displacement(remain[0], target_pt);
        double dot = nvis::inner(target_evec, last_step);
        if (dot < 0) {
            target_evec *= -1;
        }
        
        unsigned int last_id = 0;
        if (remain.size() > 1) {
            // now determine which step best aligns with target eigenvector
            std::vector< double > dots(remain.size() - 1, 0);
            for (unsigned int i = 0 ; i < remain.size() - 1 ; ++i) {
                nvis::vec2 step = metric.displacement(remain[i], target_pt);
                dots[i] = nvis::inner(step, target_evec) / nvis::norm(step);
            }
            std::vector< unsigned int > sorted_dots;
            spurt::sort(dots, sorted_dots);
            
            last_id = sorted_dots.back();
        }
        for (unsigned int i = 1 ; i <= last_id ; ++i) {
            steps.push_back(remain[i]);
        }
        steps.push_back(target_pt);
    } catch (...) {
        std::cout << "exception caught in walk()" << std::endl;
    }
}

template< typename Iterator >
bool insert(const nvis::vec2& x, Iterator& first, const Iterator& last,
            const map_metric& metric)
{
    if (first == last) {
        return false;
    }
    
    Iterator next;
    
    // locate x along the list
    next = first;
    ++next;
    
    bool found = false;
    
    while (next != last) {
        nvis::vec2 p0x = metric.displacement(*first, x);
        nvis::vec2 p1x = metric.displacement(*next, x);
        double dot = nvis::inner(p0x, p1x);
        if (dot < 0) {
            double l0 = nvis::norm(p0x);
            double l1 = nvis::norm(p1x);
            double lref = metric.distance(*first, *next);
            if (-dot / l0 / l1 < 0.9 || std::max(l0, l1) > lref) {
                return false;
            }
            
            found = true;
            
            std::cout << "insert: " << x << " is located between " << *first
                      << " and " << *next
                      << ", p0x = " << p0x << ", p1x = " << p1x << std::endl;
                      
            break;
        }
        // move to next segment
        ++first;
        ++next;
    }
    if (!found) {
        // std::cout << "no interval found for " << x << std::endl;
        return false;
    }
    
    return true;
}

void merge(std::vector< nvis::vec2 >& merged,
           const std::vector< nvis::vec2 >& fwd,
           const std::vector< nvis::vec2 >& bwd,
           const map_metric& metric)
{
    if (fwd.size() < 4 || bwd.size() < 4) {
        std::cout << "WARNING: " << fwd.size() << " vs. " << bwd.size() << std::endl;
    }
    
    typedef std::list< nvis::vec2 > List;
    List::iterator prev, next, last;
    
    List __fwd(fwd.begin(), fwd.end());
    List __bwd(bwd.begin(), bwd.end());
    nvis::vec2 start(__fwd.front());
    nvis::vec2 stop(__fwd.back());
    
    double ref_distance = 0.05 * metric.distance(fwd.front(), fwd.back());
    if (ref_distance == 0) {
        ref_distance = 0.001 * metric.diameter();
    }
    
    // remove unreliable ends from both lists
    for (List::reverse_iterator it = __fwd.rbegin() ; it != __fwd.rend() ;) {
        if (metric.distance(stop, *it) > ref_distance) {
            break;
        }
        ++it;
        if (__fwd.size() > 1) {
            __fwd.pop_back();
        } else {
            merged.clear();
            std::copy(fwd.begin(), fwd.end(), std::back_inserter(merged));
            return;
        }
        
    }
    for (List::reverse_iterator it = __bwd.rbegin() ; it != __bwd.rend() ;) {
        if (metric.distance(start, *it) > ref_distance) {
            break;
        }
        ++it;
        if (__bwd.size() > 1) {
            __bwd.pop_back();
        } else {
            merged.clear();
            std::copy(fwd.begin(), fwd.end(), std::back_inserter(merged));
            return;
        }
    }
    // append start points from reverse list to forward list
    {
        List::iterator it, previous;
        previous = __fwd.end();
        for (it = __bwd.begin() ; it != __bwd.end() ; ++it) {
            double d = metric.distance(*it, stop);
            if (d > ref_distance) {
                break;
            }
            __fwd.insert(previous, *it);
            --previous;
        }
    }
    
    // intermingle points
    List::reverse_iterator it;
    prev = __fwd.begin();
    last = __fwd.end();
    for (it = __bwd.rbegin() ; it != __bwd.rend() ; ++it) {
        const nvis::vec2& x = *it;
        if (insert(x, prev, last, metric)) {
            next = prev;
            ++next;
            __fwd.insert(next, x);
            ++prev;
        } else {
            break;
        }
    }
    merged.clear();
    std::copy(__fwd.begin(), __fwd.end(), std::back_inserter(merged));
}

bool valid_curve(const std::vector< nvis::vec2 >& curve,
                 const map_metric& metric, bool closed_loop = false)
{
    std::list< nvis::vec2 > cur_list(curve.begin(), curve.end());
    std::list< nvis::vec2 >::iterator prev, cur, next, last;
    // sanity check
    prev = cur_list.begin();
    cur = prev;
    ++cur;
    next = cur;
    ++next;
    bool draw = true;
    bool wrong = false;
    
    double avg_dot = 0;
    double avg_length = 0;
    double max_length = 0;
    avg_length += metric.distance(*prev, *cur);
    for (; next != cur_list.end() ; ++next, ++cur, ++prev) {
        nvis::vec2 a = metric.displacement(*prev, *cur);
        nvis::vec2 b = metric.displacement(*cur, *next);
        double la = nvis::norm(a);
        double lb = nvis::norm(b);
        double maxl = std::max(la, lb);
        if (maxl > max_length) {
            max_length = maxl;
        }
        double dot = nvis::inner(a, b) / (la * lb);
        if (dot < 0) {
            std::cout << "giving up at " << cur_list.front()
                      << " because dot product got negative at " << *cur << std::endl;
                      
            return false;
        }
        avg_dot += dot;
        avg_length += metric.distance(*cur, *next);
    }
    avg_dot /= (double)(cur_list.size() - 2);
    if (avg_dot < 0.9) {
        std::cout << "giving up at " << cur_list.front() << " because average dot product is too low (" << avg_dot << ")" << std::endl;
        return false;
    }
    avg_length *= 1 / (double)(cur_list.size() - 1);
    double total_distance = metric.distance(cur_list.front(), cur_list.back());
    if (closed_loop) {
        total_distance = 0.2 * metric.diameter();
    }
    
    std::cerr << "total distance = " << total_distance << '\n';
    
    if (avg_length > 0.5*total_distance) {
        std::cout << "giving up at " << cur_list.front() << " because average length is too large (" << avg_length << ")" << std::endl;
        return false;
    }
    if (max_length > 0.75*total_distance) {
        std::cout << "giving up at " << cur_list.front() << " because max length is too large (" << max_length << ")" << std::endl;
        return false;
    }
    
    return true;
}

template<typename Map>
bool transport(std::vector<separatrix_type>& separatrices,
               const Map& map, unsigned int p, const map_metric& metric,
               const std::list<spurt::fixpoint>& chain)
{
    try {
        typedef std::list< nvis::vec2 > list_type;
        
        unsigned int nbsep = separatrices.size();
        for (unsigned int n = 0 ; n < nbsep ; ++n) {
        
            list_type cur_list(separatrices[n].begin(), separatrices[n].end());
            list_type next_list, tmp_list;
            list_type::iterator prev, cur, next;
            
            for (unsigned int i = 0 ; i < p - 1 ; ++i) {
                // advance separatrix
                std::cout << "drawing separatrix after " << i + 1 << " forward iterations of the map\n";
                next_list.clear();
                
                // determine the target
                nvis::vec2 target = metric.modulo(map.map(cur_list.back(), 1));
                nvis::vec2 best;
                double d = std::numeric_limits<double>::max();
                for (std::list<spurt::fixpoint>::const_iterator it = chain.begin() ; it != chain.end() ; it++) {
                    double __d =
                }
                
                for (prev = cur_list.begin(); prev != cur_list.end() ; ++prev) {
                    next_list.push_back(metric.modulo(map.map(*prev, 1)));
                }
                
                // // check the spatial connectivity of the separatrix
                // for (list_type::iterator it = next_list.begin() ; it != next_list.end() ; ++it) {
                //  next = it;
                //  ++next;
                //  if (next == next_list.end()) break;
                //  nvis::vec2 aux = *it + metric.displacement(*it, *next);
                //  if (nvis::norm(aux-*next)>1.0e-6*metric.diameter()) {
                //      next_list.insert(next, aux);
                //      next_list.insert(next, *next - metric.displacement(*it, *next));
                //      it = next;
                //      --it;
                //  }
                // }
                
                separatrices.push_back(separatrix_type(next_list.begin(), next_list.end()));
                cur_list.swap(next_list);
            }
        }
    } catch (...) {
        std::cout << "exception caught in draw separatrix" << std::endl;
        return false;
    }
    
    return true;
}

template< typename Map >
bool process_saddle_chain(std::vector< std::vector< spurt::fixpoint > >& all_p_chains,
                          std::vector< std::vector< nvis::vec2 > >& separatrices,
                          std::vector< bool >& invalid,
                          unsigned int chain_id,
                          const Map& map,
                          const map_metric& metric)
{
    separatrices.clear();
    
    if (chain_id >= all_p_chains.size() || !all_p_chains[chain_id].size() ||
            !all_p_chains[chain_id][0].saddle) {
        return true;    // that was easy!
    }
    
    std::cout << "\n\nprocessing saddle chain of period " << all_p_chains[chain_id].size()
              << " starting at " << all_p_chains[chain_id][0].pos << std::endl;
    unsigned int p;
    bool everything_ok = true;
    try {
        // chain period
        std::vector< std::vector< spurt::fixpoint > > tmp_vec(all_p_chains.begin(), all_p_chains.end());
        std::vector< spurt::fixpoint >& saddle_chain = tmp_vec[chain_id];
        p = saddle_chain.size();
        
        const spurt::fixpoint& fp = saddle_chain[0];
        int d_id;
        std::vector< nvis::vec2 > steps;
        std::vector< nvis::vec2 > remain;
        
        bool must_escape;
        bool converging;
        bool one_chain = (p == 1);
        
        std::cerr << "Saddle chain contains " << p << " saddles\n";
        
        std::cout << "processing saddle chain #" << chain_id << " of period " << saddle_chain.size() << '\n';
        std::cout << "this saddle chain contains: " << std::endl;
        {
            nvis::vec2 x = saddle_chain[0].pos;
            for (unsigned int i = 0 ; i < p ; ++i) {
                std::cout << "saddle #" << i << ": " << saddle_chain[i].pos << " (" << x << ")" << std::endl;
                x = metric.modulo(map.map(x, 1));
            }
            std::cout << "saddle #" << p << ": " << x
                      << " (" << metric.distance(x, saddle_chain[0].pos) << ")" << std::endl;
        }
        
        Neighbor nexts[2];
        unsigned int ref_id = find_saddle_neighbors(nexts[0], nexts[1], tmp_vec, chain_id, map, metric);
        
        if (nexts[0].p_id == -1 || nexts[1].p_id == -1) {
            return false;
        }
        
        std::cout << "both saddle neighbors of saddle " << ref_id << " have been identified: " << nexts[0].p_id
                  << " and " << nexts[1].p_id << std::endl;
                  
        bool pdoubling = false;
        if (nexts[0].chain_id != chain_id) {
            // we are clearly assuming that period trippling or quadrupling or ... is not happening
            std::cout << "period doubling detected: saddle chain #" << chain_id
                      << " connected to saddle chain #"
                      << nexts[0].chain_id << std::endl;
            pdoubling = true;
            
            std::vector< spurt::fixpoint >& aux_vec = tmp_vec[nexts[0].chain_id];
            for (unsigned int i = 0 ; i < aux_vec.size() ; ++i) {
                saddle_chain.push_back(aux_vec[i]);
            }
            
            invalid[nexts[0].chain_id] = true; // this chain has been absorbed by current one
            
            nexts[0].chain_id = chain_id;
            nexts[1].chain_id = chain_id;
            nexts[0].p_id += p;
            nexts[1].p_id += p;
            p *= 2;
        }
        std::vector< std::pair< unsigned int, unsigned int > > sepids;
        std::vector< double > dirs;
        for (unsigned int i = 0 ; i < 2 ; ++i) {
            if (nexts[i].p_id < 0) {
                continue;
            }
            sepids.push_back(std::pair<unsigned int, unsigned int >(ref_id, nexts[i].p_id));
            dirs.push_back(i ? 1.0 : -1.0);
            if (pdoubling) {
                sepids.push_back(std::pair<unsigned int, unsigned int >(nexts[i].p_id, ref_id));
                dirs.push_back(0); // unknown at this point
            }
        }
        // store_primitives = true;
        for (unsigned int n = 0 ; everything_ok && n < sepids.size() ; ++n) {
        
            std::cout << "saddle chain iteration #" << n + 1 << " from " << sepids.size() << std::endl;
            
            unsigned int i = sepids[n].first;
            unsigned int j = sepids[n].second;
            std::vector< nvis::vec2 > fwd_steps, bwd_steps;
            
            std::cout << "connecting saddle #" << i << " and saddle #"
                      << j << std::endl;
                      
            std::cout << "walking forward along connection..." << std::endl;
            walk(fwd_steps, saddle_chain[i], saddle_chain[j], map, dirs[n], metric);
            
            std::cout << fwd_steps.size() << " forward steps between "
                      << fwd_steps[0] << " and " << fwd_steps.back() << std::endl;
                      
            // debugging
            for (unsigned int i = 0 ; i < fwd_steps.size() - 1 ; ++i) {
                double l = metric.distance(fwd_steps[i], fwd_steps[i+1]);
                if (l > 20) {
                    std::cout << "crazy segment between " << fwd_steps[i]
                              << " and " << fwd_steps[i+1] << std::endl;
                }
            }
            
            nvis::vec2 ev = saddle_chain[j].evec[0];
            nvis::vec2 step = metric.displacement(fwd_steps[fwd_steps.size()-2], fwd_steps.back());
            double dot = nvis::inner(step, ev);
            
            std::cout << "dot product between last step along separatrix "
                      << step << " and stable eigenvector " << ev << " is " << dot << std::endl;
                      
            std::cout << "walking backward along connection..." << std::endl;
            walk(bwd_steps, saddle_chain[j], saddle_chain[i], map, (dot > 0 ? -1 : + 1), metric, false);
            
            std::cout << bwd_steps.size() << " backward steps between "
                      << bwd_steps[0] << " and " << bwd_steps.back() << std::endl;
                      
            // if (pdoubling && n == 1) {
            //          std::cout << "connecting saddle #" << j << " and saddle #" << i << std::endl;
            // for (unsigned int l = 0 ; l < bwd_steps.size() ; ++l) {
            //      std::cout << l << ": " << bwd_steps[l] << std::endl;
            //          }
            // }
            
            std::cout << "drawing backward steps in yellow" << std::endl;
            
            // debugging
            for (unsigned int i = 0 ; i < bwd_steps.size() - 1 ; ++i) {
                double l = metric.distance(bwd_steps[i], bwd_steps[i+1]);
                if (l > 20) {
                    std::cout << "crazy segment between " << bwd_steps[i]
                              << " and " << bwd_steps[i+1] << std::endl;
                }
            }
            
            // determine proper orientation of additional separatrix
            // in case of period doubling
            if (pdoubling && !(n % 2)) {
                nvis::vec2 ev2 = saddle_chain[j].evec[1];
                if (nvis::inner(step, ev2) > 0) {
                    dirs[n+1] = -1;
                } else {
                    dirs[n+1] = 1;
                }
            }
            
            std::cout << "merging fwd and backward connections between saddles" << std::endl;
            merge(steps, fwd_steps, bwd_steps, metric);
            
            if (!valid_curve(steps, metric, one_chain)) {
            
                std::cout << "resulting curve has been deemed INVALID" << std::endl;
                everything_ok = false;
                
                continue;
            }
            
            std::cout << "resulting curve is VALID" << std::endl;
            
            std::cout << "refining valid curve through upsampling and smooth fit" << std::endl;
            std::vector< nvis::vec2 > fine;
            refine(fine, saddle_chain[i], saddle_chain[j], steps, map, metric);
            separatrices.push_back(fine);
        }
    } catch (...) {
        std::cout << "exception caught in process_saddle_chain" << std::endl;
        return false;
    }
    
    return transport(separatrices, map, p, metric);
}



} // map_analysis



#endif























































