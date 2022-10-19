#include "period.hpp"
#include <math/angle.hpp>

// #define VERBOSE_PERIOD_CPP

typedef xavier::default_metric_type metric_type;

inline bool empty(const std::pair<int, int>& p)
{
    return (p.first < 0);
}

inline bool full(const std::pair<int, int>& p)
{
    return (p.second >= 0);
}

inline void add(std::pair<int, int>& p, int i)
{
    assert(!full(p));
    if(empty(p)) {
        p.first = i;
    } else {
        p.second = i;
    }
}

bool compatible(int i, const std::pair<int, int>& p, int j,
                const std::vector<nvis::vec2>& pts,
                const metric_type& metric,
                std::ostream& os)
{
    using namespace xavier;
    
    if(empty(p)) {
        return true;
    }
    if(full(p)) {
        return false;
    } else if (i-p.first == j-i) {
        return true;
    }
    
    nvis::vec2 r = metric.displacement(pts[i], pts[p.first]);
    nvis::vec2 c = metric.displacement(pts[i], pts[j]);
#ifdef VERBOSE_PERIOD_CPP
    os << "comparing " << r << " (" << i << "->" << p.first << ") and "
       << c << " (" << i << "->" << j << "), angle=" << unsigned_angle<2>(r, c) << ")\n";
#endif
    return (to_degrees(unsigned_angle<2>(r, c)) >= 30);
}

std::ostream& operator<<(std::ostream& out, const std::pair<int, int>& p)
{
    out << "[next: " << p.first << ", prev: " << p.second << "]";
    
    return out;
}

void make_curves(std::vector<std::vector<int> >& curves,
                 std::vector<std::pair<int, int> >& edges)
{
    // connect vertices to form curves
    curves.clear();
    std::vector<bool> tagged(edges.size(), false);
    for(int i = 0 ; i < edges.size() ; ++i) {
        if (tagged[i] || empty(edges[i])) {
            continue;
        }
        
        std::list<int> c;
        c.push_back(i);
        tagged[i] = true;
        
        // forward
        c.push_back(edges[i].first);
        tagged[c.back()] = true;
        while(true) {
            int cur = c.back();
            if(empty(edges[cur])) {
                break;
            } else if(!tagged[edges[cur].first]) {
                c.push_back(edges[cur].first);
                tagged[c.back()] = true;
                continue;
            } else if(full(edges[cur]) && !tagged[edges[cur].second]) {
                c.push_back(edges[cur].second);
                tagged[c.back()] = true;
                continue;
            } else {
                break;
            }
        }
        
        // backward
        if(full(edges[i]) && !tagged[edges[i].second]) {
            c.push_front(edges[i].second);
            tagged[c.front()] = true;
            while(true) {
                int cur = c.front();
                if(empty(edges[cur])) {
                    break;
                } else if(!tagged[edges[cur].first]) {
                    c.push_front(edges[cur].first);
                    tagged[c.front()] = true;
                    continue;
                } else if(full(edges[cur]) && !tagged[edges[cur].second]) {
                    c.push_front(edges[cur].second);
                    tagged[c.front()] = true;
                    continue;
                } else {
                    break;
                }
            }
        }
        
        // did we close the loop?
        int a = c.front();
        int b = c.back();
        if(edges[a].first == b || edges[a].second == b) {
            c.push_back(a);
        }
        
        curves.push_back(std::vector<int>());
        std::copy(c.begin(), c.end(), std::back_inserter(curves.back()));
    }
}

void map_analysis::connect_by_distance(std::vector<std::vector<int> >& curves,
                                       const std::vector<nvis::vec2>& steps,
                                       const metric_type& metric)
{
    std::ostringstream os;
    
    typedef std::pair<int, int> pair_type;
    std::vector<pair_type>      v2e(steps.size(), pair_type(-1, -1));
    
    for(int i = 0 ; i < steps.size() ; ++i) {
#ifdef VERBOSE_PERIOD_CPP
        std::cerr << "\n\nprocessing vertex #" << i << " at " << metric.modulo(steps[i]) << std::endl;
        std::cerr << "this vertex currently has " << ((v2e[i].first >= 0) ? "a next neighbor " : "no next neighbor ")
                  << ((v2e[i].second >= 0) ? "and a prev neighbor" : "and no prev neighbor")
                  << std::endl;
#endif
                  
        if(is_full(v2e[i])) {
#ifdef VERBOSE_PERIOD_CPP
            std::cerr << "vertex #" << i << " is already fully connected\n";
#endif
            continue;
        }
        
        // determine closest points
        std::vector<int> neighbors;
        closest(neighbors, steps, i, steps.size(), metric);
        
        // compute local tangent vector
        nvis::vec2 tang = tangent(neighbors, steps, i, 4, metric);
#ifdef VERBOSE_PERIOD_CPP
        std::cerr << "tangent is " << tang << std::endl;
#endif
        
        std::map<double, int> pos_dist, neg_dist;
        
        for(int j = 0 ; j < neighbors.size() && (pos_dist.empty() || neg_dist.empty()) ; ++j) {
            nvis::vec2 dir = metric.displacement(steps[i], steps[neighbors[j]]);
            double dist = nvis::norm(dir);
            if(nvis::inner(dir, tang) >= 0) {
                pos_dist.insert(std::pair<double, int>(dist, neighbors[j]));
#ifdef VERBOSE_PERIOD_CPP
                std::cerr << "inserting next candidate neighbor #" << neighbors[j]
                          << " (" << metric.modulo(steps[neighbors[j]]) << ") at distance "
                          << dist << std::endl;
#endif
            } else {
                neg_dist.insert(std::pair<double, int>(dist, neighbors[j]));
#ifdef VERBOSE_PERIOD_CPP
                std::cerr << "inserting prev candidate neighbor #" << neighbors[j]
                          << " (" << metric.modulo(steps[neighbors[j]]) << ") at distance "
                          << dist << std::endl;
#endif
            }
        }
        
        if(v2e[i].first >= 0) {
            nvis::vec2 ref = metric.displacement(steps[i], steps[v2e[i].first]);
            if(nvis::inner(ref, tang) >= 0) {
                if(!neg_dist.empty()) {
                    v2e[i].second = neg_dist.begin()->second;
#ifdef VERBOSE_PERIOD_CPP
                    std::cerr << "selecting vertex " << v2e[i].second << " as prev neighbor\n";
#endif
                }
            } else if(!pos_dist.empty()) {
                v2e[i].second = pos_dist.begin()->second;
#ifdef VERBOSE_PERIOD_CPP
                std::cerr << "selecting vertex " << v2e[i].second << " as prev neighbor\n";
#endif
            }
        } else {
            if(!pos_dist.empty()) {
                v2e[i].first = pos_dist.begin()->second;
#ifdef VERBOSE_PERIOD_CPP
                std::cerr << "selecting vertex " << v2e[i].first << " as next neighbor\n";
#endif
            }
            if(!neg_dist.empty()) {
                if(v2e[i].first >= 0) {
                    v2e[i].second = neg_dist.begin()->second;
#ifdef VERBOSE_PERIOD_CPP
                    std::cerr << "selecting vertex " << v2e[i].second << " as prev neighbor\n";
#endif
                } else {
                    v2e[i].first = neg_dist.begin()->second;
#ifdef VERBOSE_PERIOD_CPP
                    std::cerr << "selecting vertex " << v2e[i].first << " as next neighbor\n";
#endif
                }
            }
        }
    }
    
    // connect vertices to form curves
    make_curves(curves, v2e);
}

void map_analysis::connect_and_check(std::vector<std::vector<int> >& curves,
                                     const std::vector<nvis::vec2>& steps,
                                     const metric_type& metric)
{
    std::ostringstream os;
    
    typedef std::pair<int, int> pair_type;
    std::vector<pair_type>      v2e(steps.size(), pair_type(-1, -1));
    
    for(int i = 0 ; i < steps.size() ; ++i) {
#ifdef VERBOSE_PERIOD_CPP
        std::cerr << "\n\nprocessing vertex #" << i << " at " << metric.modulo(steps[i]) << std::endl;
        std::cerr << "this vertex currently has " << ((v2e[i].first >= 0) ? "a next neighbor " : "no next neighbor ")
                  << ((v2e[i].second >= 0) ? "and a prev neighbor" : "and no prev neighbor")
                  << std::endl;
#endif
                  
        // determine closest points
        std::vector<int> neighbors;
        closest(neighbors, steps, i, steps.size(), metric);
        
        // compute local tangent vector
        nvis::vec2 tang = tangent(neighbors, steps, i, 4, metric);
        
#ifdef VERBOSE_PERIOD_CPP
        std::cerr << "tangent is " << tang << std::endl;
#endif
        
        std::map<double, int> pos_dist, neg_dist;
        
        for(int j = 0 ; j < neighbors.size() && (pos_dist.empty() || neg_dist.empty()) ; ++j) {
            nvis::vec2 dir = metric.displacement(steps[i], steps[neighbors[j]]);
            double dist = nvis::norm(dir);
            if(nvis::inner(dir, tang) >= 0) {
                pos_dist.insert(std::pair<double, int>(dist, neighbors[j]));
#ifdef VERBOSE_PERIOD_CPP
                std::cerr << "inserting next candidate neighbor #" << neighbors[j]
                          << " (" << metric.modulo(steps[neighbors[j]]) << ") at distance "
                          << dist << std::endl;
#endif
            } else {
                neg_dist.insert(std::pair<double, int>(dist, neighbors[j]));
#ifdef VERBOSE_PERIOD_CPP
                std::cerr << "inserting prev candidate neighbor #" << neighbors[j]
                          << " (" << metric.modulo(steps[neighbors[j]]) << ") at distance "
                          << dist << std::endl;
#endif
            }
        }
        
        if(!pos_dist.empty()) {
            v2e[i].first = pos_dist.begin()->second;
#ifdef VERBOSE_PERIOD_CPP
            std::cerr << "selecting vertex " << v2e[i].first << " as next neighbor\n";
#endif
        }
        if(!neg_dist.empty()) {
            if(v2e[i].first >= 0) {
                v2e[i].second = neg_dist.begin()->second;
#ifdef VERBOSE_PERIOD_CPP
                std::cerr << "selecting vertex " << v2e[i].second << " as prev neighbor\n";
#endif
            } else {
                v2e[i].first = neg_dist.begin()->second;
#ifdef VERBOSE_PERIOD_CPP
                std::cerr << "selecting vertex " << v2e[i].first << " as next neighbor\n";
#endif
            }
        }
    }
    
    // detect and resolve conflicts
    std::set<int> modified;
    for(int i = 0 ; i < v2e.size() ; ++i) {
        std::pair<int, int>& cur = v2e[i];
        if(cur.first >= 0) {
            int next = cur.first;
            if(v2e[next].first != i && v2e[next].second != i) {
                // incompatible connections
                cur.first = -1;
                modified.insert(i);
            }
        }
        if(cur.second >= 0) {
            int next = cur.second;
            if(v2e[next].first != i && v2e[next].second != i) {
                // incompatible connections
                cur.second = -1;
                modified.insert(i);
            }
        }
    }
    
#ifdef VERBOSE_PERIOD_CPP
    std::cerr << modified.size() << " vertices out of " << steps.size() << " have been modified\n";
#endif
    
    // connect vertices to form curves
    make_curves(curves, v2e);
}

void map_analysis::connect_by_period(std::vector<std::vector<int> >& curves,
                                     const std::vector<nvis::vec2>& steps,
                                     const metric_type& metric)
{
    std::ostringstream os;
    
    typedef std::pair<int, int> pair_type;
    std::vector<pair_type>      edges(steps.size(), pair_type(-1, -1));
    
#ifdef VERBOSE_PERIOD_CPP
    os << "\ncurrent chain is:\n";
    for (int i=0 ; i<steps.size() ; ++i) {
        os << i << ": " << metric.modulo(steps[i]) << '\n';
    }
    
    os << "\nlooping over periods:\n";
#endif
    std::set<distance_profile, Lt_dist_profile> quality;
    for(int i = 1 ; i < steps.size() ; ++i) {
        quality.insert(distance_profile(steps, i, metric));
    }
    
    typedef std::set<distance_profile, Lt_dist_profile>::const_iterator iterator_type;
    for(iterator_type it = quality.begin() ; it != quality.end() ; ++it) {
        if (it->_median > 0.05*metric.width()) {
            break;
        }
#ifdef VERBOSE_PERIOD_CPP
        os << "\nperiod " << it->_p << " has average length " << it->_mean << '\n';
#endif
        int p = it->_p;
        for(int i = 0 ; i + p < steps.size() ; ++i) {
#ifdef VERBOSE_PERIOD_CPP
            os << "vertex #" << i << " has connections " << edges[i] << '\n';
            os << "vertex #" << i << "+" << p << "=" << i+p << " has connections " << edges[i+p] << '\n';
#endif
            if(compatible(i, edges[i], i + p, steps, metric, os) &&
                    compatible(i + p, edges[i + p], i, steps, metric, os)) {
                nvis::vec2 l = metric.displacement(steps[i], steps[i+p]);
                add(edges[i], i + p);
                add(edges[i + p], i);
#ifdef VERBOSE_PERIOD_CPP
                os << "created new edges: " << edges[i] << " and " << edges[i+p] << '\n';
                
                if (nvis::norm(l) > 2) {
                    os << "\n\n\nWARNING: this edge between " << metric.modulo(steps[i])
                       << " and " << metric.modulo(steps[i+p]) << " has length " << nvis::norm(l)
                       << '\n';
                }
#endif
                
            }
        }
    }
    os << std::flush;
    
    std::cerr << os.str();
    
    // connect vertices to form curves
    make_curves(curves, edges);
    
#ifdef VERBOSE_PERIOD_CPP
    std::cerr << "exported curves:\n";
    for (int i=0 ; i<curves.size() ; ++i) {
        std::cerr << "curve #" << i << '\n';
        int last = curves[i].front();
        for (int j=0 ; j<curves[i].size() ; ++j) {
            int cur = curves[i][j];
            std::cerr << cur-last << " -> ";
            last = cur;
        }
    }
#endif
}
