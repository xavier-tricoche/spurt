#ifndef __PERIOD_ANALYSIS_HPP__
#define __PERIOD_ANALYSIS_HPP__

#include <string>
#include <list>
#include <vector>
#include <algorithm>
#include <math/fixed_vector.hpp>
#include "definitions.hpp"
#include "orbits.hpp"
#include "index.hpp"
#include <poincare/metric.hpp>
#include <data/raster.hpp>
#include "period.hpp"

#if _OPENMP
#include <omp.h>
#endif

namespace spurt {
struct map_cell_data {

    typedef std::pair<nvis::vec2, double> value_type;
    
    map_cell_data()
        : points(), unique_periods(), sum(0), depth(0),
          __min(std::numeric_limits<double>::max()),
          __max(std::numeric_limits<double>::min()) {}
    map_cell_data(const map_cell_data& data)
        : points(data.points), unique_periods(data.unique_periods), __min(data.__min),
          __max(data.__max), sum(data.sum), depth(data.depth) {}
          
    map_cell_data& operator=(const map_cell_data& cd) {
        points.clear();
        std::copy(cd.points.begin(), cd.points.end(), std::back_inserter(points));
        unique_periods.clear();
        unique_periods = cd.unique_periods;
        __min = cd.__min;
        __max = cd.__max;
        sum = cd.sum;
        depth = cd.depth;
        
        return *this;
    }
    
    void add(const nvis::vec2& x, double v) {
        points.push_back(value_type(x, v));
        if (unique_periods.find(v) == unique_periods.end()) {
            unique_periods.insert(v);
            sum += v;
            if (v<__min) {
                __min=v;
            }
            if (v>__max) {
                __max = v;
            }
        }
    }
    
    double min() const {
        return __min;
    }
    
    double max() const {
        return __max;
    }
    
    double mean() const {
        if (count()) {
            return sum/(double)count();
        } else {
            return 0;
        }
    }
    
    size_t count() const {
        return unique_periods.size();
    }
    
    double span() const {
        return __max-__min;
    }
    
    double rel_span() const {
        if (count()) {
            return span()/mean();
        } else {
            return std::numeric_limits<double>::max();
        }
    }
    
    void split(std::vector<map_cell_data>& children, const nvis::bbox2& bounds) const {
        children.clear();
        children.resize(4);
        
        nvis::vec2 half_step = bounds.size()/2.;
        const nvis::vec2& p0 = bounds.min();
        
        for (int i=0 ; i<points.size() ; ++i) {
            nvis::vec2 x = (points[i].first-p0)/half_step;
            int u = floor(x[0]);
            int v = floor(x[1]);
            int k = u + 2*v;
            children[k].add(points[i].first, points[i].second);
        }
        children[0].depth = children[1].depth = children[2].depth = children[3].depth = depth+1;
    }
    
    std::vector<value_type> points;
    std::set<double> unique_periods;
    double __min, __max, sum;
    int depth;
};

struct Lt_cell_data {
    bool operator()(const map_cell_data& c0, const map_cell_data& c1) {
        return (c0.count() < c1.count());
    }
};

struct period_convergence_predicate {
    period_convergence_predicate(int nsamples,
                                 const std::vector<double>& periods, double eps_converge, double eps_close)
        : min_nb_samples(nsamples), ref_periods(periods), tol_cvg(eps_converge),
          tol_close(eps_close) {
        std::sort(ref_periods.begin(), ref_periods.end());
    }
    
    bool split(const map_cell_data& cd) const {
        if (cd.count() < min_nb_samples || cd.rel_span() > tol_cvg) {
            return true;
        }
        std::vector<double>::const_iterator it0 =
            std::lower_bound(ref_periods.begin(), ref_periods.end(), cd.min()-tol_close);
        std::vector<double>::const_iterator it1 =
            std::upper_bound(ref_periods.begin(), ref_periods.end(), cd.max()+tol_close);
        return (it0 != it1);
    }
    
    int howmany(const map_cell_data& cd) const {
        return std::max(0, min_nb_samples - (int)cd.count());
    }
    
    int min_nb_samples;
    std::vector<double> ref_periods;
    double tol_cvg, tol_close;
};
} // spurt

namespace {

template<typename P>
double period(const nvis::vec2& x0, std::vector<nvis::vec2>& steps,
              const P& __pmap, const spurt::map_metric& metric, int n)
{
    const P* pmap = __pmap.clone();
    try {
        pmap->map(x0, steps, n);
    } catch (...) {
        return -1;
    }
    
    return map_analysis::period_x_periodic(steps, metric).first;
}
} // anonymous

namespace spurt {

template<typename MAP, typename PREDICATE>
void period_analysis(const MAP& pmap, const map_metric& metric,
                     const nvis::ivec2& res, const PREDICATE& cell_checker, int niter,
                     int max_depth,
                     std::vector<nvis::vec2>& edges,
                     std::vector<std::pair<double, std::vector<nvis::vec2> > >& orbits,
                     std::vector<std::pair<nvis::ivec3, nvis::vec3> >& quads)
{
    typedef MAP                                         pmap_type;
    typedef PREDICATE                                   split_predicate;
    typedef nvis::vec2                                  point_data;
    typedef map_cell_data                               cell_data;
    typedef std::pair<nvis::ivec3, nvis::vec3>          quad_type;
    
    const nvis::ivec2 offset[] = { nvis::ivec2(0,0),
                                   nvis::ivec2(1,0),
                                   nvis::ivec2(0,1),
                                   nvis::ivec2(1,1)
                                 };
                                 
    map_metric::bounds_type bounds = metric.bounds();
    nvis::vec2 href = bounds.size()/nvis::vec2(res[0], res[1]);
    nvis::ivec2 r(res);
    
    nvis::vec2 h = href;
    
    // initialize data structure
    std::vector<cell_data>  mesh(r[0]*r[1]);
    
    size_t nb_threads = 1;
#if _OPENMP
    nb_threads = omp_get_max_threads();
    std::cout << nb_threads << " threads available\n";
#endif
    
    // initially every cell must be considered
    std::vector<int> to_process(res[0]*res[1]);
    for (int i=0 ; i<to_process.size() ; ++i) {
        to_process[i] = i;
    }
    
    struct thread_cache {
        thread_cache() : nb_touched(0), nb_in(0), nb_skipped(0), points() {}
        
        int nb_touched, nb_in, nb_skipped;
        std::vector<std::pair<nvis::vec2, double> > points;
    };
    
    // adaptive refinement
    int actual_depth = 0;
    for (int depth=0 ; depth<=max_depth ; ++depth) {
        actual_depth = depth;
        
        std::cerr << "\n\tdepth " << depth << " out of " << max_depth << '\n';
        
        std::cerr << to_process.size() << " cells to process out of " << r[0]*r[1] << "\n";
        
        std::random_shuffle(to_process.begin(), to_process.end());
        
        thread_cache cache[nb_threads];
        std::vector<bool> flagged(mesh.size(), false);
        
        nvis::ivec2 finer = 2*r;
        
        int id_step = to_process.size()/100;
        for (int id_start = 0 ; id_start<to_process.size() ; id_start += id_step) {
            int id_end = std::min(id_start+id_step, (int)to_process.size());
            for (int i=0 ; i<nb_threads ; ++i) {
                cache[i].points.clear();
            }
            
#if 0
            std::cerr << '\n';
            for (int __h=0 ; __h<r[1] ; ++__h) {
                for (int __w=0 ; __w<r[0] ; ++__w) {
                    int n = __w + __h*r[0];
                    if (flagged[n]) {
                        std::cerr << '!';
                    } else if (mesh[n].depth < actual_depth) {
                        std::cerr << 'X';
                    } else {
                        std::cerr << std::min(mesh[n].count, 9);
                    }
                }
                std::cerr << '\n';
            }
            std::cerr << '\n';
#endif
            
            
            #pragma omp parallel
            {
                int thread_id = 0;
                
                #pragma omp for schedule(dynamic,1)
                for (int __n=id_start ; __n<id_end ; ++__n) {
#if _OPENMP
                    thread_id = omp_get_thread_num();
#endif
                    int n = to_process[__n];
                    int i = n%r[0];
                    int j = n/r[0];
                    
                    cache[thread_id].nb_in += mesh[n].count();
                    ++cache[thread_id].nb_touched;
                    
                    int m = cell_checker.howmany(mesh[n]);
                    if (!m) {
                        ++cache[thread_id].nb_skipped;
                        continue;
                    }
                    
                    std::vector< nvis::vec2 > steps;
                    for (int __m=0 ; __m<m ; ++__m) {
                        nvis::vec2 x0 = bounds.min() + nvis::vec2(i + drand48(), j + drand48())*h;
                        double q = period<MAP>(x0, steps, pmap, metric, niter);
                        if (q <= 0) {
                            flagged[n] = true;
                            break;
                        }
                        
                        for (int k=0 ; k<steps.size() ; ++k) {
                            nvis::vec2 x = metric.modulo(steps[k]);
                            cache[thread_id].points.push_back(std::pair<nvis::vec2, double>(x, q));
                        }
                    }
                    if (thread_id == 0) {
                        int skipped = 0;
                        int touched = 0;
                        int ins = 0;
                        for (int ii=0 ; ii<nb_threads ; ++ii) {
                            skipped += cache[ii].nb_skipped;
                            touched += cache[ii].nb_touched;
                            ins += cache[ii].nb_in;
                        }
                        float prec = (float)ins/(float)touched;
                        std::ostringstream os;
                        os << "\rsampled " << __n << " out of " << to_process.size()
                           << ", skipped " << skipped << " (" << 100*skipped/__n << "%)"
                           << ", " << prec << " precomputed on average           \r" << std::flush;
                        std::cout << os.str();
                    }
                }
            } // end parallel sampling
            
            // update data structure
            for (int i=0 ; i<nb_threads ; ++i) {
                for (int j=0 ; j<cache[i].points.size() ; ++j) {
                    nvis::vec2 x = cache[i].points[j].first;
                    nvis::vec2 id = (x - bounds.min()) / h;
                    int _i = floor(id[0]);
                    int _j = floor(id[1]);
                    int _n = _i + r[0]*_j;
                    
                    nvis::vec2 min = bounds.min() + nvis::vec2(_i, _j)*h;
                    nvis::bbox2 box(min, min+h);
                    mesh[_n].add(x, cache[i].points[j].second);
                }
                
                double q = std::numeric_limits<double>::min(); // invalid value
                typedef std::pair<double, std::vector<nvis::vec2> > orbit_type;
                for (int j=0 ; j<cache[i].points.size() ; ++j) {
                    double _q = cache[i].points[j].second;
                    if (_q != q) {
                        orbits.push_back(orbit_type());
                        orbits.back().first = _q;
                        q = _q;
                    }
                    orbits.back().second.push_back(cache[i].points[j].first);
                }
            }
        }
        
        if (depth < max_depth) {
            std::sort(to_process.begin(), to_process.end());
            std::set<int> active(to_process.begin(), to_process.end());
            to_process.clear();
            std::vector<cell_data>  new_mesh(4*r[0]*r[1]);
            for (int n=0 ; n<mesh.size() ; ++n) {
                int i = n%r[0];
                int j = n/r[0];
                nvis::vec2 min = bounds.min() + nvis::vec2(i,j)*h;
                nvis::bbox2 box(min, min+h);
                std::vector<cell_data> children;
                mesh[n].split(children, box);
                nvis::ivec2 base(2*i, 2*j);
                for (int k=0 ; k<4 ; ++k) {
                    nvis::ivec2 id = base + offset[k];
                    int _n = id[0] + finer[0]*id[1];
                    new_mesh[_n] = children[k];
                }
                
                bool is_active = (active.find(n) != active.end());
                
                if (!is_active || flagged[n] || !cell_checker.split(mesh[n])) {
                    for (int k=0 ; k<4 ; ++k) {
                        nvis::ivec2 id = base + offset[k];
                        int _n = id[0] + finer[0]*id[1];
                        new_mesh[_n].depth = mesh[n].depth;
                    }
                } else {
                    for (int k=0 ; k<4 ; ++k) {
                        nvis::ivec2 id = base + offset[k];
                        int _n = id[0] + finer[0]*id[1];
                        to_process.push_back(_n);
                        new_mesh[_n].depth = mesh[n].depth + 1;
                    }
                }
            }
            
            r = finer;
            h /= 2.;
            mesh.swap(new_mesh);
        }
    }
    
    std::cerr << "\nmax depth = " << actual_depth << std::endl;
    
    // export edges of finest available sampling resolution
    typedef nvis::ivec3 index_type;
    std::set<index_type, nvis::lexicographical_order> used;
    for (int n=0 ; n<mesh.size() ; ++n) {
        int i = n%r[0];
        int j = n/r[0];
        int d = 1 << (actual_depth - mesh[n].depth);
        // std::cerr << "adding " << index_type(i/d, j/d, mesh[n].depth) << " (" << d << ")" << '\n';
        used.insert(index_type(i/d, j/d, mesh[n].depth));
    }
    
    quads.clear();
    for (std::set<index_type>::const_iterator it=used.begin() ;
            it!=used.end() ; ++it) {
        int i = (*it)[0];
        int j = (*it)[1];
        int d = (*it)[2];
        int f = (1 << d);
        h = href / (double)f;
        nvis::vec2 x(h[0], 0);
        nvis::vec2 y(0, h[1]);
        nvis::vec2 p0 = bounds.min() + nvis::vec2(i,j)*h;
        edges.push_back(p0);
        edges.push_back(p0+x);
        edges.push_back(p0+x+y);
        edges.push_back(p0+y);
        
        int delta_d = actual_depth - (*it)[2];
        f = (1 << delta_d);
        i *= f;
        j *= f;
        int n = i + j*r[0];
        quads.push_back(quad_type(*it, nvis::vec3(mesh[n].min(), mesh[n].mean(), mesh[n].max())));
    }
}

}

#endif
