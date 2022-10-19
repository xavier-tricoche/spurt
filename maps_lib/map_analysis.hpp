#ifndef __MAP_ANALYSIS_HPP__
#define __MAP_ANALYSIS_HPP__

// std
#include <iostream>
#include <vector>
#include <list>
#include <map>
#include <iomanip>
#include <sstream>
#include <set>
// boost
#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <boost/rational.hpp>
// nvis
#include <math/fixed_vector.hpp>
#include <util/wall_timer.hpp>
// teem
#include <teem/nrrd.h>
// xavier
#include <data/grid.hpp>
#include <data/kdtree.hpp>
#include <data/raster_data.hpp>
#include <math/angle.hpp>
#include <misc/sort.hpp>

#if _OPENMP
#include <omp.h>
#endif

namespace {
template<typename T>
inline T _sign(T x)
{
    return (x<0) ? -1 : 1;
}
}

#define SECANT_METHOD


namespace xavier {

const double MAX_ANGLE_VARIATION = M_PI/2.;

int best_period(const std::vector<nvis::vec2>& steps, int maxp, 
                const map_metric& metric);
double safety_factor(const std::vector<nvis::vec2>& steps,
                     const map_metric& metric);
double safety_factor_bounded(const std::vector<nvis::fixed_vector<double, 2ul> >&, 
                             const xavier::map_metric&);
double safety_factor_stdmap(const std::vector<nvis::fixed_vector<double, 2ul> >&, 
                             const xavier::map_metric&);

struct orbit_data {
    static map_metric   metric;
    static int          max_period;
    
    orbit_data() : sf(0), p(0) {}
    
    void set(const std::vector<nvis::vec2>& _s) {
        std::copy(_s.begin(), _s.end(), std::back_inserter(steps));
        sf = safety_factor(steps, metric);
        p = best_period(steps, max_period, metric);
    }
    
    void set() {
        sf = safety_factor(steps, metric);
        p = best_period(steps, max_period, metric);
    }
    
    nvis::vec2 vector(int p) const {
        if (p>=steps.size()) {
            std::cerr << "requested period " << p << " while steps contains " << steps.size() << " positions\n";
            throw std::runtime_error("period exceeds available data");
        }
        return metric.displacement(steps[0], steps[p]);
    }
    
    std::vector<nvis::vec2> steps;
    double sf;
    int p;
};

struct map_analysis_param {
    typedef std::vector<nvis::vec2>                         orbit_type;
    typedef std::pair<orbit_type, int>                      tagged_orbit_type;
    typedef std::pair<nvis::vec2, int>                      tagged_vector_type;
    
    map_analysis_param()
        : upsampling_factor(1), max_depth(1), min_period(1), max_period(0),
          nb_iterations(50), record(false), verbose(false), lmin(0)
    {}
    
    nvis::ivec2                                             resolution;
    nvis::bbox2                                             bounds;
    map_metric                                              metric;
    int                                                     upsampling_factor;
    int                                                     max_depth;
    int                                                     min_period;
    int                                                     max_period;
    int                                                     nb_iterations;
    std::vector<nvis::vec2>                                 edges;
    std::vector<tagged_orbit_type>                          orbits;
    std::vector<std::pair<nvis::vec2, tagged_vector_type> > vectors;
    bool                                                    record;
    bool                                                    verbose;
    double                                                  lmin;
    std::vector<boost::rational<int> >                      valid_rationals;
};

struct rational_surface_suspected : public std::runtime_error {
    explicit rational_surface_suspected(const std::string& __arg)
        : std::runtime_error(__arg) {}
};

struct higher_precision_required : public std::runtime_error {
    explicit higher_precision_required(const std::string& __arg)
        : std::runtime_error(__arg) {}
};

struct rational_surface_found : public std::runtime_error {
    typedef std::pair<nvis::ivec2, nvis::ivec2> edge_type;
    typedef std::pair<edge_type, nvis::vec2>    edge_point;
    
    explicit rational_surface_found(const std::string& __arg)
        : std::runtime_error(__arg) {}
        
    edge_point xings[2];
    boost::rational<int> sf;
};

struct index_step_size_underflow : public std::underflow_error {
    explicit index_step_size_underflow(const std::string& __arg)
        : std::underflow_error(__arg) {}
        
    nvis::vec2 where;
};

double average_distance(const std::vector<nvis::vec2>& steps, int p,
                        const map_metric& metric);
                        
template<typename MAP>
void iterate(const nvis::vec2& x0, std::vector<nvis::vec2>& steps,
             const MAP& pmap, int niter);
             
template<typename MAP>
nvis::vec2 rhs(const nvis::vec2& x0, const MAP& pmap, int period,
               map_analysis_param& params);
               
double secant(const nvis::vec2& v0, const nvis::vec2& v1)
{
    const double threshold = 0.75*M_PI;
    const double tiny = 1.0e-3;
    
    nvis::vec2 a = v0 / nvis::norm(v0);
    nvis::vec2 b = v1 / nvis::norm(v1);
    if (acos(nvis::inner(a, b)) < threshold) {
        return 0.5;
    }
    if (nvis::norm(v0) < tiny || nvis::norm(v1) < tiny) {
        return 0.5;
    }
    
    nvis::vec2 dir = v0 - v1;
    dir /= nvis::norm(dir);
    double f0 = nvis::inner(v0, dir);
    double f1 = nvis::inner(v1, dir);
    
    return -f0/(f1-f0);
}

template<typename MAP>
double adaptive_rotation_angle(const nvis::vec2& x0, const nvis::vec2& v0,
                               const nvis::vec2& x1, const nvis::vec2& v1,
                               const MAP& pmap, unsigned int period,
                               double dx, map_analysis_param& params);
                               
template<typename MAP>
double robust_rotation_angle(const nvis::vec2& x0, const nvis::vec2& v0,
                             const nvis::vec2& x1, const nvis::vec2& v1,
                             const MAP& pmap, unsigned int period,
                             double dx, map_analysis_param& params);
                             
template<typename MAP>
int poincare_index(const grid<double, 2>& mesh, const nvis::ivec2& cell_id,
                   const raster_data<orbit_data, double, 2>& dataset,
                   const MAP& pmap, int period, double lmin,
                   map_analysis_param& params);
                   
template<typename MAP>
int __poincare_index(const grid<double, 2>& mesh,
                     const nvis::ivec2& cell_id,
                     const raster_data<orbit_data, double, 2>& dataset,
                     const MAP& pmap, int period, double lmin,
                     bool scrutinize,
                     std::vector<std::pair<int, nvis::vec2> >& sing_edges,
                     map_analysis_param& params);
                     
template<typename MAP>
int __fast_poincare_index(const grid<double, 2>& mesh,
                          const nvis::ivec2& cell_id,
                          const raster_data<orbit_data, double, 2>& dataset,
                          const MAP& pmap, int period, double lmin,
                          bool scrutinize,
                          std::vector<std::pair<int, nvis::vec2> >& sing_edges,
                          map_analysis_param& params);
                          
template<typename MAP>
void sample_raster(raster_data<orbit_data, double, 2>& dataset,
                   const grid<double, 2>& grid,
                   const MAP& pmap, map_analysis_param& params,
                   bool compute_sf = true);
                   
void period_range(std::vector<int>& ps,
                  const raster_data<orbit_data, double, 2>& dataset,
                  const nvis::ivec2& cell_id,
                  const std::vector<boost::rational<int> >& valid_rationals);
                  
void best_periods(std::vector<int>& ps,
                  const raster_data<orbit_data, double, 2>& dataset,
                  const nvis::ivec2& cell_id, const map_metric& metric,
                  int maxperiod, int nbperiods);
                  
template<typename MAP>
void
process_cell(std::vector<std::pair<int, int> >& pids,
             const raster_data<orbit_data, double, 2>& dataset,
             const grid<double, 2>& grid,
             const nvis::ivec2& cell_id, const MAP& pmap, double dx,
             map_analysis_param& params);
}

inline double xavier::average_distance(const std::vector<nvis::vec2>& steps, int p,
                                       const map_metric& metric)
{
    if (steps.size() < p) {
        return std::numeric_limits<double>::max();
    }
    double d = 0;
    int count = 0;
    for (int i=0 ; i<steps.size() ; ++i) {
        if (i+p<steps.size()) {
            d += metric.distance(steps[i], steps[i+p]);
            ++count;
        }
        if (i>=p) {
            d += metric.distance(steps[i-p], steps[i]);
            ++count;
        }
    }
    return d/(double)count;
}

inline void xavier::period_range(std::vector<int>& ps,
                                 const raster_data<orbit_data, double, 2>& dataset,
                                 const nvis::ivec2& cell_id,
                                 const std::vector<boost::rational<int> >& valid_rationals)
{

    std::vector<double> sfs;
    sfs.push_back(dataset(cell_id).sf);
    sfs.push_back(dataset(cell_id + nvis::ivec2(1,0)).sf);
    sfs.push_back(dataset(cell_id + nvis::ivec2(1,1)).sf);
    sfs.push_back(dataset(cell_id + nvis::ivec2(0,1)).sf);
    
    // insert periods associated with vertices
    std::set<int> nums;
    nums.insert(dataset(cell_id).p);
    nums.insert(dataset(cell_id + nvis::ivec2(1,0)).p);
    nums.insert(dataset(cell_id + nvis::ivec2(1,1)).p);
    nums.insert(dataset(cell_id + nvis::ivec2(0,1)).p);
    
    double min = *std::min_element(sfs.begin(), sfs.end());
    double max = *std::max_element(sfs.begin(), sfs.end());
    
    if (min == 0) {
        // integration failed at one of the vertices. give up.
        ps.clear();
        // std::cerr << "integration failed\n";
        return;
    }
    
    int N = valid_rationals.size();
    int mini = 0;
    int maxi = 0;
    double bestmin = std::numeric_limits<double>::max();
    double bestmax = std::numeric_limits<double>::max();
    for (int i=0 ; i<N ; ++i) {
        double v = value<int, double>(valid_rationals[i]);
        if (fabs(v-min) < bestmin) {
            bestmin = fabs(v-min);
            mini = i;
        }
        if (fabs(v-max) < bestmax) {
            bestmax = fabs(v-max);
            maxi = i;
        }
    }
    
    for (int i=mini ; i<=maxi ; ++i) {
        nums.insert(valid_rationals[i].numerator());
    }
    ps.clear();
    std::copy(nums.begin(), nums.end(), std::back_inserter(ps));
}

inline void xavier::best_periods(std::vector<int>& ps,
                                 const raster_data<orbit_data, double, 2>& dataset,
                                 const nvis::ivec2& cell_id, const map_metric& metric,
                                 int maxperiod, int nbperiods)
{
    typedef std::set<std::pair<int, double>, Lt_pair_second<int, double> > set_type;
    ps.clear();
    set_type periods;
    std::set<int> selected;
    nvis::ivec2 ids[4] = { cell_id, cell_id + nvis::ivec2(1,0),
                           cell_id + nvis::ivec2(2,1), cell_id + nvis::ivec2(0,1)
                         };
                         
    for (int i=0 ; i<4 ; ++i) {
        set_type candidates;
        for (int n=0 ; n<maxperiod ; ++n) {
            double d = average_distance(dataset(ids[i]).steps, n+1, metric);
            candidates.insert(std::pair<int, double>(n+1, d));
        }
        if (candidates.empty()) {
            continue;
        }
        set_type::const_iterator it = candidates.begin();
        selected.insert((it++)->first);
        for (int i=1 ; i<nbperiods && it!=candidates.end() ; ++i) {
            periods.insert(*(it++));
        }
    }
    int nbselected = selected.size();
    if (nbselected < nbperiods) {
        set_type::const_iterator it = periods.begin();
        for (int i=nbselected ; i<nbperiods && it!=periods.end() ; ++i) {
            selected.insert((it++)->first);
        }
    }
    std::copy(selected.begin(), selected.end(), std::back_inserter(ps));
}


inline double xavier::safety_factor(const std::vector<nvis::vec2>& steps,
                                    const map_metric& metric)
{
    if (steps.size()<2) {
        return 0;
    }
    
    std::ostringstream os;
    os << "steps" << std::endl;
    std::vector< double > q(steps.size()-1);
    std::vector< double > d(steps.size()-1);
    const nvis::vec2& x0 = steps[0];
    for (unsigned int i = 1 ; i < steps.size() ; ++i) {
        double dtheta = (steps[i] - x0)[0];
        q[i-1] = (double)i / dtheta * metric.width();
        d[i-1] = metric.distance(x0, steps[i]);
        os << q[i-1] << std::endl;
    }
    double min = fabs(*std::min_element(d.begin(), d.end(), Lt_fabs()));
    
    // os << std::endl;
    // os << "normalized steps" << std::endl;
    double total_q = 0;
    double total_w = 0;
    for (int i = 1 ; i < steps.size() ; ++i) {
        double w = min / fabs(d[i-1]);
        total_w += w;
        total_q += w * q[i-1];
        // os << total_q/total_w << std::endl;
    }
    
    // std::cerr << os.str();
    
    return total_q/total_w;
}

inline double xavier::safety_factor_bounded(const std::vector<nvis::vec2>& steps,
        const map_metric& metric)
{
    if (steps.size()<2) {
        return 0;
    }
    const nvis::vec2& x0 = steps[0];
    
    // compute center of gravity
    nvis::vec2 center = steps[0];
    for (int i=1 ; i<steps.size() ; ++i) {
        center += steps[i];
    }
    center /= steps.size();
    
    double last_theta = xavier::positive_angle(steps[0]-center);
    double dtheta = 0;
    
    std::vector< double > q(steps.size()-1);
    std::vector< double > d(steps.size()-1);
    for (int i=1 ; i<steps.size() ; ++i) {
        double theta = xavier::positive_angle(steps[i]-center);
        double dt = theta - last_theta;
        if (dt > M_PI) {
            dt = 2.*M_PI - dt;
        } else if (dt < -M_PI) {
            dt = 2.*M_PI + dt;
        }
        dtheta += dt;
        q[i-1] = 2.*M_PI*(double)i/dtheta;
        d[i-1] = metric.distance(x0, steps[i]);
    }
    double min = fabs(*std::min_element(d.begin(), d.end(), Lt_fabs()));
    
    double total_q = 0;
    double total_w = 0;
    for (int i = 1 ; i < steps.size() ; ++i) {
        double w = min / fabs(d[i-1]);
        total_w += w;
        total_q += w * q[i-1];
    }
    
    return total_q/total_w;
}

inline double xavier::safety_factor_stdmap(const std::vector<nvis::vec2>& steps,
        const map_metric& metric)
{
    if (steps.size()<2) {
        return 0;
    }
    const nvis::vec2& x0 = steps[0];
    
    std::vector< double > q(steps.size()-1);
    std::vector< double > d(steps.size()-1);
    for (int i=1 ; i<steps.size() ; ++i) {
        nvis::vec2 x = steps[i];
        double theta = x[0] + x[1];
        q[i-1] = (double)i/theta;
        d[i-1] = metric.distance(x0, x);
    }
    double min = fabs(*std::min_element(d.begin(), d.end(), Lt_fabs()));
    
    double total_q = 0;
    double total_w = 0;
    for (int i = 1 ; i < steps.size() ; ++i) {
        double w = min / fabs(d[i-1]);
        total_w += w;
        total_q += w * q[i-1];
    }
    
    return total_q/total_w;
}


inline int xavier::best_period(const std::vector<nvis::vec2>& steps, int maxp,
                               const map_metric& metric)
{
    std::map<double, int> dist_to_period;
    for (int p=1 ; p<=maxp ; ++p) {
        dist_to_period[average_distance(steps, p, metric)] = p;
    }
    return dist_to_period.begin()->second;
}

template<typename MAP>
void xavier::iterate(const nvis::vec2& x0, std::vector<nvis::vec2>& steps,
                     const MAP& pmap, int niter)
{
    MAP* amap = pmap.clone();
    std::vector<nvis::vec2> tmp;
    steps.clear();
    steps.push_back(x0);
    try {
        amap->map(x0, tmp, niter);
    } catch(...) {
        return;
    }
    std::copy(tmp.begin(), tmp.end(), std::back_inserter(steps));
}

template<typename MAP>
nvis::vec2 xavier::rhs(const nvis::vec2& x0, const MAP& pmap, int period,
                       const map_metric& metric)
{
    MAP* amap = pmap.clone();
    std::vector<nvis::vec2> steps;
    try {
        amap->map(x0, steps, period);
    } catch(...) {
        throw std::runtime_error("invalid result of map iteration");
    }
    if (steps.size() < period) {
        throw std::runtime_error("invalid result of map iteration");
    } else {
        nvis::vec2 v = metric.displacement(x0, steps[period-1]);
        return v;
    }
}

template<typename MAP>
nvis::vec2 xavier::rhs(const nvis::vec2& x0, const MAP& pmap, int period,
                       map_analysis_param& params)
{
    MAP* amap = pmap.clone();
    std::vector<nvis::vec2> steps;
    try {
        amap->map(x0, steps, period);
    } catch(...) {
        throw std::runtime_error("invalid result of map iteration");
    }
    if (steps.size() < period) {
        throw std::runtime_error("invalid result of map iteration");
    } else {
        nvis::vec2 v = params.metric.displacement(x0, steps[period-1]);
        if (params.record) {
            map_analysis_param::tagged_vector_type tv(v, period);
            params.vectors.push_back(std::pair<nvis::vec2,
                                     map_analysis_param::tagged_vector_type>(x0, tv));
        }
        return v;
    }
}

template<typename MAP>
void xavier::sample_raster(raster_data<orbit_data, double, 2>& dataset,
                           const grid<double, 2>& mesh, const MAP& pmap,
                           map_analysis_param& params, bool compute_sf)
{
    typedef grid<double, 2>                         grid_type;
    typedef raster_data<orbit_data, double, 2>      dataset_type;
    typedef grid_type::bounds_type                  bounds_type;
    typedef grid_type::vec_type                     vec_type;
    typedef grid_type::ivec_type                    ivec_type;
    
    const vec_type& step        = mesh.spacing();
    const bounds_type& bounds   = mesh.bounds();
    const ivec_type& resolution = mesh.dimensions();
    const map_metric& metric    = params.metric;
    
    std::cerr << "mesh spacing = " << step << '\n'
              << "mesh bounds = " << bounds << '\n'
              << "resolution = " << resolution << '\n';
              
              
    int niter = params.nb_iterations;
    int npoints = resolution[0]*resolution[1];
    
    size_t nbthreads = 1;
#if _OPENMP
    nbthreads = omp_get_max_threads();
#endif
    std::vector<map_analysis_param::tagged_orbit_type>  cached_orbits[nbthreads];
    
    #pragma omp parallel
    {
        #pragma omp for schedule(dynamic,1)
        for(int n = 0 ; n < npoints ; ++n) {
            int thread_id = 0;
#if _OPENMP
            thread_id = omp_get_thread_num();
#endif
            int j = n / resolution[0];
            int i = n % resolution[0];
            nvis::ivec2 id(i,j);
            
            if (thread_id == 0) {
                std::ostringstream os;
                os << "\rcompleted " << n << " / " << npoints << " ("
                   << 100*n/npoints << "%)      \r" << std::flush;
                std::cout << os.str();
            }
            
            nvis::vec2 x0 = bounds.min() + step * nvis::vec2(i,j);
            
            MAP* amap = pmap.clone();
            std::vector<nvis::vec2> tmp;
            try {
                amap->map(x0, tmp, niter);
            } catch(...) {
                continue;
            }
            
            if (tmp.size() < niter) {
                continue;
            }
            
            std::vector<nvis::vec2>& steps = dataset(id).steps;
            steps.push_back(x0);
            std::copy(tmp.begin(), tmp.end(), std::back_inserter(steps));
            if (compute_sf) {
                dataset(id).set();
            }
            if (params.record) {
                // display
                cached_orbits[thread_id].push_back(map_analysis_param::tagged_orbit_type());
                cached_orbits[thread_id].back().second = dataset(id).p;
                cached_orbits[thread_id].back().first = steps;
                // end display
            }
        }
    }
    std::cout << '\n';
    if (params.record) {
        // gather
        for (int i=0 ; i<nbthreads ; ++i) {
            for (int j=0 ; j<cached_orbits[i].size() ; ++j) {
                params.orbits.push_back(cached_orbits[i][j]);
            }
        }
    }
}

template<typename MAP>
double xavier::adaptive_rotation_angle(const nvis::vec2& x0, const nvis::vec2& v0,
                                       const nvis::vec2& x1, const nvis::vec2& v1,
                                       const MAP& pmap, unsigned int period,
                                       double dx, map_analysis_param& params)
{

    if (nvis::norm(v0) < 1.0e-6 || nvis::norm(v1) < 1.0e-6) {
        index_step_size_underflow e("unable to track continous vector rotation");
        e.where = 0.5*(x0 + x1);
        throw e;
    }
    
    double theta = xavier::signed_angle(v0, v1);
    
    bool verbose = params.verbose;
    
    // if (params.verbose) {
    //     std::ostringstream os;
    //     os << "adaptive angle: angle = " << theta << " between " << v0 << " at " << x0
    //        << " and " << v1 << " at " << x1 << " (d = " << nvis::norm(x0-x1) << ", dx = " << dx << ")" << std::endl;
    //     std::cerr << os.str();
    // }
    double tiny = std::max(1.0e-15, fabs(dx));
    
    if (fabs(theta) <= MAX_ANGLE_VARIATION) {
        return theta;
    } else if (params.metric.distance(x0, x1) <= tiny) {
        if (verbose) {
            std::cerr << "angle " << theta << " between " << v0 << " at " << x0 << " and " << v1 << " at " << x1
                      << " remains too large. failed\n";
        }
        
        // only bad things will come out of trying to proceed in such cases...
        index_step_size_underflow e("unable to track continous vector rotation");
        e.where = 0.5*(x0 + x1);
        throw e;
    }
    
    if (verbose) {
        std::cerr << "Distance between samples = " << params.metric.distance(x0, x1) << "\n";
    }
    
#ifndef SECANT_METHOD
    nvis::vec2 x = 0.5 * (x0 + x1);
#else
    double u = secant(v0, v1);
    if (std::isnan(u) || std::isinf(u)) {
        std::cerr << "secant method returned NaN" << std::endl;
    }
    nvis::vec2 x = (1-u)*x0 + u*x1;
#endif
    nvis::vec2 v = rhs(x, pmap, period, params);
    
    if (verbose) {
        std::cerr << x0[0] << ", " << x0[1] << ", " << v0[0] << ", " << v0[1] << '\n';
        std::cerr << x1[0] << ", " << x1[1] << ", " << v1[0] << ", " << v1[1] << '\n';
    }
    
    return adaptive_rotation_angle(x0, v0, x, v, pmap, period, dx, params) +
           adaptive_rotation_angle(x, v, x1, v1, pmap, period, dx, params);
}

template<typename MAP>
double xavier::robust_rotation_angle(const nvis::vec2& x0, const nvis::vec2& v0,
                                     const nvis::vec2& x1, const nvis::vec2& v1,
                                     const MAP& pmap, unsigned int period,
                                     double dx, map_analysis_param& params)
{
    double theta = xavier::signed_angle(v0, v1);
    
    if (fabs(theta) <= MAX_ANGLE_VARIATION) {
        return theta;
    }
    
    nvis::vec2 x = 0.5 * (x0 + x1);
    nvis::vec2 v = rhs(x, pmap, period, params);
    
    // we stop when one of two things happens:
    // 1) the norm of the vector at the halfway point is less than a
    //    prescribed epsilon used in the fixed point detection
    // 2) the rotation from left and right is not monotonic
    double norm = nvis::norm(v);
    if (norm < dx) {
        index_step_size_underflow e("unable to track continous vector rotation");
        e.where = x;
        throw e;
    }
    
    return adaptive_rotation_angle(x0, v0, x, v, pmap, period, dx, params) +
           adaptive_rotation_angle(x, v, x1, v1, pmap, period, dx, params);
}

template<typename MAP>
int xavier::__poincare_index(const grid<double, 2>& mesh,
                             const nvis::ivec2& cell_id,
                             const raster_data<orbit_data, double, 2>& dataset,
                             const MAP& pmap, int period, double lmin,
                             bool scrutinize,
                             std::vector<std::pair<int, nvis::vec2> >& sing_edges,
                             map_analysis_param& params)
{

    typedef grid<double, 2>                         grid_type;
    typedef raster_data<orbit_data, double, 2>      dataset_type;
    typedef grid_type::bounds_type                  bounds_type;
    typedef grid_type::vec_type                     vec_type;
    typedef grid_type::ivec_type                    ivec_type;
    
    const int ids[][2] = { {0,0}, {1,0}, {1,1}, {0,1} };
    
    nvis::vec2 v[5];
    nvis::vec2 x[5];
    
    sing_edges.clear();
    
    // corner values
    for (int i=0 ; i<4 ; ++i) {
        nvis::ivec2 delta(ids[i][0], ids[i][1]);
        x[i] = mesh(cell_id + delta);
        v[i] = dataset(cell_id + delta).vector(period);
    }
    // close the loop for convenience
    x[4] = x[0];
    v[4] = v[0];
    
    double dtheta[4];
    for (int i=0 ; i<4 ; ++i) {
        dtheta [i] = signed_angle(v[i], v[i+1]);
    }
    
    // compute Poincare index of the map for current period
    double theta = 0;
    for (int edge=0 ; edge<4 ; ++edge) {
        if (dtheta[edge] < MAX_ANGLE_VARIATION) {
            theta += dtheta[edge];
            continue;
        }
        
        try {
            theta += adaptive_rotation_angle(x[edge], v[edge], x[edge+1], v[edge+1],
                                             pmap, period, lmin, params);
        } catch (index_step_size_underflow& e) {
            if (scrutinize) {
                sing_edges.push_back(std::pair<int, nvis::vec2>(edge, e.where));
                if (sing_edges.size() == 2) {
                    throw rational_surface_suspected("none");
                }
                continue; // no need to further look at this edge. it is broken anyway
            } else {
                return 0;    // refined tracking failed, too.
            }
        } catch(...) {
            return 0;
        }
    }
    if (sing_edges.size()) {
        throw higher_precision_required("none");
    }
    
    // theta is a valid angle. turn it into an index
    long int idx = lrint(0.5*theta/M_PI);
    if (params.verbose) {
        std::cerr << "integral index = " << idx << '\n';
    }
    return idx;
}

template<typename MAP>
int xavier::__fast_poincare_index(const grid<double, 2>& mesh,
                                  const nvis::ivec2& cell_id,
                                  const raster_data<orbit_data, double, 2>& dataset,
                                  const MAP& pmap, int period, double lmin,
                                  bool scrutinize,
                                  std::vector<std::pair<int, nvis::vec2> >& sing_edges,
                                  map_analysis_param& params)
{

    typedef grid<double, 2>                         grid_type;
    typedef std::pair<nvis::vec2, int>              data_type;
    typedef raster_data<orbit_data, double, 2>      dataset_type;
    typedef grid_type::bounds_type                  bounds_type;
    typedef grid_type::vec_type                     vec_type;
    typedef grid_type::ivec_type                    ivec_type;
    
    const int ids[][2] = { {0,0}, {1,0}, {1,1}, {0,1} };
    
    nvis::vec2 v[5];
    nvis::vec2 x[5];
    
    sing_edges.clear();
    
    for (int i=0 ; i<4 ; ++i) {
        nvis::ivec2 delta(ids[i][0], ids[i][1]);
        x[i] = mesh(cell_id + delta);
        v[i] = dataset(cell_id + delta).vector(period);
    }
    x[5] = x[0];
    v[5] = v[0];
    
    // compute Poincare index of the map for current period
    double theta = 0;
    for (int edge=0 ; edge<4 ; ++edge) {
        try {
            theta += robust_rotation_angle(x[edge], v[edge], x[edge+1], v[edge+1],
                                           pmap, period, lmin, params);
        } catch (index_step_size_underflow& e) {
            if (scrutinize) {
                sing_edges.push_back(std::pair<int, nvis::vec2>(edge, e.where));
                if (sing_edges.size() == 2) {
                    throw rational_surface_suspected("none");
                }
                continue; // no need to further look at this edge. it is broken anyway
            } else {
                return 0;    // refined tracking failed, too.
            }
        } catch(...) {
            return 0;
        }
    }
    if (sing_edges.size()) {
        throw higher_precision_required("none");
    }
    
    // theta is a valid angle. turn it into an index
    long int idx = lrint(0.5*theta/M_PI);
    if (params.verbose) {
        std::cerr << "integral index = " << idx << '\n';
    }
    return idx;
}

template<typename MAP>
int xavier::poincare_index(const grid<double, 2>& mesh, const nvis::ivec2& cell_id,
                           const raster_data<orbit_data, double, 2>& dataset,
                           const MAP& pmap, int period, double lmin,
                           map_analysis_param& params)
{

    const int ids[][2] = { {0,0}, {1,0}, {1,1}, {0,1}, {0,0} };
    
    std::vector<std::pair<int, nvis::vec2> > sing_edges;
    
    try {
        // return __fast_poincare_index(mesh, cell_id, dataset, pmap, period, lmin,
        return __poincare_index(mesh, cell_id, dataset, pmap, period, lmin,
                                true, sing_edges, params);
    } catch (rational_surface_suspected& e) {
        if (params.verbose) {
            std::cerr << "rational surface suspected in cell " << cell_id << std::endl;
            std::cerr << "first increase the precision\n";
        }
        MAP* amap = pmap.clone();
        amap->precision(0.1*pmap.precision());
        sing_edges.clear();
        try {
            // return __fast_poincare_index(mesh, cell_id, dataset, *amap, period, lmin/10.,
            return __poincare_index(mesh, cell_id, dataset, *amap, period, lmin/10.,
                                    false, sing_edges, params);
        } catch(rational_surface_suspected& f) {
            if (params.verbose) {
                std::cerr << "rational surface suspicion remains at higher precision\n";
            }
        } catch (...) {
            if (params.verbose) {
                std::cerr << "Failed to track the orientation of the vector field. giving up.\n";
            }
            return 0;
        }
        std::vector<nvis::vec2> steps;
        const nvis::vec2 x0 = sing_edges[0].second;
        iterate<MAP>(x0, steps, pmap, params.nb_iterations);
        int period0 = best_period(steps, params.max_period, params.metric);
        double sf0 = safety_factor(steps, params.metric);
        boost::rational<int> q0 = rational_approx(sf0, params.max_period);
        
        steps.clear();
        const nvis::vec2& x1 = sing_edges[1].second;
        iterate<MAP>(x1, steps, pmap, params.nb_iterations);
        int period1 = best_period(steps, params.max_period, params.metric);
        double sf1 = safety_factor(steps, params.metric);
        boost::rational<int> q1 = rational_approx(sf1, params.max_period);
        
        if (params.verbose) {
            std::cerr << "x0 = " << x0 << ", period = " << period0
                      << ", sf = " << sf0 << ", corresponding q = " << q0 << std::endl;
            std::cerr << "x1 = " << x1 << ", period = " << period1
                      << ", sf = " << sf1 << ", corresponding q = " << q1 << std::endl;
        }
        
        if (q0 == q1 && // same rational safety factor approximation
                q0.numerator() == period && // numerator compatible with current period
                period == period0 && period == period1) { // both orbits have relevant best period
            if (params.verbose) {
                std::cerr << "a rational surface was detected for period " << period0 << std::endl;
            }
            
            int edge0 = sing_edges[0].first;
            int edge1 = sing_edges[1].first;
            nvis::ivec2 d00(ids[edge0  ][0], ids[edge0  ][1]);
            nvis::ivec2 d01(ids[edge0+1][0], ids[edge0+1][1]);
            nvis::ivec2 d10(ids[edge1  ][0], ids[edge1  ][1]);
            nvis::ivec2 d11(ids[edge1+1][0], ids[edge1+1][1]);
            
            rational_surface_found rsf("none");
            rational_surface_found::edge_type e0(cell_id + d00, cell_id + d01);
            rational_surface_found::edge_type e1(cell_id + d10, cell_id + d11);
            rsf.xings[0] = rational_surface_found::edge_point(e0, sing_edges[0].second);
            rsf.xings[1] = rational_surface_found::edge_point(e1, sing_edges[1].second);
            rsf.sf = q0;
            
            std::ostringstream os;
            os << "found a rational surface of period " << rsf.sf << std::endl;
            std::cerr << os.str();
            
            throw rsf;
        } else {
            if (params.verbose) {
                std::cerr << "conflicting values, giving up\n";
            }
        }
    } catch (higher_precision_required& e) {
        if (params.verbose) {
            std::cerr << "higher precision was requested\n";
        }
        MAP* amap = pmap.clone();
        amap->precision(0.1*pmap.precision());
        sing_edges.clear();
        try {
            // return __fast_poincare_index(mesh, cell_id, dataset, *amap, period, lmin/10.,
            return __poincare_index(mesh, cell_id, dataset, *amap, period, lmin/10.,
                                    false, sing_edges, params);
        } catch(...) {
            std::ostringstream os;
            os << "failed to track continuous vector orientation" << std::endl;
            std::cerr << os.str();
            
            return 0;
        }
    } catch(...) {
        return 0;
    }
}


template<typename MAP>
void
xavier::process_cell(std::vector<std::pair<int, int> >& pids,
                     const raster_data<orbit_data, double, 2>& dataset,
                     const grid<double, 2>& mesh,
                     const nvis::ivec2& cell_id, const MAP& pmap, double dx,
                     map_analysis_param& params)
{
    typedef grid<double, 2>                         grid_type;
    typedef raster_data<orbit_data, double, 2>      dataset_type;
    typedef grid_type::bounds_type                  bounds_type;
    typedef grid_type::vec_type                     vec_type;
    typedef grid_type::ivec_type                    ivec_type;
    
    const vec_type& step        = mesh.spacing();
    const bounds_type& bounds   = mesh.bounds();
    const ivec_type& resolution = mesh.dimensions();
    const map_metric& metric    = params.metric;
    
    const std::pair<int, int> default_result(0, 0);
    
    const int ids[][2] = { {0,0}, {1,0}, {1,1}, {0,1} };
    
    pids.clear();
    
    std::vector<int> periods;
    period_range(periods, dataset, cell_id, params.valid_rationals);
    std::set<int> available_periods(periods.begin(), periods.end());
    
    std::cerr << "available periods: " << std::endl;
    std::copy(available_periods.begin(), available_periods.end(), std::ostream_iterator<int>(std::cerr, ", "));
    std::cerr << std::endl;
    
    for (std::set<int>::const_iterator it=available_periods.begin();
            it!=available_periods.end() ; ++it) {
        if (params.verbose) {
            std::cerr << "considering period " << *it << std::endl;
        }
        
        // compute Poincare index of the map for current period
        int p = *it;
        if (p == 0) {
            return;
        }
        int pidx = poincare_index(mesh, cell_id, dataset, pmap, p, dx, params);
        if (pidx != 0) {
            pids.push_back(std::pair<int, int>(pidx, p));
        }
    }
}


#endif