#include <iostream>
#include <vector>
#include <complex>
#include <sstream>
#include <math.h>
#include <sstream>
#include <fstream>

// boost
#include <boost/format.hpp>
#include <boost/limits.hpp>
#include <boost/rational.hpp>

// nvis
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>
#include <math/bounding_box.hpp>

// map analysis
#include <math/math.hpp>
#include <poincare/metric.hpp>
#include <poincare/newton.hpp>
#include <maps/period.hpp>
#include <maps/regular_mesh.hpp>
#include <poincare/topology.hpp>
#include <poincare/fixpoints.hpp>
#include <poincare/macros.hpp>
#include <maps/sampling.hpp>
#include <maps/temporary.hpp>
#include <maps/tokatopo.hpp>

// christoph
#include <tokamak/tokamak_nimrod_parametric.hpp>
#include <tokamak/poincare_map.hpp>

// teem
#include <teem/nrrd.h>

// OpenMP
#if _OPENMP
#include <omp.h>
#endif

// meshing
#include <maps/basic_definitions.hpp>
#include <maps/triangulation.hpp>
#include <maps/adaptive_triangulation.hpp>
#include <maps/quality_control.hpp>
#include <maps/IO.hpp>

#include <maps/tokatopo.hpp>

using namespace tokatopo;

namespace {
// local functors for adaptive triangulation
struct error_measure_linear {
    typedef xavier::point_data      data_type;
    
    error_measure_linear(unsigned int period, double threshold, const xavier::map_metric& metric)
        : __period(period), __eps(threshold), __metric(metric) {}
        
    double operator()(const data_type& d,
                      const nvis::vec3& beta,
                      const data_type data[3]) const {
        nvis::vec2 ref = vector_value(d, __period, __metric);
        nvis::vec2 v[3];
        for (int i = 0 ; i < 3 ; ++i) {
            v[i] = vector_value(data[i], __period, __metric);
        }
        
        nvis::vec2 approx = beta[0] * v[0] + beta[1] * v[1] + beta[2] * v[2];
        return nvis::norm(approx - ref) / nvis::norm(ref) - __eps;
    }
    
    unsigned int __period;
    double __eps;
    xavier::map_metric __metric;
};

struct error_measure_norm {
    typedef xavier::point_data      data_type;
    
    error_measure_norm(unsigned int period, const xavier::map_metric& metric)
        : __period(period), __metric(metric) {}
        
    double operator()(const data_type& d,
                      const nvis::vec3& beta,
                      const data_type data[3]) const {
        nvis::vec2 ref = vector_value(d, __period, __metric);
        double dref = nvis::norm(ref);
        double dist[3];
        for (int i = 0 ; i < 3 ; ++i) {
            nvis::vec2 v = vector_value(data[i], __period, __metric);
            dist[i] = nvis::norm(v);
        }
        
        if (dref < *std::min_element(&dist[0], &dist[3])) {
            return 1. / dref;
        } else {
            return -1;
        }
    }
    
    unsigned int __period;
    xavier::map_metric __metric;
};

template<typename T>
inline T sign(const T& t)
{
    return (t > 0 ? 1 : -1);
}

inline double angle(const nvis::vec2& v0, const nvis::vec2& v1)
{
    if (!nvis::norm(v0) || !nvis::norm(v1)) {
        return 0;
    }
    nvis::vec2 w0 = v0 / nvis::norm(v0);
    nvis::vec2 w1 = v1 / nvis::norm(v1);
    double theta = acos(nvis::inner(w0, w1));
    double det = w0[0] * w1[1] - w0[1] * w1[0];
    return sign(det)*theta;
}

template<typename Mesh, typename Map>
struct poincare_index {

    typedef Mesh                                mesh_type;
    typedef Map                                 integrator_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::data_type       data_type;
    typedef tokatopo::interval_type             interval_type;
    typedef tokatopo::rational_type             rational_type;
    typedef std::pair<point_type, point_type>   uncertain_vec_type;
    
    poincare_index() : _log(false) {}
    
    inline int direct(const xavier::point_data d[3], unsigned int period,
                      const xavier::map_metric& metric) {
        double dtheta = 0;
        for (int i = 0 ; i < 3 ; ++i) {
            dtheta += angle(xavier::vector_value(d[i], period, metric),
                            xavier::vector_value(d[(i+1)%3], period, metric));
        }
        return (int)round(0.5*dtheta / M_PI);
    }
    
    uncertain_vec_type
    evaluate_step(const nvis::vec2& x, const integrator_type& map,
                  unsigned int period, const xavier::map_metric& metric) {
        map2d::value_type result;
        try {
            result = map.map_complete(x, period);
        } catch (...) {
            return std::make_pair(nvis::vec2(0, 0), nvis::vec2(0, 0));
        }
        
        nvis::vec2 step = metric.displacement(x, result.x);
        return uncertain_vec_type(step, result.err);
    }
    
    double smooth_rotation_angle(const nvis::vec2& x0, const nvis::vec2& x1,
                                 const integrator_type& map, unsigned int period,
                                 double dtheta, double dx,
                                 const xavier::map_metric& metric) {
        uncertain_vec_type uv0, uv1;
        uv0 = evaluate_step(x0, map, period, metric);
        uv1 = evaluate_step(x1, map, period, metric);
        const point_type& v0 = uv0.first;
        const point_type& v1 = uv1.first;
        double theta = angle(v0, v1);
        
        if (fabs(theta) < dtheta || nvis::norm(x1 - x0) <= dx) {
            return theta;
        }
        
        nvis::vec2 x = 0.5 * (x0 + x1);
        return smooth_rotation_angle(x0, x, map, period, dtheta, dx, metric) +
               smooth_rotation_angle(x, x1, map, period, dtheta, dx, metric);
    }
    
    int safe(const mesh_type& mesh, unsigned int tri, const Map& map,
             unsigned int period, double dtheta, double dx,
             const xavier::map_metric& metric) {
             
        point_type  p[3];
        data_type   d[3];
        mesh.get_triangle_info(p, d, tri);
        
        const integrator_type* pmap = map.clone();
        double delta_theta = 0;
        for (int i = 0 ; i < 3 ; ++i) {
            nvis::vec2 p0, p1;
            p0 = p[i];
            p1 = p[(i+1)%3];
            delta_theta += smooth_rotation_angle(p0, p1, *pmap, period, dtheta, dx, metric);
            if (_log) {
                std::cerr << "delta_theta along edge #" << i << " is " << delta_theta << '\n';
            }
        }
        
        return (int)round(0.5*delta_theta / M_PI);
    }
    
    bool _log;
};

struct angular_variation_priority {
    typedef nvis::vec2                  point_type;
    typedef xavier::point_data          data_type;
    typedef tokatopo::interval_type     interval_type;
    typedef tokatopo::rational_type     rational_type;
    
    angular_variation_priority(unsigned int period, double eps,
                               const xavier::map_metric& metric, double dq)
        : _p(period), _eps(eps), _dq(dq), _metric(metric) {
        for (int i = 1 ; i <= period ; ++i) {
            rational_type r(period, i);
            if (r.numerator() != period) {
                continue;
            }
            double q = xavier::value(r);
            _valid_q.push_back(interval_type(q - dq, q + dq));
        }
    }
    
    double operator()(const point_type p[3], const data_type d[3]) const {
        double dtheta = 0;
        double minq, maxq;
        minq = maxq = d[0].period();
        for (int i = 1 ; i < 3 ; ++i) {
            if (d[i].period() > maxq) {
                maxq = d[i].period();
            } else if (d[i].period() < minq) {
                minq = d[i].period();
            }
        }
        interval_type span(minq - dq, maxq + dq);
        
        bool valid = false;
        for (int i = 0 ; i < _valid_q.size() && !valid ; ++i) {
            valid = !(xavier::intersect(span, _valid_q[i]).empty());
        }
        if (!valid) {
            return -1;
        }
        
        for (int i = 0 ; i < 3 ; ++i) {
            dtheta += fabs(angle(xavier::vector_value(d[i], _p, _metric),
                                 xavier::vector_value(d[(i+1)%3], _p, _metric)));
        }
        return xavier::area(p)*(dtheta - _eps);
    }
    
    double _eps, _dq;
    unsigned int _p;
    std::vector<interval_type> _valid_q;
    const xavier::map_metric& _metric;
};

// hack
nvis::bbox2 saddle_box(nvis::vec2(84, 31.5), nvis::vec2(86, 33.5));

struct conditional_angular_variation_priority {
    typedef nvis::vec2                  point_type;
    typedef xavier::point_data          data_type;
    typedef xavier::interval<double>    interval_type;
    
    conditional_angular_variation_priority(unsigned int period,
                                           const std::vector<double>& qs,
                                           double qeps, double eps,
                                           const xavier::map_metric& metric)
        : _p(period), _ints(qs.size()), _eps(eps), _metric(metric) {
        for (int i = 0 ; i < qs.size() ; ++i) {
            _ints[i] = interval_type(qs[i] - qeps, qs[i] + qeps);
        }
    }
    
    double operator()(const point_type p[3], const data_type d[3]) const {
    
        // bool display_stuff = (_p == 1) && (saddle_box.inside(p[0]) ||
        //                                    saddle_box.inside(p[1]) ||
        //                                    saddle_box.inside(p[2]));
        
        // double min, max;
        // min = max = d[0].period();
        // for (int i = 1 ; i < 3 ; ++i) {
        //  double q = d[i].period();
        //  if (q < min) min = q;
        //  else if (q > max) max = q;
        // }
        // interval_type cell_int(min, max);
        // bool relevant = false;
        // for (int i = 0 ; i < _ints.size() && !relevant ; ++i) {
        //  relevant = !(xavier::intersect(_ints[i], cell_int).empty());
        // }
        
        // if (!relevant) {
        //  if (display_stuff) {
        //      std::cerr << "triangle: " << p[0] << ", " << p[1] << ", " << p[2]
        //                << " has been deemed irrelevant because its associated periods are "
        //                << d[0].period() << ", " << d[1].period() << ", and " << d[2].period()
        //                << '\n';
        //  }
        //  return 0;
        // }
        
        double dtheta = 0;
        double maxerror = 0;
        for (int i = 0 ; i < 3 ; ++i) {
            std::pair<nvis::vec2, nvis::vec2> vecerr0 = xavier::vector_and_error_value(d[i], _p, _metric);
            std::pair<nvis::vec2, nvis::vec2> vecerr1 = xavier::vector_and_error_value(d[(i+1)%3], _p, _metric);
            dtheta += fabs(angle(vecerr0.first, vecerr1.first));
            double err = nvis::norm(vecerr0.second);
            if (err > maxerror) {
                maxerror = err;
            }
        }
        
        // // if triangle area is smaller than half of square of side maxerror, bail
        // if (2*xavier::area(p) < maxerror*maxerror) return 0;
        
        // if (display_stuff) {
        //  std::cerr << "triangle: " << p[0] << ", " << p[1] << ", " << p[2]
        //            << " has a priority value of " << dtheta - _eps << '\n';
        // }
        return dtheta - _eps;
    }
    
    unsigned int _p;
    std::vector<interval_type> _ints;
    double _eps;
    const xavier::map_metric& _metric;
};

struct triangle_filter {
    typedef nvis::vec2                  point_type;
    typedef xavier::point_data          data_type;
    typedef tokatopo::interval_type     interval_type;
    typedef tokatopo::rational_type     rational_type;
    
    triangle_filter(unsigned int period, double eps,
                    const xavier::map_metric& metric, double dq)
        : _p(period), _eps(eps), _dq(dq), _metric(metric) {
        for (int i = 1 ; i <= period ; ++i) {
            rational_type r(period, i);
            if (r.numerator() != period) {
                continue;
            }
            double q = xavier::value(r);
            _valid_q.push_back(interval_type(q - dq, q + dq));
        }
    }
    
    bool valid(const point_type p[3], const data_type d[3]) const {
        double minq, maxq;
        minq = maxq = d[0].period();
        for (int i = 1 ; i < 3 ; ++i) {
            if (d[i].period() > maxq) {
                maxq = d[i].period();
            } else if (d[i].period() < minq) {
                minq = d[i].period();
            }
        }
        interval_type span(minq - dq, maxq + dq);
        
        bool valid = false;
        for (int i = 0 ; i < _valid_q.size() && !valid ; ++i) {
            valid = !(xavier::intersect(span, _valid_q[i]).empty());
        }
        return valid;
    }
    
    unsigned int _p;
    double _eps, _dq;
    const xavier::map_metric& _metric;
    std::vector<interval_type> _valid_q;
};

}

double approx_error;
void init(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "f",      "file",                 airTypeString,  1, 1, &file,                NULL,       "input hdf5 file name");
    hestOptAdd(&hopt, "t",      "time",                 airTypeString,  1, 1, &ts,                  NULL,       "time step string");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &outs,                NULL,       "output nrrd file name");
    hestOptAdd(&hopt, "eps",    "eps",                  airTypeDouble,  0, 1, &h,                   "1.0e-6",   "integration precision");
    hestOptAdd(&hopt, "n",      "# samples",            airTypeInt,     0, 1, &n1,                  "50000",    "number of 1D samples");
    hestOptAdd(&hopt, "m",      "# rotations",          airTypeInt,     0, 1, &m,                   "50",       "number of rotations");
    hestOptAdd(&hopt, "dq",     "safety factor step",   airTypeDouble,  0, 1, &dq,                  "0.01",     "maximum discrepancy between consecutive safety factors");
    hestOptAdd(&hopt, "mr",     "max aspect ratio",     airTypeDouble,  0, 1, &max_ratio,           "2.",       "maximum triangle aspect ratio");
    hestOptAdd(&hopt, "mt",     "max # triangles",      airTypeInt,     0, 1, &max_nb_triangles,    "500000",   "max number of triangles in adaptive sampling");
    hestOptAdd(&hopt, "ma",     "min triangle area",    airTypeDouble,  1, 1, &min_area,            NULL,       "min triangle area in adaptive sampling");
    hestOptAdd(&hopt, "Ma",     "max triangle area",    airTypeDouble,  1, 1, &max_area,            NULL,       "max triangle area in adaptive sampling");
    hestOptAdd(&hopt, "mp",     "max period",           airTypeInt,     0, 1, &max_period,          "25",       "max considered period in fixed point search");
    hestOptAdd(&hopt, "err",    "max approx error",     airTypeDouble,  0, 1, &approx_error,        "0.1",      "max relative error of linear approximation");
    // hestOptAdd(&hopt, "p",       "period",               airTypeInt,     1, 1, &period,              NULL,       "targeted period");
    
    hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
                   me, "Adaptively sample phase portrait along orbits to achieve accurate\npiecewise linear interpolation on a per-period basis",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

std::vector<rational_type> ref_values;

template<typename T>
inline bool has_zero(int period, const T& mesh, int cellid)
{
    typedef T                                               triangulation_type;
    typedef typename triangulation_type::triangle_type      triangle_type;
    typedef typename triangulation_type::index_type         index_type;
    typedef typename triangulation_type::point_type         point_type;
    typedef typename triangulation_type::data_type          data_type;
    
    const triangle_type& tri = mesh.get_triangle_vertices(cellid);
    nvis::vec2 v[3];
    for (int i = 0 ; i < 3 ; ++i) {
        const data_type& d = mesh.get_data(tri[i]);
        v[i] = xavier::vector_value(d, period, xavier::__default_metric);
    }
    
    nvis::vec2 col[2] = { v[0] - v[2], v[1] - v[2] };
    double denom, beta[3];
    denom = col[0][0] * col[1][1] - col[1][0] * col[0][1];
    if (denom == 0) {
        return false;
    }
    beta[0] = (-v[2][0] * col[1][1] + v[2][1] * col[1][0]) / denom;
    beta[1] = (-v[2][1] * col[0][0] + v[2][0] * col[0][1]) / denom;
    beta[2] = 1. - beta[0] - beta[1];
    return (beta[0] >= 0 && beta[1] >= 0 && beta[2] >= 0);
}

inline nvis::vec3 zero_coord(const nvis::vec2 v[3])
{
    nvis::vec2 col[2] = { v[0] - v[2], v[1] - v[2] };
    double denom;
    nvis::vec3 beta;
    denom = col[0][0] * col[1][1] - col[1][0] * col[0][1];
    if (denom == 0) {
        return nvis::vec3(-1, 0, 0);
    }
    beta[0] = (-v[2][0] * col[1][1] + v[2][1] * col[1][0]) / denom;
    beta[1] = (-v[2][1] * col[0][0] + v[2][0] * col[0][1]) / denom;
    beta[2] = 1. - beta[0] - beta[1];
    return beta;
}

inline double min_q_dist(double q, const std::vector<rational_type>& r)
{
    double _min = std::numeric_limits<double>::max();
    for (int i = 0 ; i < r.size() ; ++i) {
        _min = std::min(fabs(q - xavier::value(r[i])), _min);
    }
    
    return _min;
}

template<typename T, typename F>
struct Lt_point_norm {
    typedef std::pair<T, F> pair_type;
    
    bool operator()(const pair_type& p0, const pair_type& p1) {
        return p0.second < p1.second;
    }
    
};

struct p_step_func {
    typedef xavier::interval<double>    interval_type;
    
    p_step_func(unsigned int period, double eps)
        : __period(period) {
        for (int i = 1 ; i <= period ; ++i) {
            double q = (double)period / (double)i;
            qs.push_back(interval_type(q - eps, q + eps));
        }
    }
    
    int order() const {
        return 1;
    }
    
    std::string name() const {
        std::ostringstream os;
        os << __period << "-vector";
        return os.str();
    }
    
    bool is_valid(const xavier::point_data& d) const {
        bool is_inside = false;
        for (int i = 0 ; i < qs.size() && !is_inside ; ++i) {
            is_inside = qs[i].inside(d.period());
        }
        if (!is_inside) {
            return false;
        }
        
        nvis::vec2 v = xavier::vector_value(d, __period, xavier::__default_metric);
        return (v[0] != std::numeric_limits<double>::max());
    }
    
    std::string value_string(const xavier::point_data& d) const {
        bool is_inside = false;
        for (int i = 0 ; i < qs.size() && !is_inside ; ++i) {
            is_inside = qs[i].inside(d.period());
        }
        if (!is_inside) {
            return "0 0 0";
        }
        
        nvis::vec2 v = xavier::vector_value(d, __period, xavier::__default_metric);
        std::ostringstream os;
        os << v[0] << " " << v[1] << " 0";
        return os.str();
    }
    
    unsigned int __period;
    std::vector<interval_type> qs;
};

inline nvis::vec2 sanitize(const nvis::vec2& v)
{
    nvis::vec2 w(v);
    for (int i = 0 ; i < 2 ; ++i) {
        if (std::isnan(w[i]) || std::isinf(w[i])) {
            w[i] = 0.;
        }
    }
    return w;
}

struct error_func {
    error_func() {}
    
    int order() const {
        return 1;
    }
    
    std::string name() const {
        std::ostringstream os;
        os << "error_function";
        return os.str();
    }
    
    bool is_valid(const xavier::point_data& d) const {
        const nvis::vec2& err = d.error();
        return (err[0] != std::numeric_limits<double>::max());
    }
    
    std::string value_string(const xavier::point_data& d) const {
        nvis::vec2 e = sanitize(d.error());
        std::ostringstream os;
        if (e[0] == std::numeric_limits<double>::max() || nvis::norm(e) > 10.) {
            os << "0 0 0";
        } else {
            os << e[0] << " " << e[1] << " 0";
        }
        return os.str();
    }
};

void filter_orbits(std::vector<unsigned int>& ids,
                   const std::vector<xavier::orbit>& orbits,
                   double& min, double& max)
{
    ids.clear();
    if (!orbits.size()) {
        return;
    }
    double minq = min, maxq = max;
    min = std::numeric_limits<double>::max();
    max = std::numeric_limits<double>::min();
    for (int i = 0 ; i < orbits.size() ; ++i) {
        double q = orbits[i].period();
        if (q >= minq && q <= maxq) {
            ids.push_back(i);
            if (q < min) {
                min = q;
            }
            if (q > max) {
                max = q;
            }
        }
    }
}

inline void orbits_to_points(std::vector<nvis::vec2>& points,
                             std::vector<xavier::point_data>& data,
                             const std::vector<xavier::orbit>& orbits,
                             const std::vector<unsigned int>& valid_ids)
{
    points.clear();
    data.clear();
    for (unsigned int i = 0 ; i < valid_ids.size() ; ++i) {
        const xavier::orbit& orb = orbits[valid_ids[i]];
        for (unsigned int j = 0 ; j < orb.size() ; ++j) {
            points.push_back(orb[j]);
            data.push_back(xavier::point_data(valid_ids[i], j));
        }
    }
}

inline void orbits_to_points(std::vector<nvis::vec2>& points,
                             std::vector<xavier::point_data>& data,
                             const std::vector<xavier::orbit>& orbits)
{
    points.clear();
    data.clear();
    for (unsigned int i = 0 ; i < orbits.size() ; ++i) {
        const xavier::orbit& orb = orbits[i];
        for (unsigned int j = 0 ; j < orb.size() ; ++j) {
            points.push_back(orb[j]);
            data.push_back(xavier::point_data(i, j));
        }
    }
}

int main(int argc, char* argv[])
{
    typedef boost::rational<int>                rational_type;
    typedef xavier::point_data                  data_type;
    typedef std::pair<nvis::vec2, data_type>    data_point_type;
    
    init(argc, argv);
    
    tokamak_nimrod_parametric* field;
    field = new tokamak_nimrod_parametric(std::string(file), std::string(ts));
    field->periodic_coordinates(false);
    
#if 0
    {
        Nrrd* nrrd = field->to_nrrd();
        std::ostringstream os;
        os << outs << "-t=" << ts << ".nrrd";
        if (nrrdSave(os.str().c_str(), nrrd, NULL)) {
            std::cerr << biffGetDone(NRRD) << std::endl;
            exit(-1);
        }
    }
#endif
    
    bool per[2] = {true, false};
    xavier::map_metric metric(field->bounds(), per);
    
    xavier::__default_metric = metric;
    
    poincare_map pmap(field);
    pmap.precision(h);
    xavier::orbit_integrator<poincare_map> intg(pmap, m, metric);
    
    const nvis::bbox2& bounds = metric.bounds();
    nvis::vec2 diagonal = bounds.size();
    double tiny_length = 1.0e-6 * nvis::norm(diagonal);
    nvis::bbox2 inflated_bounds(bounds);
    inflated_bounds.min() -= 0.005 * diagonal;
    inflated_bounds.max() += 0.005 * diagonal;
    std::cerr << "inflated bounds are " << inflated_bounds << '\n';
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    srand48(time(0));
    std::vector< std::pair<point_type, xavier::point_data> > all_points;
    
    // reset central orbit repository
    xavier::__map_orbits.clear();
    
    // initialize base triangulation
    std::pair<nvis::vec2, xavier::point_data> boundary[4];
    std::vector<nvis::vec2> corner(4);
    corner[0] = inflated_bounds.min();
    corner[2] = inflated_bounds.max();
    corner[1] = nvis::vec2(corner[2][0], corner[0][1]);
    corner[3] = nvis::vec2(corner[0][0], corner[2][1]);
    size_t frame_orbit_id = xavier::__map_orbits.size();
    xavier::__map_orbits.push_back(xavier::orbit(corner, 0));
    for (int i = 0 ; i < 4 ; ++i) {
        boundary[i].first = corner[i];
        boundary[i].second = xavier::point_data(frame_orbit_id, i);
    }
    mesh_type base_mesh(boundary, xavier::point_locator());
    
    unsigned int res[2];
    double ratio = bounds.size()[1] / bounds.size()[0];
    res[0] = ceil(sqrt((double)n1 / ratio));
    res[1] = ceil(res[0] * ratio);
    xavier::experimental::sample_on_raster(xavier::__map_orbits,
                                           pmap, metric, res, m);
                                           
    // assign period to each orbit
    for (int i = 0 ; i < xavier::__map_orbits.size() ; ++i) {
        xavier::orbit& obt = xavier::__map_orbits[i];
        // double q = (xavier::period_x_periodic(obt.points(), metric)).first;
        double q = xavier::dist_based_x_period(obt.points(), metric);
        for (int j = 0 ; j < obt.size() ; ++j) {
            obt[j] = metric.modulo(obt[j]);
        }
        obt.period() = q;
        if (q < 0) {
            std::cerr << "Returned a negative period!\n";
            std::cerr << "points were: \n";
            for (int j = 0 ; j < obt.size() ; ++j) {
                std::cerr << j << ": " << obt[j] << '\n';
            }
            assert(false);
        }
    }
    
    {
        std::vector<nvis::vec2> points;
        std::vector<xavier::point_data> data;
        orbits_to_points(points, data, xavier::__map_orbits);
        base_mesh.insert_points(points, data);
        
        std::ostringstream os;
        os << outs << "-base_mesh.vtk";
        xavier::export_VTK(base_mesh, os.str(), "N / A", false);
    }
    base_mesh.set_tolerance(1.0e-7);
    
    xavier::map_debug::verbose_level = 1;
    
    for (int period = 1 ; period <= max_period ; ++period) {
        std::cerr << "\nprocessing period " << period << "...\n";
        
        std::cerr << "there are currently " << xavier::__map_orbits.size()
                  << " orbits on record and " << base_mesh.get_nb_triangles() << " points in base triangulation\n";
                  
        double h_p = (period == 1 ? h : 0.25 * h / (double)(period - 1));
        intg.precision(h_p);
        std::cerr << "integrator precision is currently " << 0.25*h / (double)period << '\n';
        
        mesh_type local_mesh(base_mesh); // copy mesh
        
        // determine what rationals should be considered
        std::set<rational_type> rationals;
        for (int den = 1 ; den <= period ; ++den) {
            rational_type q(period, den);
            if (q.numerator() != period) {
                continue;    // period and den are not mutually prime
            }
            rationals.insert(q);
        }
        
        std::vector<double> valid_qs;
        for (std::set<rational_type>::iterator it = rationals.begin() ; it != rationals.end() ; ++it) {
            valid_qs.push_back(xavier::value(*it));
        }
        
        std::cerr << "there are " << valid_qs.size() << " interesting rational periods:\n";
        for (int i = 0 ; i < valid_qs.size() ; ++i) {
            std::cerr << valid_qs[i] << "\t";
        }
        std::cerr << '\n';
        
        if (true) {
            std::cerr << "export subset of base mesh satisfying refinement criterion\n";
#if 0
            conditional_angular_variation_priority priority(period, valid_qs, dq, approx_error, metric);
#else
            angular_variation_priority priority(period, approx_error, metric, dq);
#endif
            std::vector<unsigned int> included_triangles;
            for (int i = 0 ; i < local_mesh.get_nb_triangles() ; ++i) {
                mesh_type::point_type p[3];
                mesh_type::data_type d[3];
                local_mesh.get_triangle_info(p, d, i);
                if (priority(p, d) > 0) {
                    included_triangles.push_back(i);
                }
            }
            
            std::ostringstream os_name, os_comment;
            os_name << outs << "-pre-mesh_for_p=" << period << ".vtk";
            p_step_func functor_all(period, dq);
            xavier::export_submesh_VTK(local_mesh, os_name.str(), "N/A", functor_all,
                                       included_triangles, true);
        }
        
        std::cerr << "building triangulation for period p = " << period << '\n';
        
#if USE_LINEAR_APPROXIMATION
        
        std::cerr << "refining triangulation around q = " << ref_values[i]
                  << " to reach piecewise linear precision\n";
        error_measure_norm priority(period, metric);
        bool ok = xavier::experimental::refine_map(local_mesh, intg, priority, max_nb_triangles);
        
#elif USE_CONDITIONAL_METHOD
        
        std::cerr << "refinement to achieve prescribed max angular variation (conditional)\n";
        conditional_angular_variation_priority priority(period, valid_qs, dq, approx_error, metric);
        bool ok = xavier::refine(local_mesh, intg, priority, max_nb_triangles);
        
#else
        
        std::cerr << "refinement to achieve prescribed max angular variation\n";
        angular_variation_priority priority(period, approx_error, metric, dq);
        bool ok = xavier::refine(local_mesh, intg, priority, max_nb_triangles);
        
#endif
        
        std::cerr << "after refinement mesh contains " << local_mesh.get_nb_triangles() << " triangles\n";
        
        {
            // export triangles matching period
            triangle_filter filter(period, approx_error, metric, dq);
            std::vector<unsigned int> included_triangles;
            for (int i = 0 ; i < local_mesh.get_nb_triangles() ; ++i) {
                mesh_type::point_type p[3];
                mesh_type::data_type d[3];
                local_mesh.get_triangle_info(p, d, i);
                if (filter.valid(p, d) > 0) {
                    included_triangles.push_back(i);
                }
            }
            
            std::ostringstream os_name;
            os_name << outs << "-" << period << "-specific_mesh.vtk";
            p_step_func functor_all(period, dq);
            xavier::export_submesh_VTK(local_mesh, os_name.str(), "N/A", functor_all, included_triangles, true);
            
            os_name.clear();
            os_name.str("");
            os_name << outs << "-" << period << "-specific_mesh-error.vtk";
            error_func functor_err;
            xavier::export_submesh_VTK(local_mesh, os_name.str(), "N/A", functor_err, included_triangles, true);
            
            // export subset of triangles matching period with nonzero index
            poincare_index<mesh_type, poincare_map> pindex;
            std::vector<unsigned int> singular_triangles;
            
            poincare_map* precise_map = pmap.clone();
            precise_map->precision(0.1*h_p);
            
            for (int i = 0 ; i < included_triangles.size() ; ++i) {
                unsigned int id = included_triangles[i];
                mesh_type::point_type p[3];
                mesh_type::data_type d[3];
                local_mesh.get_triangle_info(p, d, id);
                
                int index = pindex.direct(d, period, metric);
                if (index != 0) {
                    int real_index = pindex.safe(local_mesh, id, *precise_map, period,
                                                 0.5 * M_PI, tiny_length, metric);
                                                 
                    if (real_index != 0) {
                        singular_triangles.push_back(id);
                        std::cerr << "triangle " << p[0] << ", " << p[1] << ", " << p[2]
                                  << " contains a " << period << "-"
                                  << (real_index > 0 ? "center" : "saddle") << '\n';
                                  
                        nvis::vec2 v[3];
                        for (int j = 0 ; j < 3 ; ++j) {
                            v[j] = pindex.evaluate_step(p[j], *precise_map, period, metric).first;
                        }
                        nvis::vec3 b = zero_coord(v);
                        if (b[0] == -1) {
                            std::cerr << "invalid barycentric coordinates for zero location in singular cell!\n";
                        } else {
                            nvis::vec2 approx_zero = b[0] * p[0] + b[1] * p[1] + b[2] * p[2];
                            std::cerr << "approximate location of " << period << "-"
                                      << (real_index > 0 ? "center" : "saddle") << " is " << approx_zero << '\n';
                            double zero_norm = nvis::norm(pindex.evaluate_step(approx_zero, *precise_map, period, metric).first);
                            std::cerr << "corresponding norm of " << period << "-map is " << zero_norm << '\n';
                            if (zero_norm > 1) {
                                std::cerr << "nonsensical result detected: vertex vectors were: " << v[0] << ", "
                                          << v[1] << ", " << v[2] << '\n';
                                          
                                std::vector<nvis::vec2> pos;
                                try {
                                    precise_map->map(approx_zero, pos, 10*period);
                                } catch (...) {
                                    // std::cerr << "orbit_integrator: unable to integrate from " << x0 << std::endl;
                                }
                                if (pos.size() >= 9*period) {
                                    xavier::push_front(approx_zero, pos);
                                    double q = xavier::period_x_periodic(pos, metric).first;
                                    std::cerr << "rational period at this point is " << q << '\n';
                                } else {
                                    std::cerr << "something went wrong with the integration\n";
                                }
                            }
                        }
                        
                    }
                }
            }
            
            os_name.clear();
            os_name.str("");
            os_name << outs << "-" << period << "-specific_singular.vtk";
            xavier::export_submesh_VTK(local_mesh, os_name.str(), "N/A", functor_all, singular_triangles, true);
        }
        
        if (false) {
            mesh_type tmp_mesh(boundary, xavier::point_locator());
            tmp_mesh.set_tolerance(1.0e-8);
            std::vector<point_type> points;
            std::vector<data_type> data;
            orbits_to_points(points, data, xavier::__map_orbits);
            tmp_mesh.insert_points(points, data);
            
            std::ostringstream os_name, os_comment;
            os_name << outs << "-complete_mesh_after_" << period << "-search.vtk";
            p_step_func functor_all(period, 100);
            xavier::export_VTK(tmp_mesh, os_name.str(), "N/A", functor_all, true);
        }
    } // for all periods
}


























