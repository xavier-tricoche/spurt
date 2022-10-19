#ifndef __POINCARE_INDEX_HPP__
#define __POINCARE_INDEX_HPP__

#include <complex>
#include <math/math.hpp>
#include <math/angle.hpp>

#include "maps_lib/misc.hpp"
#include "maps_lib/definitions.hpp"
#include "math/fixed_vector.hpp"
#include "maps_lib/orbits.hpp"
#include "maps_lib/misc.hpp"

namespace map_analysis {
typedef std::pair<nvis::vec2, nvis::vec2>   step_type;

struct degenerate_point_exception : public std::exception {
    degenerate_point_exception(const nvis::vec2& x) : _x(x) {}
    
    const nvis::vec2& pos() const {
        return _x;
    }
    nvis::vec2 _x;
};

template<typename Mesh, typename Map>
struct poincare_index {

    typedef Mesh                                mesh_type;
    typedef Map                                 integrator_type;
    typedef typename mesh_type::point_type      point_type;
    typedef typename mesh_type::data_type       data_type;
    typedef typename mesh_type::triangle_type   triangle_type;
    typedef std::pair<point_type, point_type>   uncertain_vec_type;
    typedef xavier::Edge                        edge_index_type;
    
    poincare_index() {}
    
    inline int direct(const point_data d[3], unsigned int period,
                      const metric_type& metric, std::list<step_type>& steps) const {
        double dtheta = 0;
        for (int i = 0 ; i < 3 ; ++i) {
            nvis::vec2 v0 = vector_value(d[i], period, metric);
            steps.push_back(step_type(d[i].pos(), d[i].pos() + v0));
            dtheta += xavier::signed_angle(v0, vector_value(d[(i+1)%3], period, metric));
        }
        return (int)round(0.5*dtheta / M_PI);
    }
    
    inline int direct(const point_data d[3], unsigned int period, const metric_type& metric) const {
        std::list<step_type> dummy;
        return direct(d, period, metric, dummy);
    }
    
    inline std::pair<double, double>
    rotation_mag(const point_data d[3], unsigned int period,
                 const metric_type& metric, std::list<step_type>& steps) const {
        double dtheta = 0, dtheta_max = 0;
        for (int i = 0 ; i < 3 ; ++i) {
            nvis::vec2 v0 = vector_value(d[i], period, metric);
            double delta = fabs(xavier::signed_angle(v0, vector_value(d[(i+1)%3], period, metric)));
            dtheta += delta;
            dtheta_max = std::max(dtheta_max, delta);
            steps.push_back(step_type(d[i].pos(), d[i].pos() + v0));
        }
        return std::make_pair(0.5*dtheta / M_PI, 0.5*dtheta_max / M_PI);
    }
    
    inline std::pair<double, double>
    rotation_mag(const point_data d[3], unsigned int period,
                 const metric_type& metric) const {
        std::list<step_type> dummy;
        return rotation_mag(d, period, metric, dummy);
    }
    
    uncertain_vec_type
    evaluate_step(const nvis::vec2& x, const integrator_type& map,
                  unsigned int period, const metric_type& metric) const {
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
                                 const metric_type& metric,
                                 std::list<step_type>& steps,
                                 unsigned int depth,
                                 bool conservative = false) const {
        uncertain_vec_type uv0, uv1;
        uv0 = evaluate_step(x0, map, period, metric);
        uv1 = evaluate_step(x1, map, period, metric);
        const point_type& v0 = uv0.first;
        const point_type& v1 = uv1.first;
        steps.push_back(step_type(x0, x0 + v0));
        steps.push_back(step_type(x1, x1 + v1));
        double theta = xavier::signed_angle(v0, v1);
        
        if (fabs(theta) < dtheta) {
            return theta;
        } else if (nvis::norm(x1 - x0) <= dx) {
            if (conservative) {
                std::cerr << "degenerate point for period " << period << " at " << 0.5*(x0 + x1)
                          << ", theta = " << theta << ", depth = " << depth << '\n';
                // throw std::runtime_error("smooth_rotation_angle: step size underflow");
                throw degenerate_point_exception(0.5*(x0 + x1));
            }
            return theta;
        }
        
        nvis::vec2 x = 0.5 * (x0 + x1);
        return smooth_rotation_angle(x0, x, map, period, dtheta, dx, metric, steps, depth + 1, conservative) +
               smooth_rotation_angle(x, x1, map, period, dtheta, dx, metric, steps, depth + 1, conservative);
    }
    
    double smooth_rotation_angle(const nvis::vec2& x0, const nvis::vec2& x1,
                                 const integrator_type& map, unsigned int period,
                                 double dtheta, double dx,
                                 const metric_type& metric) const {
        std::list<step_type> dummy;
        return smooth_rotation_angle(x0, x1, map, period, dtheta, dx, metric, dummy, 0);
    }
    
    double safe_rotation_angle(const nvis::vec2& x0, const nvis::vec2& x1,
                               // v0 and v1 are the values associated with x0 and x1...
                               const nvis::vec2& v0, const nvis::vec2& v1,
                               // ...provided their respective flag valid? is true
                               bool valid0, bool valid1,
                               const integrator_type& map, unsigned int period,
                               double dtheta, double dx,
                               const metric_type& metric,
                               unsigned int depth,
                               std::list<step_type>& steps,
                               bool conservative = false) const {
        nvis::vec2 __v0, __v1, dummy;
        __v0 = (valid0 ? v0 : evaluate_step(x0, map, period, metric).first);
        __v1 = (valid1 ? v1 : evaluate_step(x1, map, period, metric).first);
        
        // only add steps to list upon their first computation
        if (!valid0) {
            steps.push_back(step_type(x0, x0 + __v0));
        }
        if (!valid1) {
            steps.push_back(step_type(x1, x1 + __v1));
        }
        
        double theta = xavier::signed_angle(__v0, __v1);
        
        if (fabs(theta) < dtheta) {
            return theta;
        } else if (nvis::norm(x1 - x0) <= dx) {
            if (conservative) {
                std::cerr << "degenerate point for period " << period << " at " << 0.5*(x0 + x1)
                          << ", theta = " << theta << ", depth = " << depth << '\n';
                // throw std::runtime_error("smooth_rotation_angle: step size underflow");
                throw degenerate_point_exception(0.5*(x0 + x1));
            }
            return theta;
        }
        
        nvis::vec2 x = 0.5 * (x0 + x1);
        return safe_rotation_angle(x0, x, __v0, dummy, true, false, map, period, dtheta, dx, metric, depth + 1, steps, conservative) +
               safe_rotation_angle(x, x1, dummy, __v1, false, true, map, period, dtheta, dx, metric, depth + 1, steps, conservative);
    }
    
    int safe(const mesh_type& mesh, unsigned int tri, const Map& map,
             unsigned int period, double dtheta, double dx,
             const metric_type& metric,
             const std::map<edge_index_type, double>& measured_edge_angles,
             std::list<step_type>& steps) const {
             
        typedef std::map<edge_index_type, double>::const_iterator   iterator_type;
        
        const triangle_type& ids = mesh.get_triangle_vertices(tri);
        point_type  p[3];
        data_type   d[3];
        mesh.get_triangle_info(p, d, tri);
        
        const integrator_type* pmap = map.clone();
        double delta_theta = 0;
        for (int i = 0 ; i < 3 ; ++i) {
            edge_index_type edge_id(ids[i], ids[(i+1)%3]);
            iterator_type iter = measured_edge_angles.find(edge_id);
            if (iter != measured_edge_angles.end()) {
                delta_theta += (ids[i] == iter->first.i0 ? + 1 : -1) * iter->second;
            } else {
                nvis::vec2 p0, p1;
                p0 = p[i];
                p1 = p[(i+1)%3];
                double dtheta = smooth_rotation_angle(p0, p1, *pmap, period, dtheta, dx, metric, steps, 0);
                delta_theta += dtheta;
            }
        }
        
        return (int)round(0.5*delta_theta / M_PI);
    }
    
    int safe(const mesh_type& mesh, unsigned int tri, const Map& map,
             unsigned int period, double dtheta, double dx,
             const metric_type& metric) const {
        std::list<step_type> dummy_list;
        std::map<edge_index_type, double> dummy_map;
        return safe(mesh, tri, map, period, dtheta, dx, metric, dummy_list, dummy_map);
    }
    
    
    template<typename RHS>
    double rotation_angle_linear_predictor(const nvis::vec2& x0, const nvis::vec2& x1,
                                           // v0 and v1 are the values associated with x0 and x1...
                                           const nvis::vec2& v0, const nvis::vec2& v1,
                                           // ...provided their respective flag valid? is true
                                           bool valid0, bool valid1,
                                           const RHS& rhs,
                                           double dtheta, double dx,
                                           unsigned int depth,
                                           std::list<step_type>& steps,
                                           bool conservative = false) const {
        // needed for internal consistency
        const double large_angle = 0.75 * M_PI;
        
        nvis::vec2 __v0, __v1, dummy;
        __v0 = (valid0 ? v0 : rhs(x0));
        __v1 = (valid1 ? v1 : rhs(x1));
        if (!valid0) {
            steps.push_back(step_type(x0, x0 + __v0));
        }
        if (!valid1) {
            steps.push_back(step_type(x1, x1 + __v1));
        }
        
        // easy cases first
        double theta = xavier::signed_angle(__v0, __v1);
        if (fabs(theta) < std::max(large_angle, dtheta)) {
            return theta;
        } else if (nvis::norm(x1 - x0) < dx && depth > 5) {
            if (conservative) {
                throw degenerate_point_exception(0.5*(x0 + x1));
            }
            return theta;
        }
        
        // use linear model to pick next sample
        // solve for v(x).v0 = 0
        double v0sq = nvis::inner(__v0, __v0);
        double u = v0sq / (v0sq - nvis::inner(__v0, __v1));
        nvis::vec2 x = (1 - u) * x0 + u * x1;
        return rotation_angle_linear_predictor(x0, x, __v0, __v0, true, false,
                                               rhs, dtheta, dx, depth + 1,
                                               steps, conservative) +
               rotation_angle_linear_predictor(x, x1, __v1, __v1, false, true,
                                               rhs, dtheta, dx, depth + 1,
                                               steps, conservative);
    }
    
    template<typename RHS, typename J>
    double rotation_angle_linear_predictor_with_jacobian(
        const nvis::vec2& x0, const nvis::vec2& x1,
        // v0 and v1 are the values associated with x0 and x1...
        const nvis::vec2& v0, const nvis::vec2& v1,
        // ...provided their respective flag valid? is true
        bool valid0, bool valid1,
        const RHS& rhs, const J& jacobian,
        double dtheta, unsigned int max_depth,
        unsigned int depth,
        std::list<step_type>& steps) const {
        
        // needed for internal consistency
        const double large_angle = 0.75 * M_PI;
        
        nvis::vec2 __v0, __v1, dummy;
        __v0 = (valid0 ? v0 : rhs(x0));
        __v1 = (valid1 ? v1 : rhs(x1));
        if (!valid0) {
            steps.push_back(step_type(x0, x0 + __v0));
        }
        if (!valid1) {
            steps.push_back(step_type(x1, x1 + __v1));
        }
        
        // easy cases first
        double theta = xavier::signed_angle(__v0, __v1);
        if (fabs(theta) < std::max(large_angle, dtheta)) {
            return theta;
        } else if (depth == max_depth && nvis::norm(rhs.error(x0)) > 0.1*nvis::norm(v0)) {
            // first try again at higher integration precision
            rhs.set_precision(0.1*rhs.get_precision());
            // std::cerr << "decreasing epsilon to " << rhs.get_precision() << '\n';
            rotation_angle_linear_predictor_with_jacobian(
                x0, x1, v0, v1, false, false, rhs, jacobian,
                dtheta, max_depth, depth + 1, steps);
        } else if (depth > max_depth) {
            // use local linear extrapolation to figure out the rotation direction
            nvis::mat2 j = jacobian(x0);
            nvis::vec2 dv0 = j * (x1 - x0);
            double outer_loc = v0[0] * dv0[1] - v0[1] * dv0[0];
            // verify that local rotation direction matches linear approximation along x0-x1
            double outer_edge = v0[0] * v1[1] - v0[1] * v1[0];
            if (outer_loc * outer_edge < 0) {
                throw degenerate_point_exception(0.5*(x0 + x1));
            }
            return theta;
        } else {
            // use linear model to pick next sample
            // solve for v(x).v0 = 0
            double v0sq = nvis::inner(__v0, __v0);
            double u = v0sq / (v0sq - nvis::inner(__v0, __v1));
            nvis::vec2 x = (1 - u) * x0 + u * x1;
            return rotation_angle_linear_predictor_with_jacobian(
                       x0, x, __v0, __v0, true, false, rhs,
                       jacobian, dtheta, max_depth, depth + 1,
                       steps) +
                   rotation_angle_linear_predictor_with_jacobian(
                       x, x1, __v1, __v1, false, true, rhs,
                       jacobian, dtheta, max_depth, depth + 1,
                       steps);
        }
    }
    
    template<typename RHS, typename J>
    double rotation_angle_cubic_predictor(const nvis::vec2& x0, const nvis::vec2& x1,
                                          // v0 and v1 are the values associated with x0 and x1...
                                          const nvis::vec2& v0, const nvis::vec2& v1,
                                          // J0 and J1 are the associated Jacobian matrices
                                          const nvis::mat2& J0, const nvis::mat2& J1,
                                          // ...provided their respective flag valid? is true
                                          bool valid0, bool valid1,
                                          const RHS& rhs, const J& jacobian,
                                          double dtheta, double dx,
                                          unsigned int depth,
                                          std::list<step_type>& steps,
                                          bool conservative = false) const {
        // needed for internal consistency
        const double large_angle = 0.75 * M_PI;
        
        nvis::vec2 __v0, __v1;
        __v0 = (valid0 ? v0 : rhs(x0));
        __v1 = (valid1 ? v1 : rhs(x1));
        nvis::mat2 __J0, __J1;
        __J0 = (valid0 ? J0 : jacobian(x0));
        __J1 = (valid1 ? J1 : jacobian(x1));
        if (!valid0) {
            steps.push_back(step_type(x0, x0 + __v0));
        }
        if (!valid1) {
            steps.push_back(step_type(x1, x1 + __v1));
        }
        
        // easy cases first
        double theta = xavier::signed_angle(__v0, __v1);
        if (fabs(theta) < std::max(large_angle, dtheta)) {
            return theta;
        } else if (nvis::norm(x1 - x0) < dx) {
            if (conservative) {
                throw degenerate_point_exception(0.5*(x0 + x1));
            }
            return theta;
        }
        
        // use cubic hermite model to pick next sample
        // v(x) = v0 + j0*x + (3*(v1-v0)-2*j0-j1)*x^2 + (j1-j0-2*(v1-v0))*x^3
        nvis::vec2 e = x1 - x0;
        nvis::vec2 j0 = __J0 * e;
        nvis::vec2 j1 = __J1 * e;
        // solve for v(x).v0 = 0
        double a0 = nvis::inner(__v0, __v0);
        double a1 = nvis::inner(j0, __v0);
        double a2 = nvis::inner(3.*(__v1 - __v0) - 2.*j0 - j1, __v0);
        double a3 = nvis::inner(j1 + j0 - 2.*(__v1 - __v0), __v0);
        
        std::complex<double> u[3];
        int nbroots = xavier::cubic_equation(a3, a2, a1, a0, u);
        for (int i = 0 ; i < nbroots ; ++i) {
            if (u[i].imag()) {
                continue;
            } else if (u[i].real() < 0 || u[1].real() > 1) {
                continue;
            }
            
            nvis::vec2 x = (1 - u[i].real()) * x0 + u[i].real() * x1;
            return rotation_angle_cubic_predictor(x0, x, __v0, __v0, __J0, __J0, true, false,
                                                  rhs, jacobian, dtheta, dx, depth + 1,
                                                  steps, conservative) +
                   rotation_angle_cubic_predictor(x, x1, __v1, __v1, __J1, __J1, false, true,
                                                  rhs, jacobian, dtheta, dx, depth + 1,
                                                  steps, conservative);
        }
    }
};

inline nvis::vec3 zero_bary_coord(const nvis::vec2 v[3])
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
}


#endif















































































