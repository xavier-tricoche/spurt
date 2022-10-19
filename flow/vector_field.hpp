#ifndef __FLOW_VECTOR_FIELD_HPP__
#define __FLOW_VECTOR_FIELD_HPP__

#include <cmath>
#include <string>
#include <vector>
#include <iostream>
#include <stdexcept>

// teem-related interface
#include <image/probe.hpp>
#include <image/nrrd_wrapper.hpp>

namespace spurt {

namespace {
    inline double pos_modulo(double x, double min, double max) {
        double y = std::remainder(x-min, max-min);
        if (y<0) y+=max-min;
        return min + y;
    }
}

// gage_vector_field: Vector field using Teem/gage for data reconstruction.
// This class is *NOT* thread safe because it accesses a single gageContext
// object, which itself is not thread safe.
class gage_vector_field {
    typedef gage_interface::vector_wrapper wrapper_t;
    typedef bbox3 bounds_t;
public:
    typedef wrapper_t::deriv3_t  derivative_type;
    typedef wrapper_t::point3_t  point_type;
    typedef wrapper_t::value3_t  vector_type;
    typedef double               scalar_type;

    constexpr static size_t dimension = 3;

private:
    point_type transform(const point_type& p) const {
        if (!m_is_periodic) return p;
        point_type q(p);
        for (int i=0; i<3; ++i) {
            if (m_periodic[i]) {
                if (q[i] < m_bounds.min()[i] || q[i] > m_bounds.max()[i]) {
                    q[i] = pos_modulo(q[i], m_bounds.min()[i], m_bounds.max()[i]);
                }
            }
        }
        return q;
     }

public:
    gage_vector_field() : m_wrapper(), m_name(), m_have_jacobian(false), 
        m_periodic({false, false, false}) {}
        
    gage_vector_field(const Nrrd* nin, const std::string name="unknown",
                      bool have_jac=true,
                      std::array<bool, 3> periodic = std::array<bool, 3>({false, false, false}))
        : m_wrapper(nin, spurt::gage_interface::BC_INTERP, have_jac),
          m_name(name), m_have_jacobian(have_jac), m_periodic(periodic) {
        m_wrapper.use_world();
        nrrd_traits traits(nin);
        const std::vector<double>& mins = traits.mins();
        const std::vector<double>& maxs = traits.maxs();
        m_bounds.min() = spurt::vec3(mins[1], mins[2], mins[3]);
        m_bounds.max() = spurt::vec3(maxs[1], maxs[2], maxs[3]);
        m_is_periodic = std::any_of(m_periodic.begin(), m_periodic.end(), [](bool b){return b;});
    }
    
    gage_vector_field(const gage_vector_field& other)
        : m_wrapper(other.m_wrapper), m_name(other.m_name), m_bounds(other.m_bounds),
          m_have_jacobian(other.m_have_jacobian), m_periodic(other.m_periodic), 
          m_is_periodic(other.m_is_periodic) {}
          
    const std::string& name() const { return m_name; }

    bool operator()(const point_type& x, vector_type& v) const {
        return m_wrapper.value(transform(x), v);
    }

    vector_type operator()(const point_type& x) const {
        vector_type v;
        if (!m_wrapper.value(transform(x), v)) {
            std::ostringstream os;
            point_type y = transform(x);
            os << "Probing velocity outside domain of definition at " << y << " =transform(" << x << ")";
            throw std::runtime_error(os.str());
        }
        return v;
    }

    scalar_type vorticity_scalar(const point_type& x) const {
        if (!m_have_jacobian) {
            throw std::runtime_error("Vorticity computation deactivated");
        }
        derivative_type J;
        bool ok=m_wrapper.jacobian(transform(x), J);
        if (!ok) {
            std::ostringstream os;
            os << "Probing vorticity outside domain of definition at "
                << transform(x) << " =transform(" << x << ")";
            throw std::runtime_error(os.str());
        }
        return J(1,0)-J(0,1);
    }
    scalar_type vorticity(const point_type& x) const {
        return vorticity_scalar(x);
    }

    vector_type vorticity_vector(const point_type& x) const {
        if (!m_have_jacobian) {
            throw std::runtime_error("Vorticity computation deactivated");
        }
        derivative_type J;
        bool ok=m_wrapper.jacobian(transform(x), J);
        if (!ok) {
            std::ostringstream os;
            os << "Probing vorticity outside domain of definition at "
                << transform(x) << " =transform(" << x << ")";
            throw std::runtime_error(os.str());
        }
        return vector_type(J(2,1)-J(1,2), J(0,2)-J(2,0), J(1,0)-J(0,1));
    }

    bool jacobian(const point_type& x, derivative_type& J) const {
        if (!m_have_jacobian) {
            throw std::runtime_error("Jacobian computation deactivated");
        }
       return m_wrapper.jacobian(transform(x), J);
    }

    derivative_type jacobian(const point_type& x) const {
        derivative_type J;
        if (!jacobian(transform(x), J)) {
            std::ostringstream os;
            os << "Probing Jacobian outside domain of definition at " << transform(x) << " =transform(" << x << ")";
            throw std::runtime_error(os.str());
        }
        return J;
    }

    wrapper_t           m_wrapper;
    std::string         m_name;
    bool                m_have_jacobian;
    std::array<bool, 3> m_periodic;
    bool                m_is_periodic;
    bounds_t            m_bounds;
};

} // xavier

#endif
