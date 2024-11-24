#include <flow/vector_field.hpp>


namespace spurt {


bool gage_vector_field::operator()(const gage_vector_field::point_type& x, gage_vector_field::vector_type& v) const {
    return m_wrapper.value(transform(x), v);
}

gage_vector_field::vector_type gage_vector_field::operator()(const gage_vector_field::point_type& x) const {
    vector_type v;
	// std::cout << "vector_field(" << x << ")" << std::endl;
    if (!m_wrapper.value(transform(x), v)) {
        std::ostringstream os;
		point_type y = transform(x);
        os << "Probing velocity outside domain of definition at " << y << " =transform(" << x << ")";
		throw std::runtime_error(os.str());
    }
    return v;
}

gage_vector_field::scalar_type gage_vector_field::vorticity_scalar(const gage_vector_field::point_type& x) const {
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
gage_vector_field::scalar_type gage_vector_field::vorticity(const gage_vector_field::point_type& x) const {
	return vorticity_scalar(x);
}

gage_vector_field::vector_type gage_vector_field::vorticity_vector(const gage_vector_field::point_type& x) const {
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

bool gage_vector_field::jacobian(const gage_vector_field::point_type& x, gage_vector_field::derivative_type& J) const {
    if (!m_have_jacobian) {
        throw std::runtime_error("Jacobian computation deactivated");
    }
   return m_wrapper.jacobian(transform(x), J);
}

gage_vector_field::derivative_type gage_vector_field::jacobian(const gage_vector_field::point_type& x) const {
	derivative_type J;
	if (!jacobian(transform(x), J)) {
        std::ostringstream os;
        os << "Probing Jacobian outside domain of definition at " << transform(x) << " =transform(" << x << ")";
        throw std::runtime_error(os.str());
	}
	return J;
}

}
