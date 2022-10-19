#ifndef __TOPOLOGY_DISCRETE_HPP__
#define __TOPOLOGY_DISCRETE_HPP__

#include <vector>
#include <algorithm>
#include <map>
#include <string>

#include <math/fixed_vector.hpp>

template<typename Scalar, typename Index>
struct less_than {
	bool operator()(const std::pair<Scalar, Index>& v1, const std::pair<Scalar, Index>& i2) const {
		if (i1.first < i2.first) return true;
		else if (i2.first < i1.first) return false;
		else return i1.second < i2.second;
	}
};

template<typename Triangulation, typename Scalar>
struct discrete_topology {
	typedef Triangulation triangulation_type;
	typedef Scalar scalar_type;
	typedef size_t index_type;
	typedef nvis::fixed_vector<index_type, 2> edge_type;
	typedef nvis::fixed_vector<index_type, 3> triangle_type;
	
	discrete_topology(const Triangulation& tri) {
		
	}
	
	std::vector<edge_type> m_edges;
	std::vector<triangle_type> m_triangles;
	std::vector<scalar_type> m_edge_values;
	std::vector<scalar_type> m_triangle_values;
	
};




#endif