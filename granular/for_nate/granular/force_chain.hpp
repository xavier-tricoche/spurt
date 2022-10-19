#include <vector>
#include <set>
#include <list>
#include <math/fixed_vector.hpp>
#include <math/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <granular/sparse_matrix.hpp>

typedef nvis::fvec3									vec3;
typedef nvis::fixed_vector<float, 6>				sym_mat3;
typedef xavier::mat3								mat3;

void to_mat3(mat3& out, const sym_mat3& in)
{
	out(0, 0) = in[0];
	out(0, 1) = out(1, 0) = in[1];
	out(0, 2) = out(0, 2) = in[2];
	out(1, 1) = in[3];
	out(1, 2) = out(2, 1) = in[4];
	out(2, 2) = in[5];
}

void to_sym_mat3(sym_mat3& out, const mat3& in)
{
	out[0] = in(0, 0);
	out[1] = in(0, 1);
	out[2] = in(0, 2);
	out[3] = in(1, 1);
	out[4] = in(1, 2);
	out[5] = in(2, 2);
}

extern sparse_matrix<float> forces; 	// defines interparticle forces and implicitly contact
extern std::vector<vec3> 	coords; 	// defines particle spatial locations
extern double				max_angle;	// max_angle used to control particle alignment

struct particle {

	particle() : index(-1), connect_to(-1),
			stress_direction(0), principal_stress(0), stress(0), cl(0) {}

	unsigned int 	index;
	unsigned int 	connect_to;
	vec3 			stress_direction;
	float			principal_stress;
	sym_mat3		stress;
	float 			cl;
};

struct Lt_particle {
	int operator()(const particle& pc0, const particle& pc1) const {
		return pc0.index < pc1.index;
	}
};

struct Lt_stress {
	int operator()(const particle& pc0, const particle& pc1) const {
		return pc0.principal_stress < pc1.principal_stress;
	}
};

typedef std::set<particle, Lt_particle>	chain_type;

void eigenanalysis(nvis::vec3& dir, double& val, double& cl, const sym_mat3& sigma_mat)
{
	std::vector<double> 		evals;
	std::vector<nvis::vec3>		evecs;
	mat3						mat;
	to_mat3(mat, sigma_mat);
	xavier::eigen(evals, evecs, mat);

	dir = evecs[0];
	val = evals[0];
	cl = (evals[0] - evals[1]) / (evals[0] + evals[1] + evals[2]);
}

float total_stress(const chain_type& chain)
{
	float s = 0;
	chain_type::const_iterator it;
	for (it = chain.begin() ; it != chain.end() ; ++it) {
		s += it->principal_stress;
	}
	return s;
}

float average_stress(const chain_type& chain)
{
	return total_stress(chain) / (float)chain.size();
}

vec3 average_tangent(const chain_type& chain)
{
	mat3 variance;
	chain_type::const_iterator it;
	unsigned int n = 0;
	for (it = chain.begin() ; it != chain.end() ; ++it) {
		if (it->connect_to == (unsigned int) - 1) continue;
		++n;
		vec3 t = coords[it->index] - coords[it->connect_to];
		t /= nvis::norm(t);

		for (int r = 0 ; r < 3 ; ++r) {
			for (int c = r ; c < 3 ; ++c) {
				variance(r, c) += t[r] * t[c];
			}
		}
	}

	sym_mat3 sym;
	to_sym_mat3(sym, variance);
	sym /= (float)n;
	nvis::vec3 tangent;
	double val, cl;
	eigenanalysis(tangent, val, cl, sym);
	return cl*tangent;
}

vec3 average_stress_direction(const chain_type& chain)
{
	mat3 variance;
	chain_type::const_iterator it;
	unsigned int n = 0;
	for (it = chain.begin() ; it != chain.end() ; ++it) {
		++n;
		vec3 t = it->stress_direction;

		for (int r = 0 ; r < 3 ; ++r) {
			for (int c = r ; c < 3 ; ++c) {
				variance(r, c) += t[r] * t[c];
			}
		}
	}

	sym_mat3 sym;
	to_sym_mat3(sym, variance);
	sym /= (float)n;
	nvis::vec3 tangent;
	double val, cl;
	eigenanalysis(tangent, val, cl, sym);
	return cl*tangent;
}

float total_anisotropy(const chain_type& chain)
{
	float cl = 0;
	chain_type::const_iterator it;
	unsigned int n = 0;
	for (it = chain.begin() ; it != chain.end() ; ++it) {
		++n;
		cl += it->cl;
	}

	return cl;
}

float average_anisotropy(const chain_type& chain)
{
	return total_anisotropy(chain) / (float)chain.size();
}

// helper functions
vec3 position(int i)
{
	return coords[i];
}

// return radius vector from particle i to particle j
inline vec3 radius(int i, int j)
{
	vec3 rij = coords[j] - coords[i];
	return rij / norm(rij);
}

inline float normal_force(int i, int j)
{
	return forces(i, j);
}

sym_mat3 compute_stress(int i)
{
	// compute major eigenvector of principal stress tensor
	// by applying equation (1) to all neighbors j of particle i
	// that are interaction with it with either positive (fij>0) or negative (fji>0)
	// normal force.

	// obtain neighbors
	std::vector<std::pair<unsigned int, float> >neighbors;
	// std::cerr << "particle " << i << " has... ";
	forces.row_entries(neighbors, i);
	// std::cerr << neighbors.size() << " neighbors\n";

	vec3 posi = position(i);
	mat3 sigma;
	for (int c = 0 ; c < neighbors.size() ; ++c) {

		// std::cerr << "neighbor #" << c << " is " << neighbors[c].first << '\n';

		int n = neighbors[c].first;
		vec3 ric = radius(i, n);
		ric /= norm(ric);
		// std::cerr << "looking up contact force:... ";
		float f = neighbors[c].second;
		// std::cerr << f << '\n';
		for (int row = 0 ; row < 3 ; ++row) {
			for (int col = row ; col < 3 ; ++col) {
				sigma(row, col) = 0; // initialization to zero
				sigma(row, col) +=  f * ric[row] * ric[col];
			}
		}
	}

	// std::cerr << "sigma = " << sigma << std::endl;

	sym_mat3 stress;
	to_sym_mat3(stress, sigma);

	return stress;
}

bool is_aligned(const vec3& sigma, int i, int j, bool abs = false)
{
	// determine whether particle j is within admissible angular distance
	// from the direction pointed to by sigma from particle i
	vec3 rij = radius(i, j);
	rij /= nvis::norm(rij);
	float cos_theta = nvis::inner(sigma, rij); // scalar product
	if (abs) cos_theta = fabs(cos_theta);
	return cos_theta > cos(max_angle); // max_angle = 45 degrees in paper

	// NB: test assumes that sigma is normalized
}


















































