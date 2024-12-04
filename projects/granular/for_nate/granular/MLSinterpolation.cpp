// std
#include <iostream>
#include <vector>
#include <fstream>
// kdtree++
#include <kdtree++/kdtree.hpp>
// teem
#include <teem/nrrd.h>
// nvis
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <util/timer.hpp>
// xavier
#include <math/MLS.hpp>

const double particle_radius = 0.05;
const size_t max_nb_neighbors = 100;

template<typename V>
class point
{
public:
	typedef V 						vector_type;
	typedef typename V::value_type	value_type;

	point() : __v(), __idx(-1) {}
	point(value_type x, value_type y, value_type z, unsigned int idx = -1)
			: __v(x, y, z), __idx(idx) {}
	point(const vector_type& v, unsigned int idx = -1) : __v(v), __idx(idx) {}
	point(const point& p) : __v(p.__v), __idx(p.__idx) {}

	const vector_type& pos() const {
		return __v;
	}

	vector_type& pos() {
		return __v;
	}

	unsigned int index() const {
		return __idx;
	}

	unsigned int& index() {
		return __idx;
	}

	value_type distance_to(const point& p) const {
		return norm(p.__v - __v);
	}

	value_type operator[](size_t N) const {
		return __v[N];
	}

private:
	vector_type		__v;
	unsigned int 	__idx;
};

template<typename V>
inline bool operator==(const point<V>& p0, const point<V>& p1)
{
	return (p0.pos() == p1.pos());
}

template<int N>
inline nvis::fixed_vector<float, N> to_vector(float* ptr)
{
	nvis::fixed_vector<float, N> r;
	for (int i = 0 ; i < N ; ++i) {
		r[i] = ptr[i];
	}
	return r;
}

char *coord_f, *stress_f, *output_f;
int n;
float radius;
void initialize(int argc, char* argv[])
{
	hestOpt *hopt = NULL;
	hestParm *hparm;
	airArray *mop;
	char *me;

	mop = airMopNew();
	me = argv[0];
	hparm = hestParmNew();
	airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
	hparm->elideSingleOtherType = AIR_TRUE;
	hestOptAdd(&hopt, "c", 		"coordinates", 			airTypeString, 	1, 1, &coord_f, 			NULL, 		"input coordinates file name");
	hestOptAdd(&hopt, "s", 		"stree", 				airTypeString, 	1, 1, &stress_f, 			NULL, 		"input stress value file");
	hestOptAdd(&hopt, "o", 		"output", 				airTypeString, 	1, 1, &output_f, 			NULL, 		"output stress field file");
	hestOptAdd(&hopt, "n", 		"# samples", 			airTypeInt, 	0, 1, &n, 					"100000", 	"number of random samples");
	hestOptAdd(&hopt, "r",		"support radius",		airTypeFloat,	0, 1, &radius,				"0.2",		"radius of weighting function support");

	hestParseOrDie(hopt, argc - 1, argv + 1, hparm,
	               me, "Compute smooth MLS interpolation of principal stress tensor over cylindrical domain",
	               AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
	initialize(argc, argv);

	// import particle positions
	Nrrd* nin1 = nrrdNew();
	if (nrrdLoad(nin1, coord_f, NULL)) {
		std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
		return -1;
	}
	// import particle stress tensor
	Nrrd* nin2 = nrrdNew();
	if (nrrdLoad(nin2, stress_f, NULL)) {
		std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
		return -1;
	}

	// we will be working in float precision
	typedef nvis::fixed_vector<float, 3>	vector_type;
	typedef point<vector_type>				point_type;
	typedef nvis::bounding_box<vector_type>	box_type;
	typedef nvis::fixed_vector<float, 6>	value_type;

	float *coords = (float*)nin1->data;
	float *stress = (float*)nin2->data;
	std::vector<value_type> sampled_val;
	std::vector<vector_type> sampled_pos;

	// number of valid particles in file
	int nbp = nin2->axis[1].size; // only valid positions have associated stress tensor

	// construct kd-tree
	KDTree::KDTree<3, point_type> tree;
	box_type box;

	// insert positions in tree
	for (int i = 0 ; i < nbp ; ++i) {
		vector_type v(to_vector<3>(&coords[3*i]));
		box.add(v);
		tree.insert(point_type(v, i));
		sampled_pos.push_back(v);
		sampled_val.push_back(to_vector<6>(&stress[6*i]));
	}

	// initialize CLAPACK computation space
	MLS::CLAPACK_helper helper(max_nb_neighbors, 	// max number of neighbors expected to be in radius
	                           3, 					// 3-space
	                           0, 					// 0-th order polynomial fit
	                           6, 					// N rhs for symmetric 3x3 matrix
	                           true					// use SVD to deal with degenerate cases
	                          );

	// now loop over a number of random positions
	vector_type min = box.min();
	vector_type max = box.max();
	vector_type diagonal = box.size();

	std::cout << "bounding box in input: " << min << " - " << max << std::endl;
	double cylinder_radius = max[0] + particle_radius;
	double minz = min[2] - particle_radius;
	double height = diagonal[2] + 0.1;
	std::cerr << "radius = " << radius
	          << ", minz = " << minz
	          << ", height = " << height << '\n';


	srand48(time(0));

	nvis::timer timer;
	double search_time = 0;
	double mls_time = 0;
	for (int i = 0 ; i < n ; ++i) {

		// std::cerr << "attempt #" << n + 1 << ": ";

		// select a random position within the cylinder
		double theta = 2.*M_PI * drand48();
		double dz = height * drand48();
		nvis::vec3 x(particle_radius*cos(theta), particle_radius*sin(theta), minz + dz);

		point_type p(x[0], x[1], x[2]);

		std::vector<point_type> in_cube;

		std::vector<double> coef(6*MLS::dof[1][0]);
		std::vector<nvis::vec3> pos;
		std::vector<value_type> val;

		timer.restart();
		tree.find_within_range(p, radius, std::back_inserter(in_cube));
		search_time += timer.elapsed();

		if (in_cube.empty()) {
			std::cerr << "there are no particles within distance " << radius
			          << " of " << x << std::endl;
			continue;
		}
		else {
			// std::cerr << "found " << in_cube.size() << " neighbors of " << x
			//           << " within manhattan distance " << radius << std::endl;
			for (int i = 0 ; i < in_cube.size() ; ++i) {
				const point_type& close_p = in_cube[i];
				vector_type dist = close_p.pos() - p.pos();
				double d = nvis::norm(dist);
				if (d > radius) continue;
				// std::cerr << "\tparticle #" << close_p.index() << ": " << close_p.pos() << std::flush;
				// std::cerr << " is at distance " << d << std::endl;

				pos.push_back(close_p.pos());
				val.push_back(to_vector<6>(&stress[close_p.index()*6]));
			}
		}

		if (pos.size() < 3) {
			std::cerr << "not enough particles in neighborhood (" << pos.size() << ")\n";
			continue;
		}

		timer.restart();
		int deg = MLS::MLS(coef, pos, val, x, helper, radius);
		mls_time += timer.elapsed();

		sampled_pos.push_back(x);
		value_type m;
		for (int j = 0 ; j < 6 ; ++j) {
			m[j] = coef[j];
		}
		sampled_val.push_back(m);
	}

	std::cerr << "Number of tensor values computed: " << sampled_val.size() - nbp << '\n';
	std::cerr << "Total search time: " << search_time << " s.\n";
	std::cerr << "Average search time per sample: " << search_time / (double)n << " s.\n";
	std::cerr << "Total MLS time: " << mls_time << " s.\n";
	std::cerr << "Average MLS time per sample: " << mls_time / (double)n << " s.\n";

	std::fstream vtk(output_f, std::ios::out);
	vtk << "# vtk DataFile Version 2.0\n"
	<< "Principal stress tensor computed from " << coord_f << " and " << stress_f << '\n'
	<< "ASCII\n"
	<< "DATASET POLYDATA\n"
	<< "POINTS " << sampled_pos.size() << " float\n";
	for (int i = 0 ; i < sampled_pos.size() ; ++i) {
		vtk << sampled_pos[i][0] << " "
		<< sampled_pos[i][1] << " "
		<< sampled_pos[i][2] << '\n';
	}
	vtk << "POINT_DATA " << sampled_val.size() << '\n'
	<< "TENSORS principal_stress float\n";
	for (int i = 0 ; i < sampled_val.size() ; ++i) {
		const value_type& m = sampled_val[i];
		vtk << m[0] << " " << m[1] << " " << m[2] << '\n'
		<< m[1] << " " << m[3] << " " << m[4] << '\n'
		<< m[2] << " " << m[4] << " " << m[5] << '\n';
	}
	vtk.close();

	std::cerr << "output VTK file exported\n";

	return 0;
}






















































