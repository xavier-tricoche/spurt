#include <flow/lavd.hpp>

#include <data/locator.hpp>
//Refer to reconstruction/resample.hpp

#include <sstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <cctype>
#include <regex>
#include <stdexcept>
#include <locale>
#include <iomanip>
#include <misc/option_parse.hpp>
#include <sys/stat.h>

#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <filesystem>
namespace fs = std::filesystem;

using namespace spurt::lavd;

std::ofstream log_file;
spurt::log::dual_ostream spurt::lavd::_log_(log_file, std::cout, 1, 0, 0, true);

std::array<size_t, 2> res({ 1024, 1024 });	//Resolution of the reconstructed dataset
std::array<double, 4> bnds;

std::string folderPath = "./";		//Folder where the data is stored
std::string keyword = "input";	//Unique name for data items
std::string outputName = "reshaped";
std::string outputPath = "./normalized/";	//Folder where the output data is stored
char mode = 't';					//Reconstruction methodology. t for kd tree, i for iterative (brute force)
									//b for binary (has data / doesn't have data

int max_points = 0;
int min_points = 999999999;

double gaussian_kernel(vec2 point1, vec2 point2, double h) {
	double squared_distance = (point1[0] - point2[0]) * (point1[0] - point2[0]) + (point1[1] - point2[1]) * (point1[1] - point2[1]);
	return std::exp(-squared_distance / (h*h)); //Gaussian kernel
}

double wendland_c2_kernel(vec2 point1, vec2 point2, double radius) {
	double squared_distance = (point1[0] - point2[0]) * (point1[0] - point2[0]) + (point1[1] - point2[1]) * (point1[1] - point2[1]);
	double dist = sqrt(squared_distance);

	double r = dist / radius;
	if (r > 1.0) return 0.0;
	double t = 1.0 - r;
	return t * t * t * t * (4.0 * r + 1.0);
}

Eigen::VectorXd quadratic_basis(vec2 pt) {
	Eigen::VectorXd phi(6);
	phi << pt[0] * pt[0], pt[1]* pt[1], pt[0]* pt[1], pt[0], pt[1], 1.0;
	return phi;
}

/*std::vector<int> knn(const std::vector<vec2> points, const vec2 new_point, const int k) {
	struct Neighbor {
		int index;
		double distance;

		bool operator<(const Neighbor& other) const {
			return distance < other.distance;
		}
	};

	std::vector<Neighbor> knn;
	
	for (int i = 0; i < points.size(); i++) {
		knn.push_back({i, (points[i][0] - new_point[0]) * (points[i][0] - new_point[0]) + (points[i][1] - new_point[1]) * (points[i][1] - new_point[1]) });
	}
	std::nth_element(knn.begin(), knn.begin() + k, knn.end());
	knn.resize(k);

	std::vector<int> indices;
	for (const auto& n : knn)
		indices.push_back(n.index);

	return indices;
}*/

/*
std::vector<int> knn(const std::vector<vec2> points, const vec2 new_point, const int k) {
	std::vector<int> nearest_indices(k);
	std::vector<double> nearest_distances(k);

	for (int i = 0; i < k; i++) {
		nearest_distances[i] = 10000.0;
	}

	for (int i = 0; i < points.size(); i++) {
		double dist = (points[i][0] - new_point[0]) * (points[i][0] - new_point[0]) + (points[i][1] - new_point[1]) * (points[i][1] - new_point[1]);
		int index = i;
		for (int j = 0; j < k; j++) {
			if (dist < nearest_distances[j]) {
				double temp = nearest_distances[j];
				int temp2 = nearest_indices[j];

				nearest_distances[j] = dist;
				nearest_indices[j] = index;

				dist = temp;
				index = temp2;
			}
		}
	}

	return nearest_indices;
}
*/

/*struct PointCloud {
	std::vector<vec2> pts;

	inline size_t kdtree_get_point_count() const { return pts.size(); }

	inline double kdtree_get_pt(const size_t idx, int dim) const {
		return dim == 0 ? pts[idx][0] : pts[idx][1];
	}

	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const { return false; }
};

typedef nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<double, PointCloud>,
	PointCloud,
	2 // 2D
> KDTree;

struct KNNFinder {
	PointCloud cloud;
	std::unique_ptr<KDTree> tree;

	KNNFinder(const std::vector<vec2>& points) {
		cloud.pts = points;
		tree = std::make_unique<KDTree>(2, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
		tree->buildIndex();
	}

	std::vector<int> query(const vec2& new_point, int k) const {
		std::vector<size_t> ret_indices(k);
		std::vector<double> out_dists_sqr(k);
		nanoflann::KNNResultSet<double> resultSet(k);
		resultSet.init(ret_indices.data(), out_dists_sqr.data());

		double query_pt[2] = { new_point[0], new_point[1] };
		tree->findNeighbors(resultSet, query_pt, nanoflann::SearchParameters());

		std::vector<int> indices;
		for (size_t i = 0; i < k; ++i)
			indices.push_back(static_cast<int>(ret_indices[i]));

		return indices;
	}
};

std::vector<int> knn(const std::vector<vec2>& points, const vec2& new_point, int k) {
	KNNFinder finder(points);
	return finder.query(new_point, k);
}*/

/*
double mls(std::vector<vec2> points, std::vector<double> values, vec2 new_point) {
	Eigen::MatrixXd A = Eigen::MatrixXd::Zero(6, 6);
	Eigen::VectorXd b = Eigen::VectorXd::Zero(6);

	for (int i = 0; i < points.size(); i++) {
		vec2 pt = points[i];
		double w = gaussian_kernel(pt, new_point, 0.5);
		Eigen::VectorXd phi = quadratic_basis(pt);
		A += w * phi * phi.transpose();
		b += w * values[i] * phi;
	}

	Eigen::VectorXd coeffs = A.ldlt().solve(b);
	return quadratic_basis(new_point).dot(coeffs);
}

double rbf(const std::vector<vec2>* all_points_ptr, const std::vector<double>* all_values_ptr, const vec2 new_point) {
	const std::vector<vec2> all_points = *all_points_ptr;
	const std::vector<double> all_values = *all_values_ptr;
	
	//Calculate number of relevant points
	//auto t1 = std::chrono::high_resolution_clock::now();
	std::vector<int> indices = knn(all_points, new_point, 3);
	//auto t2 = std::chrono::high_resolution_clock::now();

	size_t n = indices.size();

	std::vector<vec2> points(n);
	std::vector<double> values(n);

	for (int i = 0; i < n; i++) {
		points[i] = all_points[indices[i]];
		values[i] = all_values[indices[i]];
	}

	//auto t3 = std::chrono::high_resolution_clock::now();
	Eigen::MatrixXd A(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A(i, j) = gaussian_kernel(points[i], points[i], 1);
		}
	}

	Eigen::VectorXd b(n);
	for (int i = 0; i < n; i++) {
		b(i) = values[i];
	}
	//auto t4 = std::chrono::high_resolution_clock::now();

	//auto t5 = std::chrono::high_resolution_clock::now();
	Eigen::VectorXd w = A.colPivHouseholderQr().solve(b);
	//auto t6 = std::chrono::high_resolution_clock::now();

	double result = 0.0;
	//auto t7 = std::chrono::high_resolution_clock::now();
	for (int i = 0; i < n; i++) {
		result += w[i] * gaussian_kernel(new_point, points[i], 1);
	}
	//auto t8 = std::chrono::high_resolution_clock::now();

	//std::cout << "knn time: " << (t2 - t1).count() << "\n";
	//std::cout << "matrix generation time: " << (t4 - t3).count() << "\n";
	//std::cout << "weight generation time: " << (t6 - t5).count() << "\n";
	//std::cout << "result calculation time: " << (t8 - t7).count() << "\n";
	return result;
}
*/

/*
* Works for small COUNT (i.e. 3, 5)
* 
* 
* 

double rbf_kd_tree(const spurt::point_locator<double, double, 2>& kd_tree, const vec2 new_point) {
	const double RANGE = 0.04;
	double KERNEL = 1;
	const size_t COUNT = 3;
	//std::vector<int> indices = knn(all_points, new_point, 3);
	//printf("A\n");
	//size_t n = indices.size();

	using point = spurt::data_point<double, double, 2>;
	std::list<point> neighbors;
	kd_tree.find_k_nearest(neighbors, { new_point[0], new_point[1] }, RANGE, COUNT);
	//kd_tree.find_within_range(neighbors, { new_point[0], new_point[1] }, RANGE);
	std::list<point> op;
	double avg_dist = 0.0;
	for (const auto& neighbor : neighbors) {
		double dist = sqrt(((neighbor.coordinate()[0] - new_point[0]) * (neighbor.coordinate()[0] - new_point[0])) + ((neighbor.coordinate()[1] - new_point[1]) * (neighbor.coordinate()[1] - new_point[1])));
		if (dist < RANGE) {
			op.push_back(neighbor);
			avg_dist += dist;
		}
	}
	neighbors = op;
	//printf("B\n");
	size_t n = neighbors.size();
	//printf("%ld\n", n);
	if (n == 0) {
		//printf("NOT ENOUGH POINTS!\n");
		//exit(1);
		return 0.0;
	}

	avg_dist /= n;
	//KERNEL = avg_dist;
	
	std::vector<vec2> points(n);
	std::vector<double> values(n);
	//printf("C %ld\n", n);
	size_t i = 0;
	for (const auto& neighbor : neighbors) {
		points[i] = neighbor.coordinate();
		values[i] = neighbor.data();
		i++;
	}
	//printf("D\n");
	Eigen::MatrixXd A(n, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			A(i, j) = gaussian_kernel(points[i], points[j], KERNEL);
		}
	}
	//printf("E\n");
	double lambda = 1e-6;
	//A += lambda * Eigen::MatrixXd::Identity(n, n);
	Eigen::VectorXd b(n);
	for (int i = 0; i < n; i++) {
		b(i) = values[i];
	}
	//printf("F\n");
	//std::cout << "A.rows(): " << A.rows() << ", A.cols(): " << A.cols() << "\n";
	//std::cout << "b.size(): " << b.size() << "\n";
	//std::cout << "A.allFinite(): " << A.allFinite() << "\n";
	//std::cout << "b.allFinite(): " << b.allFinite() << "\n";
	Eigen::VectorXd w = A.colPivHouseholderQr().solve(b);
	//printf("G\n");
	double result = 0.0;
	double denom = 0.0;

	for (int i = 0; i < n; i++) {
		//printf("%lf\n", w[i]);
		double local_kernel = gaussian_kernel(new_point, points[i], KERNEL);
		result += w[i] * local_kernel;
		denom += local_kernel;
	}
	//printf("H\n");
	if (denom == 0.0) return 0.0;
	//return result / denom;
	return result;
}
*/

double rbf_kd_tree(const spurt::point_locator<double, double, 2>& kd_tree, const vec2 new_point, double radius) {
	//const double RANGE = 0.06;
	//double KERNEL = 1;
	const double RANGE = radius * 1.1;	//For genuine reconstruction - Used for all things atm!
	//const double RANGE = radius * 2.2;	//To remove holes (testing)
	//const double RANGE = radius * 24;		//For testing blur
	double KERNEL = radius;
	//const size_t COUNT = 200;
	//std::vector<int> indices = knn(all_points, new_point, 3);
	//printf("A\n");
	//size_t n = indices.size();

	using point = spurt::data_point<double, double, 2>;
	std::list<point> neighbors;
	//kd_tree.find_k_nearest(neighbors, { new_point[0], new_point[1] }, RANGE, COUNT);
	kd_tree.find_within_range(neighbors, { new_point[0], new_point[1] }, RANGE);
	std::list<point> op;
	//std::vector<double> distances;
	//double avg_dist = 0.0;
	for (const auto& neighbor : neighbors) {
		double dist = sqrt(((neighbor.coordinate()[0] - new_point[0]) * (neighbor.coordinate()[0] - new_point[0])) + ((neighbor.coordinate()[1] - new_point[1]) * (neighbor.coordinate()[1] - new_point[1])));
		if (dist <= RANGE) {
			op.push_back(neighbor);
			//avg_dist += dist;
			//distances.push_back(dist);
		}
	}
	neighbors = op;
	//printf("B\n");
	size_t n = neighbors.size();
	//if (n > 10) printf("Neighbor count for position %lf %lf: %ld\n", new_point[0], new_point[1], n);
	//All of the densest points are located around the edges. Probably where new points are seeded
	if (n == 0) {
		//printf("NOT ENOUGH POINTS!\n");
		//exit(1);
		//return std::numeric_limits<float>::max();
		return 0.0;
	}
	if (mode == 'b') {
		//printf("%ld: %lf %lf\n", n, new_point[0], new_point[1]);
		return 1.0;
	}

	//avg_dist /= n;
	//KERNEL = avg_dist;
	//KERNEL = RANGE;
	//sort(distances.begin(), distances.end());
	//double median = (n % 2 == 0 ? (distances[n / 2 - 1] + distances[n / 2]) / 2.0 : distances[n / 2]);
	//KERNEL = sqrt((median*median) / 2);
	//KERNEL = radius * 0.1;
	KERNEL = radius / 4;

	std::vector<vec2> points(n);
	std::vector<double> values(n);
	//printf("C %ld\n", n);
	size_t i = 0;
	for (const auto& neighbor : neighbors) {
		points[i] = neighbor.coordinate();
		values[i] = neighbor.data();
		i++;
	}
	//printf("D\n");
	//Eigen::MatrixXd A(n, n);
	//for (int i = 0; i < n; i++) {
	//	for (int j = 0; j < n; j++) {
	//		A(i, j) = gaussian_kernel(points[i], points[j], KERNEL);
	//	}
	//}
	//printf("E\n");
	//double lambda = 1e-6;
	//A += lambda * Eigen::MatrixXd::Identity(n, n);
	//Eigen::VectorXd b(n);
	//for (int i = 0; i < n; i++) {
	//	b(i) = values[i];
	//}
	//printf("F\n");
	//std::cout << "A.rows(): " << A.rows() << ", A.cols(): " << A.cols() << "\n";
	//std::cout << "b.size(): " << b.size() << "\n";
	//std::cout << "A.allFinite(): " << A.allFinite() << "\n";
	//std::cout << "b.allFinite(): " << b.allFinite() << "\n";
	//Eigen::VectorXd w = A.colPivHouseholderQr().solve(b);
	//printf("G\n");

	double value = 0.0;
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		double weight = gaussian_kernel(new_point, points[i], KERNEL);
		value += weight * values[i];
		sum += weight;
	}

	return value / sum;
}

/*
double bucket_rbf(size_t grid_index, nvis::vec2 point, std::vector<std::vector<vec3>> buckets, int k) {

	grid_index = 2048;

	//Grab adjacent buckets
	int x = grid_index / res[1];
	int y = grid_index % res[1];
	//printf("%ld: %d %d\n", grid_index, x, y);
	size_t neighboring_buckets[4] = { x * res[1] + y - x, x * res[1] + y - 1 - x, x * res[1] + y - res[1] - x + 1, x * res[1] + y - res[1] - x + 1 - 1 };
	if (y == 0) {
		neighboring_buckets[1] = -1;
		neighboring_buckets[3] = -1;
	}
	if (y == res[1] - 1) {
		neighboring_buckets[0] = -1;
		neighboring_buckets[2] = -1;
	}
	if (x == 0) {
		neighboring_buckets[2] = -1;
		neighboring_buckets[3] = -1;
	}
	if (x == res[0] - 1) {
		neighboring_buckets[0] = -1;
		neighboring_buckets[1] = -1;
	}

	//printf("%ld %ld %ld %ld\n", neighboring_buckets[0], neighboring_buckets[1], neighboring_buckets[2], neighboring_buckets[3]);

	//Create lists of candidate nodes
	std::vector<vec2> points;
	std::vector<double> values;

	for (int i = 0; i < 4; i++) {
		if (neighboring_buckets[i] == -1) continue;
		for (const auto& point : buckets[neighboring_buckets[i]]) {
			points.push_back(nvis::vec2(point[0], point[1]));
			values.push_back(point[2]);
		}
	}
	//knn
	std::vector<int> indices = knn(points, point, k);
	//printf("%ld: %d %d %d\n", indices.size(), indices[0], indices[1], indices[2]);
	while (points.size() < 3) {
		points.push_back(nvis::vec2(0.0, 0.0));
		values.push_back(0.0);
	}

	std::vector<vec2> points_inter(k);
	std::vector<double> values_inter(k);

	//printf("A\n");

	for (int i = 0; i < k; i++) {
		//printf("%d\n", i);
		points_inter[i] = points[indices[i]];
		//printf("%d\n", i);
		values_inter[i] = values[indices[i]];
		//printf("%d\n", i);
	}
	//printf("B\n");

	points = points_inter;
	values = values_inter;

	Eigen::MatrixXd A(k, k);
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < k; j++) {
			A(i, j) = gaussian_kernel(points[i], points[i], 1);
		}
	}

	Eigen::VectorXd b(k);
	for (int i = 0; i < k; i++) {
		b(i) = values[i];
	}

	Eigen::VectorXd w = A.colPivHouseholderQr().solve(b);

	double result = 0.0;

	for (int i = 0; i < k; i++) {
		result += w[i] * gaussian_kernel(point, points[i], 1);
	}

	return result;
}
*/

//Takes an input file, reconstructs it, and saves it
void convert_file(const char* filename) {
	//printf("Normalizing file: %s\n", filename);

	//Read in data from the file
	Nrrd* nin = spurt::nrrd_utils::readNrrd((folderPath + filename).c_str());
	spurt::nrrd_utils::nrrd_traits traits(nin);
	std::vector<double> data;
	spurt::nrrd_utils::to_vector(data, nin);

	//Sort data by meaning
	std::vector<vec2> pos;
	std::vector<double> val;
	for (int i = 0; i < data.size(); i += 3) {
		pos.push_back(nvis::vec2(data[i], data[i + 1]));
		val.push_back(data[i + 2]);													///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	}

	//Create normalized field
	const std::vector<std::string>& comments = traits.comments();
	if (comments.size() == 0) {
		//printf("Using hardcoded bounds/res\n");
		res[0] = 1024;
		res[1] = 1024;
		bnds[0] = 262.05;
		bnds[1] = 18.05;
		bnds[2] = 279.75;
		bnds[3] = 30.79;
	}
	else {
		std::regex line_regex("[ ]*\\*[ ]*(.*)[ ]*=[ ]*(.*)$");
		std::regex res_regex("(.*)[ ]*x[ ]*(.*)$");
		std::regex bnd_regex("(.*)[ ]*->[ ]*(.*)$");
		std::smatch line_match, res_match, bnd_match;
		std::string param, val_as_str;
		if (std::regex_search(comments[6], line_match, line_regex)) {
			param = line_match[1].str();
			val_as_str = line_match[2].str();
		}
		if (param == "resolution") {
			if (std::regex_search(val_as_str, res_match, res_regex)) {
				res[0] = std::stoi(res_match[1].str());
				res[1] = std::stoi(res_match[2].str());
			}
			else throw std::runtime_error("invalid resolution syntax in restart file");
		}
		//printf("%s %s\n", param.c_str(), val_as_str.c_str());

		if (std::regex_search(comments[7], line_match, line_regex)) {
			param = line_match[1].str();
			val_as_str = line_match[2].str();
		}
		if (param == "bounds") {
			if (std::regex_search(val_as_str, bnd_match, bnd_regex)) {
				std::istringstream iss1(bnd_match[1].str());
				nvis::vec2 x;
				iss1 >> x;
				bnds[0] = x[0];
				bnds[1] = x[1];
				std::istringstream iss2(bnd_match[2].str());
				iss2 >> x;
				bnds[2] = x[0];
				bnds[3] = x[1];
			}
			else throw std::runtime_error("invalid bounds syntax in restart file");
		}
	}

	double stepx = (bnds[2] - bnds[0]) / (res[0] - 1);
	double stepy = (bnds[3] - bnds[1]) / (res[1] - 1);

	double radius = sqrt((stepx * stepx) + (stepy * stepy)) / 2;
	//printf("Calculated radius: %lf\n", radius);

	//printf("%lf x %lf step size.\n", stepx, stepy);
	//exit(0);

	std::vector<vec2> normalized_pos;
	for (int i = 0; i < res[0]; i++) {
		for (int j = 0; j < res[1]; j++) {
			//y increments, then x increments
			normalized_pos.push_back(nvis::vec2(bnds[0] + stepx*i, bnds[1] + stepy*j));
			//printf("%lf %lf\n", normalized_pos.back()[0], normalized_pos.back()[1]);
		}
	}

	//Generate values for each point
	std::vector<double> normalized_vals(normalized_pos.size());
	//std::vector<double> normalized_real(normalized_pos.size(), 0);

	// ------------------------------------------------------------------------------

	//Create KD Tree
	using point = spurt::data_point<double, double, 2>;
	using locator = spurt::point_locator<double, double, 2>;

	std::vector<point> points;
	printf("Generating vector of %ld points\n", val.size());
	if (val.size() > max_points) max_points = val.size();
	if (val.size() < min_points) min_points = val.size();
	for (size_t i = 0; i < val.size(); i++) {
		points.emplace_back(point({ pos[i][0], pos[i][1] }, val[i]));
	}
	//printf("Generating kd-tree from vector\n");
	locator kd_tree(points.begin(), points.end());
	//printf("kd-tree created. Beginning normalization\n");

	// ----------------------------------------------------------------------------------

	//Bucket approach
	//printf("Generating buckets...\n");
	//std::vector<std::vector<vec3>> buckets((res[0]-1) * (res[1]-1));
	//for (size_t i = 0; i < val.size(); i++) {
	//	int x = (int)((pos[i][0] - bnds[0]) / stepx);
	//	int y = (int)((pos[1][1] - bnds[1]) / stepy);
	//	if (x < 0 || x > res[0] - 1 || y < 0 || y > res[1] - 1) {
	//		printf("Invalid bucket!\n");
	//		exit(1);
	//	}
	//	buckets[x + (res[0] - 1) * y].push_back(nvis::vec3(pos[i][0], pos[i][1], val[i]));
	//}
	//printf("Buckets generated\n");

	// ---------------------------------------------------------------------------------

	//Generate normalization timer
	spurt::ProgressDisplay progress(false);
	progress.start(normalized_vals.size(), "normalization");
	progress.set_active(true);

	#pragma omp parallel
	{
		#pragma omp for schedule(dynamic,1)
		for (size_t i = 0; i < normalized_vals.size(); i++) {
			#if _OPENMP
			const int thread = omp_get_thread_num();
			#else
			const int thread = 0;
			#endif
			if (!thread) {
				progress.update(i);
			}

			//if (mode == 'i') {
			//	normalized_vals[i] = rbf(&pos, &val, normalized_pos[i]);
			//}
			//else if (mode == 't') {
				//if (normalized_pos[i][0] < 275) continue;
				normalized_vals[i] = rbf_kd_tree(kd_tree, normalized_pos[i], radius);
				//if (normalized_vals[i] > std::numeric_limits<float>::max() / 2) {
				//	normalized_vals[i] = 0;
				//	normalized_real[i] = 1;
				//}
			//}
			//normalized_vals[i] = bucket_rbf(i, normalized_pos[i], buckets, 3);
		}
	}
	progress.end();

	//Write to NRRD
	std::vector<double> out_data(normalized_vals.size() * 3);
	for (int i = 0; i < normalized_vals.size(); i++) {
		out_data[i * 3] = normalized_pos[i][0];
		out_data[i * 3 + 1] = normalized_pos[i][1];
		out_data[i * 3 + 2] = normalized_vals[i];
		//out_data[i * 4 + 3] = normalized_real[i];
	}

	Nrrd* nrrd = nrrdNew();
	size_t size[2] = { normalized_vals.size(), 3 };
	if (nrrdWrap_nva(nrrd, out_data.data(), nrrdTypeDouble, 2, size)) {
		std::cout << "Error wrapping NRRD: " << biffGetDone(NRRD) << std::endl;
	}
	//Save to file
	std::string name = std::string(outputPath + outputName + "_") + (filename + strlen(filename) - 11);
	if (nrrdSave(name.c_str(), nrrd, NULL)) {
		std::cerr << "Error writing NRRD: " << biffGetDone(NRRD) << std::endl;
	}

	nrrdNix(nrrd);
	nrrdNuke(nin);
	std::cout << "Wrote NRRD file: " << name.c_str() << std::endl;
}

int main(int argc, const char* argv[]) {
	namespace xcl = spurt::command_line;

	xcl::option_traits
		required(true, false, "Required Options"),
		optional(false, false, "Optional Group");
	xcl::option_parser parser(argv[0],
		"Reconstruct lavd data");

	try {
		parser.use_short_symbols(true);
		parser.use_brackets(true);
		parser.add_value("input_directory", folderPath, "Input directory", optional);
		parser.add_value("output_directory", outputPath, "Output directory", optional);
		parser.add_value("input_keyword", keyword, "Keyword for input files", optional);
		parser.add_value("mode", mode, "Reconstruction method", optional);	//t for kdtree, b for binary
		parser.parse(argc, argv);
	}
	catch (std::runtime_error& err) {
		printf("An error has occurred while parsing input arguments.\n");
		exit(1);
	}

	if (mode != 't' && mode != 'b') {
		printf("Invalid mode input. Use 't' for kd tree or 'b' for iterative.\n");
		exit(1);
	}
	struct stat sb;
	if (stat(folderPath.c_str(), &sb) != 0) {
		printf("Provided input directory does not exist.\n");
		exit(1);
	}
	if (stat(outputPath.c_str(), &sb) != 0) {
		printf("Provided output directory does not exist.\n");
		exit(1);
	}

	for (const auto& entry : fs::directory_iterator(folderPath)) {
		if (entry.is_regular_file()) {
			std::string filename = entry.path().filename().string();

			//Skip files that are finished
			bool skip = false;
			char hours[5];
			for (int i = 0; i < 5; i++) {
				hours[i] = filename.c_str()[strlen(filename.c_str()) - 11 + i];
			}
			for (const auto& entry : fs::directory_iterator(outputPath)) {
				if (entry.is_regular_file()) {
					std::string outfile = entry.path().filename().string();

					if (outfile.find(hours) != std::string::npos) {
						skip = true;
						break;
					}
				}
			}
			if (skip && filename.find(".nrrd") != std::string::npos && filename.find(keyword) != std::string::npos) {
				//printf("Skipping file %s due to overlapping hour count\n", filename.c_str());
				continue;
			}

			if (filename.find(keyword) != std::string::npos && filename.find(".nrrd") != std::string::npos) {
				printf("\n");
				convert_file(filename.c_str());
				printf("\n");
			}
		}
	}
	printf("Longest file: %d\n Shortest file: %d\n", max_points, min_points);
}
