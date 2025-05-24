#include <flow/lavd.hpp>

//#include <data/locator.hpp>
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

#include <teem/nrrd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <filesystem>
namespace fs = std::filesystem;

using namespace spurt::lavd;

std::ofstream log_file;
spurt::log::dual_ostream spurt::lavd::_log_(log_file, std::cout, 1, 0, 0, true);

std::array<size_t, 2> res({ 512, 512 });
std::array<double, 4> bnds;

double gaussian_kernel(vec2 point1, vec2 point2, double h) {
	double squared_distance = (point1[0] - point2[0]) * (point1[0] - point2[0]) + (point1[1] - point2[1]) * (point1[1] - point2[1]);
	return std::exp(-squared_distance / (h*h)); //Gaussian kernel
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

double rbf(const std::vector<vec2> all_points, const std::vector<double> all_values, const vec2 new_point) {
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

void convert_file(const char* filename) {
	printf("Normalizing file: %s\n", filename);

	//Read in data from the file
	Nrrd* nin = spurt::nrrd_utils::readNrrd(filename);
	spurt::nrrd_utils::nrrd_traits traits(nin);
	std::vector<double> data;
	spurt::nrrd_utils::to_vector(data, nin);

	//Sort data by meaning
	std::vector<vec2> pos;
	std::vector<double> val;
	for (int i = 0; i < data.size(); i += 3) {
		pos.push_back(nvis::vec2(data[i], data[i + 1]));
		val.push_back(data[i + 2]);
	}

	//Create normalized field
	const std::vector<std::string>& comments = traits.comments();
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

	double stepx = (bnds[2] - bnds[0]) / (res[0] - 1);
	double stepy = (bnds[3] - bnds[1]) / (res[1] - 1);

	//printf("%lf x %lf step size.\n", stepx, stepy);
	//exit(0);

	std::vector<vec2> normalized_pos;
	for (int i = 0; i < res[0]; i++) {
		for (int j = 0; j < res[1]; j++) {
			normalized_pos.push_back(nvis::vec2(bnds[0] + stepx*i, bnds[1] + stepy*j));
			//printf("%lf %lf\n", normalized_pos.back()[0], normalized_pos.back()[1]);
		}
	}

	//Generate values for each point
	std::vector<double> normalized_vals(normalized_pos.size());

	spurt::ProgressDisplay progress(false), total_progress(false);
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
			if (!thread) progress.update(i);

			normalized_vals[i] = rbf(pos, val, normalized_pos[i]);
		}
	}
	progress.end();

	//Write to NRRD
	std::vector<double> out_data(normalized_vals.size() * 3);
	for (int i = 0; i < normalized_vals.size(); i++) {
		out_data[i * 3] = normalized_pos[i][0];
		out_data[i * 3 + 1] = normalized_pos[i][1];
		out_data[i * 3 + 2] = normalized_vals[i];
	}

	Nrrd* nrrd = nrrdNew();
	size_t size[2] = { normalized_vals.size(), 3 };
	if (nrrdWrap_nva(nrrd, out_data.data(), nrrdTypeDouble, 2, size)) {
		std::cout << "Error wrapping NRRD: " << biffGetDone(NRRD) << std::endl;
	}
	std::string name = std::string("normalized/reshaped_") + (filename + strlen(filename) - 11);
	if (nrrdSave(name.c_str(), nrrd, NULL)) {
		std::cerr << "Error writing NRRD: " << biffGetDone(NRRD) << std::endl;
	}

	nrrdNix(nrrd);
	std::cout << "Wrote NRRD file: " << name.c_str() << std::endl;
}

int main(int argc, const char* argv[]) {

	std::string folderPath = "./";
	std::string keyword = "neighbor_deletion";
	for (const auto& entry : fs::directory_iterator(folderPath)) {
		if (entry.is_regular_file()) {
			std::string filename = entry.path().filename().string();

			//Skip files that are finished
			bool skip = false;
			char hours[5];
			for (int i = 0; i < 5; i++) {
				hours[i] = filename.c_str()[strlen(filename.c_str()) - 11 + i];
			}
			std::string outputPath = "./normalized/";
			for (const auto& entry : fs::directory_iterator(outputPath)) {
				if (entry.is_regular_file()) {
					std::string outfile = entry.path().filename().string();

					if (outfile.find(hours) != std::string::npos) {
						skip = true;
						break;
					}
				}
			}
			if (skip) {
				printf("Skipping file %s due to overlapping hour count\n", filename.c_str());
				continue;
			}

			if (filename.find(keyword) != std::string::npos && filename.find(".nrrd") != std::string::npos) {
				printf("\n");
				convert_file(filename.c_str());
				printf("\n");
			}
		}
	}
}