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

std::string folderPath = "./";		//Folder where the data is stored
std::string keyword = "input";	//Unique name for data items
std::string outputName = "reshaped";
std::string outputPath = "./normalized/";	//Folder where the output data is stored
char mode = 't';					//Reconstruction methodology. t for kd tree, i for iterative (brute force)

double stepx;
double stepy;

double gaussian_kernel(vec2 point1, vec2 point2, double h) {
	double squared_distance = (point1[0] - point2[0]) * (point1[0] - point2[0]) + (point1[1] - point2[1]) * (point1[1] - point2[1]);
	return std::exp(-squared_distance / (h * h)); //Gaussian kernel
}

double bilateral_gaussian_kernel(double val1, double val2, double h) {
	double squared_distance = (val1 - val2) * (val1 - val2);
	return std::exp(-squared_distance / (h * h)); //Gaussian kernel
}

double gaussian_kernel(double dist, double h) {
	return std::exp(-(dist * dist) / (h * h));
}

double blur(const spurt::point_locator<double, double, 2>& kd_tree, const vec2 new_point, double radius) {
	const double RANGE = radius;
	const double KERNEL = radius / 4;
	const double VALUE_KERNEL = 1000;

	double value = 0.0;
	double sum = 0.0;

	using point = spurt::data_point<double, double, 2>;
	std::list<point> neighbors;
	kd_tree.find_within_range(neighbors, { new_point[0], new_point[1] }, RANGE);
	std::list<point> op;
	for (const auto& neighbor : neighbors) {
		double dist = sqrt(((neighbor.coordinate()[0] - new_point[0]) * (neighbor.coordinate()[0] - new_point[0])) + ((neighbor.coordinate()[1] - new_point[1]) * (neighbor.coordinate()[1] - new_point[1])));
		if (dist <= RANGE) {
			op.push_back(neighbor);
		}
	}
	neighbors = op;
	size_t n = neighbors.size();
	if (n == 0) {
		return 0.0;
	}

	std::vector<vec2> points(n);
	std::vector<double> values(n);
	size_t i = 0;
	double core_value = 0.0;
	for (const auto& neighbor : neighbors) {
		points[i] = neighbor.coordinate();
		values[i] = neighbor.data();
		i++;
		if (neighbor.coordinate()[0] - new_point[0] < stepx / 2 && neighbor.coordinate()[1] - new_point[1] < stepy / 2) {
			core_value = neighbor.data();
		}
	}

	for (int i = 0; i < n; i++) {
		double weight = gaussian_kernel(new_point, points[i], KERNEL);
		double bilateral_weight = bilateral_gaussian_kernel(core_value, values[i], VALUE_KERNEL);

		//Ordinary kernel - no edge saving
		//value += weight * values[i];
		//sum += weight;

		//Bilateral filter - preserves edges
		value += weight * bilateral_weight * values[i];
		sum += weight * bilateral_weight;
	}

	return value / sum;
}

//Takes an input file, reconstructs it, and saves it
void convert_file(const char* filename) {
	printf("Normalizing file: %s\n", filename);

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
		val.push_back(data[i + 2]);
	}

	std::array<size_t, 2> res({ 1024, 1024 });	//Resolution of the reconstructed dataset
	std::array<double, 4> bnds;
	res[0] = 1024;
	res[1] = 1024;
	bnds[0] = 262.05;
	bnds[1] = 18.05;
	bnds[2] = 279.75;
	bnds[3] = 30.79;
	stepx = (bnds[2] - bnds[0]) / (res[0] - 1);
	stepy = (bnds[3] - bnds[1]) / (res[1] - 1);

	std::vector<vec2> normalized_pos;
	for (int i = 0; i < res[0]; i++) {
		for (int j = 0; j < res[1]; j++) {
			//y increments, then x increments
			normalized_pos.push_back(nvis::vec2(bnds[0] + stepx * i, bnds[1] + stepy * j));
			//printf("%lf %lf\n", normalized_pos.back()[0], normalized_pos.back()[1]);
		}
	}

	//Generate values for each point
	std::vector<double> blurred_vals(val.size());

	//Create KD Tree
	using point = spurt::data_point<double, double, 2>;
	using locator = spurt::point_locator<double, double, 2>;

	std::vector<point> points;
	printf("Generating vector of %ld points\n", val.size());
	for (size_t i = 0; i < val.size(); i++) {
		points.emplace_back(point({ pos[i][0], pos[i][1] }, val[i]));
	}
	printf("Generating kd-tree from vector\n");
	locator kd_tree(points.begin(), points.end());
	printf("kd-tree created. Beginning normalization\n");

	//Generate normalization timer
	spurt::ProgressDisplay progress(false);
	progress.start(val.size(), "blur");
	progress.set_active(true);

#pragma omp parallel
	{
#pragma omp for schedule(dynamic,1)
		for (size_t i = 0; i < val.size(); i++) {
#if _OPENMP
			const int thread = omp_get_thread_num();
#else
			const int thread = 0;
#endif
			if (!thread) {
				progress.update(i);
			}

			blurred_vals[i] = blur(kd_tree, normalized_pos[i], stepx*24);
		}
	}
	progress.end();

	//Write to NRRD
	std::vector<double> out_data(blurred_vals.size() * 3);
	for (int i = 0; i < blurred_vals.size(); i++) {
		//out_data[i * 3] = normalized_pos[i][0];
		//out_data[i * 3 + 1] = normalized_pos[i][1];
		out_data[i * 3 + 2] = blurred_vals[i];
	}

	Nrrd* nrrd = nrrdNew();
	size_t size[2] = { blurred_vals.size(), 3 };
	if (nrrdWrap_nva(nrrd, out_data.data(), nrrdTypeDouble, 2, size)) {
		std::cout << "Error wrapping NRRD: " << biffGetDone(NRRD) << std::endl;
	}
	//Save to file
	std::string name = std::string(outputPath + outputName + "_") + (filename + strlen(filename) - 11);
	if (nrrdSave(name.c_str(), nrrd, NULL)) {
		std::cerr << "Error writing NRRD: " << biffGetDone(NRRD) << std::endl;
	}

	nrrdNix(nrrd);
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
		parser.add_value("mode", mode, "Reconstruction method", optional);	//t for kdtree, b for bucket, i for iterative
		parser.parse(argc, argv);
	}
	catch (std::runtime_error& err) {
		printf("An error has occurred while parsing input arguments.\n");
		exit(1);
	}

	if (mode != 't' && mode != 'i') {
		printf("Invalid mode input. Use 't' for kd tree or 'i' for iterative.\n");
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