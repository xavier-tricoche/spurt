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

std::string inputFolder;
std::string outputPath;
std::string outputName = "compressed";
std::string keyword = "reshaped";

std::array<size_t, 2> res({ 1024, 1024 });	//Resolution of the reconstructed dataset
std::array<double, 4> bnds;


void convert_file(const char* filename) {
	printf("Compressing file: %s\n", filename);
	Nrrd* nin = spurt::nrrd_utils::readNrrd((inputFolder + filename).c_str());
	spurt::nrrd_utils::nrrd_traits traits(nin);
	std::vector<double> data;
	spurt::nrrd_utils::to_vector(data, nin);

	//Sort data by meaning
	//std::vector<vec2> pos;
	std::vector<double> val;
	for (int i = 0; i < data.size(); i += 3) {
		//pos.push_back(nvis::vec2(data[i], data[i + 1]));
		val.push_back(data[i + 2]);
		//printf("%lf\n", data[i + 2]);
	}

	//Gather other information

	double stepx = 0.017302064718493507;
	double stepy = 0.012453562120307489;

	Nrrd* nrrd = nrrdNew();
	size_t size[2] = { 1024, 1024 };
	if (nrrdWrap_nva(nrrd, val.data(), nrrdTypeDouble, 2, size)) {
		std::cout << "Error wrapping NRRD: " << biffGetDone(NRRD) << std::endl;
	}
	nrrd->spaceDim = 2;

	for (unsigned int i = 0; i < 2; ++i) {
		nrrd->axis[i].spaceDirection[0] = 0.0;
		nrrd->axis[i].spaceDirection[1] = 0.0;
		//nrrd->axis[i].spaceDirection[2] = 0.0;
	}

	nrrd->axis[0].spaceDirection[0] = stepx;
	nrrd->axis[1].spaceDirection[1] = stepy;

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
		parser.add_value("input_directory", inputFolder, "Input directory", required);
		//parser.add_value("output_directory", outputPath, "Output directory", optional);
		//parser.add_value("input_keyword", keyword, "Keyword for input files", optional);
		//parser.add_value("mode", mode, "Reconstruction method", optional);	//t for kdtree, b for bucket, i for iterative
		parser.parse(argc, argv);
	}
	catch (std::runtime_error& err) {
		printf("An error has occurred while parsing input arguments.\n");
		exit(1);
	}

	outputPath = inputFolder;

	for (const auto& entry : fs::directory_iterator(inputFolder)) {
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

					if (outfile.find(hours) != std::string::npos && outfile.compare(filename) != 0) {
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