#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <math/stat.hpp>


std::string filename;
int pct_step = 10;

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;

    xcl::option_traits
            required_group(true, false, "Required Options"),
            optional_group(false, false, "Optional Group");

    xcl::option_parser parser(argv[0],
            "Analyze value distribution in NRRD file");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", filename, "Input filename", required_group);
        parser.add_value("step", pct_step, pct_step, "Percentile step", optional_group);
		parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

double percentile(const std::vector<double>& values, double pct) {
	size_t id = static_cast<size_t>(pct*static_cast<double>(values.size()-1));
	return values[id];
}


int main(int argc, const char* argv[]) {
	initialize(argc, argv);

	Nrrd* nin = spurt::nrrd_utils::readNrrd(filename);
	std::vector<double> values;
	spurt::nrrd_utils::to_vector<double>(values, nin);

	std::sort(values.begin(), values.end());
	std::pair<double, double> meanvar = spurt::meanvariance(values);
	std::cout << "mean: " << meanvar.first << '\n';
	std::cout << "variance: " << meanvar.second << '\n';
	std::cout << "percentiles:\n";
	std::cout << "  0  %: " << values.front() << '\n';
	std::cout << "  0.1%: " << percentile(values, 0.001) << '\n';
	std::cout << "  0.5%: " << percentile(values, 0.005) << '\n';
	std::cout << "  1  %: " << percentile(values, 0.01) << '\n';

    for (int pct=pct_step; pct<99; pct += pct_step) {
        std::cout << " " <<  std::setw(2) << std::setfill(' ')
                  << pct << " %: " << percentile(values, (double)pct/100.) << '\n';
    }

	// std::cout << "  5  %: " << percentile(values, 0.05) << '\n';
	// std::cout << " 10  %: " << percentile(values, 0.1) << '\n';
	// std::cout << " 20  %: " << percentile(values, 0.2) << '\n';
	// std::cout << " 30  %: " << percentile(values, 0.3) << '\n';
	// std::cout << " 40  %: " << percentile(values, 0.4) << '\n';
	// std::cout << " 50  %: " << percentile(values, 0.5) << '\n';
	// std::cout << " 60  %: " << percentile(values, 0.6) << '\n';
	// std::cout << " 70  %: " << percentile(values, 0.7) << '\n';
	// std::cout << " 80  %: " << percentile(values, 0.8) << '\n';
	// std::cout << " 90  %: " << percentile(values, 0.9) << '\n';
	// std::cout << " 95  %: " << percentile(values, 0.95) << '\n';

	std::cout << " 99  %: " << percentile(values, 0.99) << '\n';
	std::cout << " 99.9%: " << percentile(values, 0.999) << '\n';
	std::cout << "100  %: " << values.back() << '\n';

	return 0;
}
