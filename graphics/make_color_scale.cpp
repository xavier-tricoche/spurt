#include <string>
#include <iostream>
#include <fstream>

#include <misc/option_parse.hpp>
#include <graphics/colors.hpp>
#include <math/fixed_vector.hpp>


std::string name_out;
std::string scale_name = "spiral:1";
int number = 10;
int order = 1;
float r=1;

void initialize(int argc, const char* argv[])
{
    namespace xcl = xavier::command_line;

    xcl::option_traits
            required_group(true, false, "Required Options"),
            optional_group(false, false, "Optional Group");

    xcl::option_parser parser(argv[0],
            "Create a color scale");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("output", name_out, "Output base name", optional_group);
		parser.add_value("scale", scale_name, scale_name, "Color scale name: \"b2y\" (blue to yellow),"
						 "\"r2g\" (red to green), \"spiral\" (spiral, r: number of turns)", 
				         optional_group);
		parser.add_value("number", number, number, "Number of colors to create", optional_group);
		parser.add_value("order", order, order, "Scale variation order", optional_group);

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

int main(int argc, const char* argv[]) {
	initialize(argc, argv);
	
	if (name_out.empty()) name_out = "color_scale.txt";
	if (scale_name.empty()) scale_name = "spiral";
	
	size_t sep = scale_name.find(':');
	if (sep != std::string::npos) {
		r = std::stof(scale_name.substr(sep+1));
		scale_name = scale_name.substr(0, sep);
	}
	
	std::vector<nvis::fvec3> colors;
	if (scale_name == "spiral") {
		xavier::spiral_scale(colors, number, 0.2, r, 1, 1, order);
	}
	else if (scale_name == "b2y" || scale_name == "r2g") {
		nvis::fvec3 col0, col1;
		if (scale_name == "b2y") {
			col0 = xavier::blue;
			col1 = xavier::yellow;
		}
		else {
			col0 = xavier::red;
			col1 = xavier::green;
		}
		for (int i=0; i<number; ++i) {
			double u = pow((double)i/(double)(number-1), order);
			colors.push_back((1.-u)*col0 + u*col1);
		}
	}
	
	std::fstream output(name_out.c_str(), std::ios::out);
	for (int i=0; i<number; ++i) {
		output << colors[i][0] << " " << colors[i][1] << " " << colors[i][2] << '\n';
	}
	output.close();
	
	return 1;
}


