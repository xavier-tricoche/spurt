#include <data/image.hpp>
#include <math/types.hpp>
#include <iostream>
#include <misc/progress.hpp>
#include <vector>
#include <string>
#include <misc/option_parse.hpp>
#include <utils/marschner_lobb.hpp>
#include <utils/abc_function.hpp>

using namespace spurt;

int main(int argc, const char* argv[])
{
    std::string dataname, outfilename;
    svec3 res;
    size_t nsamples;

    namespace cl = command_line;

    cl::option_traits
        required_group(true, false, "Required Options"),
        positional_group(true, true, "Positional Group"),
        optional_group(false, false, "Optional Group");

    cl::option_parser parser(argv[0],
        "4D Symplectic Map Visualization Using Slab Technique");

    try {
        // parser.use_default_symbol();
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("data", dataname, "Marschner", "Name of function to reconstruct", optional_group);
        parser.add_value("res", res, {64, 64, 64}, "Resolution", optional_group)
        parser.add_value("samples", nsamples, 4000, "Number of samples", optional_group);
        parser.add_value("output", outfilename, "Output filename", required_group);

        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR(1): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    catch(std::exception& e) {
        std::cerr << "ERROR(2): " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    
    vec3 min(-1, -1, -1);
    

    
    return 0;
}