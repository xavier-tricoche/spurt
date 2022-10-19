#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <iostream>
#include <boost/filesystem.hpp>
#include <VTK/vtk_utils.hpp>
#include <format/filename.hpp>
#include <misc/strings.hpp>


namespace xcl = xavier::command_line;
typedef boost::filesystem::path path_t;


int main(int argc, const char* argv[]) {
    std::string in_nrrd_name;
    std::string out_vtk_name;
    int verbose=0;
    
    std::string exec_name=path_t(argv[0]).filename().string();
    
    xcl::option_traits 
        required_group(true, false, "Required parameters"),
        positional_group(true, true, "Positional parameters"),
        optional_group(false, false, "Optional parameters");
        
    xcl::option_parser parser(argv[0], 
        "Export contents of VTK file to a NRRD file");
        
    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("input", in_nrrd_name, "Input filename (.nrrd)", 
                         positional_group);
        parser.add_value("output", out_vtk_name, "Output filename (.vtk, .vti)",
                         optional_group);
        parser.add_value("verbose", verbose, verbose,
                         "Verbose level", optional_group);
        parser.parse(argc, argv);
    }
    catch (std::exception& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
	
	if (out_vtk_name.empty()) {
		out_vtk_name = xavier::filename::replace_extension(in_nrrd_name, "vti");
	}
	
	vtkStructuredPoints* spts;
	try {
		vtkStructuredPoints* spts = vtk_utils::load_nrrd(in_nrrd_name);
		vtk_utils::saveVTK(spts, out_vtk_name);
	}
	catch(std::runtime_error& e) {
		std::cerr << "exception caught: " << e.what() << '\n';
		exit(0);
	}
	
    return 0;
}
