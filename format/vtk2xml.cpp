// STL
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
// VTK
#include <vtkDataSet.h>
#include <vtkUnstructuredGrid.h>
#include <vtkImageData.h>
#include <vtkPolyData.h>
#include <vtkRectilinearGrid.h>

#include <VTK/vtk_utils.hpp>
#include <format/filename.hpp>



std::string me;
void usage(const std::string& message="")
{
    if (!message.empty()) {
        std::cerr << "\nERROR: " << message << '\n';
    }
    std::cerr << '\n'
              << "DESCRIPTION: Convert a file in VTK legacy format to\n"
			  << "a suitable XML VTK format\n"
              << "USAGE: " << me << " <input> <output> [options]\n"
              << "PARAMETERS:\n"
              << " <input>          Input VTK legacy filename\n"
              << " <output>         Output basename\n"
              << "OPTIONS:\n"
              << " -h | --help      Print this information\n"
              << " -v | --verbose   Activate verbose mode\n"
              << std::endl;

    exit(!message.empty());
}

int main(int argc, char* argv[])
{
	bool verbose = false;
    me = argv[0];
	if (argc < 3) usage();

    std::string in_name="", out_name="";

	int positional_count = 0;
    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage();
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else if (arg[0] == '-') {
            usage("unrecognized option: " + arg);
        }
		else {
			if (++positional_count == 1) in_name = arg;
			else if (++positional_count == 2) out_name = arg;
			else usage("unexpected positional argument: " + arg);
		}
    }

	vtkSmartPointer<vtkDataSetReader> reader = vtkSmartPointer<vtkDataSetReader>::New();
	reader->SetFileName(in_name.c_str());
	reader->Update();
	vtk_utils::saveVTK_XML(reader->GetOutput(), out_name);

	return 0;
}
