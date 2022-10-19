// STL
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <vtk/vtk_utils.hpp>
#include "dlr2vtk.hpp"

std::string me;
void usage(const std::string& message="")
{
    if (!message.empty()) {
        std::cerr << "\nERROR: " << message << '\n';
    }
    std::cerr << '\n'
              << "DESCRIPTION: Convert a DLR dataset from NetCDF format to\n"
              << "VTK format. The conversion can be restricted to the mesh\n"
              << "and / or to the simulation boundary. In the latter case,\n"
              << "the data can be segmented by connectivity, in which case\n"
              << "separate files will be generated for each connected component.\n"
              << "USAGE: " << me << " [parameters] [options]\n"
              << "PARAMETERS:\n"
              << " -g | --grid <string>       DLR grid filename\n"
              << " -o | --output <string>     Output VTK basename\n"
              << "OPTIONS:\n"
              << " -h | --help                Print this information\n"
              << " -d | --data <string>       DLR data filename\n"
              << " -a | --attribute <string>  Name of requested attribute\n"
              << " -b | --boundary            Restrict conversion to grid boundary\n"
              << " -m | --mesh                Restrict conversion to mesh information\n"
              << " -v | --verbose             Activate verbose mode\n"
              << std::endl;

    exit(!message.empty());
}

int main(int argc, char* argv[])
{
    me = argv[0];
    std::string grid_name="", val_name="", out_name="";
    std::vector<std::string> att;
    bool boundary_only = false;
    bool mesh_only = false;
    bool verbose = false;

    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage();
        } else if (arg == "-g" || arg == "--grid") {
            if (i == argc-1) {
                usage("missing grid file name");
            }
            grid_name = argv[++i];
        } else if (arg == "-d" || arg == "--data") {
            if (i == argc-1) {
                usage("missing data file name");
            }
            val_name = argv[++i];
        } else if (arg == "-a" || arg == "--attribute") {
            if (i == argc-1) {
                usage("missing attribute name");
            }
            att.push_back(argv[++i]);
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                usage("missing output file name");
            }
            out_name = argv[++i];
        } else if (arg == "-b" || arg == "--boundary") {
            boundary_only = true;
        } else if (arg == "-m" || arg == "--mesh") {
            mesh_only = true;
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else {
            usage("unrecognized argument: " + arg);
        }
    }

    // sanity check
    if (grid_name.empty()) {
        usage("missing grid file");
    }
    if (out_name.empty()) {
        usage("missing output file name");
    }
    if (val_name.empty() && att.size())
        usage("both value(s) and geometry were requested\n"
              "yet file containing values is missing");
    
    spurt::dlr2vtk_converter converter(grid_name, val_name);
    converter.set_attributes(att);
    converter.set_mesh_only(mesh_only);
    converter.set_boundary_only(boundary_only);
    converter.set_verbose(verbose);
    
    converter.import();
    VTK_SMART(vtkUnstructuredGrid) grid = converter.get_vtk_dataset();
    vtk_utils::saveVTK(grid, out_name);
    std::cout << out_name << " has been exported\n";
    return 0;
}
