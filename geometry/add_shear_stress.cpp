#include "boundary_extraction.hpp"
#include <vtk/vtk_utils.hpp>
#include <misc/option_parse.hpp>
#include <image/nrrd_wrapper.hpp>
#include <format/dlr_reader.hpp>

std::string grid_name, data_name, out_name, boundary_name;
bool verbose = false;
bool include_mesh = false;

typedef spurt::mesh mesh_type;
typedef mesh_type::face face_type;
typedef mesh_type::index_type index_type;
typedef mesh_type::vertex_type vertex_type;

void initialize(int argc, const char* argv[])
{
    namespace xcl = spurt::command_line;

    std::string cmdline = "Command line: " + std::string(argv[0]);
    for (int i=1; i<argc; i++) {
        cmdline += " " + std::string(argv[i]);
    }

    xcl::option_traits
        required(true, false, "Required Options"),
    optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
    "Extract boundary and offset surface from DLR dataset");

    verbose = false;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("grid", grid_name, "Grid filename", required);
        parser.add_value("boundary", boundary_name, "Boundary info basename", required);
        parser.add_value("data", data_name, "Data filename", required);
        parser.add_value("out", out_name, "Output filename", required);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional);

        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n"
            << e.what() << "\n"
                << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

int main(int argc, const char* argv[]) {

    initialize(argc, argv);
    std::ostringstream devnull;
    std::streambuf *cout_save;

    if (!verbose) {
        cout_save = std::cout.rdbuf();
        std::cout.rdbuf(devnull.rdbuf());
    }

    mesh_type m;

    std::cout << "Importing " << grid_name << "... " << std::flush;
    m.load_mesh(grid_name);
    std::cout << "done\n";
    std::cout << "Importing " << data_name << "... " << std::flush;
    m.load_data(data_name);
    std::cout << "done\n";
    std::cout << "Importing boundary data..." << std::flush;
    m.load_boundary(boundary_name, verbose);
    // std::cout << "done

    std::cout << "Computing shear stress values on boundary..." << std::flush;
    std::vector<mesh_type::vertex_type> newfield = m.get_shear_stress(verbose);
    std::cout << " done\n";

    std::cout << "Merging velocity and shear stress vector fields..." << std::flush;
    const std::vector<mesh_type::vertex_type>& velocity = m.get_velocity();
    for (size_t i=0; i<newfield.size(); ++i) {
        const mesh_type::vertex_type& v = newfield[i];
        if (spurt::norm(v) == 0) {
            newfield[i] = velocity[i];
        }
    }
    std::cout << " done\n";

    std::cout << "Saving to NRRD file..." << std::flush;
    size_t dims[2] = { 3, newfield.size() };
    spurt::nrrd_utils::writeNrrd(&newfield[0], out_name, nrrdTypeFloat, 2, dims);
    std::cout << " done\n";

    if (!verbose) {
        std::cout.rdbuf(cout_save);
    }
    return 0;
}
