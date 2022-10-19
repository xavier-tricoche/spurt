#include "boundary_extraction.hpp"
#include <VTK/vtk_utils.hpp>
#include <misc/option_parse.hpp>

std::string grid_name, data_name;
std::string bound_name, offset_name;
bool verbose = false;

typedef xavier::mesh mesh_type;
typedef mesh_type::face face_type;
typedef mesh_type::index_type index_type;
typedef mesh_type::vertex_type vertex_type;

void initialize(int argc, const char* argv[])
{
    namespace xcl = xavier::command_line;

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
        parser.add_value("data", data_name, "Data filename", optional);
        parser.add_value("bound", bound_name, "Boundary basename", required);
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

    mesh_type m;
    if (verbose) {
        std::cout << "Importing " << grid_name << "... " << std::flush;
    }
    m.load_mesh(grid_name);
    if (verbose) {
        std::cout << "done\n";
        std::cout << "Extracting boundary... " << std::flush;
    }
    m.extract_boundary();
    if (verbose) {
        std::cout << "done\n";
    }

    // save offset surface
    if (verbose) {
        std::cout << "Creating and saving offset geometry in VTK file... " << std::flush;
    }
    std::vector<mesh_type::face> offset_cells = m.get_offset();
    VTK_SMART(vtkPolyData) offset = vtk_utils::make_points(m.get_vertices());
    VTK_CREATE(vtkCellArray, cells2);
    size_t counter=0;
    std::cout << "There are " << offset_cells.size() << " cells on offset surface\n";
    std::for_each(offset_cells.begin(), offset_cells.end(), [&]
        (const face_type& f) {
            cells2->InsertNextCell(f.size());
            std::for_each(f.begin(), f.end(), [&](index_type i) {
               cells2->InsertCellPoint(i);
               // if (i<0 || i>=npts) {
               //     std::cerr << "cell #" << counter << ": Invalid vertex id in cell definition: " << i << '\n';
               // }
            });
            counter++;
        });
        std::cout << "cells contains " << cells2->GetNumberOfCells() << " cells\n";
    // std::cout << "\n\nbefore assigning cells\n";
    // offset->PrintSelf(std::cout, vtkIndent(0));
    offset->SetPolys(cells2);
    std::cout << "offset contains " << offset->GetNumberOfCells() << " cells\n";
    // std::cout << "\n\nafter assigning cells\n";
    // offset->PrintSelf(std::cout, vtkIndent(0));
    // offset->ComputeBounds();
    // std::cout << "\n\nafter recomputing bounds\n";
    // offset->PrintSelf(std::cout, vtkIndent(0));
    double b[6];
    offset->GetBounds(b);
    std::cout << "\n\nvertices have bounds: [" << b[0] << "," << b[1] << "," << b[2] << "," << b[3] << "," << b[4] << "," << b[5] << "]";
    vtk_utils::saveVTK(offset, bound_name + "_offset_geometry.vtp");
    if (verbose) {
        std::cout << "done\n";
    }

    // save boundary
    if (verbose) {
        std::cout << "Creating and saving boundary geometry in VTK file... " << std::flush;
    }
    VTK_SMART(vtkPolyData) boundary = vtk_utils::make_points(m.get_vertices());
    std::vector<mesh_type::face> boundary_cells = m.get_boundary();
    VTK_CREATE(vtkCellArray, cells);
    std::for_each(boundary_cells.begin(), boundary_cells.end(), [&]
        (const xavier::mesh::face& f) {
            cells->InsertNextCell(f.size());
            std::for_each(f.begin(), f.end(), [&](index_type i) {
               cells->InsertCellPoint(i);
            });
        });
    boundary->SetPolys(cells);
    if (!data_name.empty()) {
        m.load_data(data_name);
        std::vector<mesh_type::vec3_type> shearstress = m.get_shear_stress();
        vtk_utils::add_vectors(boundary, shearstress, true, "velocity", true);
    }
    else {
        vtk_utils::add_vectors(boundary, m.get_normals(), true, "normals", true);
    }
    vtk_utils::saveVTK(boundary, bound_name + "_boundary_geometry.vtp");
    if (verbose) {
        std::cout << "done\n";
    }

    m.save_boundary(bound_name);

    return 0;
}
