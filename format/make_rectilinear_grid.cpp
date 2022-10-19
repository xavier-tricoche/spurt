#include <VTK/vtk_utils.hpp>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include "teem/nrrd.h"


int verbose=0;
std::string name_x, name_y, name_z, name_data, name_out;

VTK_SMART(vtkDoubleArray) import_scalar_array(const std::string& nrrd_filename) {
    Nrrd* nin = xavier::nrrd_utils::readNrrd(nrrd_filename);
    size_t n = nin->axis[0].size;
    VTK_CREATE(vtkDoubleArray, array);
    array->SetNumberOfTuples(n);
    array->SetNumberOfComponents(1);
    xavier::nrrd_utils::nrrd_data_wrapper<double> wrapper(nin);
    for (size_t i=0; i<n; ++i) {
        array->SetTuple1(i, wrapper[i]);
    }
    return array;
}

int main(int argc, const char* argv[]) {

    std::string name_id, name_coef, name_rhs;
    namespace xcl = xavier::command_line;
    int verbose=0;

    xcl::option_traits
        required_group(true, false, "Required parameters"),
        positional_group(true, true, "Positional parameters"),
        optional_group(false, false, "Optional parameters");

    xcl::option_parser parser(argv[0],
        "Create a rectilinear grid from several nrrd files");

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("x", name_x, "X coordinates filename (.nrrd)",
                         positional_group);
        parser.add_value("y", name_y, "Y coordinates filename (.nrrd)",
                         positional_group);
        parser.add_value("z", name_z, "Z coordinates filename (.nrrd)",
                         positional_group);
        parser.add_value("data", name_data, "data filename (.nrrd)",
                         positional_group);
        parser.add_value("output", name_out, "Output filename (.vtr)",
                         positional_group);
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
    
    VTK_SMART(vtkDoubleArray) xs = import_scalar_array(name_x);
    VTK_SMART(vtkDoubleArray) ys = import_scalar_array(name_y);
    VTK_SMART(vtkDoubleArray) zs = import_scalar_array(name_z);
    
    VTK_CREATE(vtkRectilinearGrid, grid);
    grid->SetDimensions(xs->GetNumberOfTuples(), 
                        ys->GetNumberOfTuples(), 
                        zs->GetNumberOfTuples());
    grid->SetXCoordinates(xs);
    grid->SetYCoordinates(ys);
    grid->SetZCoordinates(zs);
    
    VTK_SMART(vtkDoubleArray) data = import_scalar_array(name_data);
    data->SetName(xavier::filename::remove_extension(name_data).c_str());
    grid->GetPointData()->SetScalars(data);
    
    VTK_CREATE(vtkXMLRectilinearGridWriter, writer);
    writer->SetFileName(name_out.c_str());
    writer->SetInputData(grid);
    writer->Write();
    
    exit(0);
}
