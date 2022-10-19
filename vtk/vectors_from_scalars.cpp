#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include <misc/option_parse.hpp>
#include <vtkDataSet.h>
#include <vtkSmartPointer.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkDataSetWriter.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>


std::string input_name, output_name;

// void initialize(int argc, char* argv[]) {
//     namespace xcl = xavier::command_line;
//
//     xcl::option_traits
//             required_group(true, false, "Required Options"),
//             optional_group(false, false, "Optional Group");
//     xcl::option_parser parser(argv[0],
//             "Merge 3 scalar fields to form a vector field");
//
//     try {
//         parser.use_short_symbols(false);
//         parser.use_brackets(true);
//         parser.add_value("input", input_name, "Input filename", required_group);
//         parser.add_value("output", output_name, "Output filename", optional_group);
//         parser.add_tuple<3>("scalars", scalar_names, "Scalar names", optional_group);
//
//         parser.parse(argc, const_cast<const char**>(argv));
//     }
//     catch(std::runtime_error& e) {
//         std::cerr << "ERROR: " << argv[0] << " threw exception: "
//                   << e.what() << "\n"
//                   << "Command line options entered so far:\n"
//                   << parser.print_self(false, true, false) << "\n\n\n";
//         exit(1);
//     }
// }

int main(int argc, char* argv[]) {
    
    vtkSmartPointer<vtkXMLUnstructuredGridReader> reader = vtkSmartPointer<vtkXMLUnstructuredGridReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();
	
	vtkUnstructuredGrid* data = vtkUnstructuredGrid::SafeDownCast(reader->GetOutput());
	
	size_t npts = data->GetNumberOfPoints();
	
	vtkSmartPointer<vtkFloatArray> vecs = vtkSmartPointer<vtkFloatArray>::New();
	vecs->SetNumberOfComponents(3);
	vecs->SetNumberOfTuples(npts);
	
	data->GetPointData()->SetActiveAttribute("x_array", vtkDataSetAttributes::SCALARS);
	vtkDataArray* scl1 = reader->GetOutput()->GetPointData()->GetScalars();
	data->GetPointData()->SetActiveAttribute("y_array", vtkDataSetAttributes::SCALARS);
	vtkDataArray* scl2 = reader->GetOutput()->GetPointData()->GetScalars();
	data->GetPointData()->SetActiveAttribute("z_array", vtkDataSetAttributes::SCALARS);
	vtkDataArray* scl3 = reader->GetOutput()->GetPointData()->GetScalars();
	
	for (size_t i=0; i<npts; ++i) {
		double x = scl1->GetTuple1(i);
		double y = scl2->GetTuple1(i);
		double z = scl3->GetTuple1(i);
		vecs->SetTuple3(i, x, y, z);
	}
	
	data->GetPointData()->SetVectors(vecs);
	
	vtkSmartPointer<vtkDataSetWriter> writer = vtkSmartPointer<vtkDataSetWriter>::New();
	writer->SetInputData(data);
	writer->SetFileName(argv[2]);
	writer->Write();
    
    return 0;
}
