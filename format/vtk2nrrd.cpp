#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <iostream>
#include <boost/filesystem.hpp>
#include <VTK/vtk_utils.hpp>
#include <format/filename.hpp>
#include <misc/strings.hpp>

#include <vtkImageIterator.h>

namespace xcl = xavier::command_line;
typedef boost::filesystem::path path_t;

template<typename T>
void save_nrrd_image(vtkImageData* img, size_t* dims, int* ext, const std::string& name) {
	typedef T data_type;
	data_type* data = (data_type*)(calloc(dims[0]*dims[1]*dims[2]*dims[3], sizeof(data_type)));
	data_type* out_ptr = data;
	vtkImageIterator<data_type> iter(img, ext);
	while (!iter.IsAtEnd()) {
		for (auto span_iter=iter.BeginSpan(); span_iter!=iter.EndSpan();
		     *out_ptr++ = *span_iter++) {};
    	iter.NextSpan();
	}

	int idx = xavier::nrrd_utils::nrrd_value_traits_from_type<data_type>::index;
	xavier::nrrd_utils::writeNrrd(data, name, idx, (dims[0]>1 ? 4 : 3), (dims[0]>1 ? dims : &dims[1]));
}

int main(int argc, const char* argv[]) {
    std::string input_name;
    std::string output_name;
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
        parser.add_value("input", input_name, "Input filename (.vtk, .vti)",
                         positional_group);
        parser.add_value("output", output_name, "Output filename (.nrrd)",
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

	vtkDataSet* ds;
	try {
		ds = vtk_utils::readVTK(input_name);
	}
	catch(std::runtime_error& e) {
		std::cerr << "exception caught: " << e.what() << '\n';
		exit(0);
	}

    vtkImageData* img = vtkImageData::SafeDownCast(ds);
	if (!img) {
		std::cerr << "input dataset type is incompatible with NRRD format\n";
		exit(-1);
	}

	size_t res[4];
	size_t nvals = img->GetNumberOfScalarComponents();
	int* dims = img->GetDimensions();
	res[0] = nvals;
	res[1] = dims[0];
	res[2] = dims[1];
	res[3] = dims[2];

	int* ext = img->GetExtent();
	int type = img->GetScalarType();

	switch (type) {
		case VTK_CHAR:
		case VTK_SIGNED_CHAR:    save_nrrd_image<char>(img, res, ext, output_name); break;
		case VTK_UNSIGNED_CHAR:  save_nrrd_image<unsigned char>(img, res, ext, output_name); break;
		case VTK_SHORT:          save_nrrd_image<short>(img, res, ext, output_name); break;
		case VTK_UNSIGNED_SHORT: save_nrrd_image<unsigned short>(img, res, ext, output_name); break;
		case VTK_INT:            save_nrrd_image<int>(img, res, ext, output_name); break;
		case VTK_UNSIGNED_INT:   save_nrrd_image<unsigned int>(img, res, ext, output_name); break;
		case VTK_LONG:           save_nrrd_image<long>(img, res, ext, output_name); break;
		case VTK_UNSIGNED_LONG:  save_nrrd_image<unsigned long>(img, res, ext, output_name); break;
		case VTK_FLOAT:          save_nrrd_image<float>(img, res, ext, output_name); break;
		case VTK_DOUBLE:         save_nrrd_image<double>(img, res, ext, output_name); break;
		default: {
			std::cerr << "Unrecognized / unsupported VTK data type\n";
			exit(-1);
		}
	}
    return 0;
}
