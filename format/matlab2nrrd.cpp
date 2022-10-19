#include <matio.h>
#include <image/nrrd_wrapper.hpp>
#include <misc/option_parse.hpp>
#include <iostream>
#include <boost/filesystem.hpp>

namespace xcl = xavier::command_line;
typedef boost::filesystem::path path_t;

int main(int argc, const char* argv[]) {
    std::string input_name;
    std::string output_basename;
    int verbose=0;
    
    std::string exec_name=path_t(argv[0]).filename().string();
    
    xcl::option_traits 
        required_group(true, false, "Required parameters"),
        positional_group(true, true, "Positional parameters"),
        optional_group(false, false, "Optional parameters");
        
    xcl::option_parser parser(argv[0], 
        "Export contents of Matlab file to one or several NRRD files");
        
    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("input", input_name, "Input filename", 
                         positional_group);
        parser.add_value("output", output_basename, "Output filename",
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
    
    path_t p(output_basename);
    output_basename=path_t(output_basename).replace_extension().string();
    
    if (verbose) std::cout << "Reading " << input_name << '\n';
    
    mat_t *mat_file=Mat_Open(input_name.c_str(), MAT_ACC_RDONLY);
    if (mat_file==NULL) {
        std::cerr << argv[0] << ": Unable to open " << input_name << '\n';
        exit(1);
    }
    
    matvar_t *variable;
    while ((variable=Mat_VarReadNext(mat_file))!=NULL) {
        std::string name;
        int type;
        int dim;
        size_t *sizes;
        void *data;
        name=variable->name;
        dim=variable->rank;
        sizes=variable->dims;
        data=variable->data;
        if (verbose) std::cout << "\nCurrent variable is " << name << '\n';
        switch(variable->data_type) {
            case MAT_T_INT8:   type=nrrdTypeShort;  break;
            case MAT_T_UINT8:  type=nrrdTypeUShort; break;
            case MAT_T_INT16:  type=nrrdTypeInt;    break;
            case MAT_T_UINT16: type=nrrdTypeUInt;   break;
            case MAT_T_INT32:  type=nrrdTypeLLong;  break;
            case MAT_T_UINT32: type=nrrdTypeULLong; break;
            case MAT_T_SINGLE: type=nrrdTypeFloat;  break;
            case MAT_T_DOUBLE: type=nrrdTypeDouble; break;
            default: {
                std::cerr << "Data type: " << variable->data_type 
                          << " not supported. Skipping...\n";
                continue;
            }
        }
        
        std::string filename=output_basename+'_'+name+".nrrd";
        Nrrd* nout=nrrdNew();
        if (nrrdWrap_nva(nout, data, type, dim, sizes)) {
            std::cerr << xavier::nrrd_utils::error_msg(exec_name, "NRRD wrapping failed") 
                      << '\n';
            exit(1);
        }
        nrrdCommentAdd(nout, 
            ("variable "+name+" extracted by "+exec_name+
             " from "+input_name).c_str());
        if (nrrdSave(filename.c_str(), nout, NULL)) {
            std::cerr << xavier::nrrd_utils::error_msg(exec_name, "NRRD export failed");
            exit(1);
        }
        std::cout << name << " saved to " << filename << '\n';
    }
    return 0;
}
