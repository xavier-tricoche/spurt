#include <array>
#include <iomanip>

#include <misc/progress.hpp>

#include <image/nrrd_wrapper.hpp>
#include <image/nrrd_matrix_manip.hpp>
#include <misc/option_parse.hpp>
#include <format/filename.hpp>
#include <misc/strings.hpp>

bool show_progress=false;
bool verbose=false;
bool ok=true;
bool is_sym=false;
std::vector<std::string> input_file_names;
std::vector<std::string> output_file_names;
std::vector<int> dims;
std::string opname="none";
std::string tname="implicit";

using namespace spurt;

template<typename T> 
void transposer(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input, 
                const std::vector<int>& dims) {
        assert(input.size()==1 && dims.size()==2);
        output.resize(1);
        output[0]=spurt::nrrd_manip::transpose<T>(input[0], dims[0], dims[1], show_progress);
};

template<typename T> 
void multiplier(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                const std::vector<int>& dims) {
    assert(input.size()==2 && dims.size()==3);
    output.resize(1);
    output[0]=spurt::nrrd_manip::product<T>(input[0], input[1], dims[0], dims[1], dims[2], show_progress);
}

template<typename T> 
void addition(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
              const std::vector<int>& dims) {
    assert(input.size()==2 && dims.size()==2);
    output.resize(1);
    output[0]=spurt::nrrd_manip::sum<T>(input[0], input[1], dims[0], dims[1], show_progress);
}

template<typename T> 
void subtraction(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                 const std::vector<int>& dims) {
    assert(input.size()==2 && dims.size()==2);
    output.resize(1);
    output[0]=spurt::nrrd_manip::subtraction<T>(input[0], input[1], dims[0], dims[1], show_progress);
}

template<typename T> 
void SVDecomposer(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                  const std::vector<int>& dims) {
    assert(input.size()==1 && dims.size()==2);
    output.resize(3);
    spurt::nrrd_manip::SVD<T>(output[0], output[1], output[2], input[0], dims[0], dims[1], show_progress);
}

template<typename T>
void Eigendecomposer(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                     const std::vector<int>& dims) {
    assert(input.size()==1 && 
           (dims.size()==1 || (dims.size()==2 && dims[0]==dims[1]))); 
    output.resize(2);
    spurt::nrrd_manip::Eigendecomposition<T>(output[0], output[1], input[0], dims[0], is_sym, show_progress);
}

template<typename T> 
void transpose_multiplier(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                          const std::vector<int>& dims) {
    assert(input.size()==2 && dims.size()==3);
    output.resize(1);
    output[0]=spurt::nrrd_manip::trans_product<T>(input[0], input[1], dims[0], dims[1], dims[2], show_progress);
}

template<typename T> 
void CG(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
        const std::vector<int>& dims) {
    assert(input.size()==1 && 
           (dims.size()==1 || (dims.size()==2 && dims[0]==dims[1])));
    output.resize(1);
    output[0]=spurt::nrrd_manip::cauchy_green<T>(input[0], dims[0], show_progress);
}

template<typename T> 
void FTLE(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
        const std::vector<int>& dims) {
    assert(input.size()==1 && 
           (dims.size()==1 || (dims.size()==2 && dims[0]==dims[1])));
    output.resize(1);
    output[0]=spurt::nrrd_manip::FTLE<T>(input[0], dims[0], show_progress);
}

template<typename T>
void inverter(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
              const std::vector<int>& dims) {
    assert(input.size()==1 && (dims.size()==1 || (dims.size()==2 && dims[0]==dims[1])));
    output.resize(1);
    // output[0]=spurt::nrrd_manip::invert<T>(input[0], dims[0], show_progress);
    output[0]=spurt::nrrd_manip::inverse<T>(input[0], dims[0], show_progress);
}

template<typename T>
void determinant(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
              const std::vector<int>& dims) {
    assert(input.size()==1 && (dims.size()==1 || (dims.size()==2 && dims[0]==dims[1])));
    output.resize(1);
    // output[0]=spurt::nrrd_manip::invert<T>(input[0], dims[0], show_progress);
    output[0]=spurt::nrrd_manip::determinant<T>(input[0], dims[0], show_progress);
}

template<typename T>
void trace(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
           const std::vector<int>& dims) {
    assert(input.size()==1 && (dims.size()==1 || (dims.size()==2 && dims[0]==dims[1])));
    output.resize(1);
    // output[0]=spurt::nrrd_manip::invert<T>(input[0], dims[0], show_progress);
    output[0]=spurt::nrrd_manip::trace<T>(input[0], dims[0], show_progress);
}

template<typename T>
void select_operator(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input, 
                     const std::vector<int>& dims, 
                     const std::string& operator_name) {
     if (operator_name=="transpose") return transposer<T>(output, input, dims);
     else if (operator_name=="product") return multiplier<T>(output, input, dims);
     else if (operator_name=="sum") return addition<T>(output, input, dims);
     else if (operator_name=="subtraction") return subtraction<T>(output, input, dims);
     else if (operator_name=="invert") return inverter<T>(output, input, dims);
     else if (operator_name=="determinant") return determinant<T>(output, input, dims);
     else if (operator_name=="trace") return trace<T>(output, input, dims);
     else if (operator_name=="SVD") return SVDecomposer<T>(output, input, dims);
     else if (operator_name=="eigen") return Eigendecomposer<T>(output, input, dims);
     else if (operator_name=="transpose+product") return transpose_multiplier<T>(output, input, dims);
     else if (operator_name=="CG") return CG<T>(output, input, dims);
     else if (operator_name=="FTLE") return FTLE<T>(output, input, dims);
     else throw std::runtime_error("unrecognized operator name: "+operator_name);
}


void select_type(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input, 
                 const std::vector<int>& dims, const std::string& operator_name,
                 const std::string& type_name) {
    int typeID;
    if (type_name=="float") typeID=nrrdTypeFloat;
    else if (type_name=="double") typeID=nrrdTypeDouble;
    else if (type_name=="unknown" || type_name=="implicit") {
        // pick highest precision available among arguments
        typeID=input[0]->type;
        for (int i=1; i<input.size(); ++i) {
            typeID=std::max(typeID, input[i]->type);
        }
        if (typeID<nrrdTypeFloat) typeID=nrrdTypeFloat;
    }
    
    switch(typeID) {
/*  case nrrdTypeChar:   return select_operator<char>           (output, input, dims, operator_name);
    case nrrdTypeUChar:  return select_operator<unsigned char>  (output, input, dims, operator_name);
    case nrrdTypeShort:  return select_operator<short>          (output, input, dims, operator_name);
    case nrrdTypeUShort: return select_operator<unsigned short> (output, input, dims, operator_name);
    case nrrdTypeInt:    return select_operator<int>            (output, input, dims, operator_name);
    case nrrdTypeUInt:   return select_operator<unsigned int>   (output, input, dims, operator_name);
    case nrrdTypeLLong:  return select_operator<long>           (output, input, dims, operator_name);
    case nrrdTypeULLong: return select_operator<unsigned long>  (output, input, dims, operator_name);
*/  case nrrdTypeFloat:  return select_operator<float>          (output, input, dims, operator_name);
    case nrrdTypeDouble: return select_operator<double>         (output, input, dims, operator_name);
    default: throw std::runtime_error("unrecognized type ID");
    }
}

int main(int argc, const char* argv[]) {
    namespace xcl = spurt::command_line;
    
    xcl::option_traits
        positional_group(true, true, "Positional Group"),
        required_group(true, false, "Required Options"), 
        optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
        "Apply linear algebra operators to matrices contained in NRRD files\n\n"
        "Supported operators:\n"
        "+/-/* : Add/subtract/multiply matrices (2 Nrrds)\n"
        "transpose: Transpose matrices (1 Nrrd)\n"
        "SVD: Singular value decomposition (1 Nrrd)\n"
        "invert: Invert square matrices (1 Nrrd)\n"
        "determinant: Self-explanatory (1 Nrrd)\n"
        "trace: Self-explanatory (1 Nrrd)\n"
        "transpose-product: Transpose first then multiply with second (2 Nrrds)\n"
        "Eigen: Eigenvalues and eigenvectors (1 Nrrd)\n"
        "FTLE: Finite-time Lyapunov exponent from flow map gradient (1 Nrrd)\n"
        "Cauchy-Green: Cauchy-Green tensor from flow map gradient (1 Nrrd)\n");

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("operation", opname, opname, "Operator to apply", positional_group);
        parser.add_sequence("input", input_file_names, "Input filename(s)", required_group);
        parser.add_sequence("output", output_file_names, "Output filename(s)", required_group);
        parser.add_sequence("dims", dims, "Matrix dimensions", required_group);
        parser.add_value("type", tname, tname, "Type to use for computation", optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional_group);
        parser.add_flag("sym", is_sym, "Matrices are symmetric", optional_group);
        
        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n" 
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
    
    assert(dims.size()>=1);
    if (verbose) show_progress=true;
    
    std::cout << "input filenames=";
    std::copy(input_file_names.begin(), input_file_names.end(), 
              std::ostream_iterator<std::string>(std::cout, ", "));
    std::cout << "\noutput filenames=";
    std::copy(output_file_names.begin(), output_file_names.end(), 
              std::ostream_iterator<std::string>(std::cout, ", "));
    std::cout << "\noperator name=" << opname << '\n';
    std::cout << "type name=" << tname << '\n';
    
    std::vector<Nrrd*> input, output;
    input.resize(input_file_names.size());
    for (int i=0; i<input.size(); ++i) {
        input[i]=spurt::nrrd_utils::readNrrd(input_file_names[i]);
    }
    
    tname=spurt::lower_case(tname);
    opname=spurt::lower_case(opname);
    
    timer _timer;
    
    // determine selected operator
    try {
        if (opname=="transpose" || opname=="t") {
            opname="t";
            assert(input_file_names.size()==1 &&
                   output_file_names.size()==1);
            select_type(output, input, dims, "transpose", tname);
        }
        else if (opname=="product" || opname=="x" || opname=="*" ||
                 opname=="mult" || opname=="multiply") {
            opname="x";
            assert(input_file_names.size()==2 &&
                   output_file_names.size()==1 &&
                   dims.size()==3);
            select_type(output, input, dims, "product", tname);
        }
        else if (opname=="sum" || opname=="+" || 
                 opname=="add" || opname=="addition") {
            opname="+";
            assert(input_file_names.size()==2 &&
                   output_file_names.size()==1 &&
                   dims.size()==2);
            select_type(output, input, dims, "sum", tname);
        }
        else if (opname=="subtract" || opname=="-" || 
                 opname=="sub" || opname=="subtraction") {
            opname="-";
            assert(input_file_names.size()==2 &&
                   output_file_names.size()==1 &&
                   dims.size()==2);
            select_type(output, input, dims, "subtraction", tname);
        }
        else if (opname=="svd") {
            assert(input_file_names.size()==1 &&
                   (output_file_names.size()==1 || 
                    output_file_names.size()==3));
            select_type(output, input, dims, "SVD", tname);
        }
        else if (opname=="invert" || opname=="/") {
            opname="/";
            assert(input_file_names.size()==1 &&
                   output_file_names.size()==1 &&
                   (dims.size()==1 || (dims.size()==2 && dims[0]==dims[1])));
            select_type(output, input, dims, "invert", tname);
        }
        else if (opname=="determinant" || opname=="det") {
            opname="det";
            assert(input_file_names.size()==1 &&
                   output_file_names.size()==1 &&
                   (dims.size()==1 || (dims.size()==2 && dims[0]==dims[1])));
            select_type(output, input, dims, "determinant", tname);
        }
        else if (opname=="trace" || opname=="tr") {
            opname="trace";
            assert(input_file_names.size()==1 &&
                   output_file_names.size()==1 &&
                   (dims.size()==1 || (dims.size()==2 && dims[0]==dims[1])));
            select_type(output, input, dims, "trace", tname);
        }
        else if (opname=="transprod" || opname=="transx" || 
                 opname=="tprod" || opname=="tx") {
            opname="tx";
            assert(input_file_names.size()==2 &&
                   output_file_names.size()==1 &&
                   dims.size()==3);
            select_type(output, input, dims, "transpose+product", tname);
        }
        else if (opname=="eigen" || opname=="ev") {
            opname="eigen";
            assert(input_file_names.size()==1 && 
                   (output_file_names.size()==1 || output_file_names.size()==2));
            select_type(output, input, dims, "eigen", tname);
        }
        else if (opname=="ftle") {
            opname="FTLE";
            assert(input_file_names.size()==1 && 
                   output_file_names.size()==1);
            select_type(output, input, dims, "FTLE", tname);
        }
        else if (opname=="cg" || opname=="cauchygreen" || opname=="cauchy-green") {
            opname="CG";
            assert(input_file_names.size()==1 && 
                   output_file_names.size()==1);
            select_type(output, input, dims, "CG", tname);
        }
        else {
            std::cout << "ERROR: unrecognized operator string: " << opname << '\n';
            ok=false;
        }
    }
    catch(std::exception& e) {
        std::cout << "Exception caught: " << e.what() << '\n';
        ok=false;
    }
    
    std::cout << "Overall processing time: " << _timer.elapsed() << " s.\n";
    
    if (ok) {
        if (output.size()==3 && output_file_names.size()==1 && opname=="svd") {
            std::string basename=spurt::filename::remove_extension(output_file_names[0]);
            output_file_names.resize(3);
            output_file_names[0]=basename+"-sinvals.nrrd";
            output_file_names[1]=basename+"-leftvec.nrrd";
            output_file_names[2]=basename+"-rightvec.nrrd";
        }
        else if (output.size()==2 && output_file_names.size()==1 && opname=="eigen") {
            std::string basename=spurt::filename::remove_extension(output_file_names[0]);
            output_file_names.resize(2);
            output_file_names[0]=basename+"-eigenvals.nrrd";
            output_file_names[1]=basename+"-eigenvecs.nrrd";
        }
        for (int i=0; i<output.size(); ++i) {
            if (nrrdSave(output_file_names[i].c_str(), output[i], NULL)) {
                std::cout << "error while exporting nrrd file\n";
                break;
            }
        }
    }
    
    for (int i=0; i<input.size(); ++i) {
        nrrdNuke(input[i]);
    }
    
    return 0;
}
