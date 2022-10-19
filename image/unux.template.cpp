#include <array>
#include <iomanip>

#include <util/timer.hpp>

#include <image/nrrd_wrapper.hpp>
#include <image/nrrd_manip.hpp>
#include <misc/option_parse.hpp>
#include <format/filename.hpp>

struct transpose_wrapper {
    template<typename T, int Nrows, int Ncols, int, int> 
    static void apply(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input) {
        assert(input.size()==1);
        xavier::nrrd_matrix_transpose<T, Nrows, Ncols> trans_op;
        output.resize(0);
        output[0]=trans_op(input[0]);
    }
};

struct product_wrapper {
    template<typename T, int Nrows, int Ncols, int Ncols2, int, int> 
    static void apply(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input) {
        assert(input.size()==2);
        xavier::nrrd_matrix_product<T, Nrows, Ncols, Ncols2> prod_op;
        output.resize(0);
        output[0]=prod_op(input[0], input[1]);
    }
};

struct svd_wrapper {
    template<typename T, int Nrows, int Ncols, int, int> 
    static void apply(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input) {
        assert(input.size()==1);
        xavier::nrrd_matrix_svd<T, Nrows, Ncols> svd_op;
        output.resize(3);
        svd_op(input[0], output[0], output[1], output[2]);
    }
};

struct trans_prod_wrapper {
    template<typename T, int Nrows, int Ncols, int Ncols2, int> 
    static void apply(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input) {
        assert(input.size()==2);
        xavier::nrrd_matrix_transpose_product<T, Nrows, Ncols, Ncols2> transx_op;
        output.resize(1);
        output[0]=transx_op(input[0], input[1]);
    }
};

struct invert_wrapper {
    template<typename T, int N, int, int, int> 
    static void apply(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input) {
        assert(input.size()==1);
        xavier::nrrd_matrix_invert<T, N> inv_op;
        output.resize(1);
        output[0]=inv_op(input[0]);
    }
};

struct dummy_operator {
    template<typename, int, int, int, int> 
    static void apply(std::vector<Nrrd*>&, const std::vector<Nrrd*>&) {
        throw std::runtime_error("No operator defined with 4 size parameters");
    }
};

template<typename T, typename Functor, int M, int N, int P, int Q>
struct operator_wrapper {
    void operator()(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                    Functor op) {
         return op.template apply<T, M, N, P, Q>(output, input);               
    }
};

template<typename T, typename Functor, int M, int N, int P>
struct operator_wrapper<T, Functor, M, N, P, -1> {
    void operator()(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                    Functor op) {
        // return op.template apply<T, M, N, P>(output, input);
    }  
};

template<typename T, typename Functor, int M, int N>
struct operator_wrapper<T, Functor, M, N, -1, -1> {
    void operator()(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                    Functor op) {
        return op.template apply<T, M, N>(output, input);
    }
};

template<typename T, typename Functor, int M>
struct operator_wrapper<T, Functor, M, -1, -1, -1> {
    void operator()(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                    Functor op) {
        return op.template apply<T, M>(output, input);
    }
};

template<typename Functor, int M, int N=-1, int P=-1, int Q=-1>
void do_something_sized_untyped(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                                Functor op, int typeID) {
    if (typeID==nrrdTypeUnknown) {
        // pick highest precision available among arguments
        typeID=input[0]->type;
        for (int i=1; i<input.size(); ++i) {
            typeID=std::max(typeID, input[i]->type);
        }
    }
    switch(typeID) {
    case nrrdTypeChar: return operator_wrapper<char, Functor, M, N, P, Q>()(output, input, op);
    case nrrdTypeUChar: return operator_wrapper<unsigned char, Functor, M, N, P, Q>()(output, input, op);
    case nrrdTypeShort: return operator_wrapper<short, Functor, M, N, P, Q>()(output, input, op);
    case nrrdTypeUShort: return operator_wrapper<unsigned short, Functor, M, N, P, Q>()(output, input, op);
    case nrrdTypeInt: return operator_wrapper<int, Functor, M, N, P, Q>()(output, input, op);
    case nrrdTypeUInt: return operator_wrapper<unsigned int, Functor, M, N, P, Q>()(output, input, op);
    case nrrdTypeLLong: return operator_wrapper<long, Functor, M, N, P, Q>()(output, input, op);
    case nrrdTypeULLong: return operator_wrapper<unsigned long, Functor, M, N, P, Q>()(output, input, op);
    case nrrdTypeFloat: return operator_wrapper<float, Functor, M, N, P, Q>()(output, input, op);
    case nrrdTypeDouble: return operator_wrapper<double, Functor, M, N, P, Q>()(output, input, op);
    default: throw std::runtime_error("unrecognized type ID");
    }
}

template<typename Functor, int M, int N, int P, int Q>
void do_something_sized4(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                         const std::vector<int>& sizes, Functor op, int typeID) {
    if (sizes.size()>3) {
        throw std::runtime_error("5 size parameters not supported");
    }
    do_something_sized_untyped<Functor, M, N, P, Q>(output, input, op, typeID);
}

template<typename Functor, int M, int N, int P>
void do_something_sized3(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                         const std::vector<int>& sizes, Functor op, int typeID) {
    if (sizes.size()==3) do_something_sized_untyped<Functor, M, N, P>(output, input, op, typeID);
    // else if (sizes[3]==1) do_something_sized4<Functor, M, N, P, 1>(output, input, sizes, op, typeID);
    // else if (sizes[3]==2) do_something_sized4<Functor, M, N, P, 2>(output, input, sizes, op, typeID);
    // else if (sizes[3]==3) do_something_sized4<Functor, M, N, P, 3>(output, input, sizes, op, typeID);
    // else if (sizes[3]==4) do_something_sized4<Functor, M, N, P, 4>(output, input, sizes, op, typeID);
    else throw std::runtime_error("Unsupported 4th size parameter");
}

template<typename Functor, int M, int N>
void do_something_sized2(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                         const std::vector<int>& sizes, Functor op, int typeID) {
    if (sizes.size()==2) do_something_sized_untyped<Functor, M, N>(output, input, op, typeID);
    else if (sizes[2]==1) do_something_sized3<Functor, M, N, 1>(output, input, sizes, op, typeID);
    else if (sizes[2]==2) do_something_sized3<Functor, M, N, 2>(output, input, sizes, op, typeID);
    else if (sizes[2]==3) do_something_sized3<Functor, M, N, 3>(output, input, sizes, op, typeID);
    else if (sizes[2]==4) do_something_sized3<Functor, M, N, 4>(output, input, sizes, op, typeID);
    else throw std::runtime_error("Unsupported 3rd size parameter");
}

template<typename Functor, int M>
void do_something_sized1(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                         const std::vector<int>& sizes, Functor op, int typeID) {
    if (sizes.size()==1) do_something_sized_untyped<Functor, M>(output, input, op, typeID);
    else if (sizes[1]==1) do_something_sized2<Functor, M, 1>(output, input, sizes, op, typeID); 
    else if (sizes[1]==2) do_something_sized2<Functor, M, 2>(output, input, sizes, op, typeID); 
    else if (sizes[1]==3) do_something_sized2<Functor, M, 3>(output, input, sizes, op, typeID);   
    else if (sizes[1]==4) do_something_sized2<Functor, M, 4>(output, input, sizes, op, typeID);    
    else throw std::runtime_error("Unsupported 2nd size parameter");                   
}

template<typename Functor>
void do_something_unsized(std::vector<Nrrd*>& output, const std::vector<Nrrd*>& input,
                          const std::vector<int>& sizes, Functor op, int typeID) {
     if      (sizes[0]==1) do_something_sized1<Functor, 1>(output, input, sizes, op, typeID);
     else if (sizes[0]==2) do_something_sized1<Functor, 2>(output, input, sizes, op, typeID);
     else if (sizes[0]==3) do_something_sized1<Functor, 3>(output, input, sizes, op, typeID);
     else if (sizes[0]==4) do_something_sized1<Functor, 4>(output, input, sizes, op, typeID);
     else throw std::runtime_error("Unsupported first size parameter");
}


int main(int argc, char* argv[]) {
    namespace xcl = xavier::command_line;
    std::vector<std::string> input_file_names;
    std::vector<std::string> output_file_names;
    std::vector<int> dims;
    bool verbose=false;
    int typeID=nrrdTypeUnknown;
    bool ok=true;
    
    std::string operation[][2] = {
        {"transpose", "t"},
        {"multiply", "x"},
        {"invert", "/"},
        {"SVD", "svd"},
    };
    
    std::string opname="none";
    std::string tname="none";
    
    xcl::option_traits
        positional_group(true, true, "Positional Group"),
        required_group(true, false, "Required Options"), 
        optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
        "Apply linear algebra operators to matrices contained in NRRD files");

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("operation", opname, opname, "Operator to apply", positional_group);
        parser.add_sequence("input", input_file_names, "Input filename(s)", required_group);
        parser.add_sequence("output", output_file_names, "Output filename(s)", required_group);
        parser.add_sequence("dims", dims, "Matrix dimensions", required_group);
        parser.add_value("type", typeID, typeID, "Type to use for computation", optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional_group);
        
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
    
    std::vector<Nrrd*> input, output;
    input.resize(input_file_names.size());
    for (int i=0; i<input.size(); ++i) {
        input[i]=xavier::readNrrd(input_file_names[i]);
    }
    
    // determine selected operator
    try {
    if (opname=="transpose" || opname=="t" || opname=="T") {
        assert(input_file_names.size()==1 &&
               output_file_names.size()==1);
        do_something_unsized(input, output, dims, transpose_wrapper(), typeID);
    }
    else if (opname=="product" || opname=="x" || opname=="mult" || opname=="multiply") {
        assert(input_file_names.size()==2 &&
               output_file_names.size()==1 &&
               dims.size()==3);
        do_something_unsized(input, output, dims, product_wrapper(), typeID);
    }
    else if (opname=="SVD" || opname=="svd") {
        assert(input_file_names.size()==1 &&
               (output_file_names.size()==1 || 
                output_file_names.size()==3));
        do_something_unsized(input, output, dims, svd_wrapper(), typeID);
    }
    else if (opname=="invert" || opname=="/") {
        assert(input_file_names.size()==1 &&
               output_file_names.size()==1);
        do_something_unsized(input, output, dims, invert_wrapper(), typeID);
    }
    else if (opname=="transprod" || opname=="transx" || opname=="tprod" || opname=="Tprod" || opname=="tx" || opname=="Tx") {
        assert(input_file_names.size()==2 &&
               output_file_names.size()==1 &&
               dims.size()==3);
        do_something_unsized(input, output, dims, trans_prod_wrapper(), typeID);
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
    
    if (ok) {
        if (output.size()==3 && output_file_names.size()==1) {
            std::string basename=xavier::get_basename(output_file_names[0]);
            output_file_names.resize(3);
            output_file_names[0]=basename+"-sinvals.nrrd";
            output_file_names[1]=basename+"-leftvec.nrrd";
            output_file_names[2]=basename+"-rightvec.nrrd";
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
