#include <image/nrrd_wrapper.hpp>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <misc/option_parse.hpp>

std::string cmdline;
std::string input, output;
double length;
bool verbose;

void initialize(int argc, const char* argv[]) {
    namespace xcl = spurt::command_line;
    verbose = false;

    cmdline = "Command line: " + std::string(argv[0]);
    for (int i=1; i<argc; i++) {
        cmdline += " " + std::string(argv[i]);
    }

    xcl::option_traits
        required(true, false, "Required Options"),
    optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
    "Compute FTLE from flow map");

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", input, "Input flow map filename", required);
        parser.add_value("output", output, "Output base name", required);
        parser.add_value("T", length, "Integration time", required);
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

    Nrrd* nin = spurt::nrrd_utils::readNrrd(input);
    assert(nin->dim == 4);
    size_t res[3] = { nin->axis[1].size, nin->axis[2].size, nin->axis[3].size };
    size_t N = res[0]*res[1]*res[2];
    size_t shift[3] = { 1, res[0], res[0]*res[1] };
    double spc[3] = { nin->axis[1].spacing, nin->axis[2].spacing, nin->axis[3].spacing };

    spurt::nrrd_utils::nrrd_data_wrapper<double> wrapper(nin);

    typedef Eigen::Matrix<double, 3, 3> mat_t;
    typedef Eigen::Matrix<double, 3, 1> vec_t;
    typedef Eigen::SelfAdjointEigenSolver<mat_t> solver_t;

    double* ftle = (double*)calloc(res[0]*res[1]*res[2], sizeof(double));

    for (size_t n=0; n<N; n++) {
        size_t c[3] = { n % res[0], (n / res[0]) % res[1], n / (res[0]*res[1]) };

        mat_t J, C;
        vec_t f1, f2, f;
        solver_t solver;

        for (int dim=0; dim<3; ++dim) {
            size_t s = shift[dim];
            if (c[dim]>0 && c[dim]<res[dim]-1) {
                f1[0] = wrapper[3*(n+s)  ];
                f1[1] = wrapper[3*(n+s)+1];
                f1[2] = wrapper[3*(n+s)+2];
                f2[0] = wrapper[3*(n-s)  ];
                f2[1] = wrapper[3*(n-s)+1];
                f2[2] = wrapper[3*(n-s)+2];
                f = f2-f1;
                f *= 1./(2*spc[dim]);
            }
            else if (c[dim]==0) {
                f2[0] = wrapper[3*(n+s)  ];
                f2[1] = wrapper[3*(n+s)+1];
                f2[2] = wrapper[3*(n+s)+2];
                f1[0] = wrapper[3*n  ];
                f1[1] = wrapper[3*n+1];
                f1[2] = wrapper[3*n+2];
                f = f2-f1;
                f *= 1./spc[dim];
            }
            else {
                f1[0] = wrapper[3*(n-s)  ];
                f1[1] = wrapper[3*(n-s)+1];
                f1[2] = wrapper[3*(n-s)+2];
                f2[0] = wrapper[3*n  ];
                f2[1] = wrapper[3*n+1];
                f2[2] = wrapper[3*n+2];
                f = f2-f1;
                f *= 1./spc[dim];
            }

            J(dim, 0) = f[0];
            J(dim, 1) = f[1];
            J(dim, 2) = f[2];

            C = J.transpose() * J;
            solver.compute(C);
            double lambda = solver.eigenvalues()[2];
            if (lambda > 0) {
                ftle[n] = 1/length * log(sqrt(lambda));
            }
        }
    }

    Nrrd *nout = nrrdNew();

    if (nrrdWrap_nva(nout, ftle, nrrdTypeDouble, 3, res)) {
        throw std::runtime_error(spurt::nrrd_utils::error_msg("writeNrrd: error while wrapping"));
    }
    nout->axis[0].spacing = spc[0];
    nout->axis[1].spacing = spc[1];
    nout->axis[2].spacing = spc[2];
    if (nrrdSave(output.c_str(), nout, NULL)) {
        throw std::runtime_error(spurt::nrrd_utils::error_msg("writeNrrd: error while saving"));
    }

    nrrdNix(nout); // only deletes the structure not the data
}
