#include "vtk_utils.hpp"
#include <format/filename.hpp>
#include <misc/option_parse.hpp>
#include <misc/strings.hpp>
#include <iostream>
#include <vector>
#include <string>

int main(int argc, const char* argv[]) {
    namespace xcl = spurt::command_line;

    std::string ofname, ifname;
    std::string expr;
    nvis::bbox3 bounds(nvis::vec3(0), nvis::vec3(1,1,0));
    nvis::ivec3 res(10, 10, 1);
    nvis::vec4 bounds_as_vec(0,0,1,1);
    bool verbose=false;
    double length=100;
    double eps=0.001;
    bool enhance=false;
    std::vector<std::string> exprs;
    int factor=1;

    xcl::option_traits
        required_group(true, false, "Required Options"),
        optional_group(false, false, "Optional parameters");

    xcl::option_parser parser(argv[0], "Test VTK's GPU-based LIC computation");

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("output", ofname, "Output filename",
                         required_group);
        parser.add_value("rhs", expr, "RHS expression (use ';' to separate " \
                         "multiple expressions OR @ prefix to indicate " \
                         "filename)", required_group);
        parser.add_value("length", length, length, "Integration length",
                         optional_group);
        parser.add_value("eps", eps, eps, "Integration step size",
                         optional_group);
        parser.add_value("factor", factor, factor, "LIC magnification factor",
                         optional_group);
        parser.add_value("enhance", enhance, enhance, "Enhanced LIC",
                         optional_group);
        parser.add_tuple<4>("bounds", bounds_as_vec, bounds_as_vec,
                            "Bounds (xmin, ymin, xmax, ymax)", optional_group);
        parser.add_tuple<3>("res", res, res, "Image resolution",
                            optional_group);
        parser.add_value("verbose", verbose, verbose, "verbose level",
                         optional_group);

        parser.parse(argc, argv);
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n"
                  << e.what() << '\n'
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        return 1;
    }

    res[2]=1;

    if (nvis::norm(bounds_as_vec)!=0) {
        bounds.min()[0]=bounds_as_vec[0];
        bounds.min()[1]=bounds_as_vec[1];
        bounds.max()[0]=bounds_as_vec[2];
        bounds.max()[1]=bounds_as_vec[3];
    }
    if (verbose) {
        std::cout << "selected bounds=\n" << bounds << '\n';
    }

    if (expr[0]=='@') {
        ifname=expr.substr(1);
    }
    else {
        spurt::tokenize(exprs, expr, ";");

        if (ifname.empty() && exprs.size()!=1 && exprs.size()!=2) {
            std::cerr << "ERROR: " << argv[0] << ": invalid number of terms: "
                << exprs.size() << '\n';
            return 1;
        }

        if (exprs.size()==2) {
            exprs.push_back("0");
        }

        if (verbose) {
            std::cout << "There are " << exprs.size()
                << " expressions in input:\n";
            for (int i=0; i<exprs.size(); ++i) {
                std::cout << "\texpr #" << i << ": " << exprs[i] << '\n';
            }
        }
    }

    ofname = spurt::filename::remove_extension(ofname) + ".png";
    if (verbose) {
        std::cout << "name of exported image will be "
                  << ofname << '\n';
    }

    vtkSmartPointer<vtkImageData> rhs_img;
    if (ifname.empty()) {
        if (verbose) {
            std::cout << "converting expression(s) to an image...\n";
        }
        rhs_img=vtk_utils::expr_to_image<float>(exprs, bounds, res);
    }
    else {
        rhs_img=
            vtk_utils::
            load_image_of_given_format<vtkStructuredPointsReader>(ifname);
        if (verbose) {
            nvis::ivec3 dims;
            if (rhs_img.GetPointer()) {
                rhs_img->GetDimensions(static_cast<int*>(dims.begin()));
                std::cout << "image dimensions: " << dims << '\n';
            }
        }
    }

    if (verbose) {
        std::cout << "\timage created\n";
    }

    if (verbose && ifname.empty()) {
        std::cout << "Saving rhs image to file\n";
        vtk_utils::save_image(rhs_img, "vf_img.vtk");
    }

    if (verbose) {
        std::cout << "computing LIC image...\n";
    }

    vtkSmartPointer<vtkImageData> lic=vtk_utils::do_lic(rhs_img, length, eps,
                                                        enhance, factor);

    if (verbose) {
        std::cout << "\tLIC image computed\n";
    }
    vtk_utils::save_image(lic, ofname);
    if (verbose) {
        std::cout << "LIC image exported successfully to " << ofname << '\n';
    }

    return 0;
}
