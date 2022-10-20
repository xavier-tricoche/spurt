#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include <vtk/vtk_utils.hpp>
#include <vtk/vtk_camera_helper.hpp>
#include <misc/option_parse.hpp>
#include <math/fixed_vector.hpp>

#include <vtkColorTransferFunction.h>
#include <vtkDataSet.h>

typedef nvis::fixed_vector<double, 3> color_type;

std::string lines_name, scalars_name, img_name, cam_in, cam_out;
bool verbose=false;
color_type bg_col(0, 0, 0);
nvis::ivec2 res(800, 800);
float radius=0.01;

void initialize(int argc, char* argv[]) {
    namespace xcl = spurt::command_line;

    xcl::option_traits
            required_group(true, false, "Required Options"),
            optional_group(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Visualize a set of lines with per-line scalar attributes");

    try {
        parser.use_short_symbols(false);
        parser.use_brackets(true);
        parser.add_value("line", lines_name, "Lines filename", required_group);
        parser.add_value("scalar", scalars_name, "Scalars filename", optional_group);
        parser.add_value("save", img_name, "Snapshot filename", optional_group);
        parser.add_value("radius", radius, "Tube radius", optional_group);
        parser.add_value("camin", cam_in, "Camera input filename", optional_group);
        parser.add_value("camout", cam_out, "Camera output filename", optional_group);
        parser.add_value("bg", bg_col, bg_col, "Background color", optional_group);
        parser.add_tuple<2>("res", res, res, "Image resolution", optional_group);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional_group);

        parser.parse(argc, const_cast<const char**>(argv));
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception: "
                  << e.what() << "\n"
                  << "Command line options entered so far:\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

void read_color(std::vector< std::pair<double, color_type> >& cps, const std::string filename) {
    std::fstream input(filename.c_str(), std::ios::in);
    if (!input) {
        throw std::runtime_error("Unable to open " + filename + ": no such file or directory");
    }
    cps.clear();
    std::pair<double, color_type> point;
    while (!input.eof()) {
        input >> point.first >> point.second[0] >> point.second[1] >> point.second[2];
        cps.push_back(point);
    }
    input.close();
    if (verbose) {
        std::cout << cps.size() << " color control points successfully imported\n";
    }
}

int main(int argc, char* argv[]) {
    initialize(argc, argv);

    vtkSmartPointer<vtkPolyDataReader> lines_reader = vtkSmartPointer<vtkPolyDataReader>::New();
    lines_reader->SetFileName(lines_name.c_str());
    lines_reader->Update();

    double min, max;
    vtkSmartPointer<vtkPolyData> poly;

    vtkSmartPointer<vtkDataSetReader> scalars_reader = vtkSmartPointer<vtkDataSetReader>::New();
    if (!scalars_name.empty()) {
        scalars_reader->SetFileName(scalars_name.c_str());
        scalars_reader->Update();
        vtkSmartPointer<vtkDataArray> scalars = scalars_reader->GetOutput()->GetPointData()->GetScalars();

        poly = lines_reader->GetOutput();
        vtkSmartPointer<vtkCellArray> lines = poly->GetLines();
        lines->InitTraversal();
        vtkIdType npts;
        const vtkIdType* ids;
        size_t ncells = lines->GetNumberOfCells();

        std::vector<float> data(poly->GetNumberOfPoints(), 0);
        for (size_t i=0; i<ncells; ++i) {
            lines->GetNextCell(npts, ids);
            double val = scalars->GetTuple1(i);
            for (size_t n=0; n<npts; ++n) {
                data[ids[n]] = val;
            }
        }

        min = *std::min_element(data.begin(), data.end());
        max = *std::max_element(data.begin(), data.end());

        vtk_utils::add_scalars(poly, data);
    }

    vtkSmartPointer<vtkTubeFilter> tubes = vtkSmartPointer<vtkTubeFilter>::New();
    tubes->SetInputData(poly);
    tubes->SetRadius(radius);
    tubes->SetNumberOfSides(6);
    if (!scalars_name.empty()) {
        tubes->SetVaryRadiusToVaryRadiusByScalar();
    }

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(tubes->GetOutputPort());

    vtkSmartPointer<vtkScalarBarActor> color_bar;
    if (!scalars_name.empty()) {
        vtkSmartPointer<vtkColorTransferFunction> ctf = vtkSmartPointer<vtkColorTransferFunction>::New();
        ctf->AddRGBPoint(min, 0.2, 0.2, 0.2);
        ctf->AddRGBPoint(max, 1, 0, 0);
        mapper->SetLookupTable(ctf);
        vtk_utils::colorbar_param param;
        param.pos[1] = 0.2;
        param.title = scalars_name;
        color_bar = vtk_utils::colorbar(ctf, param);
    }

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    renderer->AddActor(actor);
    if (!scalars_name.empty()) {
        renderer->AddActor(color_bar);
    }
    renderer->SetBackground(bg_col[0], bg_col[1], bg_col[2]);

    vtkSmartPointer<vtkRenderWindow> window = vtkSmartPointer<vtkRenderWindow>::New();
    window->AddRenderer(renderer);

    vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    interactor->SetRenderWindow(window);

    renderer->ResetCamera();
    window->SetSize(res[0], res[1]);
    interactor->Initialize();
    window->Render();
    interactor->Start();

    return 0;
}
