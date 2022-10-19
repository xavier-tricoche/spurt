#include <list>
#include <vector>
#include <string>
#include <exception>

#include <geometry/tetrahedrization.hpp>
#include <VTK/vtk_utils.hpp>
#include <misc/option_parse.hpp>
#include <misc/sort.hpp>
#include <misc/progress.hpp>

bool verbose;
bool refine;
std::string name_in, name_out;

typedef xavier::fixed_sorted_vector<size_t, 2> edge_type;
typedef nvis::fixed_vector<double, 3> pos_type;

void initialize(int argc, char* argv[]) {
    namespace xcl = xavier::command_line;

    xcl::option_traits
            required(true, false, "Required Options"),
            optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
            "Tetrahedrize and refine mesh");

    verbose = false;
    refine = false;
    name_out = "refined.vtu";

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("input", name_in, "Input filename", required);
        parser.add_value("output", name_out, name_out, "Output filename", optional);
        parser.add_value("refine", refine, refine, "Refine mesh", optional);
        parser.add_value("verbose", verbose, verbose, "Verbose output", optional);

        parser.parse(argc, const_cast<const char**>(argv));
    }
    catch(std::runtime_error& e) {
        std::cerr << "ERROR: " << argv[0] << " threw exception:\n"
                  << e.what() << "\n"
                  << parser.print_self(false, true, false) << "\n\n\n";
        exit(1);
    }
}

VTK_SMART(vtkUnstructuredGrid) tetrahedrize(VTK_SMART(vtkDataSet) dataset) {
    size_t ncells = dataset->GetNumberOfCells();

    std::cout << "Tetrahedrization...\n";

    VTK_CREATE(vtkCellArray, tets);

    xavier::ProgressDisplay progress(true);
    progress.begin(ncells, "Splitting tets");
    for (size_t i=0; i<ncells; ++i) {
        progress.update(i);
        VTK_CREATE(vtkGenericCell, cell);
        dataset->GetCell(i, cell);
        int celltype = cell->GetCellType();
        std::list<xavier::tetrahedron> newtets;
        if (celltype == VTK_TETRA) {
            // simply copy cell
            tets->InsertNextCell(cell);
            continue;
        }

        VTK_SMART(vtkIdList) id_list = cell->GetPointIds();
        if (celltype == VTK_VOXEL || celltype == VTK_HEXAHEDRON) {
            xavier::hexahedron hexa;
            for (int k=0; k<8; ++k) {
                hexa[k] = id_list->GetId(k);
            }
            xavier::decompose_hexahedron(newtets, hexa);
        }
        else if (celltype == VTK_WEDGE) {
            xavier::prism pr;
            for (int k=0; k<6; ++k) {
                pr[k] = id_list->GetId(k);
            }
            xavier::decompose_prism(newtets, pr);
        }
        else if (celltype == VTK_PYRAMID) {
            xavier::pyramid py;
            for (int k=0; k<5; ++k) {
                py[k] = id_list->GetId(k);
            }
            xavier::decompose_pyramid(newtets, py);
        }
        else {
            std::cerr << "Unrecognized / unsupported cell type: " << celltype << '\n';
            continue;
        }

        for (auto it=newtets.begin(); it!=newtets.end(); ++it) {
            tets->InsertNextCell(4);
            tets->InsertCellPoint((*it)[0]);
            tets->InsertCellPoint((*it)[1]);
            tets->InsertCellPoint((*it)[2]);
            tets->InsertCellPoint((*it)[3]);
        }
    }
    progress.end();

    VTK_CREATE(vtkUnstructuredGrid, tetmesh);
    tetmesh->SetCells(VTK_TETRA, tets);
    if (vtkPointSet::SafeDownCast(dataset) != nullptr) {
        tetmesh->SetPoints(vtkPointSet::SafeDownCast(dataset)->GetPoints());
    }
    else {
        VTK_CREATE(vtkPoints, pts);
        size_t npts = dataset->GetNumberOfPoints();
        pts->SetNumberOfPoints(npts);
        for (size_t i=0; i<npts; ++i) {
            double d[3];
            dataset->GetPoint(i, d);
            pts->SetPoint(i, d);
        }
        tetmesh->SetPoints(pts);
    }
    for (size_t i=0; i<dataset->GetPointData()->GetNumberOfArrays(); ++i) {
        tetmesh->GetPointData()->AddArray(dataset->GetPointData()->GetArray(i));
    }
    return tetmesh;
}

VTK_SMART(vtkUnstructuredGrid) refine_tetmesh(VTK_SMART(vtkUnstructuredGrid) tetmesh) {
    size_t ncells = tetmesh->GetNumberOfCells();
    size_t npts = tetmesh->GetNumberOfPoints();

    VTK_CREATE(vtkCellArray, tets);
    std::vector<pos_type> all_new_pos;
    std::map<edge_type, size_t> edge_to_index;
    xavier::ProgressDisplay progress(true);
    progress.begin(ncells, "Subdividing tetrahedra");
    for (size_t i=0; i<ncells; ++i) {
        progress.update(i);
        VTK_CREATE(vtkGenericCell, cell);
        tetmesh->GetCell(i, cell);
        int celltype = cell->GetCellType();
        VTK_SMART(vtkIdList) id_list = cell->GetPointIds();
        std::list<xavier::tetrahedron> newtets;
        std::vector<pos_type> old_pos(4);
        for (int k=0; k<4; ++k) {
            tetmesh->GetPoint(id_list->GetId(k), &old_pos[k][0]);
        }

        std::vector<edge_type> edges(6);
        edges[0] = edge_type(id_list->GetId(0), id_list->GetId(1));
        edges[1] = edge_type(id_list->GetId(1), id_list->GetId(2));
        edges[2] = edge_type(id_list->GetId(2), id_list->GetId(0));
        edges[3] = edge_type(id_list->GetId(0), id_list->GetId(3));
        edges[4] = edge_type(id_list->GetId(1), id_list->GetId(3));
        edges[5] = edge_type(id_list->GetId(2), id_list->GetId(3));
        for (int k=0; k<6; ++k) {
            const edge_type& e = edges[k];
            if (edge_to_index.find(e) == edge_to_index.end()) {
                pos_type p = 0.5*(old_pos[e[0]] + old_pos[e[1]]);
                all_new_pos.push_back(p);
                edge_to_index[e] = all_new_pos.size()-1;
            }
        }

        tets->InsertNextCell(4);
        tets->InsertCellPoint(id_list->GetId(0));
        tets->InsertCellPoint(npts+edge_to_index[edges[0]]);
        tets->InsertCellPoint(npts+edge_to_index[edges[2]]);
        tets->InsertCellPoint(npts+edge_to_index[edges[3]]);

        tets->InsertNextCell(4);
        tets->InsertCellPoint(id_list->GetId(1));
        tets->InsertCellPoint(npts+edge_to_index[edges[1]]);
        tets->InsertCellPoint(npts+edge_to_index[edges[0]]);
        tets->InsertCellPoint(npts+edge_to_index[edges[4]]);

        tets->InsertNextCell(4);
        tets->InsertCellPoint(id_list->GetId(2));
        tets->InsertCellPoint(npts+edge_to_index[edges[2]]);
        tets->InsertCellPoint(npts+edge_to_index[edges[1]]);
        tets->InsertCellPoint(npts+edge_to_index[edges[5]]);

        tets->InsertNextCell(4);
        tets->InsertCellPoint(npts+edge_to_index[edges[3]]);
        tets->InsertCellPoint(npts+edge_to_index[edges[4]]);
        tets->InsertCellPoint(npts+edge_to_index[edges[5]]);
        tets->InsertCellPoint(id_list->GetId(3));

        // split octahedron into 2 pyramids
        xavier::pyramid py1, py2;
        py1[0] = npts+edge_to_index[edges[3]];
        py1[1] = npts+edge_to_index[edges[4]];
        py1[2] = npts+edge_to_index[edges[1]];
        py1[3] = npts+edge_to_index[edges[2]];
        py1[4] = npts+edge_to_index[edges[5]];

        py2[0] = npts+edge_to_index[edges[3]];
        py2[1] = npts+edge_to_index[edges[4]];
        py2[2] = npts+edge_to_index[edges[1]];
        py2[3] = npts+edge_to_index[edges[2]];
        py2[4] = npts+edge_to_index[edges[0]];

        std::list<xavier::tetrahedron> _tets;
        xavier::decompose_pyramid(_tets, py1);
        xavier::decompose_pyramid(_tets, py2);
        for (auto it=_tets.begin(); it!=_tets.end(); ++it) {
            const xavier::tetrahedron& t = *it;
            tets->InsertNextCell(4);
            tets->InsertCellPoint(t[0]);
            tets->InsertCellPoint(t[1]);
            tets->InsertCellPoint(t[2]);
            tets->InsertCellPoint(t[3]);
        }
    }
    progress.end();

    VTK_CREATE(vtkUnstructuredGrid, finer);
    finer->SetCells(VTK_TETRA, tets);

    VTK_CREATE(vtkDoubleArray, coords);
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(npts+all_new_pos.size());
    for (size_t i=0; i<npts; ++i) {
        double c[3];
        tetmesh->GetPoint(i, c);
        coords->SetTuple3(i, c[0], c[1], c[2]);
    }
    for (size_t i=0; i<all_new_pos.size(); ++i) {
        const pos_type& p=all_new_pos[i];
        coords->SetTuple3(npts+i, p[0], p[1], p[2]);
    }
    VTK_CREATE(vtkPoints, points);
    points->SetData(coords);
    finer->SetPoints(points);
    return finer;
}

int main(int argc, char* argv[]) {
    initialize(argc, argv);

    VTK_SMART(vtkDataSet) input = vtk_utils::readVTK(name_in);
    VTK_SMART(vtkUnstructuredGrid) tets = tetrahedrize(input);
    VTK_SMART(vtkUnstructuredGrid) finetets = refine_tetmesh(tets);

    vtk_utils::saveVTK(finetets, name_out);
    return 0;
}
