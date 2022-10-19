#ifdef _OPENMP
#include <omp.h>
#endif
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <geometry/boundary_extraction.hpp>
#include <format/filename.hpp>
#include <teem/nrrd.h>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <boost/math/constants/constants.hpp>
#include <math.h>

#include <format/filename.hpp>
#include <VTK/vtk_utils.hpp>
#include <misc/option_parse.hpp>

#include <vtkCellIterator.h>

#include <Eigen/Eigen>

typedef nvis::bounding_box<nvis::fixed_vector<float, 3> > fbox3;

/*
void load_VTK(celltree::mesh** amesh, celltree::variable** scalars, celltree::variable** vectors, const std::string& filename) {

    VTK_SMART(vtkDataSet) vtk = vtk_utils::readVTK(filename);

    size_t ncells = vtk->GetNumberOfCells();
    celltree::mesh::cell* cells = new celltree::mesh::cell[ncells];

    size_t npts = vtk->GetNumberOfPoints();
    std::vector<float> coords;
    std::vector<unsigned int> ids;

    float* points = new float[npts*3];
    double c[3];
    for (size_t pid=0; pid<npts; ++pid) {
        vtk->GetPoint(pid, c);
        std::copy(c, c+3, &points[3*pid]);
    }

    vtkCellIterator* iter = vtk->NewCellIterator();
    celltree::mesh::cell* cit = cells;
    size_t counter = 0;
    for (iter->InitTraversal(); !iter->IsDoneWithTraversal(); iter->GoToNextCell())
    {
        vtkIdList *pointIds = iter->GetPointIds();
        int nptsincell = iter->GetNumberOfPoints();
        for (int i=0; i<nptsincell; ++i) {
            ids.push_back(pointIds->GetId(i));
        }
        celltree::cell_kind type;
        switch (iter->GetCellType()) {
            VTK_TETRA: type = celltree::TETRAHEDRON; break;
            VTK_HEXAHEDRON: type = celltree::HEXAHEDRON; break;
            VTK_PYRAMID: type = celltree::PYRAMID; break;
            VTK_WEDGE: type = celltree::PRISM; break;
        }
        (*cit++) = celltree::mesh::cell(type, counter);
        counter += nptsincell;
    }
    unsigned int* indices = new unsigned int[counter];
    std::copy(ids.begin(), ids.end(), indices);

    *amesh = new celltree::mesh(npts, ncells, counter, points, cells, indices);

    *scalars = nullptr;
    if (vtk->GetPointData()->GetScalars() != NULL) {
        float* data = new float[npts];
        for (size_t i=0; i<npts; ++i) {
            data[i] = *vtk->GetPointData()->GetScalars()->GetTuple(i);
        }
        *scalars = new celltree::variable(1, npts, data);
    }

    *vectors = nullptr;
    if (vtk->GetPointData()->GetVectors() != NULL) {
        float* data = new float[3*npts];
        double vec[3];
        for (size_t i=0; i<npts; ++i) {
            vtk->GetPointData()->GetVectors()->GetTuple(i, vec);
            std::copy(vec, vec+3, &data[3*i]);
        }
        *vectors = new celltree::variable(3, npts*3, data);
    }
}
*/

VTK_SMART(vtkUnstructuredGrid) load_data(const std::string& filename) {



    celltree::dataset* ds;
    try {
        ds = celltree::dataset::create(filename);
    }
    catch( std::exception& e) {
        std::cout << "exception caught while importing " << filename << ": " << e.what() << '\n';
        return false;
    }
    *m = ds->read_mesh();
    std::cout << "Mesh read\n";
    *s = ds->read_scalar_variable(0, "pressure");
    std::cout << "Pressure read\n";
    *v = ds->read_vector_variable(0, "velocity");
    std::cout << "Velocity read\n";
    return true;
}

typedef Eigen::Matrix<unsigned int, 3, 1> uint3;
typedef Eigen::Matrix<double, 3, 1> vec3;
typedef Eigen::Matrix<float, 3, 1> fvec3;
typedef Eigen::Matrix<int, 3, 1> int3;
typedef std::pair< vec3, vec3 > box3;
typedef Eigen::Matrix<double, 6, 1> vec6;

// parameters
size_t npts;
uint3 resolution;
double radius;
bool verbose;
std::string in_name, in_boundary, out_name, geom_name;

void initialize(int argc, char* argv[])
{
    namespace xcl = xavier::command_line;

    xcl::option_traits
        required(true, false, "Required Options"),
    optional(false, false, "Optional Group");
    xcl::option_parser parser(argv[0],
    "Decimate an unstructured grid and remesh remaining points. Boundary is preserved.");
    verbose = false;
    in_parallel = true;

    try {
        parser.use_short_symbols(true);
        parser.use_brackets(true);
        parser.add_value("volume", in_name, "Input unstructured grid filename", required);
        parser.add_value("boundary", in_boundary, "File containing the boundary geometry", optional)
        parser.add_value("output", out_name, "Output base name", required);
        parser.add_value("number", npts, "Number of points to retain", required);
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

typedef xavier::mesh mesh_type;
typedef xavier::mesh_type::face face_type;
typedef xavier::mesh_type::cell_kind cell_kind;
typedef xavier::mesh_type::vertex_type vertex_type;

unsigned int vtk_celltype_to_dlr_celltype(unsigned char vtkid) {
    if (vtkid == VTK_TRIANGLE) return cell_kind::TRIANGLE;
    else if (vtkid == VTK_QUAD) return cell_kind::QUADRILATERAL;
    else if (vtkid == VTK_HEXAHEDRON) return cell_kind::HEXAHEDRON;
    else if (vtkid == VTK_TETRAHEDRON) return cell_kind::TETRAHEDRON;
    else if (vtkid == VTK_WEDGE) return cell_kind::PRISM;
    else if (vtkid == VTK_PYRAMID) return cell_kind::PYRAMID;
}

int main(int argc, char* argv[]) {
    const int dim = 3;

    initialize(argc, argv);

    std::vector<long int> boundary_ids;
    std::vector<face_type> boundary_faces;
    mesh_type amesh;

    if (xavier::filename::extension(in_name) == "grid") {
        amesh.load(in_name);
    }
    else {
        VTK_SMART(vtkDataSet) _data = vtk_utils::readVTK(in_name);
        VTK_SMART(vtkUnstructuredGrid) data = vtkUnstructuredGrid::SafeDownCast(_data);
        if (!data) {
            throw std::runtime_error("Invalid grid type. Only unstructured meshes are supported.")
        }
        VTK_SMART(vtkUnsignedCharArray) cell_types = data->GetCellTypesArray();
        VTK_SMART(vtkDataArrayTemplate) offsets = data->GetOffsetsArray();
        VTK_SMART(vtkDataArray) indices = data->GetConnectivityArray();
        VTK_SMART(vtkPoints) points = data->GetPoints();
        std::vector<std::pair<cell_kind, long>> cell_offsets(offsets->size());
        std::vector<vertex_type> vertices(points->GetNumberOfPoints());
        std::vector<long> cell_indices(indices->GetNumberOfTuples());
        for (size_t i=0; i<cell_offsets.size()-1; ++i) {
            int _type = cell_types->GetTuple1(i);
            int _offset = offsets->GetTuple1(i);
            cell_offsets[i].first = vtk_celltype_to_dlr_celltype(_type);
            cell_offsets[i].second = _offset;
        }
        // last entry points 1 past the last valid cell index
        cell_offsets.back().first = 0; // dummy cell type
        cell_offsets.back().second = indices->GetNumberOfTuples();

        for (size_t i=0; i<cell_indices.size(); ++i) {
            cell_indices[i] = indices->GetTuple1(i);
        }
        for (size_t i=0; i<vertices.size(); ++i) {
            double c[3];
            points->GetPoint(i, c);
            for (int j=0; j<3; ++j) vertices[i][j] = c[j];
        }
        amesh = mesh_type(vertices, cell_indices, cell_offsets);
    }

    if (!in_boundary) {
        boundary_faces = amesh.get_boundary();
    }
    else {
        VTK_SMART(vtkDataSet) _boundary = vtk_utils::readVTK(in_boundary);
        boundary_faces.resize(_boundary->GetNumberOfCells());
        for (size_t i=0; i<boundary_faces.size(); ++i) {
            VTK_SMART(vtkCell) c = _boundary->GetCell(i);
            VTK_SMART(vtkIdList) ids = c->GetPointIds();
            int nids = ids->GetNumberOfIds();
            if (nids == 3) boundary_faces[i] = face_type(ids->GetId(0), ids->GetId(1), ids->GetId(2));
            else if (nids == 4) {
                boundary_faces[i] = face_type(ids->GetId(0), ids->GetId(1), ids->GetId(2), ids->GetId(3));
            }
        }
    }

    std::set<long int> unique_ids;
    for (auto it=boundary_faces.begin(); it!=boundary_faces.end(); ++it) {
        for (auto jt=it->begin(); jt!=it.end(); ++jt) {
            unique_ids.insert(*jt);
        }
    }
    std::copy(unique_ids.begin(), unique_ids.end(), boundary_ids.begin());
    std::cout << "There are " << boundary_ids.size() << " vertices on the boundary and "
              << boundary_faces.size() << " cells\n";
              
    return 0;
}
