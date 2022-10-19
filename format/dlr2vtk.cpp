// STL
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
// Boost
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
// VTK
#include <vtkCellArray.h>
#include <vtkDataSet.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkUnstructuredGrid.h>
// nvis
#include <math/bounding_box.hpp>
#include <math/fixed_vector.hpp>

#include "DLRreader.hpp"
#include <VTK/vtk_utils.hpp>

typedef nvis::fvec3                     point_type;
typedef nvis::bounding_box<point_type>  box_type;
typedef xavier::DLRreader::cell_type    cell_type;
typedef std::pair<cell_type, long>      cell_entry;

typedef std::pair< std::vector<double>, std::string > pair_type;
std::vector< pair_type > all_vector_fields;
std::vector< pair_type > all_scalar_fields;

// helper functions
struct cell_array_helper {
    cell_array_helper(const std::vector<long>& ids,
                      const std::vector<cell_entry>& types)
        : __ids(ids), __types(types) {

        std::cerr << ids.size() << " vertex ids in input\n";
        std::cerr << types.size() << " cell types in input\n";

        std::set<long> all_ids;
        for (size_t i=0 ; i<ids.size() ; ++i) all_ids.insert(ids[i]);
        std::cout << "there are " << all_ids.size() << " unique indices in mesh\n";
        std::cout << "there are " << types.size()-1 << " cells in mesh\n";
    }

    long start(const long& cell_id) const {
        return __types[cell_id].second;
    }

    long end(const long& cell_id) const {
        return __types[cell_id+1].second;
    }

    cell_type type(const long& cell_id) const {
        return __types[cell_id].first;
    }

    long vertex(const long& loc_id, const long& cell_id) const {
        return __ids[this->start(cell_id) + loc_id];
    }

    long vertex(const long& id) const {
        return __ids[id];
    }

    long size(const long& cell_id) const {
        return this->end(cell_id) - this->start(cell_id);
    }

    long nb_cells() const {
        return __types.size()-1;
    }

    const std::vector<long>&       __ids;
    const std::vector<cell_entry>& __types;
};

// Compute connected components of a mesh using boost's
// implementation of graph's DFS algorithm
void split_cc(std::vector<std::set<long> >& components,
              const cell_array_helper& helper)
{
    using namespace boost;
    typedef adjacency_list<vecS, vecS, undirectedS> graph_type;

    graph_type G;
    for (long cell_id=0 ; cell_id<helper.nb_cells(); ++cell_id) {
        for (long id=helper.start(cell_id) ; id<helper.end(cell_id) ; ++id) {
            long next = id+1;
            if (next >= helper.end(cell_id)) next = helper.start(cell_id);
            long vert_id0 = helper.vertex(id);
            long vert_id1 = helper.vertex(next);

            add_edge(vert_id0, vert_id1, G);
        }
    }

    std::vector<int> cc(num_vertices(G));
    int num = connected_components(G, &cc[0]);

    std::vector<std::set<long> > tmp_cc(num);
    for (long cell_id=0 ; cell_id<helper.nb_cells() ; ++cell_id) {
        // component id of first vertex is component id of entire cell
        long vert_id = helper.vertex(0, cell_id);
        size_t cc_id = cc[vert_id];
        tmp_cc[cc_id].insert(cell_id);
    }

    // Note: using vecS as vertex container is fast but boost
    // internally considers a contiguous set of indices as composing
    // the graph. In our case, those indices correspond to the vertices
    // that do not belong to the boundary. These vertices are considered
    // singletons in the computation of connected components. Hence most
    // of the CC found are trivial.
    components.clear();
    for (size_t i=0 ; i<tmp_cc.size() ; ++i) {
        if (!tmp_cc[i].empty()) {
            components.push_back(std::set<long>());
            components.back().swap(tmp_cc[i]);
        }
    }
}

// Map indices of used points to contiguous range
void reindex(std::map<long, long>& old2new,
             const std::set<long>& cc,
             const cell_array_helper& helper)
{
    old2new.clear();
    long count = 0;
    std::map<long, long>::iterator m_it;
    std::set<long>::const_iterator s_it;
    if (!cc.empty()) {
        for (s_it=cc.begin() ; s_it!=cc.end() ; ++s_it) {
            long cell_id = *s_it;
            for (long loc_id=0 ; loc_id<helper.size(cell_id) ; ++loc_id) {
                long vert_id = helper.vertex(loc_id, cell_id);
                m_it = old2new.find(vert_id);
                if (m_it == old2new.end()) {
                    old2new[vert_id] = count++;
                }
            }
        }
    } else {
        for (long cell_id=0 ; cell_id<helper.nb_cells() ; ++cell_id) {
            for (long loc_id=0 ; loc_id<helper.size(cell_id) ; ++loc_id) {
                long vert_id = helper.vertex(loc_id, cell_id);
                m_it = old2new.find(vert_id);
                if (m_it == old2new.end()) {
                    old2new[vert_id] = count++;
                }
            }
        }
    }
}

vtkPoints*
vtk_points(const std::vector<point_type>& pos)
{
    vtkFloatArray* __coords = vtkFloatArray::New();
    __coords->SetNumberOfComponents(3);
    __coords->SetNumberOfTuples(pos.size());
    for (size_t i=0 ; i<pos.size() ; ++i) {
        __coords->SetTuple(i, &(pos[i][0]));
    }
    vtkPoints* __points = vtkPoints::New();
    __points->SetData(__coords);
    __coords->Delete();
    return __points;
}

VTK_SMART(vtkUnstructuredGrid)
vtk_cellarray(const cell_array_helper& helper)
{
    std::cout << "in vtk_cellarray...\n";
    using namespace xavier;
    VTK_CREATE(vtkCellArray, __cells);
    __cells->SetNumberOfCells(helper.nb_cells());
    for (long cell_id=0 ; cell_id<helper.nb_cells() ; ++cell_id) {
        __cells->InsertNextCell(helper.size(cell_id));
        for (long loc_id=0 ; loc_id<helper.size(cell_id) ; ++loc_id) {
            __cells->InsertCellPoint(helper.vertex(loc_id, cell_id));
        }
    }

    std::cout << "cell points inserted\n";
    std::cout << "there are " << __cells->GetNumberOfCells() << " cells"
    << " and " << __cells->GetNumberOfConnectivityIds() << " connectivity ids\n";

    VTK_CREATE(vtkUnsignedCharArray, __types);
    VTK_CREATE(vtkIdTypeArray, __locations);
    __types->SetNumberOfComponents(1);
    __types->SetNumberOfTuples(helper.nb_cells());
    __locations->SetNumberOfComponents(1);
    __locations->SetNumberOfTuples(helper.nb_cells());
    std::cout << "arrays initialized\n";
    for (long cell_id=0 ; cell_id<helper.nb_cells() ; ++cell_id) {
        unsigned char type_name;
        switch(helper.type(cell_id)) {
            case DLRreader::TRIANGLE:
                type_name = VTK_TRIANGLE;
                break;
            case DLRreader::QUADRILATERAL:
                type_name = VTK_QUAD;
                break;
            case DLRreader::TETRAHEDRON:
                type_name = VTK_TETRA;
                break;
            case DLRreader::HEXAHEDRON:
                type_name = VTK_HEXAHEDRON;
                break;
            case DLRreader::PRISM:
                type_name = VTK_WEDGE;
                break;
            case DLRreader::PYRAMID:
                type_name = VTK_PYRAMID;
                break;
            default:
                throw std::runtime_error("invalid cell type");
        }
        __types->SetValue(cell_id, type_name);
        __locations->SetValue(cell_id, helper.start(cell_id));
    }
    std::cout << "cell types defined\n";

    VTK_CREATE(vtkUnstructuredGrid, grid);
    grid->SetCells(__types, __locations, __cells);
    std::cout << "grid cells set\n";
    // std::cout << "grid is now:\n";
    // grid->PrintSelf(std::cout, vtkIndent(0));
    return grid;
}

std::string me;
void usage(const std::string& message="")
{
    if (!message.empty()) {
        std::cerr << "\nERROR: " << message << '\n';
    }
    std::cerr << '\n'
              << "DESCRIPTION: Convert a DLR dataset from NetCDF format to\n"
              << "VTK format. The conversion can be restricted to the mesh\n"
              << "and / or to the simulation boundary. In the latter case,\n"
              << "the data can be segmented by connectivity, in which case\n"
              << "separate files will be generated for each connected component.\n"
              << "USAGE: " << me << " [parameters] [options]\n"
              << "PARAMETERS:\n"
              << " -g | --grid <string>       DLR grid filename\n"
              << " -o | --output <string>     Output VTK basename\n"
              << "OPTIONS:\n"
              << " -h | --help                Print this information\n"
              << " -d | --data <string>       DLR data filename\n"
              << " -a | --attribute <string>  Name of requested attribute\n"
              << " -b | --boundary            Restrict conversion to grid boundary\n"
              << "    | --velname <string>    Name of velocity components (x_<name>, y_<name>, z_<name>)\n"
              << " -s | --split               Split connected components into different files\n"
              << " -r | --reindex             Reindex used vertices to reduce file size\n"
              << " -m | --mesh                Restrict conversion to mesh information\n"
              << " -v | --verbose             Activate verbose mode\n"
              << std::endl;

    exit(!message.empty());
}

void transform(std::vector<point_type>& new_points,
               std::vector<long>& new_ids,
               std::vector<std::pair<cell_type, long> >& new_types,
               std::vector< pair_type >& new_scalars,
               std::vector< pair_type >& new_vectors,
               const std::set<long>& cc,
               const std::map<long, long>& old2new,
               const std::vector<point_type>& points,
               const cell_array_helper& helper)
{
    std::map<long, long>::const_iterator m_it;
    if (!old2new.empty()) {
        new_points.resize(old2new.size());
        for (m_it=old2new.begin() ; m_it!=old2new.end() ; ++m_it) {
            new_points[m_it->second] = points[m_it->first];
        }
        if (!all_vector_fields.empty()) {
            new_vectors.resize(all_vector_fields.size());
            for (int i=0; i<new_vectors.size(); ++i) {
                new_vectors[i].first.resize(3*old2new.size());
                for (m_it=old2new.begin() ; m_it!=old2new.end() ; ++m_it) {
                    new_vectors[i].first[3*m_it->second  ] =
                        all_vector_fields[i].first[3*m_it->first  ];
                    new_vectors[i].first[3*m_it->second+1] =
                        all_vector_fields[i].first[3*m_it->first+1];
                    new_vectors[i].first[3*m_it->second+2] =
                        all_vector_fields[i].first[3*m_it->first+2];
                }
                new_vectors[i].second = all_vector_fields[i].second;
            }
        }
        if (!all_scalar_fields.empty()) {
            new_scalars.resize(all_scalar_fields.size());
            for (int i=0; i<new_scalars.size(); ++i) {
                new_scalars[i].first.resize(old2new.size());
                for (m_it=old2new.begin() ; m_it!=old2new.end() ; ++m_it) {
                    new_scalars[i].first[m_it->second] = all_scalar_fields[i].first[m_it->first];
                }
                new_scalars[i].second = all_scalar_fields[i].second;
            }
        }
    }

    new_ids.clear();
    new_types.clear();
    std::set<long>::const_iterator s_it;
    for (s_it=cc.begin() ; s_it!=cc.end() ; ++s_it) {
        size_t old_cell_id = *s_it;
        cell_type type = helper.type(old_cell_id);
        new_types.push_back(cell_entry(type, new_ids.size()));
        if (!old2new.empty()) {
            for (long loc_id=0; loc_id<helper.size(old_cell_id) ; ++loc_id) {
                long vert_id = helper.vertex(loc_id, old_cell_id);
                m_it = old2new.find(vert_id);
                new_ids.push_back(m_it->second);
            }
        } else {
            for (long loc_id=0 ; loc_id<helper.size(old_cell_id) ; ++loc_id) {
                new_ids.push_back(helper.vertex(loc_id, old_cell_id));
            }
        }
    }
    new_types.push_back(cell_entry(cell_type::NONE, new_ids.size()));

}

int main(int argc, char* argv[])
{
    me = argv[0];
    std::string grid_filename, data_filename, out_filename, comp_name="velocity";
    std::vector<std::string> att;
    bool boundary_only = false;
    bool mesh_only = false;
    bool verbose = false;
    bool split = false;
    bool re = false;

    for (int i=1 ; i<argc ; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            usage();
        } else if (arg == "-g" || arg == "--grid") {
            if (i == argc-1) {
                usage("missing grid file name");
            }
            grid_filename = argv[++i];
        } else if (arg == "-d" || arg == "--data") {
            if (i == argc-1) {
                usage("missing data file name");
            }
            data_filename = argv[++i];
        } else if (arg == "-a" || arg == "--attribute") {
            if (i == argc-1) {
                usage("missing attribute name");
            }
            att.push_back(argv[++i]);
        } else if (arg == "--velname") {
            if (i == argc-1) {
                usage("missing velocity component name");
            }
            comp_name = argv[++i];
        } else if (arg == "-o" || arg == "--output") {
            if (i == argc-1) {
                usage("missing output file name");
            }
            out_filename = argv[++i];
        } else if (arg == "-b" || arg == "--boundary") {
            boundary_only = true;
        } else if (arg == "-m" || arg == "--mesh") {
            mesh_only = true;
        } else if (arg == "-s" || arg == "--split") {
            split = true;
        } else if (arg == "-r" || arg == "--reindex") {
            re = true;
        } else if (arg == "-v" || arg == "--verbose") {
            verbose = true;
        } else {
            usage("unrecognized argument: " + arg);
        }
    }

    // sanity check
    if (grid_filename.empty()) {
        usage("missing grid file");
    }
    if (out_filename.empty()) {
        usage("missing output file name");
    }
    if (data_filename.empty() && att.size())
        usage("both value(s) and geometry were requested\n"
              "yet file containing values is missing");
    if (!boundary_only && (split || re)) {
        std::cout << "WARNING: split and reindex option only valid\n"
                  << "with boundary restriction\n";
        split = false;
        re = false;
    }

    xavier::DLRreader reader(grid_filename, data_filename);
    std::vector<point_type> points;
    std::vector<long> cell_indices;
    std::vector<cell_entry> cell_types;
    reader.read_mesh(boundary_only, points, cell_indices, cell_types, verbose);

    cell_array_helper helper(cell_indices, cell_types);
    std::cout << "vertices range: " << *helper.__ids.begin() << " -> " << *helper.__ids.rbegin() << '\n';

    if (!mesh_only) {
        for (int i=0; i<att.size() ; ++i) {
            if (att[i] == "velocity" || att[i] == "vel" || att[i] == "v") {
                std::vector<double> val0, val1, val2;
                reader.read_data_from_file(data_filename, "x_" + comp_name, val0, verbose);
                reader.read_data_from_file(data_filename, "y_" + comp_name, val1, verbose);
                reader.read_data_from_file(data_filename, "z_" + comp_name, val2, verbose);
                std::cout << "val0.size()=" << val0.size() << '\n'
                << "val1.size()=" << val1.size() << '\n'
                << "val2.size()=" << val2.size() << '\n';
                assert(val0.size() == val1.size() && val1.size() == val2.size());
                all_vector_fields.push_back(pair_type());
                all_vector_fields.back().second = "velocity";
                std::vector<double>& values = all_vector_fields.back().first;
                values.resize(3*val0.size());
                long k=0;
                for (long n=0; n<val0.size(); ++n) {
                    values[k++] = val0[n];
                    values[k++] = val1[n];
                    values[k++] = val2[n];
                }

                std::cout << "x_velocity range: "
                    << *std::min_element(val0.begin(), val0.end()) << " - "
                    << *std::max_element(val0.begin(), val0.end()) << '\n';

                std::cout << "y_velocity range: "
                    << *std::min_element(val1.begin(), val1.end()) << " - "
                    << *std::max_element(val1.begin(), val1.end()) << '\n';

                std::cout << "z_velocity range: "
                    << *std::min_element(val2.begin(), val2.end()) << " - "
                    << *std::max_element(val2.begin(), val2.end()) << '\n';

            }
            else {
                all_scalar_fields.push_back(pair_type());
                std::vector<double>& values = all_scalar_fields.back().first;
                reader.read_data(att[i], values, verbose);
                all_scalar_fields.back().second = att[i];
            }
        }
    }

    box_type bounds;
    for (size_t i=0 ; i<points.size() ; ++i) {
        bounds.add(points[i]);
    }
    if (true || verbose) {
        std::cout << "there are " << points.size() << " points in input mesh\n"
                  << "associated bounding box: " << bounds << '\n';
    }

    // compute connected components, if requested
    std::vector<std::set<long> > ccs;
    if (split) {
        if (verbose) {
            std::cout << "computing connected components..." << std::flush;
        }
        split_cc(ccs, helper);
        if (verbose) {
            std::cout << " done: " << ccs.size() << " connected components found\n";
            std::cout << "cc0 has size " << ccs[0].size() << '\n';
            std::cout << "cc1 has size " << ccs[1].size() << '\n';
        }
    } else {
        std::vector<long> _ids(cell_types.size()-1);
        for (size_t i=0 ; i<_ids.size() ; ++i) {
            _ids[i] = i;
        }
        ccs.push_back(std::set<long>(_ids.begin(), _ids.end()));
    }

    // reindex used vertices, if requested
    std::vector<std::map<long, long> > index_maps(std::max(ccs.size(), (size_t)1));
    if (re) {
        for (size_t i=0 ; i<ccs.size() ; ++i) {
            if (verbose) {
                std::cout << "reindexing component " << i+1 << " from "
                    << ccs.size() << "... ";
            }
            reindex(index_maps[i], ccs[i], helper);
            if (verbose) {
                std::cout << "done\n";
            }
        }
    }

    // export individual submeshes
     for (size_t i=0 ; i<ccs.size() ; ++i) {
        std::cout << "processing connected component #" << i << '\n';
        std::string name = out_filename;
        // remove existing vtk extension
        std::string ext = xavier::filename::extension(name);
        name = xavier::filename::remove_extension(name);
        if (ccs.size()>1) {
            std::ostringstream os;
            os << name << "_" << i+1 << "_of_" << ccs.size() << '.' << ext;
            name = os.str();
        }
        else name += '.' + ext;

        std::cout << "name = " << name << '\n';

        VTK_SMART(vtkUnstructuredGrid) grid;
        if (split || re) {
            std::vector<point_type> new_points;
            std::vector<long> new_ids;
            std::vector<std::pair<cell_type, long> > new_types;
            std::vector<pair_type> new_vectors;
            std::vector<pair_type> new_scalars;
            transform(new_points, new_ids, new_types, new_scalars, new_vectors,
                      ccs[i], index_maps[i], points, helper);
            grid = vtk_cellarray(cell_array_helper(new_ids, new_types));
            if (re) grid->SetPoints(vtk_points(new_points));
            else grid->SetPoints(vtk_points(points));

            for (int i=0; i<new_vectors.size(); ++i)
                vtk_utils::add_vectors_from_numbers<std::vector<double>, vtkUnstructuredGrid*,3>(grid, new_vectors[i].first, true, new_vectors[i].second);
            for (int i=0; i<new_scalars.size(); ++i)
                vtk_utils::add_scalars(grid, new_scalars[i].first, true, new_scalars[i].second);
        }
        else {
            std::cout << "about to create grid...\n";
            grid = vtk_cellarray(helper);
            std::cout << "grid created. About to add points...\n";
            grid->SetPoints(vtk_points(points));
            std::cout << "points added. About to add vectors...\n";
            for (int i=0; i<all_vector_fields.size(); ++i) {
                std::cout << "adding vector field #" << i << " of size " << all_vector_fields[i].first.size() << std::endl;
                vtk_utils::add_vectors_from_numbers<std::vector<double>, vtkUnstructuredGrid*,3>(grid, all_vector_fields[i].first, true, all_vector_fields[i].second);
            }
            std::cout << "vectors added. About to add scalars...\n";
            for (int i=0; i<all_scalar_fields.size(); ++i) {
                std::cout << "adding scalar field #" << i << "of size " << all_scalar_fields[i].first.size() << std::endl;
                vtk_utils::add_scalars(grid, all_scalar_fields[i].first, true, all_scalar_fields[i].second);
            }
            std::cout << "scalars added\n";
        }

        double* __bounds = grid->GetBounds();
        box_type vtk_bounds(point_type(__bounds[0], __bounds[2], __bounds[4]),
                            point_type(__bounds[1], __bounds[3], __bounds[5]));

        if (verbose) {
            std::cout << "after conversion, exported grid has following properties:\n";
            // grid->PrintSelf(std::cout, vtkIndent(0));
            std::cout << '\n'
                      << "# points:     " << grid->GetNumberOfPoints() << '\n'
                      << "# cells:      " << grid->GetCells()->GetNumberOfCells() << '\n'
                      << "bounding box: " << vtk_bounds << '\n';
            if (grid->GetPointData() != nullptr) {
                if (grid->GetPointData()->GetVectors() != nullptr)
                    std::cout << "number of vectors: " << grid->GetPointData()->GetVectors()->GetNumberOfTuples() << '\n';
                if (grid->GetPointData()->GetScalars() != nullptr)
                    std::cout << "number of scalars: " << grid->GetPointData()->GetScalars()->GetNumberOfTuples() << '\n';
            }

            // grid->PrintSelf(std::cout, vtkIndent(0));
        }

        vtk_utils::saveVTK(grid, name);
        std::cout << name << " has been exported\n";
    }

    return 0;
}
