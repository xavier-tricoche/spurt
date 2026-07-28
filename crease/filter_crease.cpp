#include <vtk/vtk_utils.hpp>
#include <vtk/filter_ccs.hpp>
#include <vtk/vtk_mesh.hpp>
#include <misc/cxxopts.hpp>
#include <string>
#include <math/types.hpp>
#include <vtkTransformFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkAOSDataArrayTemplate.h>
#include <vtkCleanPolyData.h>
#include <math/stat.hpp>
#include <data/locator.hpp>
#include <format/filename.hpp>
using namespace spurt;

constexpr double minus_infinity = std::numeric_limits<double>::lowest();
constexpr double plus_infinity = std::numeric_limits<double>::max();
constexpr double epsilon = std::numeric_limits<double>::min();
constexpr int very_large = std::numeric_limits<int>::max();

template<typename T = double>
inline spurt::small_vector<T, 3> cross_product_from_three(const spurt::small_vector<T, 3>& a, 
                                                const spurt::small_vector<T, 3>& b, 
                                                const spurt::small_vector<T, 3>& c) {
    spurt::small_vector<T, 3> ab = b-a;
    spurt::small_vector<T, 3> ac = c-a;
    return spurt::cross(ab, ac);
}

template<typename T = double>
inline spurt::small_vector<T, 3> normal_from_three(const spurt::small_vector<T, 3>& a,
                                         const spurt::small_vector<T, 3>& b, 
                                         const spurt::small_vector<T, 3>& c) {
    auto n = cross_product_from_three(a, b, c);
    T len = spurt::norm(n);
    if (len == 0.) {
        // std::cout << "Zero normal: input points were " << a << ", " << b << ", " << c << '\n';
        return n;
    }
    else return n / len;
}

template<typename T = double>
inline spurt::small_vector<T, 3> normal_from_ids(vtkIdType p1, vtkIdType p2, vtkIdType p3, 
                                                  VTK_SMART(vtkPolyData) polydata) {
    auto coords = polydata->GetPoints();
    double t1[3];
    double t2[3];
    double t3[3];
    coords->GetPoint(p1, t1);
    coords->GetPoint(p2, t2);
    coords->GetPoint(p3, t3);
    spurt::small_vector<T, 3> a(t1[0], t1[1], t1[2]), b(t2[0], t2[1], t2[2]), c(t3[0], t3[1], t3[2]);
    auto n = normal_from_three(a, b, c);
    return n;
}

inline bool check_orientation(const std::array<vtkIdType, 2>& edge, vtkIdList* pts) {
    auto n = pts->GetNumberOfIds();
    if (n < 2) return false;
    for (vtkIdType eid=0; eid<n; eid++) {
        vtkIdType p0 = pts->GetId(eid);
        vtkIdType p1 = pts->GetId((eid + 1) % n);
        if (p0 == edge[0] && p1 == edge[1]) {
            // oriented incorrectly (e.g., in same direction)
            return false;
        }
        if (p0 == edge[1] && p1 == edge[0]) {
            // opposite orientation: correct
            return true;
        }
    }
    return false;
}

VTK_SMART(vtkPolyData) translate(VTK_SMART(vtkPolyData) input, const vec3& t)
{
    if (spurt::norm(t) == 0) return input;
    VTK_CREATE(vtkTransform, translate);
    translate->Identity();
    translate->Translate(t[0], t[1], t[2]);
    VTK_CREATE(vtkTransformFilter, filter);
    filter->SetTransform(translate);
    filter->SetInputData(input);
    filter->Update();
    return vtkPolyData::SafeDownCast(filter->GetOutput());
}

std::vector<double> array2vector(VTK_SMART(vtkAbstractArray) _array) {
    VTK_SMART(vtkDataArray) array = vtkDataArray::SafeDownCast(_array);
    std::vector<double> values;
    long ntuples = array->GetNumberOfTuples();
    int ncomps = array->GetNumberOfComponents();
    std::cout << "there are " << ntuples << " tuples and " << ncomps << " comps\n";
    for (long i=0; i<ntuples; ++i) {
        double* ptr = array->GetTuple(i);
        if (ncomps == 1) {
            values.push_back(*ptr);
        }
        else {
            values.push_back(sqrt(std::inner_product(ptr, ptr+ncomps, ptr, 0)));
        }
    }
    std::cout << "values contains " << values.size() << " entries\n";
    return values;
}

void array2vector2(std::vector<double>& values, VTK_SMART(vtkAbstractArray) _array) {
    VTK_SMART(vtkDataArray) array = vtkDataArray::SafeDownCast(_array);
    long ntuples = array->GetNumberOfTuples();
    int ncomps = array->GetNumberOfComponents();
    std::cout << "there are " << ntuples << " tuples and " << ncomps << " comps\n";
    values.resize(ntuples);
    for (long i=0; i<ntuples; ++i) {
        double* ptr = array->GetTuple(i);
        if (ncomps == 1) {
            values[i] = (*ptr);
        }
        else {
            values[i] = sqrt(std::inner_product(ptr, ptr+ncomps, ptr, 0));
        }
    }
    std::cout << "values contains " << values.size() << " entries\n";
}

VTK_SMART(vtkPolyData) filter_cells(VTK_SMART(vtkPolyData) input, 
                                    const std::vector<bool>& selected)
{
    std::cout << "entering filter_cells\n";
    VTK_SMART(vtkCellArray) vertices = input->GetVerts();
    VTK_SMART(vtkCellArray) lines = input->GetLines();
    VTK_SMART(vtkCellArray) polys = input->GetPolys();
    int nverts = input->GetNumberOfVerts();
    int nlines = input->GetNumberOfLines();
    int npolys = input->GetNumberOfPolys();
    assert(nverts + nlines + npolys == selected.size());

    VTK_CREATE(vtkCellArray, newverts);
    VTK_CREATE(vtkCellArray, newlines);
    VTK_CREATE(vtkCellArray, newpolys);

    VTK_SMART(vtkCellData) celldata = input->GetCellData();
    int natts = celldata->GetNumberOfArrays();
    std::cout << "there are " << natts << " cell attributes\n";
    std::vector< VTK_SMART(vtkDataArray) > out_attributes(natts);
    for (int i=0; i<natts; ++i) {
        VTK_SMART(vtkDataArray) data = celldata->GetArray(i);
        std::cout << "attribute #" << i << " is called " << data->GetName() << '\n';
        std::cout << "this data set contains " << data->GetNumberOfTuples() << " tuples and each tuple contains " << data->GetNumberOfComponents() << " coefficients\n";
        out_attributes[i] = vtkDataArray::CreateDataArray(data->GetDataType());
        out_attributes[i]->SetName(data->GetName());
        out_attributes[i]->SetNumberOfComponents(data->GetNumberOfComponents());
        std::cout << "copy attribute " << i << " is called " << out_attributes[i]->GetName() << '\n';
        std::cout << "It contains " << out_attributes[i]->GetNumberOfTuples() << " tuples and each tuple contains " << out_attributes[i]->GetNumberOfComponents() << " coefficients\n";
    }

    // Note: what about strips?
    for (int what=0; what<3; ++what)
    {
        VTK_SMART(vtkCellArray) old_cells;
        VTK_SMART(vtkCellArray) new_cells;
        int offset = 0;
        if (what == 0) { // Verts
            old_cells = vertices;
            new_cells = newverts;
        }
        else if (what == 1) { // Lines
            old_cells = lines;
            offset = nverts;
            new_cells = newlines;
        } 
        else { // Polys
            old_cells = polys;
            offset = nverts + nlines;
            new_cells = newpolys;
        }
        VTK_CREATE(vtkIdList, alist);
        int cellid = offset;
        for (auto it=old_cells->NewIterator(); 
             !it->IsDoneWithTraversal(); it->GoToNextCell(), ++cellid) {
            if (selected[cellid]) {
                it->GetCurrentCell(alist);
                new_cells->InsertNextCell(alist);
                for (int i=0; i<natts; ++i) {
                    VTK_SMART(vtkDataArray) data = celldata->GetArray(i);
                    out_attributes[i]->InsertNextTuple(data->GetTuple(cellid));
                }
            }
        }
    }

    VTK_CREATE(vtkPolyData, output);
    output->DeepCopy(input);
    output->SetPolys(newpolys);
    output->SetLines(newlines);
    output->SetVerts(newverts);
    for (int i=0; i<natts; ++i) {
        output->GetCellData()->RemoveArray(i);
    }
    for (int i=0; i<natts; ++i) {
        output->GetCellData()->AddArray(out_attributes[i]);
    }
    std::cout << "leaving filter_cells: output contains " 
              << output->GetNumberOfPoints() << " points and " 
              << output->GetNumberOfCells() << " cells\n";
    
    return output;
}

VTK_SMART(vtkPolyData) filter_points(VTK_SMART(vtkPolyData) input, 
                                     const std::vector<bool>& selected)
{
    std::cout << "entering filter_points\n";
    assert(selected.size() == input->GetNumberOfPoints());
    std::vector<bool> selected_cells(input->GetNumberOfCells(), false);
    VTK_CREATE(vtkIdList, alist);
    int nselected = 0;
    for (int cellid=0; cellid<selected_cells.size(); ++cellid) 
    {
        input->GetCellPoints(cellid, alist);
        bool included = true;
        for (const vtkIdType& id : *alist) {
            if (!selected[id]) {
                included = false;
                break;
            }
        }
        selected_cells[cellid] = included;
        if (selected_cells[cellid]) ++nselected;
    }
    std::cout << "leaving filter_points\n";
    std::cout << "There were " << nselected << " selections\n";
    return filter_cells(input, selected_cells);
}

template<typename T>
VTK_SMART(vtkPolyData) filter_by_value(VTK_SMART(vtkPolyData) input,
                                       const std::string& array_name, 
                                       T min, T max, bool cellwise)
{ 
    if (min > max) return input;

    std::cout << "\n\nFiltering a dataset with " 
        << input->GetNumberOfPoints() << " points and " 
        << input->GetNumberOfCells() << " cells by " 
        << array_name << " " << (cellwise ? "cellwise" : "pointwise") << "\n\n";
    std::cout << "value range is " << min << " -> " << max << '\n';

    VTK_SMART(vtkDataArray) values;
    std::cout << "getting data array... " << std::flush;
    if (cellwise) {
        values = input->GetCellData()->GetArray(array_name.c_str());
    }
    else {
        values = input->GetPointData()->GetArray(array_name.c_str());
    }
    std::cout << "done\n";
    std::vector<bool> selected(values->GetNumberOfTuples(), false);
    std::cout << "checking individual tuples... " << std::flush;
    int nselected = 0;
    for (int i=0; i<selected.size(); ++i)
    {
        T v = values->GetTuple1(i);
        selected[i] = (v >= min) && (v <= max);
        if (selected[i]) ++nselected;
    }
    std::cout << "done\n";
    std::cout << "There were " << nselected << " selections\n";
    if (cellwise) return filter_cells(input, selected);
    else return filter_points(input, selected);
}

VTK_SMART(vtkPolyData) shrink(VTK_SMART(vtkPolyData) input)
{
    std::map<long, long> old2new;
    std::vector<long> new2old;

    VTK_SMART(vtkCellArray) vertices = input->GetVerts();
    VTK_SMART(vtkCellArray) lines = input->GetLines();
    VTK_SMART(vtkCellArray) polys = input->GetPolys();
    VTK_CREATE(vtkCellArray, newverts);
    VTK_CREATE(vtkCellArray, newlines);
    VTK_CREATE(vtkCellArray, newpolys);
    VTK_CREATE(vtkPolyData, shrunk);

    // loop over all cells of all cell kinds and reindex vertices
    for (int what=0; what<3; ++what)
    {
        VTK_SMART(vtkCellArray) old_cells;
        VTK_SMART(vtkCellArray) new_cells;
        if (what == 0) {
            old_cells = vertices;
            new_cells = newverts;
        }
        else if (what == 1) {
            old_cells = lines;
            new_cells = newlines;
        }
        else if (what == 2) {
            old_cells = polys;
            new_cells = newpolys;
        }

        VTK_CREATE(vtkIdList, alist);
        for (auto it=old_cells->NewIterator(); 
             !it->IsDoneWithTraversal(); it->GoToNextCell()) {
            it->GetCurrentCell(alist);
            long nids = alist->GetNumberOfIds();
            for (long n=0; n<nids; ++n) {
                long id = alist->GetId(n);
                auto iter = old2new.find(id);
                if (iter == old2new.end()) {
                    old2new[id] = new2old.size();
                    alist->SetId(n, new2old.size());
                    new2old.push_back(id);
                }
                else {
                    alist->SetId(n, iter->second);
                }
            }
            new_cells->InsertNextCell(alist);
        }
    }
    shrunk->SetVerts(newverts);
    shrunk->SetLines(newlines);
    shrunk->SetPolys(newpolys);

    // Update points array
    VTK_SMART(vtkDataArray) oldpos = input->GetPoints()->GetData();
    VTK_CREATE(vtkFloatArray, newpos);
    newpos->SetNumberOfComponents(3);
    newpos->SetNumberOfTuples(new2old.size());
    for (long i=0; i<new2old.size(); ++i) {
        auto ptr = oldpos->GetTuple3(new2old[i]);
        newpos->SetTuple3(i, ptr[0], ptr[1], ptr[2]);
    }
    VTK_CREATE(vtkPoints, newpoints);
    newpoints->SetData(newpos);
    shrunk->SetPoints(newpoints);
    auto nnew = shrunk->GetNumberOfPoints();
    auto nold = input->GetNumberOfPoints();
    std::cout << "shrunk dataset contains " << nnew << " points (from " << nold << ")\n";
    std::cout << "ratio: " << (double)nnew/(double)nold*100. << "%\n";

    // Update point attributes
    VTK_SMART(vtkPointData) pointdata = input->GetPointData();
    int natts = pointdata->GetNumberOfArrays();
    std::cout << "there are " << natts << " points attributes\n";
    VTK_SMART(vtkPointData) newpointdata = shrunk->GetPointData();
    for (int n=0; n<natts; ++n) {
        VTK_SMART(vtkDataArray) data = pointdata->GetArray(n);
        std::cout << "attribute #" << n << " is called " << data->GetName() << '\n';
        std::cout << "this data set contains " << data->GetNumberOfTuples() << " tuples and each tuple contains " << data->GetNumberOfComponents() << " coefficients\n";
        VTK_SMART(vtkDataArray) att = vtkDataArray::CreateDataArray(data->GetDataType());
        att->SetName(data->GetName());
        att->SetNumberOfComponents(data->GetNumberOfComponents());
        std::cout << "copy attribute " << n << " is called " << att->GetName() << '\n';
        std::cout << "It contains " << att->GetNumberOfTuples() << " tuples and each tuple contains " << att->GetNumberOfComponents() << " coefficients\n";
        for (long i=0; i<new2old.size(); ++i) {
            att->InsertNextTuple(data->GetTuple(new2old[i]));
        }
        newpointdata->AddArray(att);
    }

    // Pass along cell attributes (unchanged)
    VTK_SMART(vtkCellData) celldata = input->GetCellData();
    std::cout << "there are " << natts << " cells attributes\n";
    natts = celldata->GetNumberOfArrays();
    for (int n=0; n<natts; ++n) {
        VTK_SMART(vtkAbstractArray) att = celldata->GetAbstractArray(n);
        std::cout << "attribute #" << n << " is called " << att->GetName() << '\n';
        std::cout << "this data set contains " << att->GetNumberOfTuples() << " tuples and each tuple contains " << att->GetNumberOfComponents() << " coefficients\n";
        shrunk->GetCellData()->AddArray(att);
    }
    return shrunk;
}

template<typename T=double>
std::ostream &operator<<(std::ostream &os, const typename std::vector<T>::iterator &i) {
    os << &i;
    return os;
}

/*
   Extract connected components of input mesh while treating sharp and non-
   manifold edges as boundary edges. Individual connected components are 
   reindexed such that no vertex index is shared across CC's. In other words, 
   sharps and non-manifold edges are duplicated to create cuts. This function 
   returns a new polydata object since it modifies both positions and cells.
*/
template<typename T=double>
std::array<VTK_SMART(vtkPolyData), 4> feature_aware_ccs(VTK_SMART(vtkPolyData) input, double angle_threshold=45.)
{
    typedef T                              val_t;
    typedef vtkAOSDataArrayTemplate<val_t> valarray_t;
    typedef vtkCellArray                   cellarray_t;
    typedef vtkIntArray                    intarray_t;
    typedef spurt::small_vector<val_t, 3>  vec3_t;
    typedef vtkIdType                      idx_t;
    typedef vtkIdList                      idxlist_t;

    constexpr val_t _deg_to_rad = 3.14159265358979323846 / 180.;
    val_t threshold = std::cos(angle_threshold * _deg_to_rad);
    std::cout << "threshold=" << threshold << '\n';

    std::cout << "initially, polygons contain " << input->GetPolys()->GetNumberOfCells() << " cells ";
    std::cout << "and " << input->GetPoints()->GetNumberOfPoints() << " points\n";

    // first, clean up the mesh
    spurt::ProgressDisplay progress;

    std::cout << "Cleaning up input mesh...\n";
    VTK_SMART(vtkPolyData) cleaned = vtk_utils::clean_polydata_mesh(input);
    auto polygons = cleaned->GetPolys();
    auto coordinates = cleaned->GetPoints();
    std::cout << "done\n";
    std::cout << "After cleaning,  polygons contain " << polygons->GetNumberOfCells() << " cells ";
    std::cout << "and " << coordinates->GetNumberOfPoints() << " points\n";

    std::vector<vec3_t> normals(polygons->GetNumberOfCells(), {0., 0., 0.});
    VTK_CREATE(idxlist_t, acell_point_ids);
    acell_point_ids->Reserve(20);

    progress.begin(polygons->GetNumberOfCells(), "Computing cell normals");
    size_t nzeros = 0;
    std::set<idx_t> zero_norm_cells;
    VTK_CREATE(vtkDoubleArray, normals_array);
    normals_array->SetName("MyNormals");
    normals_array->SetNumberOfComponents(3);
    for (idx_t i=0; i<polygons->GetNumberOfCells(); i++)
    {
        polygons->GetCellAtId(i, acell_point_ids);
        if (acell_point_ids->GetNumberOfIds() < 3) continue;
        normals[i] = normal_from_ids(acell_point_ids->GetId(0), 
                                     acell_point_ids->GetId(1), 
                                     acell_point_ids->GetId(2), cleaned);
        if (spurt::norm(normals[i]) < 0.5) {
            // std::cout << "Warning: suspicious normal length for cell " << i << ": " << spurt::norm(normals[i]) << '\n';
            ++nzeros;
            zero_norm_cells.insert(i);
        }
        normals_array->InsertNextTuple(normals[i].data());
        progress.update(i);
    }
    progress.end();
    cleaned->GetCellData()->AddArray(normals_array);
    std::cout << "there were " << nzeros << " zero normals out of " << polygons->GetNumberOfCells() << '\n';
    std::cout << "There are " << normals_array->GetNumberOfTuples() << " normal values in array\n";
    {
        VTK_CREATE(vtkXMLPolyDataWriter, writer2);
        writer2->SetFileName("cleaned_with_normals.vtp");
        writer2->SetInputData(cleaned);
        writer2->Write();
    }

    std::vector<std::string> point_data_array_names;
    for (int i=0; i<cleaned->GetPointData()->GetNumberOfArrays(); ++i) {
        point_data_array_names.push_back(cleaned->GetPointData()->GetArray(i)->GetName());
        std::cout << "Point data array " << i << " is " << point_data_array_names.back() << '\n';
    }
    std::vector<std::string> cell_data_array_names;
    for (int i=0; i<cleaned->GetCellData()->GetNumberOfArrays(); ++i) {
        cell_data_array_names.push_back(cleaned->GetCellData()->GetArray(i)->GetName());
        std::cout << "Cell data array " << i << " is " << cell_data_array_names.back() << '\n';
    }

    // compute connected components while watching for sharp angles
    std::map<idx_t, idx_t> point_id_new_to_old;
    std::map<idx_t, idx_t> cell_id_new_to_old;
    std::vector<idx_t> poly_to_cc(polygons->GetNumberOfCells(), -1);
    VTK_CREATE(vtkPolyData, new_dataset);
    VTK_CREATE(vtkPoints, new_coordinates);
    VTK_CREATE(cellarray_t, new_polygons);
    VTK_CREATE(intarray_t, region_ids); // CC ID
    VTK_CREATE(intarray_t, sizes); // CC size
    VTK_CREATE(valarray_t, laciness); // CC laciness
    VTK_CREATE(cellarray_t, boundary_edges);
    VTK_CREATE(cellarray_t, feature_edges);
    VTK_CREATE(valarray_t, sharpness);
    region_ids->SetName("CCIDs");
    region_ids->SetNumberOfComponents(1);
    sizes->SetName("CCsizes");
    sizes->SetNumberOfComponents(1);
    laciness->SetName("CC Laciness");
    laciness->SetNumberOfComponents(1);
    sharpness->SetName("Edge Sharpness");
    sharpness->SetNumberOfComponents(1);
    
    std::map<idx_t, idx_t> vertex_index_remap; // vertex reindexing per CC
    idx_t cc_id = 0;
    progress.begin(polygons->GetNumberOfCells(), "Computing connected components");
    for (idx_t i=0; i<polygons->GetNumberOfCells(); i++)
    {
        progress.update(i);
        if (poly_to_cc[i] != -1) continue; // cell has been processed already
        if (zero_norm_cells.find(i) != zero_norm_cells.end()) continue; // skip zero normal cells

        // We are starting a new CC...
        size_t nedges = 0;
        size_t nvertices = 0;
        size_t nboundary_edges = 0;
        std::vector<idx_t> new_cell_ids;
        std::list<idx_t> queue;
        poly_to_cc[i] = cc_id;
        queue.push_back(i);
        // Clunky
        VTK_CREATE(idxlist_t, new_cell_point_ids);
        VTK_CREATE(idxlist_t, neighbor_cell_point_ids);
        VTK_CREATE(idxlist_t, neighbor_cell_ids);
        VTK_CREATE(idxlist_t, cell_point_ids);
        neighbor_cell_ids->Reserve(20);
        cell_point_ids->Reserve(20);
        new_cell_point_ids->Reserve(20);
        neighbor_cell_point_ids->Reserve(20);
        neighbor_cell_ids->Reserve(20);
        // Clunky
        while (!queue.empty()) {
            idx_t cell_id = queue.front();
            queue.pop_front();
            // Clunky
            cell_point_ids->Reset();
            new_cell_point_ids->Reset();
            // Clunky
            polygons->GetCellAtId(cell_id, cell_point_ids);
            auto cell_normal = normals[cell_id];
            // we discard hanging lines and points
            if (cell_point_ids->GetNumberOfIds() < 3) continue;
            // add current cell to our new mesh
            for (idx_t j=0; j<cell_point_ids->GetNumberOfIds(); j++) {
                idx_t a = cell_point_ids->GetId(j);
                idx_t b = cell_point_ids->GetId((j+1)%cell_point_ids->GetNumberOfIds());
                auto it = vertex_index_remap.find(a);
                if (it == vertex_index_remap.end()) {
                    // first time we see this vertex as part of this CC
                    auto p = coordinates->GetPoint(a);  
                    idx_t new_id = new_coordinates->GetNumberOfPoints();
                    point_id_new_to_old[new_id] = a;
                    new_coordinates->InsertNextPoint(p[0], p[1], p[2]);
                    vertex_index_remap[a] = new_id;
                    new_cell_point_ids->InsertNextId(new_id);
                    ++nvertices;
                }
                else {
                    new_cell_point_ids->InsertNextId(it->second);
                }
                std::array<vtkIdType, 2> edge{a, b};
                ++nedges;
                neighbor_cell_ids->Reset();
                cleaned->GetCellEdgeNeighbors(cell_id, a, b, neighbor_cell_ids);
                if (neighbor_cell_ids->GetNumberOfIds() == 0) {
                    // boundary edge
                    ++nboundary_edges;
                    boundary_edges->InsertNextCell(2, edge.data());
                    continue;
                }
                else if (neighbor_cell_ids->GetNumberOfIds() == 1) {
                    // regular edge: check for sharpness
                    idx_t neighbor_cell_id = neighbor_cell_ids->GetId(0);
                    if (neighbor_cell_id >= poly_to_cc.size()) {
                        std::cout << "Invalid neighbor id: " << neighbor_cell_id << " > " 
                                    << polygons->GetNumberOfCells() << std::endl;
                        ++nboundary_edges;
                        // boundary_edges->InsertNextCell(2, edge.data());
                        continue;
                    }
                    if (poly_to_cc[neighbor_cell_id] != -1) continue; // cell has been processed already
                    neighbor_cell_point_ids->Reset();
                    polygons->GetCellAtId(neighbor_cell_id, neighbor_cell_point_ids);
                    vec3_t neighbor_normal = normals[neighbor_cell_id];
                    if (!check_orientation(edge, neighbor_cell_point_ids)) {
                        neighbor_normal = -neighbor_normal;
                    }
                    val_t dot = spurt::inner(cell_normal, neighbor_normal);
                    if (dot > threshold) {
                        // add to cc
                        poly_to_cc[neighbor_cell_id] = cc_id;
                        queue.push_back(neighbor_cell_id);
                    }
                    else {
                        ++nboundary_edges;
                        // boundary_edges->InsertNextCell(2, edge.data());
                        feature_edges->InsertNextCell(2, edge.data());
                        sharpness->InsertNextTuple1(dot);
                    }
                }
                else if (neighbor_cell_ids->GetNumberOfIds() > 1) {
                    // non-manifold edge: find the  eighbor with the largest dot product
                    val_t maxdot = -1;
                    idx_t best_neigh = -1;
                    for (idx_t k=0; k<neighbor_cell_ids->GetNumberOfIds(); k++) {
                        idx_t neighbor_cell_id = neighbor_cell_ids->GetId(k);
                        if (neighbor_cell_id >= poly_to_cc.size()) {
                            std::cout << "Invalid neighbor id: " << neighbor_cell_id << " > " 
                                      << polygons->GetNumberOfCells()
                                      << std::endl;
                            boundary_edges->InsertNextCell(2, edge.data());
                            continue;
                        }
                        if (poly_to_cc[neighbor_cell_id] != -1) continue; // cell has been processed already
                        neighbor_cell_point_ids->Reset();
                        polygons->GetCellAtId(neighbor_cell_id, neighbor_cell_point_ids);
                        vec3_t neighbor_normal = normals[neighbor_cell_id];
                        if (!check_orientation(edge, neighbor_cell_point_ids)) {
                            neighbor_normal = -neighbor_normal;
                        }
                        val_t dot = spurt::inner(cell_normal, neighbor_normal);
                        if (dot > threshold && dot > maxdot) {
                            maxdot = dot;
                            best_neigh = neighbor_cell_id;
                        }
                    }
                    if (best_neigh != -1) {
                        // add to cc
                        poly_to_cc[best_neigh] = cc_id;
                        queue.push_back(best_neigh);
                    }
                    else {
                        ++nboundary_edges;
                        // boundary_edges->InsertNextCell(2, edge.data());
                        feature_edges->InsertNextCell(2, edge.data());
                        sharpness->InsertNextTuple1(maxdot);
                    }
                }
            }
            // cell has been processed. add it to our new mesh
            auto new_id = new_polygons->InsertNextCell(new_cell_point_ids);
            cell_id_new_to_old[new_id] = i;
            new_cell_ids.push_back(new_id);
            region_ids->InsertNextTuple1(cc_id);
        }
        cc_id++;
        vertex_index_remap.clear();
        val_t _laciness = static_cast<val_t>(nboundary_edges) / static_cast<val_t>(nedges);
        for (auto id : new_cell_ids) {
            laciness->InsertNextTuple1(_laciness);
            sizes->InsertNextTuple1(new_cell_ids.size());
        }
    }
    progress.end();
    new_dataset->SetPoints(new_coordinates);
    new_dataset->SetPolys(new_polygons);
    new_dataset->GetCellData()->AddArray(region_ids);
    new_dataset->GetCellData()->AddArray(laciness);
    new_dataset->GetCellData()->AddArray(sizes);

    VTK_CREATE(vtkPolyData, boundary_dataset);
    boundary_dataset->SetPoints(coordinates);
    boundary_dataset->SetLines(boundary_edges);

    VTK_CREATE(vtkPolyData, feature_dataset);
    feature_dataset->SetPoints(coordinates);
    feature_dataset->SetLines(feature_edges);
    feature_dataset->GetCellData()->AddArray(sharpness);

    // copy initial point-wise and cell-wise attributes
    auto npoints = new_coordinates->GetNumberOfPoints();
    for (const auto& name : point_data_array_names) {
        auto array = cleaned->GetPointData()->GetArray(name.c_str());
        VTK_CREATE(vtkDoubleArray, newarray);
        newarray->SetName(name.c_str());
        newarray->SetNumberOfComponents(array->GetNumberOfComponents());
        newarray->SetNumberOfTuples(npoints);
        for (int j=0; j<npoints; ++j) {
            auto k = point_id_new_to_old[j];
            newarray->SetTuple(j, array->GetTuple(k));
        }
        new_dataset->GetPointData()->AddArray(newarray);
    }
    auto ncells = new_polygons->GetNumberOfCells();
    for (const auto& name : cell_data_array_names) {
        auto array = cleaned->GetCellData()->GetArray(name.c_str());
        VTK_CREATE(vtkDoubleArray, newarray);
        newarray->SetName(name.c_str());
        newarray->SetNumberOfComponents(array->GetNumberOfComponents());
        newarray->SetNumberOfTuples(ncells);
        for (int j=0; j<ncells; ++j) {
            auto k = cell_id_new_to_old[j];
            newarray->SetTuple(j, array->GetTuple(k));
        }
        new_dataset->GetCellData()->AddArray(newarray);
    }

    return { cleaned, new_dataset, boundary_dataset, feature_dataset };
}

int main(int argc, const char* argv[]) 
{
    cxxopts::Options options("filter_crease", "Manipulate crease surface");
    options.add_options()
        ("i,input", "Input filename", cxxopts::value<std::string>())
        ("o,output", "Output filename", cxxopts::value<std::string>())
        // no default values
        ("s,size", "Lower bound on CC size", cxxopts::value<int>())
        ("n,number", "Max number of CCs", cxxopts::value<int>())
        ("value", "Min scalar value", cxxopts::value<std::string>())
        ("strength", "Min ridge strength (>0)", cxxopts::value<std::string>())
        ("angle", "Angle threshold for feature edges", cxxopts::value<double>())
        ("namev", "value name", cxxopts::value<std::string>()->default_value("values"))
        ("names", "strength name", cxxopts::value<std::string>()->default_value("ridge_strength"))
        ("pre", "Apply value filters before CC filters", cxxopts::value<bool>())
        ("t,translate", "Translate mesh", cxxopts::value<std::array<double,3>>())
        ("shrink", "Shrink dataset by removing unused vertices", cxxopts::value<bool>())
        ("stats", "Compute dataset stats and exit", cxxopts::value<bool>())
        ("v,verbose", "Verbose output", cxxopts::value<bool>())
        ("h,help", "Print usage information");
    
    auto result = options.parse(argc, argv);

    if (result.count("help") || !result.count("input") || (!result.count("output") && !result.count("stats"))) {
        std::cout << options.help() << '\n';
        exit(0);
    }

    std::string input = result["input"].as<std::string>();
    std::string output;
    if (result.count("output")) output = result["output"].as<std::string>();

    int minsize = 0;
    int maxsize = very_large;
    int minnb = 0;
    int maxnb = very_large;
    double maxangle = 45.;
    double minstr = minus_infinity;
    double maxstr = plus_infinity;
    double minval = minus_infinity;
    double maxval = plus_infinity;

    if (result.count("size")) minsize = result["size"].as<int>();
    if (result.count("number")) maxnb = result["number"].as<int>();;
    if (result.count("angle")) maxangle = result["angle"].as<double>();
    vec3 t(0,0,0);
    if (result.count("translate")) {
        std::array<double, 3> _t = result["translate"].as<std::array<double, 3>>();
        t[0] = _t[0];
        t[1] = _t[1];
        t[2] = _t[2];
    }

    bool verbose=false;
    if (result.count("verbose")) verbose = true;

    bool doshrink = false;
    if (result.count("shrink")) doshrink = true;

    if (verbose) std::cout << "loading dataset... " << std::flush; 

    VTK_CREATE(vtkXMLPolyDataReader, reader);
    reader->SetFileName(input.c_str());
    reader->Update();
    VTK_SMART(vtkPolyData) data = reader->GetOutput();
    VTK_SMART(vtkPolyData) boundary_edges;
    VTK_SMART(vtkPolyData) feature_edges;
    VTK_SMART(vtkPolyData) cleaned;
    if (verbose) std::cout << "done.\n";

    VTK_CREATE(vtkPolyDataNormals, normals);
    normals->SetInputData(data);
    normals->ConsistencyOn();
    normals->ComputeCellNormalsOn();
    normals->ComputePointNormalsOff();
    normals->AutoOrientNormalsOn();
    normals->Update();
    data = normals->GetOutput();

    if (result.count("stats")) {
        int natts = data->GetPointData()->GetNumberOfArrays();
        std::cout << "There are " << natts << " data arrays in input\n";
        for (int n=0; n<natts; ++n) {
            VTK_SMART(vtkAbstractArray) arr = data->GetPointData()->GetAbstractArray(n);
            // std::vector<double> vals = array2vector(arr);
            std::vector<double> vals;
            array2vector2(vals, arr);
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            double _mode = spurt::mode(vals.begin(), vals.end(), 1024);
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            double _min = spurt::min(vals.begin(), vals.end());
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            double _max = spurt::max(vals.begin(), vals.end());
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            double _median = spurt::median(vals.begin(), vals.end());
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            std::pair<double, double> _meanvar = spurt::meanvariance(vals.begin(), vals.end());
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            std::vector<double> per = spurt::percentiles(vals.begin(), vals.end(), 21);
            std::cout << "vals.begin()=" << vals.begin() << ", vals.end()=" << vals.end() << '\n'; 
            std::cout << "Stats for " << arr->GetName() << ":\n";
            std::cout << "min: " << _min << "\nmax: " << _max << "\nmean: " << _meanvar.first << "\nvariance: " << _meanvar.second << "\nmedian: " << _median << "\nmode: " << _mode << '\n'; 
            std::cout << "percentiles:\n";
            std::copy(per.begin(), per.end(), std::ostream_iterator<double>(std::cout, ", "));
            std::cout << '\n';
        }
        return 0;
    }

    if (result.count("value")) {
        std::string valasstr = result["value"].as<std::string>();
        if (valasstr.back() == '%') {
            int pct = std::stoi(valasstr.substr(0, valasstr.size()-1));
            std::vector<double> vals = array2vector(data->GetPointData()->GetAbstractArray(result["namev"].as<std::string>().c_str()));
            std::sort(vals.begin(), vals.end());
            minval = vals[std::floor(pct*vals.size()/100.)];
        }
        else minval = std::stof(result["value"].as<std::string>());
    }
    if (result.count("strength")) {
        std::string strasstr = result["strength"].as<std::string>();
        if (strasstr.back() == '%') {
            int pct = 100-std::stoi(strasstr.substr(0, strasstr.size()-1));
            std::vector<double> vals = array2vector(data->GetPointData()->GetAbstractArray(result["names"].as<std::string>().c_str()));
            std::sort(vals.begin(), vals.end());
            maxstr = vals[std::floor(pct*vals.size()/100.)];
        }
        else maxstr = std::stof(result["strength"].as<std::string>());
    }

    std::cout << "thresholds are set to: value: " << minval << ", ridge strength: " << maxstr << '\n';

    if (verbose)
        std::cout << "Initially, there are " 
                  << data->GetNumberOfCells() << " cells\n";

    if (result["pre"].as<bool>()) {
        if (verbose) std::cout << "pre is TRUE\n";
        // first filter by value ...
        data = filter_by_value<double>(data, result["namev"].as<std::string>(), minval, maxval, false);

        if (verbose)
            std::cout << "After value filtering (>" << minval << "), there are " << data->GetNumberOfCells() << " cells\n";

        // ... then filter by ridge strength ...
        data = filter_by_value<double>(data, result["names"].as<std::string>(), minstr, maxstr, false);

        if (verbose) 
            std::cout << "Aftering strength filtering (<" << minstr << "), there are " << data->GetNumberOfCells() << " cells\n";

        // spurt::compute_cc_sizes(data);
        auto r = feature_aware_ccs(data, maxangle);
        cleaned = r[0];
        data = r[1];
        boundary_edges = r[2];
        feature_edges = r[3];

        if (verbose)
            std::cout << "There are initially " << data->GetCellData()->GetArray("CCIDs")->GetMaxNorm() << " connected components\n";

        // ... then filter by CC size ...
        data = filter_by_value<int>(data, "CCsizes", minsize, maxsize, true);

        if (verbose) 
            std::cout << "After filtering by CC size (>" << minsize << "), there are " << data->GetNumberOfCells() << " cells and " << data->GetCellData()->GetArray("CCIDs")->GetMaxNorm() << " connected components\n";

        // ... then filter by number of CCs
        data = filter_by_value<int>(data, "CCIDs", minnb, maxnb-1, true);

        if (verbose) 
            std::cout << "After filtering by number of CCs (<" << maxnb << "), there are " << data->GetNumberOfCells() << " cells\n";
    }
    else {
        if (verbose) std::cout << "pre is FALSE\n";
        if (verbose) std::cout << "computing connected components and their sizes\n";
        // spurt::compute_cc_sizes(data);
        auto r = feature_aware_ccs(data, maxangle);
        cleaned = r[0];
        data = r[1];
        boundary_edges = r[2];
        feature_edges = r[3];

        if (verbose) {
            std::cout << "There are initially " << data->GetCellData()->GetArray("CCIDs")->GetMaxNorm() << " connected components\n";
            std::cout << "filtering by size...\n";
            std::cout << "filtering size range is " << minsize << " to " << maxsize << '\n';
        }

        // First filter by CC size ...
        if (result.count("size")) data = filter_by_value<int>(data, "CCsizes", minsize, maxsize, true);

        if (verbose)
            std::cout << "After filtering by CC size (>" << minsize << "), there are " << data->GetNumberOfCells() << " cells and " << data->GetCellData()->GetArray("CCIDs")->GetMaxNorm() << " connected components\n";

        // ... then filter by number of CCs ...
        if (verbose) std::cout << "filtering by number of connected components, between " << minnb << " and " << maxnb << '\n';
        if (result.count("number")) data = filter_by_value<int>(data, "CCIDs", minnb, maxnb-1, true);

        if (verbose) 
            std::cout << "After filtering by number of CCs (<" << maxnb << "), there are " << data->GetNumberOfCells() << " cells\n";

        // ... then filter by value ...
        if (result.count("value")) data = filter_by_value<double>(data, result["namev"].as<std::string>(), minval, maxval, false);

        if (verbose)
            std::cout << "After value filtering (>" << minval << "), there are " << data->GetNumberOfCells() << " cells\n";

        // ... then filter by ridge strength.
        data = filter_by_value<double>(data, result["names"].as<std::string>(), minstr, maxstr, false);

        if (verbose)
            std::cout << "After strength filtering (<" << minstr << "), there are " << data->GetNumberOfCells() << " cells\n";
    }

    if (doshrink) {
        if (verbose) std::cout << "shrinking dataset...\n";
        data = shrink(data);
    }

    if (verbose) std::cout << "\ncomputing laciness...\n";
    // spurt::compute_laciness(data);
    VTK_CREATE(vtkXMLPolyDataWriter, writer);
    writer->SetCompressorTypeToLZ4();
    writer->SetCompressionLevel(5);

    auto basename = spurt::filename::remove_extension(output);
    std::cout << "exporting filtered crease mesh in " << basename << ".vtp... " << std::flush;
    writer->SetFileName((basename + ".vtp").c_str());
    writer->SetInputData(data);
    writer->Write();
    std::cout << "done.\n";
    std::cout << "exporting " << basename << "_boundaries.vtp... " << std::flush;
    writer->SetFileName((basename + "_boundaries.vtp").c_str());
    writer->SetInputData(boundary_edges);
    writer->Write();
    std::cout << "done.\n";
    std::cout << "exporting " << basename << "_features.vtp... " << std::flush;
    writer->SetFileName((basename + "_features.vtp").c_str());
    writer->SetInputData(feature_edges);
    writer->Write();
    std::cout << "done.\n";
    std::cout << "exporting " << basename << "_cleaned.vtp... " << std::flush;
    writer->SetFileName((basename + "_cleaned.vtp").c_str());
    writer->SetInputData(cleaned);
    writer->Write();
    std::cout << "done.\n";
    return 0;
}

