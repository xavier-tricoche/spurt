#pragma once

#include <vtk/vtk_utils.hpp>
#include <vtkConnectivityFilter.h>
#include <vtkCellArrayIterator.h>
#include <vtkCellDataToPointData.h>
#include <vtkExtractEdges.h>
#include <vtkFeatureEdges.h>
#include <vtkLogger.h>
#include <vtkDoubleArray.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <map>
#include <vector>
#include <algorithm>
#include <set>
#include <array>
#include <iterator>

#include <misc/progress.hpp>

namespace spurt {

template<typename T = double>
inline std::array<T, 3> _cross_product(const T* a, const T* b, const T* c) {
    T ab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
    T ac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
    return {ab[1]*ac[2] - ab[2]*ac[1], ab[2]*ac[0] - ab[0]*ac[2], ab[0]*ac[1] - ab[1]*ac[0]};
}

template<typename T = double>
inline std::array<T, 3> _normal(const T* a, const T* b, const T* c) {
    auto n = _cross_product(a, b, c);
    T len = std::sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    if (len == 0.) return {0., 0., 0.};
    else return {n[0]/len, n[1]/len, n[2]/len};
}

template<typename T = double>
inline std::array<T, 3> _normal_from_ids(vtkIdType p1, vtkIdType p2, vtkIdType p3, VTK_SMART(vtkPolyData) polydata) {
    auto coords = polydata->GetPoints();
    auto a = coords->GetPoint(p1);
    auto b = coords->GetPoint(p2);
    auto c = coords->GetPoint(p3);
    return _normal(a, b, c);
}

template<typename T = double>
inline T _dot_product(const std::array<T, 3>& a, const std::array<T, 3>& b) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

inline bool _check_orientation(const std::array<vtkIdType, 2>& edge, vtkIdList* pts) {
    auto n = pts->GetNumberOfIds();
    if (n < 2) return false;
    for (vtkIdType eid=0; eid<n; eid++) {
        vtkIdType p0 = pts->GetId(eid);
        vtkIdType p1 = pts->GetId((eid + 1) % n);
        if (p0 == edge[0] && p1 == edge[1]) {
            // oriented correctly
            return true;
        }
        if (p0 == edge[1] && p1 == edge[0]) {
            // opposite orientation
            return false;
        }
    }
    return false;
}

inline std::array<double, 3> _get_edge_neighbor_normal(vtkIdType cellid, vtkIdType p1, vtkIdType p2, VTK_SMART(vtkPolyData) polydata) {
    auto polys = polydata->GetPolys();
    auto coords = polydata->GetPoints();
    vtkIdList* pts;
    polys->GetCellAtId(cellid, pts);
    auto q0 = coords->GetPoint(pts->GetId(0));
    auto q1 = coords->GetPoint(pts->GetId(1));
    auto q2 = coords->GetPoint(pts->GetId(2));

    if (_check_orientation(std::array<vtkIdType, 2>{p1, p2}, pts)) 
       return _normal(q0, q1, q2);
    else return _normal(q1, q0, q2);
}

inline void compute_cc_sizes(VTK_SMART(vtkPolyData) inout, bool only_faces=true)
{
    if (only_faces) {
        inout->SetVerts(nullptr);
        inout->SetLines(nullptr);
    }

    std::cout << "Running connectivity filter... " << std::flush;
    VTK_CREATE(vtkConnectivityFilter, connec);
    connec->SetExtractionModeToAllRegions();
    connec->ColorRegionsOn();
    connec->SetInputData(inout);
    connec->Update();
    std::cout << "done.\n";

    VTK_SMART(vtkPolyData) pd = dynamic_cast<vtkPolyData*>(connec->GetOutput());
    // std::cout << "output polydata is \n"; 
    // pd->PrintSelf(std::cout, vtkIndent(1)); 
    std::cout << "trying to access RegionId attribute... " << std::flush;
    auto regionids = pd->GetCellData()->GetArray("RegionId");
    std::cout << "done\n";
    inout->GetCellData()->AddArray(regionids);

    // std::cout << "Region ids is " << regionids << '\n';
    std::map<int, std::vector<int> > region_to_cells;
    
    spurt::ProgressDisplay progress;
    progress.start(regionids->GetNumberOfTuples(), "Computing CC sizes");
    for (int i=0; i<regionids->GetNumberOfTuples(); ++i) {
        int id = regionids->GetTuple1(i);
        auto it = region_to_cells.find(id);
        if (it == region_to_cells.end()) {
            region_to_cells[id] = std::vector<int>(1, id);
        }
        else {
            it->second.push_back(id);
        }
        progress.update(i);
    }
    progress.end();
    std::cout << "There are " << region_to_cells.size() << " regions\n";

    std::map<int, int> id2size;
    std::vector<int> sizes;
    std::map<int, std::vector<int>> size2ids;
    progress.start(region_to_cells.size(), "Mapping cc's to cells");
    int count=0;
    for (auto iter = region_to_cells.begin(); iter!=region_to_cells.end(); ++iter) 
    {
        int rid = iter->first;
        int ncells = iter->second.size();
        id2size[rid] = ncells;
        sizes.push_back(ncells);
        auto it = size2ids.find(ncells);
        if (it == size2ids.end()) {
            size2ids[ncells] = std::vector<int>(1, rid);
        }
        else it->second.push_back(rid);
        progress.update(++count);
    }
    progress.end();
    std::sort(sizes.begin(), sizes.end(), [&](int a, int b) { return a > b; });
    std::cout << "10 largest regions are:\n";
    for (int i=0; i<10; ++i) {
        std::cout << sizes[i] << '\n';
    }

    std::cout << "10 largest regions with ids are:\n";
    int nshown=0;
    for (auto iter=size2ids.rbegin(); iter!=size2ids.rend() && nshown<10; ++iter, ++nshown) {
        std::cout << iter->first << ": ";
        std::copy(iter->second.begin(), iter->second.end(),
                  std::ostream_iterator<int>(std::cout, ", "));
    }

    // compute argsort of region ids by decreasing size
    std::vector<int> sorted_ids(region_to_cells.size());
    std::iota(sorted_ids.begin(), sorted_ids.end(), 0);
    std::sort(sorted_ids.begin(), sorted_ids.end(), [&]
        (int a, int b) {
            return (region_to_cells.at(a).size() > region_to_cells.at(b).size());
        }
    );

    std::map<int, int> old2new;
    // renumber regions by decreasing order of size
    for (int i=0; i<sorted_ids.size(); ++i) {
        old2new[sorted_ids[i]] = i;
    }
    std::vector<int> newregionids(regionids->GetNumberOfTuples());
    std::vector<int> regionsizes(regionids->GetNumberOfTuples());
    progress.start(newregionids.size(), "Renumbering regions");
    for (int cellid=0; cellid<newregionids.size(); ++cellid) {
        int old_region_id = regionids->GetTuple1(cellid);
        newregionids[cellid] = old2new[old_region_id];
        regionsizes[cellid] = region_to_cells.at(old_region_id).size();
        progress.update(cellid);
    }
    progress.end();
    vtk_utils::add_scalars(inout, newregionids, false, "CCIDs", false);
    vtk_utils::add_scalars(inout, regionsizes, false, "CCsizes", false);
    VTK_CREATE(vtkCellDataToPointData, cell2point);
    cell2point->SetInputData(inout);
    cell2point->PassCellDataOn();
    cell2point->ProcessAllArraysOn();
    cell2point->Update();
    inout = cell2point->GetPolyDataOutput();
    for (int i=0; i<inout->GetPointData()->GetNumberOfArrays(); ++i) {
        std::cout << "Point data array " << i << " is " << inout->GetPointData()->GetArray(i)->GetName() << '\n';
    }
    for (int i=0; i<inout->GetCellData()->GetNumberOfArrays(); ++i) {
        std::cout << "Cell data array " << i << " is " << inout->GetCellData()->GetArray(i)->GetName() << '\n';
    }
}

inline void prune_mesh(VTK_SMART(vtkPolyData) inout, int minsize=-1, 
                       int maxnbclusters=-1, bool only_faces = true)
{
    compute_cc_sizes(inout, only_faces);

    // Note: cell ids run across verts, lines, strips, and polys,
    // presumably in that order (though unclear about polys vs. strips)
    int nverts = inout->GetNumberOfVerts();
    int nlines = inout->GetNumberOfLines();
    int npolys = inout->GetNumberOfPolys();

    auto ids = inout->GetCellData()->GetArray("CCIDs");
    auto sizes = inout->GetCellData()->GetArray("CCsizes");
    std::map<int, int> id2size;
    std::map<int, int> id2size_alt;
    spurt::ProgressDisplay progress;
    progress.start(ids->GetNumberOfTuples(), "Pruning mesh");
    for (int i=0; i<ids->GetNumberOfTuples(); ++i) {
        int id = ids->GetTuple1(i);
        auto it = id2size.find(id);
        if (it == id2size.end()) {
            id2size[id] = 1;
        }
        else ++it->second;
        int val = sizes->GetTuple1(i);
        auto it_alt = id2size_alt.find(id);
        if (it_alt == id2size_alt.end()) {
            id2size_alt[id] = val;
        }
        else
        {
            assert(it_alt->second == val);
        }
        progress.update(i);
    }
    progress.end();
    std::cout << "first few regions in output\n";
    int count=0;
    for (auto it=id2size.begin(); it!=id2size.end() && count<10; ++it, ++count) {
        std::cout << "Region " << it->first << " has size " << it->second;
        std::cout << "Precomputed value for region " << it->first << " was " << id2size_alt.at(it->first) << '\n';
    }

    VTK_CREATE(vtkCellArray, newtris);
    VTK_CREATE(vtkCellArray, newlines);
    VTK_CREATE(vtkCellArray, newverts); // only needed if no filtering
    auto oldtris = inout->GetPolys();
    auto oldlines = inout->GetLines();
    auto oldverts = inout->GetVerts();
    VTK_CREATE(vtkIdList, alist);
    VTK_SMART(vtkCellArrayIterator) cellit;
    VTK_SMART(vtkCellArray) newcells;
    std::vector<int> newids;
    std::vector<int> newsizes;
    
    auto regionsizes = inout->GetCellData()->GetArray("CCsizes");
    auto regionids = inout->GetCellData()->GetArray("CCIDs");
    for (int c=0; c<3; ++c) {
        int offset = 0;
        if (c==0) {
            if (minsize > 1 || maxnbclusters > 0 || only_faces) continue;
            cellit = oldverts->NewIterator();
            newcells = newverts;
        }
        else if (c==1) {
            if (only_faces) continue;
            cellit = oldlines->NewIterator();
            offset = nverts;
            newcells = newlines;
        }
        else {
            cellit = oldtris->NewIterator();
            offset = nverts + nlines;
            newcells = newtris;
        }

        int cellid = offset;
        for (; !cellit->IsDoneWithTraversal(); cellit->GoToNextCell(), ++cellid)
        {
            if (maxnbclusters > 0 && 
                regionids->GetTuple1(cellid) >= maxnbclusters) 
                continue;
            if (minsize > 0 && 
                regionsizes->GetTuple1(cellid) < minsize) 
                continue;
            
            cellit->GetCurrentCell(alist);
            newcells->InsertNextCell(alist);
            newids.push_back(int(regionids->GetTuple1(cellid)));
            newsizes.push_back(int(regionsizes->GetTuple1(cellid)));
        }
    }
    inout->GetCellData()->RemoveArray("CCIDs");
    inout->GetCellData()->RemoveArray("CCsizes");
    if (!only_faces) {
        inout->SetVerts(newverts);
        inout->SetLines(newlines);
    }
    inout->SetPolys(newtris);
    vtk_utils::add_scalars(inout, newids, false, "CCIDs", false);
    vtk_utils::add_scalars(inout, newsizes, false, "CCsizes", false);
    std::cout << "On output there are " 
        << inout->GetNumberOfVerts() << " verts, "
        << inout->GetNumberOfLines() << " lines, and " 
        << inout->GetNumberOfPolys() << " triangles\n";
}

// struct EdgeInfo {
//     vtkIdType p1, p2; // the id of the vertices 
//     vtkIdType cellid; // the id of the cell that contains this edge
//     double angle; // the angle between the normals of the two cells sharing this edge
//     vtkIdType best_neigh; // the id of the best neighbor cell
//     int atype; // 0: boundary, 1: nonmanifold, 2: sharp, 3: line 

//     EdgeInfo(vtkIdType _p1, vtkIdType _p2, vtkIdType _cellid, double _angle, int _atype, vtkIdType _best_neigh=vtkIdType(-1))
//         : p1(_p1), p2(_p2), cellid(_cellid), angle(_angle), best_neigh(_best_neigh), atype(_atype), best_neigh(_best_neigh) {}
// };

/*
   Extract connected components of input mesh while treating sharp and non-
   manifold edges as boundary edges. Individual connected components are 
   reindexed such that no vertex index is shared across CC's. In other words, 
   sharps and non-manifold edges are duplicated to create cuts. This function 
   returns a new polydata object since it modifies both positions and cells.
*/
inline VTK_SMART(vtkPolyData) compute_smooth_ccs(VTK_SMART(vtkPolyData) input, double angle_threshold=45.)
{
    typedef double                   val_t;
    typedef vtkDoubleArray           valarray_t;
    typedef vtkCellArray             cellarray_t;
    typedef vtkIntArray              intarray_t;
    typedef std::array<val_t, 3>     vec3_t;
    typedef vtkIdType                idx_t;
    typedef std::array<idx_t, 2>     edge_t;
    typedef vtkIdList                idxlist_t;

    constexpr val_t _deg_to_rad = 3.14159265358979323846 / 180.;
    val_t threshold = std::cos(angle_threshold * _deg_to_rad);

    auto polygons = input->GetPolys();
    std::cout << "polygons contain " << polygons->GetNumberOfCells() << " cells\n";
    input->BuildLinks();
    auto coordinates = input->GetPoints();
    VTK_CREATE(vtkIdList, neighbor_cell_ids);
    VTK_CREATE(vtkIdList, cell_point_ids);
    neighbor_cell_ids->Allocate(20);
    cell_point_ids->Allocate(20);

    // first, compute cell normals
    std::vector<vec3_t> normals(polygons->GetNumberOfCells(), {0., 0., 0.});
    for (idx_t i=0; i<polygons->GetNumberOfCells(); i++)
    {
        polygons->GetCellAtId(i, cell_point_ids);
        if (cell_point_ids->GetNumberOfIds() < 3) continue;
        normals[i] = _normal_from_ids(cell_point_ids->GetId(0), 
                                      cell_point_ids->GetId(1), 
                                      cell_point_ids->GetId(2), input);
    }

    // compute connected components while watching for sharp angles
    std::vector<idx_t> poly_to_cc(polygons->GetNumberOfCells(), -1);
    VTK_CREATE(vtkPolyData, new_dataset);
    VTK_CREATE(vtkPoints, new_coordinates);
    VTK_CREATE(cellarray_t, new_polygons);
    VTK_CREATE(intarray_t, region_ids); // CC ID
    VTK_CREATE(intarray_t, sizes); // CC size
    VTK_CREATE(valarray_t, laciness); // CC laciness
    region_ids->SetName("CCIDs");
    region_ids->SetNumberOfComponents(1);
    sizes->SetName("CCsizes");
    sizes->SetNumberOfComponents(1);
    laciness->SetName("CC Laciness");
    laciness->SetNumberOfComponents(1);
    
    std::map<idx_t, idx_t> vertex_index_remap; // vertex reindexing per CC
    idx_t cc_id = 0;
    for (idx_t i=0; i<polygons->GetNumberOfCells(); i++)
    {
        if (poly_to_cc[i] != -1) continue; // cell has been processed already

        // We are starting a new CC...
        size_t nedges = 0;
        size_t nvertices = 0;
        size_t nboundary_edges = 0;
        std::vector<idx_t> new_cell_ids;
        std::list<idx_t> queue;
        poly_to_cc[i] = cc_id;
        queue.push_back(i);
        while (!queue.empty()) {
            idx_t cell_id = queue.front();
            queue.pop_front();
            polygons->GetCellAtId(cell_id, cell_point_ids);
            auto cell_normal = normals[cell_id];
            VTK_CREATE(idxlist_t, new_cell_point_ids);
            new_cell_point_ids->Allocate(20);
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
                    new_coordinates->InsertNextPoint(p[0], p[1], p[2]);
                    vertex_index_remap[a] = new_id;
                    new_cell_point_ids->InsertNextId(new_id);
                    ++nvertices;
                }
                else {
                    new_cell_point_ids->InsertNextId(it->second);
                }
                edge_t edge{a, b};
                ++nedges;
                input->GetCellEdgeNeighbors(cell_id, a, b, neighbor_cell_ids);
                if (neighbor_cell_ids->GetNumberOfIds() == 0) {
                    // boundary edge
                    ++nboundary_edges;
                    continue;
                }
                else if (neighbor_cell_ids->GetNumberOfIds() == 1) {
                    // regular edge: check for sharpness
                    idx_t neighbor_cell_id = neighbor_cell_ids->GetId(0);
                    if (poly_to_cc[neighbor_cell_id] != -1) continue; // cell has been processed already
                    VTK_CREATE(vtkIdList, neighbor_cell_point_ids);
                    neighbor_cell_point_ids->Allocate(20);
                    polygons->GetCellAtId(neighbor_cell_id, neighbor_cell_point_ids);
                    vec3_t neighbor_normal = normals[neighbor_cell_id];
                    if (!_check_orientation(edge, neighbor_cell_point_ids)) {
                        neighbor_normal = {-neighbor_normal[0], -neighbor_normal[1], -neighbor_normal[2]};
                    }
                    val_t dot = _dot_product(cell_normal, neighbor_normal);
                    if (dot > threshold) {
                        // add to cc
                        poly_to_cc[neighbor_cell_id] = cc_id;
                        queue.push_back(neighbor_cell_id);
                    }
                    else ++nboundary_edges;
                }
                else if (neighbor_cell_ids->GetNumberOfIds() > 1) {
                    // non-manifold edge: find the  eighbor with the largest dot product
                    val_t maxdot = -1;
                    idx_t best_neigh = -1;
                    for (idx_t k=0; k<neighbor_cell_ids->GetNumberOfIds(); k++) {
                        idx_t neighbor_cell_id = neighbor_cell_ids->GetId(k);
                        if (poly_to_cc[neighbor_cell_id] != -1) continue; // cell has been processed already
                        // TODO: process neighbor
                        VTK_CREATE(vtkIdList, neighbor_cell_point_ids);
                        neighbor_cell_point_ids->Allocate(20);
                        polygons->GetCellAtId(neighbor_cell_id, neighbor_cell_point_ids);
                        vec3_t neighbor_normal = normals[neighbor_cell_id];
                        if (!_check_orientation(edge, neighbor_cell_point_ids)) {
                            neighbor_normal = {-neighbor_normal[0], -neighbor_normal[1], -neighbor_normal[2]};
                        }
                        val_t dot = _dot_product(cell_normal, neighbor_normal);
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
                    else ++nboundary_edges;
                }
            }
            // cell has been processed. add it to our new mesh
            auto new_id = new_polygons->InsertNextCell(new_cell_point_ids);
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
    new_dataset->SetPoints(new_coordinates);
    new_dataset->SetPolys(new_polygons);
    new_dataset->GetCellData()->AddArray(region_ids);
    new_dataset->GetCellData()->AddArray(laciness);
    new_dataset->GetCellData()->AddArray(sizes);
    return new_dataset;
}

// inline void compute_feature_edges(std::vector<EdgeInfo>& fedges, 
//     VTK_SMART(vtkPolyData) input, double angle_threshold=45.0)
// {
//     double threshold = std::cos(angle_threshold * M_PI / 180.0);
//     typedef std::array<vtkIdType, 2> edge_t;
//     vtkIdType nneigh, npts, cellid;
//     const vtkIdType* pts;
//     vtkIdList* neighbors;
//     vtkIdType* neigh_pts;
//     input->BuildLinks();
//     auto polys = input->GetPolys();
//     auto coords = input->GetPoints();
//     for (cellid = 0, polys->InitTraversal(); polys->GetNextCell(npts, pts);
//        cellid++) {
//         if (npts < 2) {
//             continue;
//         }
//         else if (npts == 2) {
//             fedges.push_back(EdgeInfo{pts[0], pts[1], cellid, 0.0, 3});
//             continue;
//         }
//         else if (npts >= 3) {
//             auto q0 = coords->GetPoint(pts[0]);
//             auto q1 = coords->GetPoint(pts[1]);
//             auto q2 = coords->GetPoint(pts[2]);
//             auto n = _normal(q1, q0, q2);
//             for (vtkIdType i = 0; i < npts; i++) {
//                 const vtkIdType& p1 = pts[i];
//                 const vtkIdType& p2 = pts[(i + 1) % npts];
//                 input->GetCellEdgeNeighbors(cellid, p1, p2, neighbors);
//                 nneigh = neighbors->GetNumberOfIds();
//                 if (nneigh == 0) {
//                     fedges.push_back(EdgeInfo{p1, p2, cellid, 0.0, 0});
//                 }
//                 else if (nneigh > 1) {
//                     // only record if current cell is min of all cells sharing that edge
//                     neighbors->Sort();
//                     if (cellid < neighbors->GetId(0)) {
//                         std::vector<double> cosines;
//                         for (vtkIdType j = 0; j < nneigh; j++) {
//                             auto neigh_id = neighbors->GetId(j);
//                             auto neigh_n = _get_edge_neighbor_normal(neigh_id, p1, p2, input);
//                             cosines.push_back(_dot(n, neigh_n));
//                         }
//                         // best neighbor is one that maximizes cosine 
//                         auto max_it = std::max_element(cosines.begin(), cosines.end());
//                         auto max_idx = std::distance(cosines.begin(), max_it);
//                         auto best_neigh_id = neighbors->GetId(max_idx);
//                         fedges.push_back(EdgeInfo{p1, p2, cellid, 0.0, 1, best_neigh_id});
//                     }
//                 }
//                 else if (cellid < neighbors->GetId(0)) {
//                     vtkIdType n_npts;
//                     polys->GetCellAtId(neighbors->GetId(0), n_npts, neigh_pts);
//                     std::vector<vtkIdType> neigh_pts_vec(std::begin(neigh_pts), std::begin(neigh_pts) + n_npts);
//                     _reorient(edge_t{p1, p2}, neigh_pts_vec);
//                     auto nq0 = coords->GetPoint(neigh_pts_vec[0]);
//                     auto nq1 = coords->GetPoint(neigh_pts_vec[1]);
//                     auto nq2 = coords->GetPoint(neigh_pts_vec[2]);
//                     auto nn = _normal(nq1, nq0, nq2);
//                     double dot = _dot_product(n, nn);
//                     if (dot < threshold) {
//                         fedges.push_back(EdgeInfo{p1, p2, cellid, std::acos(dot), 2});
//                     }
//                 }
//             }
//         }
//     }
// }

// inline void cut_mesh_by_edges(VTK_SMART(vtkPolyData) inout, const std::vector<EdgeInfo>& fedges)
// {
//     std::vector<std::array<double, 3>> newpoints;
//     VTK_CREATE(vtkCellArray, newpolys);
//     VTK_CREATE(vtkDoubleArray, newcoords);
//     std::map<vtkIdType, vtkIdType> old2new;

//     auto polys = inout->GetPolys();
//     polys->BuildLinks();
//     for ( auto& edge : fedges ) {
//         if (edge.atype == 0) continue;
//         else if (edge.atype == 1) {
//             // nonmanifold edge
//             vtkIdList* neighs;
//             vtkIdType* pts;
//             vtkIdType npts;
//             polys->GetCellEdgeNeighbors(edge.cellid, edge.p1, edge.p2, neighs);
//             // find the cells that share this edge
//             // TODO: handle type 1 edges
//         }
//         else {
//             // TODO: handle type 2 edges
//         }
//     }
// }

// inline void compute_laciness(VTK_SMART(vtkPolyData) inout)
// {
//     spurt::ProgressDisplay progress;
//     auto regionids = inout->GetCellData()->GetArray("CCIDs");
//     // group cells by region id
//     std::map<int, std::list<int>> id2cells;
//     progress.start(inout->GetNumberOfCells(), "Computing region to cell mapping");
//     for (int i=0; i<inout->GetNumberOfCells(); ++i) {
//         auto id = regionids->GetTuple1(i);
//         auto it = id2cells.find(int(id));
//         if (it == id2cells.end()) {
//             id2cells[int(id)] = std::list<int>();
//         }
//         id2cells[int(id)].push_back(i);
//         progress.update(i);
//     }
//     progress.end();

//     VTK_CREATE(vtkFloatArray, laciness);
//     laciness->SetNumberOfComponents(1);
//     laciness->SetNumberOfTuples(inout->GetNumberOfCells());
//     laciness->SetName("Laciness");
//     VTK_CREATE(vtkExtractEdges, all_edges);
//     all_edges->SetInputData(inout);
//     all_edges->Update();
//     auto _edges = all_edges->GetOutput();
//     std::map<int, std::list<int>> ccid_to_edges;
//     progress.start(_edges->GetNumberOfCells(), "Computing edge to region mapping");
//     for (int i=0; i<_edges->GetNumberOfCells(); ++i) {
//         auto id = _edges->GetCellData()->GetArray("CCIDs")->GetTuple1(i);
//         auto it = ccid_to_edges.find(int(id));
//         if (it == ccid_to_edges.end()) {
//             ccid_to_edges[int(id)] = std::list<int>();
//         }
//         ccid_to_edges[int(id)].push_back(i);
//         progress.update(i);
//     }
//     progress.end();
//     VTK_CREATE(vtkFeatureEdges, feature_edges);
//     feature_edges->SetInputData(inout);
//     feature_edges->BoundaryEdgesOn();
//     feature_edges->FeatureEdgesOff();
//     feature_edges->ManifoldEdgesOff();
//     feature_edges->NonManifoldEdgesOn();
//     feature_edges->Update();
//     auto _feature_edges = feature_edges->GetOutput();
//     std::map<int, std::list<int>> ccid_to_feature_edges;
//     progress.start(_feature_edges->GetNumberOfCells(), "Computing feature edge to region mapping");
//     for (int i=0; i<_feature_edges->GetNumberOfCells(); ++i) {
//         auto id = _feature_edges->GetCellData()->GetArray("CCIDs")->GetTuple1(i);
//         auto it = ccid_to_feature_edges.find(int(id));
//         if (it == ccid_to_feature_edges.end()) {
//             ccid_to_feature_edges[int(id)] = std::list<int>();
//         }
//         ccid_to_feature_edges[int(id)].push_back(i);
//         progress.update(i);
//     }
//     progress.end();
//     std::map<int, float> ccid_to_laciness;
//     progress.start(ccid_to_edges.size(), "Computing laciness");
//     int count = 0;
//     for (const auto& [id, edges] : ccid_to_edges) {
//         auto feature_edges = ccid_to_feature_edges.find(id);
//         if (feature_edges == ccid_to_feature_edges.end()) {
//             ccid_to_laciness[id] = 0.0f;
//         } else {
//             ccid_to_laciness[id] = float(feature_edges->second.size()) / float(edges.size());
//         }
//         progress.update(count++);
//     }
//     progress.end();
//     for (int i=0; i<inout->GetNumberOfCells(); ++i) {
//         auto id = inout->GetCellData()->GetArray("CCIDs")->GetTuple1(i);
//         laciness->SetTuple1(i, ccid_to_laciness[int(id)]);
//     }

//     inout->GetCellData()->AddArray(laciness);
// }

} // namespace spurt