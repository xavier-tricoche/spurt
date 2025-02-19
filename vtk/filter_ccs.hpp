#pragma once

#include <vtk/vtk_utils.hpp>
#include <vtkConnectivityFilter.h>
#include <vtkCellArrayIterator.h>
#include <vtkCellDataToPointData.h>
#include <map>
#include <vector>
#include <algorithm>
#include <set>

namespace spurt {

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

    // std::cout << "Region ids is " << regionids << '\n';
    std::map<int, std::vector<int> > region_to_cells;
    for (int i=0; i<regionids->GetNumberOfTuples(); ++i) {
        int id = regionids->GetTuple1(i);
        auto it = region_to_cells.find(id);
        if (it == region_to_cells.end()) {
            region_to_cells[id] = std::vector<int>(1, id);
        }
        else {
            it->second.push_back(id);
        }
    }

    std::map<int, int> id2size;
    std::vector<int> sizes;
    std::map<int, std::vector<int>> size2ids;
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
    }
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
    for (int cellid=0; cellid<newregionids.size(); ++cellid) {
        int old_region_id = regionids->GetTuple1(cellid);
        newregionids[cellid] = old2new[old_region_id];
        regionsizes[cellid] = region_to_cells.at(old_region_id).size();
    }
    vtk_utils::add_scalars(inout, newregionids, false, "CCIDs", false);
    vtk_utils::add_scalars(inout, regionsizes, false, "CCsizes", false);
    VTK_CREATE(vtkCellDataToPointData, cell2point);
    cell2point->SetInputData(inout);
    cell2point->PassCellDataOn();
    cell2point->ProcessAllArraysOn();
    cell2point->Update();
    inout = cell2point->GetPolyDataOutput();
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
    }
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

} // namespace spurt