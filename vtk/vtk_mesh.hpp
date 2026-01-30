#pragma once

#include <map>
#include <vector>
#include <algorithm>
#include <set>
#include <array>
#include <iterator>

#include <vtkCellArray.h>
#include <vtkCellArrayIterator.h>
#include <vtkCellDataToPointData.h>
#include <vtkConnectivityFilter.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>
#include <vtkExtractEdges.h>
#include <vtkFeatureEdges.h>
#include <vtkIdTypeArray.h>
#include <vtkIntArray.h>
#include <vtkLogger.h>
#include <vtkPolyData.h>

#include <vtk/vtk_utils.hpp>
#include <misc/progress.hpp>
#include <math/small_vector.hpp>

namespace {

template<typename T, size_t N>
struct epsilon_lexico {
    typedef T value_t;
    typedef spurt::small_vector<value_t, N> vec_t;
    const size_t size=N;
    value_t _M_eps;
    
    epsilon_lexico(value_t eps=1.0e-9) : _M_eps(eps) {}
    
    bool operator()(const vec_t& a, const vec_t& b) const {
        for (size_t i = 0; i < size; ++i) {
            if (a[i] + _M_eps < b[i]) return true;
            else if (b[i] + _M_eps < a[i]) return false;
        }
        return false;
    }
};

} // anonymous namespace

namespace vtk_utils {

void unique_vertices(VTK_SMART(vtkPoints) newpoints, std::vector<vtkIdType>& point_remap, 
    std::map<vtkIdType, vtkIdType>& new_to_old_points, 
    VTK_SMART(vtkPoints) oldpoints, double eps=1.0e-9, vtkIdType binsize=20) {
    typedef vtkIdType idx_t;
    typedef double value_t;
    typedef spurt::small_vector<value_t, 3> point_t;
    typedef epsilon_lexico<value_t, 3> epsilon_lexico_t;
    typedef std::map<point_t, idx_t, epsilon_lexico_t> point_map_t;

    idx_t npoints = oldpoints->GetNumberOfPoints();
    std::vector<idx_t> initial_remap(npoints, -1);
    oldpoints->ComputeBounds();
    double bnds[6];
    oldpoints->GetBounds(bnds);
    double dx = bnds[1]-bnds[0];
    double dy = bnds[3]-bnds[2];
    double dz = bnds[5]-bnds[4];
    double volume = dx*dy*dz; 
    double mu = pow((float)npoints/(volume*binsize), 1./3.);
    double maxdim = std::max({dx, dy, dz});
    spurt::small_vector<long, 3> resolution;
    resolution[0] = std::max<long>((long)(mu*dx), 1);
    resolution[1] = std::max<long>((long)(mu*dy), 1);
    resolution[2] = std::max<long>((long)(mu*dz), 1);
    std::vector< std::list<idx_t> > grid(resolution[0]*resolution[1]*resolution[2]);
    std::vector<std::mutex> grid_mtx(resolution[0]*resolution[1]*resolution[2]);
    idx_t neighbors[3][3][3];
    for (int i=-1; i<2; ++i) { 
        for (int j=-1; j<2; ++j) {
            for (int k=-1; k<2; ++k) {
                neighbors[i+1][j+1][k+1] = i + resolution[0]*(j + resolution[1]*k);
            }
        }
    }
    spurt::ProgressDisplay progress;
    std::mutex update_mutex;
    progress.begin(npoints, "Assign vertices to bins");
    std::atomic<idx_t> counter = 0;
    tbb::parallel_for(tbb::blocked_range<idx_t>(0, npoints),
        [&](tbb::blocked_range<idx_t> r)
        {
            point_t pt;
            for (idx_t n = r.begin(); n < r.end(); ++n)
            {
                ++counter;
                std::unique_lock<std::mutex> lock(update_mutex, std::defer_lock);
                if (lock.try_lock())
                {
                    progress.update(counter);
                }
                oldpoints->GetPoint(n, const_cast<double*>(pt.data()));
                double u = (pt[0]-bnds[0])*resolution[0]/dx;
                double v = (pt[1]-bnds[2])*resolution[1]/dy;
                double w = (pt[2]-bnds[4])*resolution[2]/dz;
                idx_t i = floor(u);
                idx_t j = floor(v);
                idx_t k = floor(w);
                idx_t idx = i + resolution[0]*(j + resolution[1]*k);
                double du = u-(double)i;
                double dv = v-(double)j;
                double dw = w-(double)k;
                bool valid[3][3][3] = {true};
                if (i==0 || du>eps) {
                    for (int j=-1; j<2; ++j) {
                        for (int k=-1; k<2; ++k) {
                            valid[0][j+1][k+1] = false;
                        }
                    }
                }
                else if (i==resolution[0]-1 || du<1-eps) {
                    for (int j=-1; j<2; ++j) {
                        for (int k=-1; k<2; ++k) {
                            valid[2][j+1][k+1] = false;
                        }
                    }
                }
                if (j==0 || dv<1-eps) {
                    for (int i=-1; i<2; ++i) {
                        for (int k=-1; k<2; ++k) {
                            valid[i+1][0][k+1] = false;
                        }
                    }
                }
                else if (j==resolution[1]-1 || dv<1-eps) {
                    for (int i=-1; i<2; ++i) {
                        for (int k=-1; k<2; ++k) {
                            valid[i+1][2][k+1] = false;
                        }
                    }
                }
                if (k==0 || dw>eps) {
                    for (int i=-1; i<2; ++i) {
                        for (int j=-1; j<2; ++j) {
                            valid[i+1][j+1][0] = false;
                        }
                    }
                }
                else if (k==resolution[2]-1 || dw<1-eps) {
                    for (int i=-1; i<2; ++i) {
                        for (int j=-1; j<2; ++j) {
                            valid[i+1][j+1][2] = false;
                        }
                    }
                }
                for (int i=-1; i<2; ++i) {
                    for (int j=-1; j<2; ++j) {
                        for (int k=-1; k<2; ++k) {
                            if (valid[i+1][j+1][k+1]) {
                                idx_t idx2 = idx+neighbors[i+1][j+1][k+1];
                                std::scoped_lock lock(grid_mtx[idx]);
                                grid[idx].push_back(n);
                            }
                        }
                    }
                }
            }
        });
    progress.end();

    progress.begin(npoints, "Consolidate vertices per bin");
    counter = 0;
    std::mutex remap_mtx;
    tbb::parallel_for(tbb::blocked_range<idx_t>(0, grid.size()),
        [&](tbb::blocked_range<idx_t> r)
        {
            point_t pt;
            for (idx_t n = r.begin(); n < r.end(); ++n)
            {
                ++counter;
                if (grid[n].empty()) continue;
                std::unique_lock<std::mutex> lock(grid_mtx[n], std::defer_lock);
                if (lock.try_lock()) {
                    progress.update(counter);
                }
                epsilon_lexico_t eps_order(eps);
                point_map_t _map(eps_order);
                for (auto it=grid[n].begin(); it!=grid[n].end(); ++it) {
                    idx_t idx = *it;
                    oldpoints->GetPoint(idx, const_cast<double*>(pt.data()));
                    auto _it = _map.find(pt);
                    if (_it != _map.end()) {
                        std::scoped_lock lock(remap_mtx);
                        initial_remap[idx] = _it->second;
                    }
                    else {
                        std::scoped_lock lock(remap_mtx);
                        _map[pt] = idx;
                        initial_remap[idx] = idx;
                    }
                }
            }
        });
    progress.end();

    progress.begin(grid.size(), "Reindexing");
    point_remap.resize(npoints);
    for (idx_t n=0; n<npoints; ++n) {
        progress.update(n);
        if (initial_remap[n] == n) {
            idx_t newid = newpoints->InsertNextPoint(oldpoints->GetPoint(n));
            point_remap[n] = newid;
            new_to_old_points[newid] = n;
        }
        else {
            point_remap[n] = point_remap[initial_remap[n]];
        }
    }
    progress.end();
}

VTK_SMART(vtkPolyData) clean_polydata_mesh(VTK_SMART(vtkPolyData) input, double eps=1.0e-9) {
    typedef double value_t;
    typedef spurt::small_vector<value_t, 3> point_t;
    typedef epsilon_lexico<value_t, 3> epsilon_lexico_t;
    typedef vtkIdType idx_t;

    VTK_SMART(vtkPoints) old_points = input->GetPoints();
    VTK_SMART(vtkCellArray) old_cells = input->GetPolys();
    
    // a global spatial sorting structure for unique points up to epsilon L1-distance
    epsilon_lexico_t eps_order(eps);
    std::map<point_t, idx_t, epsilon_lexico_t> point_map(eps_order);
    // the new, consolidated indices of the new points
    std::vector<idx_t> point_remap(old_points->GetNumberOfPoints(), -1);
    // mapping from new point index to old point index
    std::map<idx_t, idx_t> new_to_old_points;
    // mapping from new cell index to old cell index
    std::map<idx_t, idx_t> new_to_old_cells;
    point_t pt;
    size_t pt_id = 0;

    // Compute unique points
    spurt::ProgressDisplay progress;
    VTK_CREATE(vtkPoints, new_points);
    progress.begin(old_points->GetNumberOfPoints(), "Find unique vertices");
    for (int i=0; i<old_points->GetNumberOfPoints(); ++i) {
        progress.update(i);
        old_points->GetPoint(i, const_cast<double*>(pt.data()));
        auto it = point_map.find(pt);
        if (it != point_map.end()) {
            point_remap[i] = it->second;
        } else {
            point_remap[i] = pt_id;
            point_map[pt] = pt_id;
            new_points->InsertNextPoint(pt.data());
            new_to_old_points[pt_id] = i;
            pt_id++;
        }
    }
    progress.end();
    std::cout << "Point uniquefication reduced point count from " << old_points->GetNumberOfPoints() 
        << " to " << new_points->GetNumberOfPoints() << std::endl;

    // VTK_CREATE(vtkPoints, newnewpoints);
    // std::vector<idx_t> newpoint_remap;
    // std::map<idx_t, idx_t> new_to_new_to_old_points;
    // unique_vertices(newnewpoints, newpoint_remap, new_to_new_to_old_points, old_points);
    // std::cout << "for comparison, alternative method yielded " << newnewpoints->GetNumberOfPoints() << " unique points" << std::endl;

    // Compute new cells
    VTK_CREATE(vtkCellArray, new_cells);
    idx_t npts = 0;
    idx_t pts[20];
    old_cells->InitTraversal();
    progress.begin(old_cells->GetNumberOfCells(), "Find valid cells");
    for (idx_t i=0; i<old_cells->GetNumberOfCells(); ++i) {
        progress.update(i);
        old_cells->GetCellAtId(i, npts, pts);
        std::vector<idx_t> new_pts(npts);
        for (idx_t j = 0; j < npts; ++j) {
            new_pts[j] = point_remap[pts[j]];
        }
        std::sort(new_pts.begin(), new_pts.end());
        if (std::unique(new_pts.begin(), new_pts.end()) != new_pts.end()) {
            continue;
        }
        new_cells->InsertNextCell(npts, new_pts.data());
        new_to_old_cells[new_cells->GetNumberOfCells() - 1] = i;
    }
    progress.end();

    VTK_CREATE(vtkPolyData, cleaned);
    cleaned->SetPoints(new_points);
    cleaned->SetPolys(new_cells);

    // Copy point data
    progress.begin(input->GetPointData()->GetNumberOfArrays(), "Copy point data");
    for (idx_t i=0; i<input->GetPointData()->GetNumberOfArrays(); ++i) {
        progress.update(i);
        auto old_array = input->GetPointData()->GetArray(i);
        VTK_SMART(vtkDataArray) new_array = vtkDataArray::CreateDataArray(old_array->GetDataType());
        new_array->SetName(old_array->GetName());
        new_array->SetNumberOfTuples(new_points->GetNumberOfPoints());
        for (idx_t j=0; j<new_points->GetNumberOfPoints(); ++j) {
            new_array->SetTuple(j, old_array->GetTuple(new_to_old_points[j]));
        }
        cleaned->GetPointData()->AddArray(new_array);
    }
    progress.end();

    // Copy cell data
    progress.begin(input->GetCellData()->GetNumberOfArrays(), "Copy cell data");
    for (idx_t i=0; i<input->GetCellData()->GetNumberOfArrays(); ++i) {
        progress.update(i);
        auto old_array = input->GetCellData()->GetArray(i);
        VTK_SMART(vtkDataArray) new_array = vtkDataArray::CreateDataArray(old_array->GetDataType());
        new_array->SetName(old_array->GetName());
        new_array->SetNumberOfTuples(new_cells->GetNumberOfCells());
        for (idx_t j=0; j<new_cells->GetNumberOfCells(); ++j) {
            new_array->SetTuple(j, old_array->GetTuple(new_to_old_cells[j]));
        }
        cleaned->GetCellData()->AddArray(new_array);
    }
    progress.end();
    std::cout << "Compute neighborhood data structure... " << std::flush;
    cleaned->BuildLinks();
    std::cout << "done" << std::endl;
    return cleaned;
}

} // namespace vtk_utils