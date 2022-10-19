#include <VTK/vtk_utils.hpp>


int main(int argc, char* argv[]) {
    vtkDataSet* dataset = vtk_utils::readVTK(argv[1]);
    vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast(dataset);
    vtk_utils::add_jacobian(grid);
    vtk_utils::add_vorticity(grid);
    vtk_utils::add_lambda2(grid);
    
    grid->PrintSelf(std::cout, vtkIndent(0));
    
    vtk_utils::saveVTK(grid, argv[2]);
    dataset->Delete();
    
    return 0;
}