#ifndef __BOUNDARY_AWARE_RECT_GRID_H__
#define __BOUNDARY_AWARE_RECT_GRID_H__


#include <vtkRectilinearGrid.h>
#include <vtkBitArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkUnsignedCharArray.h>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <unsupported/Eigen/AutoDiff>
#include <unsupported/Eigen/Splines>
#include <vector>

class boundaryAwareRectGrid : public vtkRectilinearGrid
{
public:
	enum class INOUT { OUT = 0, IN = 1 };

private:
    vtkSmartPointer<vtkFloatArray> u_val_float_array = nullptr;
	vtkSmartPointer<vtkUnsignedCharArray> cutcell_code_array = nullptr;

	vtkSmartPointer<vtkDoubleArray> CP_array = nullptr;
	int bspline_degree = -1;
	int ndim = -1;

	//Eigen related variables
    Eigen::RowVectorXd x_knot;
	Eigen::RowVectorXd y_knot;
	Eigen::RowVectorXd z_knot;

	std::vector<Eigen::RowVectorXd> CP_Eigen_vec;

	bool prepare_knot_array();

	template<int DEG>
	vtkIdType bspline_interp_impl_2d(double query_coord[3], double* return_val);
	template<int DEG>
	vtkIdType bspline_interp_impl_3d(double query_coord[3], double* return_val);
	template<int DEG>
	vtkIdType bspline_deriv_impl_2d(double query_coord[3], std::vector<double>& return_val, int order);
	template<int DEG>
	vtkIdType bspline_deriv_impl_3d(double query_coord[3], std::vector<double>& return_val, int order);

	//double* N_ret variable is added specific for flow attachment fix. (for 3D only)
	//The avg normalized normal vector of the boudnary triangles is returned in N_ret if it's not null
	INOUT query_pt_inout_3d(double x[3], vtkIdType ci = -2, double* pcoords = nullptr, double* N_ret = nullptr);
	INOUT query_pt_inout_2d(double x[3], vtkIdType ci = -2, double* pcoords = nullptr);

public:

    static boundaryAwareRectGrid* New();
    vtkTypeMacro(boundaryAwareRectGrid, vtkRectilinearGrid); //defines GetClassName, IsTypeOf and SafeDownCast

    void SetArray(); //parse the u-array and pt-label array

    // Override the definition of FindCell to use boundary data
    vtkIdType FindCell(double x[3], vtkCell* cell, vtkIdType cellId, double tol2, int& subId, double pcoords[3], double* weights) override;

	// If use_bdry_aware==true and project_on_bdry==true, the interp_val (of 3 components) is assumed to be a vector,
	// and any point that is considered "out" will have its vector value projected onto the boundary surface.
	// If use_bdry_aware==true and project_on_bdry==false, points that are "out" are set to zero.
	vtkIdType BsplineInterpolate(double x[3], double * ret_val, bool use_bdry_aware= true, bool project_on_bdry = false);

	vtkIdType BsplineAllDerivatives(double x[3], std::vector<double>& ret_val, int order, bool use_bdry_aware= true);

	int get_degree() const;

    void DeepCopy(vtkDataObject * src) override;
    void ShallowCopy(vtkDataObject * src) override;
};


#endif
