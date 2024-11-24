#include "boundaryAwareRectGrid.h"
#include "triangle_list.h"

#include <vtkPointData.h>
#include <vtkCellData.h>

#include <vector>
#include <string>
#include <cmath>
#include <math.h>

#ifndef STANDALONE_TEST
    #include <nvis-math/fixed_vector.hpp>
    #include <nvis-math/dopri5.hpp>
#endif
#include "point.h"

#include <stdexcept>

vtkStandardNewMacro(boundaryAwareRectGrid);

bool boundaryAwareRectGrid::prepare_knot_array() {

	const int degree = this->bspline_degree;
	const int * rgrid_dim = this->GetDimensions();
	int knot_res[3];
	for (int i = 0; i <ndim; i++)
		knot_res[i] = rgrid_dim[i] - degree + 1;

	Eigen::RowVectorXd x_knot_no_mul = Eigen::RowVectorXd(knot_res[0]);
	Eigen::RowVectorXd y_knot_no_mul = Eigen::RowVectorXd(knot_res[1]);

	//populate x_knot
	for (int i = 0; i < knot_res[0]; i++) {
		x_knot_no_mul(i) = this->GetXCoordinates()->GetTuple1(i);
	}
	//populate y_knot
	for (int i = 0; i < knot_res[1]; i++) {
		y_knot_no_mul(i) = this->GetYCoordinates()->GetTuple1(i);
	}
	x_knot = Knots_with_multiplicity(x_knot_no_mul, degree);
	y_knot = Knots_with_multiplicity(y_knot_no_mul, degree);

	if (ndim == 3)
	{
		Eigen::RowVectorXd z_knot_no_mul = Eigen::RowVectorXd(knot_res[2]);

		for (int i = 0; i < knot_res[2]; i++) {
			z_knot_no_mul(i) = this->GetZCoordinates()->GetTuple1(i);    ////DANA can die if this is the problem
		}

		z_knot = Knots_with_multiplicity(z_knot_no_mul, degree);
	}

	//perpare ctrl points
	// construct control points matrix
	vtkSmartPointer<vtkDataArray> array;
	array = this->GetPointData()->GetArray("CP");
	int ncomp = array->GetNumberOfComponents();

	int CP_size;
	if (ndim == 2)
		CP_size = rgrid_dim[0] * rgrid_dim[1];
	else if (ndim ==3)
		CP_size = rgrid_dim[0] * rgrid_dim[1]* rgrid_dim[2];

	 CP_Eigen_vec= std::vector<Eigen::RowVectorXd>(ncomp, Eigen::RowVectorXd(CP_size));

    for (int i = 0; i < ncomp; i++)

		for (int k = 0; k < CP_size; k++)
		{
			CP_Eigen_vec[i](k) = array->GetComponent(k, i);
		}

    return true;
}

void boundaryAwareRectGrid::SetArray()
{
	bool has_error = false;

	// cutcell code
	vtkSmartPointer<vtkDataArray> cutcell_code_array_untyped = this->GetCellData()->GetArray("cutcell_code");
	if (cutcell_code_array_untyped == nullptr) {
		cerr << "cut-cell code array not found" << endl;
		has_error = true;
	}
	else if (!(cutcell_code_array = vtkUnsignedCharArray::SafeDownCast(cutcell_code_array_untyped))) {
		cerr << "cut-cell code type mismatch" << endl;
		has_error = true;
	}

	// u-val (edge intersection) array
	vtkSmartPointer<vtkDataArray> u_val_array_untyped = this->GetPointData()->GetArray("u_values_vector");
	if (u_val_array_untyped == nullptr) {
		cerr << "u-val (edge intersection) array not found" << endl;
		has_error = true;
	}
	else if (!(u_val_float_array = vtkFloatArray::SafeDownCast(u_val_array_untyped))) {
		cerr << "u-val (edge intersection) array type mismatch" << endl;
		has_error = true;
	}

	// control points
	vtkSmartPointer<vtkDataArray> CP_aray_untyped = this->GetPointData()->GetArray("CP");
	if (CP_aray_untyped == nullptr) {
		cerr << "Control Point array not found " << endl;
		has_error = true;
	}
	else if (!(CP_array = vtkDoubleArray::SafeDownCast(CP_aray_untyped))) {
		cerr << "Control Point array type mismatch" << endl;
		has_error = true;
	}

	int* dim = this->GetDimensions();
	ndim = dim[2] == 1 ? 2 : 3;

	//bspline degree
	auto bspline_degree_array = this->GetFieldData()->GetArray("bspline_degree");
	if (!bspline_degree_array) {
		this->bspline_degree = -1;
	}
	else {
		this->bspline_degree = bspline_degree_array->GetComponent(0, 0);

		if (!prepare_knot_array()) {
			cerr << "Prepare knot array failed" << endl;
			has_error = true;
		}
	}

	if (has_error)
		exit(1);
}

int boundaryAwareRectGrid::get_degree() const
{
	return this->bspline_degree;
}



  /******************************************************************************************************************************/

void get_marching_triangles(int code, vtkSmartPointer<vtkIdList> grid_vid,
	vtkSmartPointer<vtkFloatArray> u_valFloatArray, std::vector<Point3>& tricoords_ls)
{
	const int* T_list = triTable[code]; //list of the edges the marching triangles lie on

	int edge_count;
	for (edge_count = 0; edge_count < 16; edge_count++) {
		if (T_list[edge_count] == -1)
			break;
	}

	tricoords_ls.resize(edge_count);
	for (int i = 0; i < edge_count; i++)
	{
		int edge_i = T_list[i];
		//find edge corresponding u
		int vi = cube_edge_to_grid_edge_idx[edge_i][0];
		int edge_dir = cube_edge_to_grid_edge_idx[edge_i][1];
		int rect_pt_index = grid_vid->GetId(vi);
		double u = u_valFloatArray->GetComponent(rect_pt_index, edge_dir);

        //find correspoinging vertices of the edge
		Point3& pt1 = unit_cube_coords.at(cube_edge2vertex[edge_i][0]);
		Point3& pt2 = unit_cube_coords.at(cube_edge2vertex[edge_i][1]);
		Point3 interp_pt = interpolate_linear(pt1, pt2, u);

		tricoords_ls[i] = interp_pt;
	}
}

void get_marching_lines(int code, vtkSmartPointer<vtkIdList> grid_vid,
	vtkSmartPointer<vtkFloatArray> u_valFloatArray, std::vector<Point3>& tricoords_ls)
{
	const int* T_list = marching_square_table[code]; //list of the edges the marching triangles lie on
	vtkIdType grid_vid_ordered[4];
	int order[4] = { 0,1,3,2 };
	for (int i = 0; i < 4; i++) {
		grid_vid_ordered[i] = grid_vid->GetId(order[i]);
	}

	int edge_count;
	for (edge_count = 0; edge_count < 4; edge_count++) {
		if (T_list[edge_count] == -1)
			break;
	}
	assert(edge_count == 2 || edge_count == 4);

	tricoords_ls.resize(edge_count);
	for (int i = 0; i < edge_count; i++)
	{
		int edge_i = T_list[i];
		//find edge corresponding u
		int vi = square_edge_to_grid_edge_idx[edge_i][0];
		int edge_dir = square_edge_to_grid_edge_idx[edge_i][1];
		int rect_pt_index = grid_vid_ordered[vi];
		double u = u_valFloatArray->GetComponent(rect_pt_index, edge_dir);

		assert(u != -1);

		//find correspoinging vertices of the edge
		Point3& pt1 = unit_square_coords.at(square_edge2vert[edge_i][0]);
		Point3& pt2 = unit_square_coords.at(square_edge2vert[edge_i][1]);
		Point3 interp_pt = interpolate_linear(pt1, pt2, u);

		tricoords_ls[i] = interp_pt;
	}
}

boundaryAwareRectGrid::INOUT boundaryAwareRectGrid::query_pt_inout_3d(double x[3], vtkIdType ci, double* pcoords, double* N_ret)
{
	//cell_id and pcoords can be passed in from another FindCell call to save computation
	// but if not, find the cell and the parametric coord
	double pcoords_backup[3];
	if (ci == -2) {
		int loc[3];
		if (this->ComputeStructuredCoordinates(x, loc, pcoords_backup) == 0) {
			return INOUT::OUT;
		}
		ci = this->ComputeCellId(loc); //get the cell id
		pcoords = pcoords_backup;
	}

	assert(ci >= 0);

	unsigned char code_chr;
	cutcell_code_array->GetTypedTuple(ci, &code_chr);
	int code = code_chr;
	//code = cutcell_code_array->GetComponent(ci, 0);

	if (code == 0)
		return INOUT::OUT;
	else if (code == 255)
		return INOUT::IN;

	vtkNew<vtkIdList> rcell_vids;
	this->GetCellPoints(ci, rcell_vids);

	//determine inside/outside through marching triangles
	std::vector<Point3> triangle_list;
	get_marching_triangles(code, rcell_vids, u_val_float_array, triangle_list);
	assert(triangle_list.size() > 0);

	Point3 pcoords_pt(pcoords);

	int num_tri = triangle_list.size() / 3;
	std::vector<double> dir_list(num_tri);
	double dir;
	int ii = 0;
	double N_avg[3] = { 0 }; //average normal
	for (int i = 0; i < triangle_list.size(); i += 3)
	{
		//dot product the normal of the triangle
		Point3& p1 = triangle_list[i];
		Point3& p2 = triangle_list[i + 1];
		Point3& p3 = triangle_list[i + 2];
		Point3 v1 = p2 - p1;
		Point3 v2 = p3 - p1;
		Point3 N = crossProduct(v1, v2);
		N = N / sqrt(dotProduct(N, N));  //normalize the normal vector
		N_avg[0] += N[0]; N_avg[1] += N[1]; N_avg[2] += N[2];
		Point3 v_query = pcoords_pt - p1;
		double  v_query_norm = sqrt(dotProduct(v_query, v_query));
		v_query = v_query / v_query_norm;
		dir = -1 * dotProduct(N, v_query);
		dir_list[ii++] = dir;
	}

	if (N_ret) { //save the normalized normal
		for (int j = 0; j < 3; j++)
			N_ret[j] += N_avg[j] / double(ii-1);
	}

	if (!all_negative(dir_list))
		return boundaryAwareRectGrid::INOUT::IN;
	else
		return boundaryAwareRectGrid::INOUT::OUT;
}

boundaryAwareRectGrid::INOUT boundaryAwareRectGrid::query_pt_inout_2d(double x[3], vtkIdType ci, double* pcoords)
{
	//cell_id and pcoords can be passed in from another FindCell call to save computation
	// but if not, run the findcell call
	double pcoords_backup[3];
	if (ci == -2) {
		int loc[3];
		if (this->ComputeStructuredCoordinates(x, loc, pcoords_backup) == 0) {
			return INOUT::OUT;
		}
		ci = this->ComputeCellId(loc); //get the cell id
		pcoords = pcoords_backup;
	}


	unsigned char code_chr;
	cutcell_code_array->GetTypedTuple(ci, &code_chr);
	int code = code_chr;
	//code = cutcell_code_array->GetComponent(ci, 0);

	if (code == 0)
		return INOUT::OUT;
	else if (code == 15)
		return INOUT::IN;

	vtkNew<vtkIdList> rcell_vids;
	this->GetCellPoints(ci, rcell_vids);

	//determine inside/outside through marching triangles
	std::vector<Point3> triangle_list;
	triangle_list.reserve(4);
	get_marching_lines(code, rcell_vids, u_val_float_array, triangle_list);

	Point3 pcoords_pt(pcoords[0], pcoords[1], pcoords[2]);

	int num_tri = triangle_list.size() / 2;
	std::vector<double> dir_list(num_tri);
	int ii = 0;
	double mat[3][3] = { 1 };
	for (int i = 0; i < triangle_list.size(); i += 2)
	{
		//dot product the normal of the triangle
		Point3 p1 = triangle_list[i];
		Point3 p2 = triangle_list[i + 1];

		mat[0][0] = pcoords_pt[0]; mat[1][0] = p1[0]; mat[2][0] = p2[0];
		mat[0][1] = pcoords_pt[1]; mat[1][1] = p1[1]; mat[2][1] = p2[1];
		mat[0][2] = 1;            mat[1][2] = 1;    mat[2][2] = 1;
		double dir = determinant_3x3(mat) / 2;
		dir_list[ii++] = dir;
	}

	if (all_negative(dir_list))
		return boundaryAwareRectGrid::INOUT::IN;
	else
		return boundaryAwareRectGrid::INOUT::OUT;

}


//**********************************************************************
vtkIdType boundaryAwareRectGrid::FindCell(double x[3], vtkCell* cell, vtkIdType cellId, double tol2, int& subId, double pcoords[3], double* weights)
{
	assert(false);
	vtkIdType ci = vtkRectilinearGrid::FindCell(x, cell, cellId, tol2, subId, pcoords, weights);
	if (ci == -1)
		return -1;

	INOUT is_inside;
	if (ndim == 3) {
		is_inside = query_pt_inout_3d(x, ci, pcoords);
		if (is_inside == boundaryAwareRectGrid::INOUT::OUT) {
			weights[0] = weights[1] = weights[2] = weights[3] = weights[4] = weights[5] = weights[6] = weights[7] = 0;
		}
	}
	else {
		is_inside = query_pt_inout_2d(x, ci, pcoords);
		if (is_inside == boundaryAwareRectGrid::INOUT::OUT) {
			weights[0] = weights[1] = weights[2] = weights[3] = 0;
		}
	}

	return ci;
}


// 2D implementation for bspline interpolating
template<int DEG>
vtkIdType boundaryAwareRectGrid::bspline_interp_impl_2d(double query_coord[3], double* return_val)
{
	const int NDIM = 2;
	using Spline1D = Eigen::Spline<double, 1, DEG>;
	using BasisVectorType = typename Eigen::SplineTraits<Spline1D>::BasisVectorType;

	const int CP_count_x = this->GetDimensions()[0];

	const BasisVectorType nonzero_basis_x = Spline1D::BasisFunctions(query_coord[0], DEG, x_knot);
	const BasisVectorType nonzero_basis_y = Spline1D::BasisFunctions(query_coord[1], DEG, y_knot);

	int span_start_pt_i[2];
	span_start_pt_i[0] = Spline1D::Span(query_coord[0], DEG, x_knot) - DEG;
	span_start_pt_i[1] = Spline1D::Span(query_coord[1], DEG, y_knot) - DEG;

    Eigen::RowVectorXd CP_Y;
	for (int c = 0; c < this->CP_array->GetNumberOfComponents(); c++)
	{
        CP_Y = Eigen::RowVectorXd::Zero(DEG + 1);
		for (int i = 0; i < (DEG + 1); i++)
		{
			for (int j = 0; j < DEG + 1; j++) {
				double v = CP_array->GetTypedComponent((span_start_pt_i[1] + i) * CP_count_x + span_start_pt_i[0] + j, c);
				CP_Y(i) += nonzero_basis_x(j) * v;
			}
		}

        return_val[c] = CP_Y.dot(nonzero_basis_y.matrix());
		// double ret_v = 0;
		// for (int i = 0; i < DEG + 1; i++) {
		// 	ret_v += CP_Y_arr[i] * nonzero_basis_y(i);
		// }
		// return_val[c] = ret_v;
	}
	return 0;
}

// 2D implementation for bspline derivatives
template<int DEG>
vtkIdType boundaryAwareRectGrid::bspline_deriv_impl_2d(double query_coord[3], std::vector<double>& return_val, int deriv_order)
{
	const int NDIM = 2;
	using Spline1D = Eigen::Spline<double, 1, DEG>;
	using BasisVectorType = typename Eigen::SplineTraits<Spline1D>::BasisVectorType;
    using MatrixDerivType = typename Spline1D::BasisDerivativeType;

	const int CP_count_x = this->GetDimensions()[0];

    // count the number of unique derivatives (including 0th):
    // (1 + 2 + 3 + ... + deriv_order+1)
    // (deriv_order+1)*(deriv_order+2)/2
    int ncomp = this->CP_array->GetNumberOfComponents();
    int nderiv = (deriv_order+1)*(deriv_order+2)/2;
    int nvalues = ncomp*nderiv;
    return_val.resize(nvalues);

    // basis function derivatives
	const MatrixDerivType nonzero_deriv_basis_x = Spline1D::BasisFunctionDerivatives(query_coord[0], deriv_order, DEG, x_knot);
	const MatrixDerivType nonzero_deriv_basis_y = Spline1D::BasisFunctionDerivatives(query_coord[1], deriv_order, DEG, y_knot);

	int span_start_pt_i[2];
	span_start_pt_i[0] = Spline1D::Span(query_coord[0], DEG, x_knot) - DEG;
	span_start_pt_i[1] = Spline1D::Span(query_coord[1], DEG, y_knot) - DEG;

	std::vector<Eigen::RowVectorXd> CP_Y(nderiv);

	for (int c = 0; c < this->CP_array->GetNumberOfComponents(); c++)
	{
        for (int i=0; i<CP_Y.size(); ++i) {
            CP_Y[i] = Eigen::RowVectorXd::Zero(DEG + 1);
        }
		for (int i = 0; i < (DEG + 1); i++)
		{
			for (int j = 0; j < DEG + 1; j++)
            {
				double v = CP_array->GetTypedComponent((span_start_pt_i[1] + i) * CP_count_x + span_start_pt_i[0] + j, c);

                int derivid = 0;
                for (int order=0; order<=deriv_order; ++order) {
                    for (int xorder=order; xorder>=0; --xorder, ++derivid) {
                        CP_Y[derivid](i) += v * nonzero_deriv_basis_x(xorder,j);
                    }
                }
            }
		}
        int derivid=0;
        int offset = c*nderiv;
        for (int order=0; order<=deriv_order; ++order) {
            for (int xorder=order; xorder>=0; --xorder, ++derivid) {
                return_val[offset+derivid] = (CP_Y[derivid]*nonzero_deriv_basis_y.row(order-xorder).matrix())[0];
            }
        }
	}
	return 0;
}

// 3D implementation
template<int DEG>
vtkIdType boundaryAwareRectGrid::bspline_interp_impl_3d(double query_coord[3], double * return_val) {
	const int NDIM = 3;
	using Spline1D = Eigen::Spline<double, 1, DEG>;
	using VectorDegP1Type = typename Eigen::Matrix<double, DEG + 1, 1>;
    using MatrixDegP1Type = typename Eigen::Matrix<double, DEG + 1, DEG + 1>;

	const int * rgrid_dim = this->GetDimensions();
	const int CP_count_x = rgrid_dim[0];
	const int CP_count_y = rgrid_dim[1];

	const VectorDegP1Type basis_val_x = Spline1D::BasisFunctions(query_coord[0], DEG, x_knot);
	const VectorDegP1Type basis_val_y = Spline1D::BasisFunctions(query_coord[1], DEG, y_knot);
	const VectorDegP1Type basis_val_z = Spline1D::BasisFunctions(query_coord[2], DEG, z_knot);

	int span_start_idx[3];
	span_start_idx[0] = Spline1D::Span(query_coord[0], DEG, x_knot) - DEG;
	span_start_idx[1] = Spline1D::Span(query_coord[1], DEG, y_knot) - DEG;
	span_start_idx[2] = Spline1D::Span(query_coord[2], DEG, z_knot) - DEG;

	for (int c = 0; c < this->CP_array->GetNumberOfComponents(); c++)
	{
		MatrixDegP1Type CP_Y_mat = MatrixDegP1Type::Zero(DEG+1, DEG+1);

		for (int d_z = 0; d_z < (DEG + 1); d_z++) // second dim
		{
			for (int d_y = 0; d_y < (DEG + 1); d_y++) //first dim
			{
				int idx = ((span_start_idx[2] + d_z) * CP_count_y  + (span_start_idx[1] + d_y)) * CP_count_x + span_start_idx[0];

				for (int d_x = 0; d_x < DEG + 1; d_x++) {
					double v = CP_array->GetTypedComponent(idx + d_x, c);
					CP_Y_mat(d_z, d_y) += v * basis_val_x(d_x);
				}
			}
		}

		VectorDegP1Type CP_Z_vec = CP_Y_mat * basis_val_y;
	    return_val[c] = basis_val_z.transpose() * CP_Z_vec;
	}
	return 0;
}

template<int DEG>
vtkIdType boundaryAwareRectGrid::bspline_deriv_impl_3d(double query_coord[3], std::vector<double>& return_val, int deriv_order) {
    // std::cout << "Entering bspline_deriv_impl_3d with p=[" << query_coord[0] << ", " << query_coord[1] << ", " << query_coord[2] << "], deriv order=" << deriv_order << "\n";

	const int NDIM = 3;
	using Spline1D = Eigen::Spline<double, 1, DEG>;
	using VectorDegP1Type = typename Eigen::Matrix<double, DEG+1, 1>;
    using MatrixDegP1Type = typename Eigen::Matrix<double, DEG+1, DEG+1>;
    using MatrixDerivType = typename Spline1D::BasisDerivativeType;

	const int * rgrid_dim = this->GetDimensions();
	const int CP_count_x = rgrid_dim[0];
	const int CP_count_y = rgrid_dim[1];
    // std::cout << "CP_count_x = " << CP_count_x << ", CP_count_y = " << CP_count_y << '\n';

    // count number of derivatives to be computed
    int nderiv = 0;
    for (int order=0; order<=deriv_order+1; ++order) {
        nderiv += order*(order+1)/2;
    }
    // std::cout << "nderiv=" << nderiv << '\n';
    int ncomp = this->CP_array->GetNumberOfComponents();
    int nvalues = ncomp * nderiv;
    return_val.resize(nvalues);
    // std::cout << "nvalues=" << nvalues << '\n';

    // basis function and their derivatives
	const MatrixDerivType nonzero_deriv_basis_x = Spline1D::BasisFunctionDerivatives(query_coord[0], deriv_order, DEG, x_knot);
	const MatrixDerivType nonzero_deriv_basis_y = Spline1D::BasisFunctionDerivatives(query_coord[1], deriv_order, DEG, y_knot);
	const MatrixDerivType nonzero_deriv_basis_z = Spline1D::BasisFunctionDerivatives(query_coord[2], deriv_order, DEG, z_knot);
    // std::cout << "nonzero_deriv_basis_x=\n" << nonzero_deriv_basis_x << '\n';
    // std::cout << "nonzero_deriv_basis_y=\n" << nonzero_deriv_basis_y << '\n';
    // std::cout << "nonzero_deriv_basis_z=\n" << nonzero_deriv_basis_z << '\n';


	int span_start_idx[3];
	span_start_idx[0] = Spline1D::Span(query_coord[0], DEG, x_knot) - DEG;
	span_start_idx[1] = Spline1D::Span(query_coord[1], DEG, y_knot) - DEG;
	span_start_idx[2] = Spline1D::Span(query_coord[2], DEG, z_knot) - DEG;
    // std::cout << "span_start_idx=" << span_start_idx[0] << " " << span_start_idx[1] << " " << span_start_idx[2] << '\n';

    std::vector<MatrixDegP1Type> CP_Y_mat(nderiv);

	for (int c = 0; c < this->CP_array->GetNumberOfComponents(); c++)
	{
        for (int i=0; i<nderiv; ++i) {
            CP_Y_mat[i] = MatrixDegP1Type::Zero(DEG+1,DEG+1);
        }

		for (int d_z = 0; d_z < DEG+1; d_z++) // second dim
		{
			for (int d_y = 0; d_y < DEG+1; d_y++) //first dim
			{
				int idx = ((span_start_idx[2] + d_z) * CP_count_y  + (span_start_idx[1] + d_y)) * CP_count_x + span_start_idx[0];

				for (int d_x = 0; d_x < DEG+1; d_x++) {
					double v = CP_array->GetTypedComponent(idx + d_x, c);

                    int derivid=0;
                    for (int order=0; order<=deriv_order; ++order) {
                        for (int xorder=order; xorder>=0; --xorder) {
                            for (int yorder=order-xorder; yorder>=0; --yorder, ++derivid) {
                                CP_Y_mat[derivid](d_z, d_y) +=
                                v * nonzero_deriv_basis_x(xorder, d_x);
                            }
                        }
                    }
				}
			}
		}

        int derivid=0;
        std::vector<VectorDegP1Type> CP_Z(nderiv);
        for (int order=0; order<=deriv_order; ++order) {
            for (int xorder=order; xorder>=0; --xorder) {
                for (int yorder=order-xorder; yorder>=0; --yorder, ++derivid) {
                    VectorDegP1Type CP_Z = CP_Y_mat[derivid]*VectorDegP1Type(nonzero_deriv_basis_y.row(yorder));
                    return_val[nderiv*c + derivid] = nonzero_deriv_basis_z.row(order-xorder-yorder).matrix()*CP_Z;
                }
            }
        }
	}
	return 0;
}

vtkIdType boundaryAwareRectGrid::BsplineInterpolate(double query_coord[3], double* return_val, bool use_bdry_aware,
	bool project_onto_closest_bdry_surface)
{
	vtkIdType ret_id;
    assert(this->bspline_degree != -1);
	double normal_vec[3];

	if (use_bdry_aware)
	{
		if (ndim == 3 && project_onto_closest_bdry_surface) {
			//if out, we want to project interpolated vector onto triangle plane.
			// first we get the plane normal
			INOUT inout = query_pt_inout_3d(query_coord, -2, nullptr, normal_vec);
			project_onto_closest_bdry_surface = inout == INOUT::OUT;
			assert(abs(normal_vec[0]) + abs(normal_vec[1]) + abs(normal_vec[2]) > 0);
		}
		else if ((ndim == 2 && query_pt_inout_2d(query_coord) == INOUT::OUT) ||
			(ndim == 3 && query_pt_inout_3d(query_coord) == INOUT::OUT))
		{
			for (int i = 0; i < this->CP_array->GetNumberOfComponents(); i++)
				return_val[i] = 0;
			return -1;
		}
	}
	else { //not using boundary-aware...
		//simply check for bspline domain
		double bounds[6] = { x_knot[bspline_degree], x_knot[x_knot.size() - bspline_degree - 1] ,
							  y_knot[bspline_degree], y_knot[y_knot.size() - bspline_degree - 1],
							  0, 0};
		if (ndim == 3) {
			bounds[4] = z_knot[bspline_degree];
			bounds[5] = z_knot[z_knot.size() - bspline_degree] - 1;
		}

		for (int d = 0; d < ndim; d++) {
			if (query_coord[d] < bounds[d * 2] || query_coord[d] > bounds[d * 2 + 1])
				return -1;
		}
	}

    if (ndim == 2) {
        switch (bspline_degree) {
		case 1:
			ret_id = bspline_interp_impl_2d<1>(query_coord, return_val); break;
        case 2:
			ret_id = bspline_interp_impl_2d<2>(query_coord, return_val); break;
        case 3:
			ret_id = bspline_interp_impl_2d<3>(query_coord, return_val); break;
        case 4:
			ret_id = bspline_interp_impl_2d<4>(query_coord, return_val); break;
        case 5:
			ret_id = bspline_interp_impl_2d<5>(query_coord, return_val); break;
        case 6:
			ret_id = bspline_interp_impl_2d<6>(query_coord, return_val); break;
        default:
            cerr << "Unsupported b-spline degree" << endl;
			exit(-1);
        }
    }
    else { //3D
        switch (bspline_degree) {
		case 1:
			ret_id = bspline_interp_impl_3d<1>(query_coord, return_val); break;
        case 2:
			ret_id = bspline_interp_impl_3d<2>(query_coord, return_val); break;
        case 3:
			ret_id = bspline_interp_impl_3d<3>(query_coord, return_val); break;
        case 4:
			ret_id = bspline_interp_impl_3d<4>(query_coord, return_val); break;
        case 5:
			ret_id = bspline_interp_impl_3d<5>(query_coord, return_val); break;
        case 6:
			ret_id = bspline_interp_impl_3d<6>(query_coord, return_val); break;
        default:
            cerr << "Unsupported b-spline degree" << endl;
			exit(-1);
        }
    }

	if ( ndim == 3 && project_onto_closest_bdry_surface) {
		// v_ret = v - (dot(N,v) * N);
		double Ndotv = normal_vec[0] * return_val[0] +
						normal_vec[1] * return_val[1] +
						normal_vec[2] * return_val[2];
		for (int d = 0; d < 3; d++)
			return_val[d] -= normal_vec[d] * Ndotv;
		assert(ret_id >= 0);
	}

	return ret_id;

}

vtkIdType boundaryAwareRectGrid::BsplineAllDerivatives(double query_coord[3], std::vector<double>& return_val, int order, bool use_bdry_aware)
{
	vtkIdType ret_id;
    assert(this->bspline_degree != -1);
	double normal_vec[3];

	if (use_bdry_aware)
	{
        if ((ndim == 2 && query_pt_inout_2d(query_coord) == INOUT::OUT) ||
			(ndim == 3 && query_pt_inout_3d(query_coord) == INOUT::OUT))
		{
			return_val.clear();
			return -1;
		}
	}
	else { //not using boundary-aware...
		//simply check for bspline domain
		double bounds[6] = { x_knot[bspline_degree], x_knot[x_knot.size() - bspline_degree - 1] ,
							  y_knot[bspline_degree], y_knot[y_knot.size() - bspline_degree - 1],
							  0, 0};
		if (ndim == 3) {
			bounds[4] = z_knot[bspline_degree];
			bounds[5] = z_knot[z_knot.size() - bspline_degree] - 1;
		}

		for (int d = 0; d < ndim; d++) {
			if (query_coord[d] < bounds[d * 2] || query_coord[d] > bounds[d * 2 + 1])
				return -1;
		}
	}

    if (ndim == 2) {
        switch (bspline_degree) {
		case 1:
			ret_id = bspline_deriv_impl_2d<1>(query_coord, return_val, order); break;
        case 2:
			ret_id = bspline_deriv_impl_2d<2>(query_coord, return_val, order); break;
        case 3:
			ret_id = bspline_deriv_impl_2d<3>(query_coord, return_val, order); break;
        case 4:
			ret_id = bspline_deriv_impl_2d<4>(query_coord, return_val, order); break;
        case 5:
			ret_id = bspline_deriv_impl_2d<5>(query_coord, return_val, order); break;
        case 6:
			ret_id = bspline_deriv_impl_2d<6>(query_coord, return_val, order); break;
        default:
            cerr << "Unsupported b-spline degree" << endl;
			exit(-1);
        }
    }
    else { //3D
        switch (bspline_degree) {
		case 1:
			ret_id = bspline_deriv_impl_3d<1>(query_coord, return_val, order); break;
        case 2:
			ret_id = bspline_deriv_impl_3d<2>(query_coord, return_val, order); break;
        case 3:
			ret_id = bspline_deriv_impl_3d<3>(query_coord, return_val, order); break;
        case 4:
			ret_id = bspline_deriv_impl_3d<4>(query_coord, return_val, order); break;
        case 5:
			ret_id = bspline_deriv_impl_3d<5>(query_coord, return_val, order); break;
        case 6:
			ret_id = bspline_deriv_impl_3d<6>(query_coord, return_val, order); break;
        default:
            cerr << "Unsupported b-spline degree" << endl;
			exit(-1);
        }
    }

	return ret_id;
}

void boundaryAwareRectGrid::DeepCopy( vtkDataObject * src )
{
	 vtkRectilinearGrid::DeepCopy(src);
     this->SetArray();
}

void boundaryAwareRectGrid::ShallowCopy( vtkDataObject * src)
{
	 vtkRectilinearGrid::ShallowCopy(src);
     this->SetArray();
}
