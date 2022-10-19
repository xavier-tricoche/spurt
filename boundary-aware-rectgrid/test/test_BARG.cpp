/*
* Reads a normal vtkXMLRectilinearGrid or our boundaryAwareRectilinearGrid, do some queries inside the grid, and draw the query points.
* Queries can be random points or an image
*/
//Data used for the testing 
//ICE 
//"C:\xavier_data\recon\v_ICE_equal_axis_May_31.vtr"
//"C:\Data\ICE\ice_30deg_vel_pres.vtu"
// Delta
//"C:\xavier_data\recon\v_EDelta_modified_equal_axis_June_15.vtr"
//"C:\Data\EDelta\delta65.pval.unsteady_i=700_t=7.0100e-01.modified_velocity.vtu"
//cylinder 
//"C:\xavier_data\recon\flow_past_cyn_p_4_June_21_sol.vtr"
//"C:\Data\flow_past_cylinder\flow.vtu"

#include "boundaryAwareRectGrid.h"
#include <vtkSmartPointer.h>
#include <vtkProperty.h>
#include <vtkCellTreeLocator.h>
#include <vtkDataSetMapper.h>
#include <vtkXMLRectilinearGridReader.h>
#include <vtkXMLUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkScalarBarActor.h>
#include <vtkColorTransferFunction.h>
#include <vtkLookupTable.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkXMLRectilinearGridWriter.h>
#include <vtkSphereSource.h>
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkPointSource.h>
#include <vtkBoundedPointSource.h>
#include <vtkGlyph3D.h>
#include <vtkImplicitPlaneWidget2.h>
#include <vtkImplicitPlaneRepresentation.h>
#include <vtkClipPolyData.h>
#include <vtkPlane.h>
#include <vtkPlaneSource.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkImageDataGeometryFilter.h>
#include <vtkImageData.h>
#include <chrono>
#include <math.h>
#include <Eigen/Core>

using namespace std;
using namespace std::chrono;

/* Generate random double in a range */
double fRand(double fMin, double fMax)
{
	double f = (double)rand() / RAND_MAX;
	return fMin + f * (fMax - fMin);
}


//Read the rectilinear grid file with or without boundary-aware
template<typename RGridType>
vtkSmartPointer<RGridType> read_rect_grid(const char* inputFilename);

template<typename GridType>
vtkSmartPointer<vtkPolyData> query_image_poly_data(int n, vtkSmartPointer<GridType> grid);


template<>
vtkSmartPointer<vtkRectilinearGrid> read_rect_grid(const char* inputFilename)
{ //read vtkRectilinearGrid
	vtkNew<vtkXMLRectilinearGridReader> reader;
	reader->SetFileName(inputFilename);
	reader->Update();
	auto rect_grid = reader->GetOutput();
	return rect_grid;
}

template<>
vtkSmartPointer<boundaryAwareRectGrid> read_rect_grid(const char* inputFilename)
{ // Read vtkRectilinearGrid, then convert it to boundaryAwareRectGrid
	vtkSmartPointer<vtkRectilinearGrid> rect_grid = read_rect_grid<vtkRectilinearGrid>(inputFilename);
    if (!rect_grid) return nullptr;
	vtkNew<boundaryAwareRectGrid> bdry_aware_rgrid;
	bdry_aware_rgrid->ShallowCopy(rect_grid);
	
	return bdry_aware_rgrid;
}
//reader for the unstr grid
vtkSmartPointer<vtkUnstructuredGrid> read_unstr_grid(const char* inputFilename)
{ //read vtkRectilinearGrid
    vtkNew<vtkXMLUnstructuredGridReader> un_reader;
    un_reader->SetFileName(inputFilename);
    un_reader->Update();
    auto unstr = un_reader->GetOutput();
    return unstr;
}

class vtkIPWCallback : public vtkCommand
{
public:
	static vtkIPWCallback *New()
	{
		return new vtkIPWCallback;
	}
	virtual void Execute(vtkObject *caller, unsigned long, void*)
	{
		vtkImplicitPlaneWidget2 *planeWidget =
			reinterpret_cast<vtkImplicitPlaneWidget2*>(caller);
		vtkImplicitPlaneRepresentation *rep =
			reinterpret_cast<vtkImplicitPlaneRepresentation*>(planeWidget->GetRepresentation());
		rep->GetPlane(this->Plane);
	}
	vtkIPWCallback() :Plane(0), Actor(0) {}
	vtkPlane *Plane;
	vtkActor *Actor;

};

vtkSmartPointer<vtkPolyData> rand_point_poly_data(int n, double* bound) {
    //generate n random points within specified bounds

    vtkSmartPointer<vtkBoundedPointSource> pointSource =
        vtkSmartPointer<vtkBoundedPointSource>::New();
    pointSource->SetBounds(bound);
    pointSource->SetNumberOfPoints(n);
    pointSource->ProduceCellOutputOn();
    pointSource->Update();
    pointSource->GetOutput()->Print(cout);

    return pointSource->GetOutput();
}

template<>
vtkSmartPointer<vtkPolyData> query_image_poly_data(int n, vtkSmartPointer<vtkRectilinearGrid> grid) 
{
    // n is the goal number of total points.
    double* bounds = grid->GetBounds();
    double xlen = bounds[1] - bounds[0];
    double ylen = bounds[3] - bounds[2];
    double zlen = bounds[5] - bounds[4];

    if (grid->GetDataDimension() == 2) {
        double k = xlen / ylen;
        double nx = sqrt(double(n) / k);
        double ny = k * nx;

        vtkNew<vtkPlaneSource> plane_source;
        plane_source->SetResolution(int(nx), int(ny));
        plane_source->SetOrigin(bounds[0], bounds[2], bounds[4]);
        plane_source->SetPoint1(bounds[0], bounds[3], bounds[4]);
        plane_source->SetPoint2(bounds[1], bounds[2], bounds[4]);

        plane_source->Update();
        return plane_source->GetOutput();
    }
    else {
        double k1 = xlen / ylen;
        double k2 = xlen / zlen;
        double nx = cbrt(double(n) / k1 / k2);
        double ny = k1 * nx;
        double nz = k2 * nx;

        vtkNew<vtkImageData> image_data;
        image_data->SetDimensions(nx, ny, nz);
        image_data->SetOrigin(bounds[0], bounds[2], bounds[4]);
        image_data->SetSpacing(
            (bounds[1] - bounds[0]) / nx,
            (bounds[3] - bounds[2]) / ny,
            (bounds[5] - bounds[4]) / nz
        );

        vtkNew<vtkImageDataGeometryFilter> imageDataGeometryFilter;
        imageDataGeometryFilter->SetInputData(image_data);
        imageDataGeometryFilter->Update();

        vtkNew<vtkVertexGlyphFilter> vert_filter;
        vert_filter->SetInputConnection(imageDataGeometryFilter->GetOutputPort());
        vert_filter->Update();
        return vert_filter->GetOutput();
    }


}
template<>
vtkSmartPointer<vtkPolyData> query_image_poly_data(int n, vtkSmartPointer<vtkUnstructuredGrid> grid)
{
    // n is the goal number of total points.
    int dimens = 3;
    double* bounds = grid->GetBounds();
    double xlen = bounds[1] - bounds[0];
    double ylen = bounds[3] - bounds[2];
    double zlen = bounds[5] - bounds[4];
    if (zlen == 0) dimens = 2;
    if (dimens == 2) {
        double k = xlen / ylen;
        double nx = sqrt(double(n) / k);
        double ny = k * nx;

        vtkNew<vtkPlaneSource> plane_source;
        plane_source->SetResolution(int(nx), int(ny));
        plane_source->SetOrigin(bounds[0], bounds[2], bounds[4]);
        plane_source->SetPoint1(bounds[0], bounds[3], bounds[4]);
        plane_source->SetPoint2(bounds[1], bounds[2], bounds[4]);

        plane_source->Update();
        return plane_source->GetOutput();
    }
    else {
        double k1 = xlen / ylen;
        double k2 = xlen / zlen;
        double nx = cbrt(double(n) / k1 / k2);
        double ny = k1 * nx;
        double nz = k2 * nx;

        vtkNew<vtkImageData> image_data;
        image_data->SetDimensions(nx, ny, nz);
        image_data->SetOrigin(bounds[0], bounds[2], bounds[4]);
        image_data->SetSpacing(
            (bounds[1] - bounds[0]) / nx,
            (bounds[3] - bounds[2]) / ny,
            (bounds[5] - bounds[4]) / nz
        );

        vtkNew<vtkImageDataGeometryFilter> imageDataGeometryFilter;
        imageDataGeometryFilter->SetInputData(image_data);
        imageDataGeometryFilter->Update();

        vtkNew<vtkVertexGlyphFilter> vert_filter;
        vert_filter->SetInputConnection(imageDataGeometryFilter->GetOutputPort());
        vert_filter->Update();
        return vert_filter->GetOutput();
    }
}

void print_usage() {
    std::cout << "Usage: test_BARG"
        << " Filename.vtr [random/image]" << std::endl;
}


bool ends_with(std::string str, std::string piece)
{
    size_t n = str.rfind(piece.c_str());
    size_t len = str.length();
    return n == (len - piece.size());
}
//**********************************************************************************

enum class Interp_Type { linear, bspline, unstructure};

int main(int argc, char* argv[])
{
    if (argc == 1) {
        print_usage();
        return EXIT_FAILURE;
    }

    char* inputFilename = argv[1];
    int  ndim = 3;
    const bool use_boundary_aware_interp = true;
    Interp_Type interpolation_type = Interp_Type::bspline;
    vtkSmartPointer<boundaryAwareRectGrid> rgrid;
    vtkSmartPointer<vtkUnstructuredGrid> unstr_grid;
    if (ends_with(inputFilename, "vtu"))
    {
        unstr_grid = read_unstr_grid(inputFilename);
        interpolation_type = Interp_Type::unstructure;

        if (!unstr_grid) {
            cerr << "Input file not found" << endl;
            exit(1);
        }

    }
    else if (ends_with(inputFilename, "vtr"))
    {
        rgrid = read_rect_grid<boundaryAwareRectGrid>(inputFilename);
        if (!rgrid) {
            cerr << "Input file not found" << endl;
            exit(1);
        }
    }
    char query_source = 'i';

    //const bool use_boundary_aware_interp = true;
    //Interp_Type interpolation_type = Interp_Type::bspline; 
    if (interpolation_type != Interp_Type::unstructure)
    {
        auto rgrid_array = rgrid->GetPointData()->GetArray(0);
        if (interpolation_type == Interp_Type::bspline && rgrid->get_degree() < 0) {
            cout << "Bspline interpolation is chosen, but input dataset does not contain bspline data." << endl;
            cout << "Linear interpolation is used instead." << endl;
            interpolation_type = Interp_Type::linear;
        }

        if (interpolation_type == Interp_Type::bspline) {
            cout << "interpolation type is bspline of degree " << rgrid->get_degree() << " and use_boundar_aware = " << use_boundary_aware_interp << " \n";
        }
        else {
            cout << "interpolation type is linear " << "and use_boundar_aware =" << use_boundary_aware_interp << endl;
        }
    }
    else
    {
        cout << "interpolation type is unstructured" << endl;

    }

    vtkNew<vtkDoubleArray> interpolated_pt_array;
    vtkSmartPointer<vtkPolyData> pt_polydata;
    double* rng;
    //double grid_len;
    double* grid_bounds;
    //vtkSmartPointer<vtkDoubleArray> rgrid_array;
    //vtkSmartPointer<vtkDoubleArray> ugrid_array;


    if (interpolation_type != Interp_Type::unstructure)
    {
        //range and extent of the data 
        auto rgrid_array = rgrid->GetPointData()->GetArray(0);
        double* rng = rgrid_array->GetRange();
        double rgrid_len = rgrid->GetLength();
        double* grid_bounds = rgrid->GetBounds();
        const int* rgrid_res = rgrid->GetDimensions();
        int ndim = rgrid->GetDimensions()[2] == 1 ? 2 : 3;
        
        
        //vtkSmartPointer<vtkPolyData> pt_polydata;
#ifdef NDEBUG
        int poly_points_num = 1000000;
#else
        int poly_points_num = 50000;
#endif
        if (query_source == 'r')
            pt_polydata = rand_point_poly_data(poly_points_num, grid_bounds);
        else
            pt_polydata = query_image_poly_data<vtkRectilinearGrid>(poly_points_num, rgrid);


        //pt_polydata->Print(cout);

        //vtkNew<vtkDoubleArray> interpolated_pt_array;
        interpolated_pt_array->SetNumberOfComponents(1);
        interpolated_pt_array->SetNumberOfTuples(pt_polydata->GetNumberOfPoints());
        auto start = high_resolution_clock::now();
        const int ncomp = rgrid_array->GetNumberOfComponents();
        double interpolated_val[3] = { 0 };
        for (int i = 0; i < pt_polydata->GetNumberOfPoints(); i++)
        {
            double pt[3];
            pt_polydata->GetPoint(i, pt);
            if (interpolation_type == Interp_Type::bspline)
            {
                auto ci = rgrid->BsplineInterpolate(pt, interpolated_val, use_boundary_aware_interp);
                if (ci == -1) { //cell not found
                    //do nothing
                }
            }
            else if (interpolation_type == Interp_Type::linear)
            {
                int subid;
                double pcoords[3];
                double weights[8];
                vtkIdType ci;
                if (use_boundary_aware_interp)
                    ci = rgrid->vtkRectilinearGrid::FindCell(pt, nullptr, 0, 0, subid, pcoords, weights);
                else
                    ci = rgrid->vtkRectilinearGrid::FindCell(pt, nullptr, 0, 0, subid, pcoords, weights);

                if (ci >= 0) {
                    vtkNew<vtkIdList> rcell_vids;
                    rgrid->GetCellPoints(ci, rcell_vids);

                    for (int k = 0; k < ncomp; k++) {
                        interpolated_val[k] = 0;
                        for (int j = 0; j < rcell_vids->GetNumberOfIds(); j++) {
                            double vi = rcell_vids->GetId(j);
                            interpolated_val[k] += rgrid_array->GetComponent(vi, k) * weights[j];
                        }
                    }
                }
            }
            //reduce to 1 component for visualization purpose
            double val_1comp = 0;
            for (int k = 0; k < ncomp; k++) {
                val_1comp += pow(interpolated_val[k], 2);
            }
            val_1comp = sqrt(val_1comp);
            interpolated_pt_array->SetTuple1(i, val_1comp);
        }

        pt_polydata->GetPointData()->SetScalars(interpolated_pt_array);

        //print solution duration 
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        std::cout << "Interpolation in rectilinear grid " << poly_points_num << " points took " << duration.count() / 1000000.0 << " seconds" << endl;
    }














    else //if the interpolation is unstr
    {
    //range and extent of the data 
    auto ugrid_array = unstr_grid->GetPointData()->GetArray(1);
    rng = ugrid_array->GetRange();
    double ugrid_len = unstr_grid->GetLength();
    grid_bounds = unstr_grid->GetBounds();
    double zlen = grid_bounds[5] - grid_bounds[4];
    if (zlen == 0) 
        ndim = 2;
    
    //vtkSmartPointer<vtkPolyData> pt_polydata;
#ifdef NDEBUG
    int poly_points_num = 1000000;
#else
    int poly_points_num = 50000;
#endif
    if (query_source == 'r')
        pt_polydata = rand_point_poly_data(poly_points_num, grid_bounds);
    else
        pt_polydata = query_image_poly_data<vtkUnstructuredGrid>(poly_points_num, unstr_grid);


    //pt_polydata->Print(cout);

    //vtkNew<vtkDoubleArray> interpolated_pt_array;
    interpolated_pt_array->SetNumberOfComponents(1);
    interpolated_pt_array->SetNumberOfTuples(pt_polydata->GetNumberOfPoints());
    vtkSmartPointer<vtkCellTreeLocator> m_locator = vtkSmartPointer<vtkCellTreeLocator>::New();
    m_locator->SetDataSet(unstr_grid);
    auto start_time_cell_tree = high_resolution_clock::now();
    m_locator->BuildLocator();
    auto stop_time_cell_tree = high_resolution_clock::now();
    auto cell_tree_duration = duration_cast<microseconds>(stop_time_cell_tree - start_time_cell_tree);
    std::cout << "cell tree duration" << cell_tree_duration.count() / 1000000.0 << " seconds" << endl;
     auto start = high_resolution_clock::now();
    const int ncomp = ugrid_array->GetNumberOfComponents();
    double interpolated_val[3] = { 0 };
    for (int i = 0; i < pt_polydata->GetNumberOfPoints(); i++)
    {
        double pt[3];
        pt_polydata->GetPoint(i, pt);
        int subid;
        double pcoords[3];
        double weights[8];
        double tol = 0;
        vtkGenericCell* g_cell;
        vtkIdType ci;
        //ci =unstr_grid->FindCell(pt, nullptr, 0, 0, subid, pcoords, weights);
        ci = m_locator->vtkCellTreeLocator::FindCell(pt);
        //not sure how to get the weights 
        if (ci >= 0) {
            vtkNew<vtkIdList> rcell_vids;
            unstr_grid->GetCellPoints(ci, rcell_vids);

            for (int k = 0; k < ncomp; k++) {
                interpolated_val[k] = 0;
                for (int j = 0; j < rcell_vids->GetNumberOfIds(); j++) {
                    double vi = rcell_vids->GetId(j);
                    interpolated_val[k] += ugrid_array->GetComponent(vi, k) * weights[j];
                }
            }
        }
        //reduce to 1 component for visualization purpose
        double val_1comp = 0;
        for (int k = 0; k < ncomp; k++) {
            val_1comp += pow(interpolated_val[k], 2);
        }
        val_1comp = sqrt(val_1comp);
        interpolated_pt_array->SetTuple1(i, val_1comp);
    }

    pt_polydata->GetPointData()->SetScalars(interpolated_pt_array);

    //print solution duration 
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "Interpolation in unstructured grid " << poly_points_num << " points took " << duration.count() / 1000000.0 << " seconds" << endl;
    }
    // Write the interpolated points out to file
	vtkNew< vtkXMLPolyDataWriter> poly_writer;
	poly_writer->SetFileName("poly_data_out.vtp");
	poly_writer->SetInputData(pt_polydata);
	poly_writer->Write();

	// Visualize
#ifdef NDEBUG  //skip rendering if not in debug mode, for timing.
    if(0)
#endif
    {
        bool use_clipper = ndim == 3;

        //colormap
        rng = interpolated_pt_array->GetRange();
        cout << "color array range = [" << rng[0] << " " << rng[1] << "]" << endl;

        vtkNew<vtkColorTransferFunction> ctf;
        ctf->AddRGBPoint(rng[0], 0, 0, 0);
        double w = 0.99999;
		ctf->AddRGBPoint(rng[0]*w + rng[1]*(1.0-w), 1, 0, 0);
        ctf->AddRGBPoint((rng[0] + rng[1]) / 2.0, 1, 1, 1);
        ctf->AddRGBPoint(rng[1], 0, 0, 1);
        ctf->Build();

        //clipper 
        vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        plane->SetNormal(1, 0, 0);
        plane->SetOrigin(grid_bounds[0], grid_bounds[2], grid_bounds[4]);
        vtkSmartPointer<vtkClipPolyData> clipper = vtkSmartPointer<vtkClipPolyData>::New();
        clipper->SetClipFunction(plane);
        clipper->SetInputData(pt_polydata);

        //Mapper
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputConnection(clipper->GetOutputPort());
        //mapper->SetScalarRange(pt_polydata->GetScalarRange());
        mapper->SetLookupTable(ctf);
        
        vtkNew<vtkActor> pointActor;
        pointActor->SetMapper(mapper);
        pointActor->GetProperty()->SetPointSize(10.0);
        
        vtkNew<vtkScalarBarActor> scalarBar;
        scalarBar->SetLookupTable(mapper->GetLookupTable());
        scalarBar->SetNumberOfLabels(5);

        vtkNew<vtkRenderer> renderer;
        vtkNew<vtkRenderWindow> renderWindow;
        renderWindow->AddRenderer(renderer);
        renderWindow->SetSize(1000, 700);
        vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;

        vtkNew<vtkInteractorStyleTrackballCamera> style;
        renderWindowInteractor->SetInteractorStyle(style);
        renderWindowInteractor->SetRenderWindow(renderWindow);

        renderer->AddActor(pointActor);
        renderer->AddActor(scalarBar);
        renderer->SetBackground(.4, .4, .4); // Background gray green

        renderWindow->Render();


        // the callback for clipper plane
        vtkSmartPointer<vtkIPWCallback> clipper_callback = vtkSmartPointer<vtkIPWCallback>::New();
        clipper_callback->Plane = plane;
        clipper_callback->Actor = pointActor;

        vtkSmartPointer<vtkImplicitPlaneRepresentation> rep = vtkSmartPointer<vtkImplicitPlaneRepresentation>::New();
        rep->SetPlaceFactor(1.1); // This must be set prior to placing the widget
        rep->PlaceWidget(pt_polydata->GetBounds());
        rep->SetNormal(plane->GetNormal());
        rep->OutlineTranslationOff();

        vtkSmartPointer<vtkImplicitPlaneWidget2> planeWidget = vtkSmartPointer<vtkImplicitPlaneWidget2>::New();
        planeWidget->SetInteractor(renderWindowInteractor);
        planeWidget->SetRepresentation(rep);
        planeWidget->AddObserver(vtkCommand::InteractionEvent, clipper_callback);

        renderWindowInteractor->Initialize();
        renderWindow->Render();
        if (use_clipper)
            planeWidget->On();

        // Begin mouse interaction 
        renderWindowInteractor->Start();
    }


	return EXIT_SUCCESS;
}



