#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>

#include <math/fixed_vector.hpp>
#include <teem/nrrd.h>
#include <netcdf.h>
#include <stdexcept>
#include <sstream>
#include <fstream>

#include "nek5000reader.hpp"
#include "ncio.hpp"
#include "vtk/vtk_utils.hpp"

char *metafile, *outfile;
int timestep;
bool do_grid;


void initialize(int argc, char* argv[])
{
	hestOpt *hopt = NULL;
	hestParm *hparm;
	airArray *mop;
	char *me;

	mop = airMopNew();
	me = argv[0];
	hparm = hestParmNew();
	airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
	hparm->elideSingleOtherType = AIR_TRUE;
	hestOptAdd(&hopt, "i", "metafile",  airTypeString, 1, 1, &metafile, NULL, "metafile name");
	hestOptAdd(&hopt, "o", "output",    airTypeString, 1, 1, &outfile,  NULL, "output file name");
	hestOptAdd(&hopt, "t", "time step", airTypeInt,    1, 1, &timestep, NULL, "time step");

	hestParseOrDie(hopt, argc - 1, const_cast<const char**>(argv + 1), hparm,
	               me, "Convert Nek5000 data file to VTK binary format",
	               AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

struct vec_equal {
	bool operator()(const nvis::fvec3& a, const nvis::fvec3& b) const {
		return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
	}
};

int main(int argc, char* argv[])
{
	initialize(argc, argv);

    spurt::nek5000reader reader(metafile);

	std::vector<unsigned int> hexind; // hexahedra vertices' indices
	std::vector<float>        hexpts; // vertices' coordinates
    std::vector<unsigned int> rindex; // indices of unique vertex ids
    std::vector<nvis::fvec3>  vals;   // vector values

    reader.read_grid(hexind, hexpts, rindex);
    reader.read_vectors(vals, timestep, rindex);
    
    VTK_CREATE(vtkFloatArray, coords);
    coords->SetNumberOfComponents(3);
    coords->SetNumberOfTuples(hexpts.size()/3);
    vtkIdType _count=0;
    float x[3];
    for (size_t i=0; i<hexpts.size(); i+=3, ++_count) {
        x[0] = hexpts[i];
        x[1] = hexpts[i+1];
        x[2] = hexpts[i+2];
        coords->SetTypedTuple(_count, x);
    }
    VTK_PTR(vtkPoints, points);
    points->SetData(coords);
    
    VTK_CREATE(vtkCellArray, cells);
    for (int i=0; i<hexind.size(); i+=8) {
        cells->InsertNextCell(8);
        cells->InsertCellPoint(hexind[i]);
        cells->InsertCellPoint(hexind[i+1]);
        cells->InsertCellPoint(hexind[i+2]);
        cells->InsertCellPoint(hexind[i+3]);
        cells->InsertCellPoint(hexind[i+4]);
        cells->InsertCellPoint(hexind[i+5]);
        cells->InsertCellPoint(hexind[i+6]);
        cells->InsertCellPoint(hexind[i+7]);
    }
    
    VTK_CREATE(vtkUnstructuredGrid, grid);
    grid->SetPoints(points);
    grid->SetCells(VTK_HEXAHEDRON, cells);
    vtk_utils::add_vectors(grid, vals, true, "velocity");
    
    VTK_CREATE(vtkXMLUnstructuredGridWriter, writer);
    writer->SetFileName(outfile);
    writer->SetInputData(grid);
    writer->SetCompressorTypeToZLib();
    writer->Write();
    
    /*
    
    size_t npts = rindex.size();

	std::fstream file;
	file.open(outfile, std::ios::out);

	file << "# vtk DataFile Version 2.0\n"
	<< "Converted from Nek5000 file " << metafile << ", timestep = " << timestep << '\n'
	<< "BINARY\n"
	<< "DATASET UNSTRUCTURED_GRID\n";

	// float *hexpts_ptr = (float*)calloc(hexpts.size() * 3, sizeof(float));
	// std::copy(hexpts.begin(), hexpts.end(), reinterpret_cast<nvis::fvec3 *>(hexpts_ptr));

	file << "POINTS " << hexpts.size() << " float\n";
	file.write(reinterpret_cast<char*>(&hexpts.front()), sizeof(float)*hexpts.size());
	file << '\n';
	// delete[] hexpts_ptr;

	size_t ncells = hexind.size() / 8;

	int* hexind_ptr = (int*)calloc(9 * ncells, sizeof(int));
	file << "CELLS " << ncells << " " << 9*ncells << '\n';
	for (int i = 0 ; i < ncells ; ++i) {
		hexind_ptr[9*i] = 8;
		for (int j = 0 ; j < 8 ; ++j) {
			hexind_ptr[9*i+j+1] = hexind[8*i+j];
		}
	}
	file.write(reinterpret_cast<char*>(hexind_ptr), sizeof(int)*9*ncells);
	delete[] hexind_ptr;

	file << '\n'
	<< "CELL_TYPES " << ncells << '\n';
	int* types = (int*)calloc(ncells, sizeof(int));
	std::fill(&types[0], &types[ncells], 12);
	file.write(reinterpret_cast<char*>(types), sizeof(int)*ncells);
	delete[] types;

	// float *vecf = (float*)calloc(3*npts, sizeof(float));
	//     std::copy(vals.begin(), vals.end(), reinterpret_cast<nvis::fvec3*>(vecf));
	file << "POINT_DATA " << npts << '\n';
	file << "VECTORS velocity float\n";
	file.write(reinterpret_cast<char*>(&vals.front()), sizeof(float)*3*npts);

	file.close();
    */

	return 0;
}

//---------------------------------------------------------------------------
