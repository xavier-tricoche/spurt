#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>

#include <math/fixed_vector.hpp>
#include <teem/nrrd.h>
#include <netcdf.h>
#include <stdexcept>
#include <sstream>

#include "avtNek5000FileFormat.hpp"
#include "ncio.hpp"

using namespace std;

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
	hestOptAdd(&hopt, "i", 		"metafile", 		airTypeString, 	1, 1, &metafile, 	NULL, 		"metafile name");
	hestOptAdd(&hopt, "o", 		"basename", 		airTypeString, 	1, 1, &outfile, 	NULL, 		"basename for output file");
	hestOptAdd(&hopt, "t",		"time step",		airTypeInt,		1, 1, &timestep,	NULL,		"time step");
	hestOptAdd(&hopt, "g", 		"convert grid",		airTypeBool,	1, 1, &do_grid,		NULL,		"convert grid vs. data file");

	hestParseOrDie(hopt, argc - 1, const_cast<const char **>(argv + 1), hparm,
	               me, "Convert Nek5000 data file to NetCDF format",
	               AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline void check_nc(int status)
{
	if (status != NC_NOERR)
		throw std::runtime_error(nc_strerror(status));
}

bool ncdefvarexist(int ncid, const char *name,
                   nc_type xtype, int ndims, const int *dimidsp, int *varidp)
{
	int status = nc_def_var(ncid, name, xtype, ndims, dimidsp, varidp);

	if (status == NC_ENAMEINUSE) {

		check_nc(nc_inq_varid(ncid, name, varidp));

		int newNdims;
		int newDimids[20];

		check_nc(nc_inq_varndims(ncid, *varidp, &newNdims));

		if (newNdims != ndims || newNdims > 20)
			throw std::runtime_error(string(name) + " number of dimensions differ");

		check_nc(nc_inq_vardimid(ncid, *varidp, newDimids));

		for (int i = 0;i < ndims;i++)
			if (newDimids[i] != dimidsp[i])
				throw std::runtime_error(string(name) + " dimensions differ");
	}
	else
		check_nc(status);

	return status == NC_NOERR;

}

bool ncdefdimexist(int fileid, const char*name, size_t length, int * id)
{
	int status = nc_def_dim(fileid, name, length, id);

	if (status == NC_ENAMEINUSE) {
		check_nc(nc_inq_dimid(fileid, name, id));
		size_t newlen;
		check_nc(nc_inq_dimlen(fileid, *id, &newlen));
		if (length != newlen)
			throw std::runtime_error(string(name) + " lengths differ!");
	}
	else
		check_nc(status);

	return status == NC_NOERR;
}

//---------------------------------------------------------------------------

void save_grid(const std::string& filename,
               const std::vector<float>& pos,
               const std::vector<unsigned int>& ids)
{
	size_t nPoints = pos.size() / 3;
	size_t nCells = ids.size() / 8;

	int fileId;
	int status = nc_create(filename.c_str(), NC_WRITE, &fileId);

	if (status == NC_EEXIST) {
		check_nc(nc_open(filename.c_str(), NC_WRITE, &fileId));
		check_nc(nc_redef(fileId));
	}
	else check_nc(status);

	int nPointsId;
	int pointIds[3];

	ncdefdimexist(fileId, "no_of_points", nPoints, &nPointsId);

	for (int i = 0 ; i < 3 ; ++i) {
		ncdefvarexist(fileId,
		              (string("points_") + char(i + 'x') + "c").c_str(),
		              NC_FLOAT, 1, &nPointsId, pointIds + i);
	}

	//define cell dimensions
	int varid = -1;
	{
		int dimids[2];
		ncdefdimexist(fileId, "no_of_hexaeders", nCells, &dimids[0]);
		ncdefdimexist(fileId, "points_per_hexaeder", 8, &dimids[1]);
		ncdefvarexist(fileId, "points_of_hexaeders", NC_LONG, 2, dimids, &varid);
	}

	check_nc(nc_enddef(fileId));

	//now writing variables
	//writing points, dimension by dimension
	for (int i = 0 ; i < 3 ; ++i) {
		if (i == 2) {
			for (size_t j = 0 ; j < nPoints ; ++j) {
				float nullfloat = 0;
				check_nc(nc_put_var1_float(fileId,
				                           pointIds[i],
				                           &j, &nullfloat));
			}
		}
		else {
			for (size_t j = 0 ; j < nPoints ; ++j) {
				check_nc(nc_put_var1_float(fileId,
				                           pointIds[i],
				                           &j, &(pos[3*j+i])));
			}
		}
	}

	size_t start[2] = {0, 0}, count[2] = {1, 8};
	long cell[8];
	for (size_t j = 0 ; j < nCells ; j++) {
		std::copy(&ids[8*j], &ids[8*(j+1)], cell);
		check_nc(nc_put_vara_long(fileId, varid,
		                          start, count, cell));
		start[0]++;
	}

	check_nc(nc_close(fileId));
}

struct vec_equal {
	bool operator()(const nvis::fvec3& a, const nvis::fvec3& b) const {
		return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
	}
};

int main(int argc, char* argv[])
{
	initialize(argc, argv);

	avtNek5000FileFormat fmt(metafile);

	if (do_grid) {
		std::vector<unsigned int> hexind;
		std::vector<float>        hexpts;

		fmt.GetMesh(hexpts, hexind);

		int npts = hexpts.size() / 3;

		nvis::fvec3 *pts = (nvis::fvec3*) & hexpts.front();

		std::vector<nvis::fvec3> spts(pts, pts + npts);

		std::sort(spts.begin(), spts.end(), nvis::lexicographical_order());
		spts.resize(std::unique(spts.begin(), spts.end(), vec_equal()) - spts.begin());

		std::cout << spts.size() << " unique points = " << (100*spts.size()) / npts << "%%\n";

		std::cout << "Reindexing..." << std::flush;
		std::vector<unsigned int> indmap(npts);
		std::vector<unsigned int> rindex(spts.size());

		for (unsigned int i = 0; i < npts; ++i) {
			indmap[i] = std::lower_bound(spts.begin(), spts.end(), pts[i], nvis::lexicographical_order()) - spts.begin();
			rindex[indmap[i]] = i;
		}
		rindex.resize(spts.size());

		for (std::vector<unsigned int>::iterator ii = hexind.begin(); ii != hexind.end(); ++ii)
			*ii = indmap[*ii];
		std::cout << " done\n";

		// save reindexing information
		std::cout << "saving reindexing information... " << std::flush;
		std::string idxname = std::string(outfile) + ".idx";
		try {
			ncio::put(idxname.c_str(), rindex, "indices");
		}
		catch (std::runtime_error& e) {
			std::cerr << "exception caught while trying to save indexing information in " << idxname << std::endl;
			std::cerr << e.what() << std::endl;
			return -1;
		}
		std::cout << "done\n";

		hexpts = std::vector<float>();
		hexpts.resize(spts.size()*3);
		std::copy(spts.begin(), spts.end(), (nvis::fvec3*)&hexpts[0]);

		std::vector<unsigned int> hexfind(hexind.size());
		std::copy(hexind.begin(), hexind.end(), hexfind.begin());
		hexind = std::vector<unsigned int>();

		std::string filename = std::string(outfile) + ".grid";
		save_grid(filename, hexpts, hexfind);
	}
	else {
		std::vector<float> vec;
		fmt.GetVectorVar(timestep, vec);

		std::cout << "load reindexing information... " << std::flush;
		std::string idxname = std::string(outfile) + ".idx";
		std::vector<unsigned int> rindex;
		try {
			ncio::get(idxname.c_str(), rindex, "indices");
		}
		catch (std::runtime_error& e) {
			std::cerr << "exception caught while trying to read index file: " << idxname << std::endl;
			std::cerr << e.what() << std::endl;
			return -1;
		}
		std::cout << "done\n";

		size_t npts = rindex.size();
		std::vector<nvis::fvec3> vecd(npts);
		std::cout << "reindexing vector field... " << std::flush;
		for (unsigned int i = 0; i < npts / 3; ++i) {
			for (unsigned int j = 0; j < 3; ++j)
				vecd[i][j] = vec[3*rindex[i] + j];
		}
		std::cout << "done\n";

		std::ostringstream os;
		os << outfile << "-" << timestep << ".val";
		ncio::put(os.str().c_str(), vecd, "velocity");
	}

	return 0;
}

//---------------------------------------------------------------------------
