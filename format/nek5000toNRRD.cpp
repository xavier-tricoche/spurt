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

#include "avtNek5000FileFormat.hpp"
#include "ncio.hpp"

#include <boost/lexical_cast.hpp>

#include <image/nrrd_wrapper.hpp>

using namespace std;

char *metafile, *outfile;
int timestep;
int do_grid, do_vals;

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
    hestOptAdd(&hopt, "i", "metafile",       airTypeString, 1, 1, &metafile, NULL, "metafile name");
    hestOptAdd(&hopt, "o", "output prefix",  airTypeString, 1, 1, &outfile,  NULL, "prefix of output file name");
    hestOptAdd(&hopt, "t", "time step",      airTypeInt,    1, 1, &timestep, NULL, "time step");
    hestOptAdd(&hopt, "g", "convert grid",   airTypeInt,    0, 0, &do_grid,  NULL, "convert grid");
    hestOptAdd(&hopt, "v", "convert values", airTypeInt,    0, 0, &do_vals,  NULL, "convert values");

    hestParseOrDie(hopt, argc - 1, const_cast<const char**>(argv) + 1, hparm,
                   me, "Convert Nek5000 data file to NRRD format",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

struct vec_equal {
    bool operator()(const nvis::fvec3& a, const nvis::fvec3& b) const {
        return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
    }
};

//---------------------------------------------------------------------------

void save_grid(size_t nPoints, float *pos, size_t nCells, int *ids)
{
    std::vector<size_t> dims(2);
    std::vector<double> spc(2);
    spc[0] = spc[1] = airNaN();

    dims[0] = 3;
    dims[1] = nPoints;
    std::string name = std::string(outfile) + "-points.nrrd";
    xavier::nrrd_utils::writeNrrdFromContainers(pos, name, dims, spc);

    dims[0] = 8;
    dims[1] = nCells;
    name = std::string(outfile) + "-cells.nrrd";
    xavier::nrrd_utils::writeNrrdFromContainers(ids, name, dims, spc);
}


int main(int argc, char* argv[])
{
    initialize(argc, argv);

    avtNek5000FileFormat fmt(metafile);

    if (!do_grid && !do_vals) {
        std::cerr << "Neither grid nor values need to be exported. Done.\n";
        return -1;
    }

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

        hexpts.clear();
        float *_pts = (float*)calloc(3 * spts.size(), sizeof(float));
        for (int i = 0 ; i < spts.size() ; ++i) {
            _pts[3*i  ] = spts[i][0];
            _pts[3*i+1] = spts[i][1];
            _pts[3*i+2] = spts[i][2];
        }
        size_t nPoints = spts.size();
        spts.clear();

        int *cellids = (int*)calloc(hexind.size(), sizeof(int));
        for (int i = 0 ; i < hexind.size() ; ++i) {
            cellids[i] = hexind[i];
        }
        size_t nCells = hexind.size() / 8;
        hexind.clear();

        save_grid(nPoints, _pts, nCells, cellids);
        delete[] _pts;
        delete[] cellids;
    }

    if (do_vals) {
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
        float *vals = (float*)calloc(3 * npts, sizeof(float));
        std::cout << "reindexing vector field... " << std::flush;
        for (unsigned int i = 0; i < npts; ++i) {
            for (unsigned int j = 0; j < 3; ++j)
                vals[3*i+j] = vec[3*rindex[i] + j];
        }
        std::cout << "done\n";

        std::vector<size_t> dims(2);
        dims[0] = 3;
        dims[1] = npts;
        std::vector<double> spc(2);
        spc[0] = spc[1] = airNaN();

        std::string name = std::string(outfile) + "-vals-t=" + boost::lexical_cast<std::string>(timestep) + ".nrrd";
        xavier::nrrd_utils::writeNrrdFromContainers(vals, name, dims, spc);
        delete[] vals;
    }

    return 0;
}

//---------------------------------------------------------------------------
