#include <array>
#include <iostream>
#include <assert.h>

#include <teem/nrrd.h>
#include <teem/unrrdu.h>

#include <math/fixed_vector.hpp>
#include "double_point.hpp"
#include <data/raster.hpp>
#include <image/nrrd_wrapper.hpp>
#include <format/filename.hpp>

char* name_out;
size_t nsamples[3];
double rel;

void initialize(int argc, char* argv[], hestOpt* hopt)
{
    hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1,  1,  &name_out,      NULL,   "output name");
    hestOptAdd(&hopt, "s",      "sz0 sz1 sz2",      airTypeSize_t,  3,  3,  nsamples,       NULL,   "number of samples per axis");
    hestOptAdd(&hopt, "r",      "dpl rel distance", airTypeDouble,  0,  1,  &rel,           "0.5",  "relative distance between single point loads in procedural double point load model. Ignored if an input file is selected");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Generate a double point load data set in NRRD format",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    using namespace spurt;
    
    hestOpt *hopt;
    initialize(argc, argv, hopt);

    DoublePointLoad *dpl = new DoublePointLoad(50, 50, 20, rel);
    dpl->set_check_inside(false);
    std::cerr << "generating a synthetic double point load tensor field" << std::endl;
    nvis::ivec3 res(nsamples[0], nsamples[1], nsamples[2]);
    std::cerr << "Resolution = " << res << std::endl;
    raster_grid<3> sampling_grid(res, dpl->bounds());
    int npoints = sampling_grid.size();

    float *data = (float*)calloc(7 * npoints, sizeof(float));
    for (int i = 0 ; i < npoints ; ++i) {
        nvis::ivec3 c = sampling_grid.coordinates(i);
        nvis::vec3 x = sampling_grid(c);
        nvis::fixed_vector<double, 7> t = (*dpl)(x);
        for (int j = 0 ; j < 7 ; ++j) {
            data[i*7+j] = t[j];
        }
    }

    spurt::nrrd_utils::nrrd_params<float, 3> params;
    nvis::vec3 st = sampling_grid.spacing();
    params.spacings()[0] = std::numeric_limits<float>::quiet_NaN();
    params.sizes()[0] = 7;
    params.mins()[0] = std::numeric_limits<float>::quiet_NaN();
    for (int i = 0 ; i < 3 ; ++i) {
        params.spacings()[i+1] = st[i];
        params.sizes()[i+1] = nsamples[i];
        params.mins()[i+1] = dpl->bounds().min()[i];
    }
    params.centers()[0] = nrrdCenterUnknown;
    for (int i=1; i<4; ++i) params.centers()[i] = nrrdCenterNode;
    spurt::nrrd_utils::writeNrrdFromParams(data, name_out, params);

    delete[] data;

    return 0;
}