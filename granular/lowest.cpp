#include <iostream>
#include <vector>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <data/raster.hpp>
#include <image/nrrd_wrapper.hpp>
#include <fstream>

#include <teem/hest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

char* name_in, *name_out;
int n;

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",                airTypeString,  1, 1, &name_in,         NULL,       "input TXT file");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &name_out,            NULL,       "output NRRD file");
    hestOptAdd(&hopt, "n",      "npart",                airTypeInt,         1, 1, &n,           NULL,       "number of particles");
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute bounding box of particle assembly over time",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    std::fstream input(name_in, std::ios::in);
    std::vector<double> bounds;
    while (!input.eof()) {
        nvis::bounding_box<nvis::vec2> box;
        double x, y, z, t;
        int id;
        for (int i=0 ; i<n ; ++i) {
            input >> id >> x >> y >> z >> t;
            nvis::vec2 p(x,y);
            box.add(p);
        }
        bounds.push_back(box.min()[0]);
        bounds.push_back(box.min()[1]);
        bounds.push_back(box.max()[0]);
        bounds.push_back(box.max()[1]);
        std::cerr << box << std::endl;
    }
    
    double* data = (double*)calloc(bounds.size(), sizeof(double));
    for (int i=0 ; i<bounds.size() ; ++i) {
        data[i] = bounds[i];
    }
    std::vector<size_t> dims(2);
    dims[0] = 4;
    dims[1] = bounds.size()/4;
    spurt::writeNrrdFromContainers(data, std::string(name_out), /*nrrdTypeDouble,*/ dims);
    
    return 0;
}
