#include <iostream>
#include <vector>
#include <string>
#include <math/fixed_vector.hpp>
#include <data/raster.hpp>
#include <image/nrrd_wrapper.hpp>
#include <sstream>
#include <iterator>

#include <teem/hest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

char* in;
double t, relaxation, _gamma, frequency;
int range[2];

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
    hestOptAdd(&hopt, "i",      "input",        airTypeString,  1, 1, &in,              NULL,       "input NRRD file");
    // hestOptAdd(&hopt, "t",       "time step",    airTypeDouble,  0, 1, &t,               "0",        "simulation time (s)");
    //  hestOptAdd(&hopt, "r",      "relaxation",   airTypeDouble,  0, 1, &relaxation,      "0",        "relaxation time (s)");
    //  hestOptAdd(&hopt, "f",      "frequency",    airTypeDouble,  0, 1, &frequency,       "7.5",      "tap frequency (Hz)");
    //  hestOptAdd(&hopt, "g",      "gamma",        airTypeDouble,  0, 1, &_gamma,          "0",        "tap gamma coefficient");
    hestOptAdd(&hopt, "r",  "valid range",  airTypeInt,     2, 2, &range,           "0 1",      "valid Y range (index space)");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute wave location in 1D density plot",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    Nrrd* nin = spurt::readNrrd(in);
    
    std::vector<float> data;
    spurt::to_vector<float>(data, nin);
    
    std::cout << std::distance(data.begin(), std::min_element(data.begin()+range[0], data.begin()+range[1]+1));
    
    exit(0);
}
