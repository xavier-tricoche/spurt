#include <iostream>
#include <vector>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <data/raster.hpp>
#include <image/nrrd_wrapper.hpp>
#include <sstream>
#include <kdtree/kdtree.hpp>

#include <teem/hest.h>

#ifdef _OPENMP
#include <omp.h>
#endif

double t, pour, relaxation, _gamma, frequency;

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
    hestOptAdd(&hopt, "t",  "time step",    airTypeDouble,  0, 1, &t,               "0",        "simulation time (s)");
    hestOptAdd(&hopt, "p",      "pour",   airTypeDouble,  0, 1, &pour,              "0",            "pouring time (s)");
    hestOptAdd(&hopt, "r",  "relaxation",   airTypeDouble,  0, 1, &relaxation,      "0",        "relaxation time (s)");
    hestOptAdd(&hopt, "f",  "frequency",    airTypeDouble,  0, 1, &frequency,       "7.5",      "tap frequency (Hz)");
    hestOptAdd(&hopt, "g",  "gamma",        airTypeDouble,  0, 1, &_gamma,          "0",        "tap gamma coefficient");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute bump height in tapping simulation",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    double height = 0;
    if (t > pour) {
        double bump_time = 0.5/frequency;
        double tap_time = relaxation + bump_time;
        double loc_t = fmod(t-pour, tap_time);
        
        std::cerr << "t=" << t << ", relaxation=" << relaxation << ", bump_time=" << bump_time
                  << ", tap_time=" << tap_time << ", loc_t=" << loc_t << std::endl;
                  
        if (loc_t < bump_time) {
            double omega = 2.*M_PI*frequency;
            double amplitude = 9.81*_gamma/(omega*omega);
            height = amplitude*sin(omega*loc_t);
        }
    }
    
    std::cout << height << std::flush;
    
    exit(0);
}
