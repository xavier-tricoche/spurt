#include <teem/nrrd.h>
#include <iostream>
#include <math/fixed_vector.hpp>

char* input;
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
    hestOptAdd(&hopt, "i",  "input",    airTypeString,  1, 1, &input,   NULL,   "input file name");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Read point coordinates from file",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    Nrrd* nin = nrrdNew();
    if (nrrdLoad(nin, input, NULL)) {
        std::cerr << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    
    int dim = nin->axis[0].size;
    int npts = nin->axis[1].size;
    float* coords = (float*)nin->data;
    
    for (int i=0 ; i<std::min(npts, 100) ; ++i) {
        std::cout << "position #" << i << " has coordinates (";
        for (int j=0 ; j<dim-1 ; ++j) {
            std::cout << coords[dim*i+j] << ", ";
        }
        std::cout << coords[dim*(i+1)-1] << ")\n";
    }
    
    return 0;
}
