#ifdef _OPENMP
#include <omp.h>
#endif

#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>
#include <teem/nrrd.h>
#include <vector>
#include <image/nrrd_wrapper.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <math.h>
#include <reconstruction/functions.hpp>


// parameters
unsigned int     res[3];
char             *name;
char             *function;
float             _min[3], _max[3];

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
    hestOptAdd(&hopt, "f",    "function",        airTypeString,     1,     1,    &function,        NULL,            "function name (\"franke\", \"sinc\", \"marlobb\")");
    hestOptAdd(&hopt, "r",    "resolution",    airTypeInt,     3,     3,    &res,            "128 128 128",    "size of output volume");
    hestOptAdd(&hopt, "min","minimum",        airTypeFloat,     3,     3,    &_min,            "-1 -1 -1",        "minimum of bounding volume");
    hestOptAdd(&hopt, "max","maximum",        airTypeFloat,     3,     3,    &_max,            "1 1 1",        "maximum of bounding volume");
    hestOptAdd(&hopt, "o",    "output",        airTypeString,     1,     1,    &name,            NULL,            "output name");

    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   me, "Generate procedurally defined 2D/3D scalar dataset",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[]) {
    using namespace reconstruction;
    
    initialize(argc, argv);
    
    int fun = 0;
    if (!strcmp(function, "franke") || !strcmp(function, "Franke")) {
        res[2] = 1;
        fun = 1;
    }
    else if (!strcmp(function, "sinc") || !strcmp(function, "Sinc") || !strcmp(function, "SINC")) {
        fun = 2;
    }
    else if (!strcmp(function, "Marschner") || !strcmp(function, "MarschnerLobb") || 
            !strcmp(function, "marschner") || !strcmp(function, "Marschner-Lobb") ||
            !strcmp(function, "ml") || !strcmp(function, "marschnerlobb") ||
            !strcmp(function, "marschner-lobb") || !strcmp(function, "marlobb")) {
        fun = 3;
    }
    
    std::cerr << "fun = " << fun << std::endl;
    
    float *data = (float*)calloc(res[0]*res[1]*res[2], sizeof(float));
    nvis::vec3 lo(_min[0], _min[1], _min[2]);
    nvis::vec3 hi(_max[0], _max[1], _max[2]);
    nvis::vec3 span = hi - lo;
    nvis::vec3 size(res[0]-1, res[1]-1, res[2]-1);
    nvis::vec3 step = span/size;
    std::cerr << "lo = " << lo << ", hi = " << hi << std::endl;
    
    if (res[2] <= 1) {
        _min[2] = _max[2] = 0.;
    }
    
    for (int k=0 ; k<res[2] ; ++k) {
        for (int j=0 ; j<res[1] ; ++j) {
            for (int i=0 ; i<res[0] ; ++i) {
                int idx = i + res[0]*(j + res[1]*k);
                double x = _min[0] + i*step[0];
                double y = _min[1] + j*step[1];
                double z = _min[2] + k*step[2];
                switch (fun) {
                    case 1: {
                        data[idx] = franke(x,y);
                        break;
                    }
                    case 2: {
                        data[idx] = sinc(x)*sin(y)*sin(z);
                        break;
                    }
                    case 3: {
                        data[idx] = MarschnerLobb(x,y,z);
                        break;
                    }
                    default: {
                        std::cerr << "requested function is unknown\n";
                        exit(-1);
                    }
                }
            }
        }
    }
    
    std::vector<size_t> sz(res[2]>1 ? 3 : 2);
    for (int i=0 ; i<sz.size() ; ++i) sz[i] = size[i]+1;
    std::vector<double> spc(res[2]>1 ? 3 : 2);
    for (int i=0 ; i<spc.size() ; ++i) spc[i] = step[i];
    spurt::writeNrrd(data, name, nrrdTypeFloat, sz, spc);
    
    return 0;
}
