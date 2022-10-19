#include <teem/nrrd.h>
#include <iostream>

int dim, npts;
char* output;
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
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1, 1, &output,      NULL,       "output file name");
    hestOptAdd(&hopt, "d",      "dimension",        airTypeInt,     0, 1, &dim,         "2",        "space dimension");
    hestOptAdd(&hopt, "n",      "# points",         airTypeInt,     0, 1, &npts,        "100",      "number of randomly selected points");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Generate a random point set in arbitrary dimensions (bounding box is [-1, 1]^N)",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    srand48(time(0));
    
    float* coord = (float*)calloc(npts*dim, sizeof(float));
    for (int i=0 ; i<npts*dim ; ++i) {
        coord[i] = -1. + 2.*drand48();
    }
    
    Nrrd* nout = nrrdNew();
    size_t size[2] = {dim, npts};
    nrrdWrap_nva(nout, coord, nrrdTypeFloat, 2, size);
    nrrdSave(output, nout, NULL);
}
