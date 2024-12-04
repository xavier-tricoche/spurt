#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <teem/nrrd.h>
#include <image/nrrd_wrapper.hpp>

char* in_name, *out_name;
void init(int argc, char* argv[])
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
    hestOptAdd(&hopt, "i",      "input file",       airTypeString,  1, 1, &in_name,     NULL,   "input file name (NRRD)");
    hestOptAdd(&hopt, "o",      "output name",      airTypeString,  1, 1, &out_name,    NULL,   "output file name (TXT)");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Convert NRRD particle dataset to Voro++ compatible text format",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    init(argc, argv);
    
    Nrrd* nin = spurt::readNrrd(in_name);
    std::vector<float> data;
    spurt::to_vector<float>(data, nin);
    
    std::fstream output(out_name, std::ios::out);
    if (!output) {
        std::cerr << argv[0] << ": unable to open " << out_name;
        exit(-1);
    }
    nrrdNuke(nin);
    
    int npart = data.size()/3;
    std::cerr << npart << " particles\n";
    for (int i=0 ; i<npart ; ++i) {
        output << i << " " << data[3*i] << " " << data[3*i+1] << " " << data[3*i+2] << '\n';
    }
    output.close();
    
    return 0;
}