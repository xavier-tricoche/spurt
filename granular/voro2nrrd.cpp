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
    hestOptAdd(&hopt, "i",      "input file",       airTypeString,  1, 1, &in_name,     NULL,   "input file name (TXT)");
    hestOptAdd(&hopt, "o",      "output name",      airTypeString,  1, 1, &out_name,    NULL,   "output file name (NRRD)");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Convert Voro++ output format to NRRD",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    init(argc, argv);
    
    std::fstream input(in_name, std::ios::in);
    if (!input) {
        std::cerr << argv[0] << ": unable to open " << in_name;
        exit(-1);
    }
    
    float x, y, z, v;
    int i;
    
    std::map<int, float> density;
    
    while (!input.eof()) {
        input >> i >> x >> y >> z >> v;
        density[i] = v;
    }
    
    int npart = density.back().first();
    std::cerr << "max index is " << npart << '\n';
    
    Nrrd* nin = nrrdNew();
    size_t sz[] = {2, npart+1};
    
    
    xavier::readNrrd(in_name);
    std::vector<float> data;
    xavier::to_vector<float>(data, nin);
    nrrdNuke(nin);
    
    int npart = data.size()/3;
    std::cerr << npart << " particles\n";
    for (int i=0 ; i<npart ; ++i) {
        output << i << " " << data[3*i] << " " << data[3*i+1] << " " << data[3*i+2] << '\n';
    }
    output.close();
    
    return 0;
}