#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <image/nrrd_wrapper.hpp>
#include <graphics/colors.hpp>

char* in_name, *out_name;
int nlayers;
float thickness;
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
    hestOptAdd(&hopt, "i",      "input file",       airTypeString,  1, 1, &in_name,     NULL,   "input file name (nrrd)");
    hestOptAdd(&hopt, "o",      "output base",      airTypeString,  1, 1, &out_name,    NULL,   "output file name (nrrd)");
    hestOptAdd(&hopt, "n",      "# layers",         airTypeInt,     1, 1, &nlayers,     NULL,   "number of layers");
    hestOptAdd(&hopt, "t",      "thickness",        airTypeFloat,   1, 1, &thickness,   NULL,   "layer thickness");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Extract layers in particle column and assign colors to them",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    init(argc, argv);
    
    Nrrd* nin = spurt::readNrrd(in_name);
    
    std::vector<float> coords;
    spurt::to_vector<float>(coords, nin);
    
    int npart = coords.size()/3;
    
    std::vector<float> y(npart);
    for (int i=0 ; i<npart ; ++i) {
        y[i] = coords[3*i+1];
    }
    
    float min = *std::min_element(y.begin(), y.end());
    float max = *std::max_element(y.begin(), y.end());
    float eps = (max - min - nlayers*thickness)/(float)(nlayers-1);
    float slice = thickness + eps;
    
    typedef std::pair<int, int> layer_type;
    std::vector<layer_type> layers;
    
    for (int i=0 ; i<y.size() ; ++i) {
        int l = floor((y[i]-min)/slice);
        float u = y[i] - l*slice;
        if (u<thickness) {
            layers.push_back(std::make_pair(i, l));
        }
    }
    
    int* data = (int*)calloc(2*layers.size(), sizeof(int));
    for (int i=0 ; i<layers.size() ; ++i) {
        data[2*i  ] = layers[i].first;
        data[2*i+1] = layers[i].second;
    }
    
    size_t size[2] = {2, layers.size()};
    
    Nrrd* nout = nrrdNew();
    if (nrrdWrap_nva(nout, data, nrrdTypeInt, 2, size) ||
            nrrdSave(out_name, nout, NULL)) {
        std::cerr << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    std::cerr << "exported:\t" << out_name << std::endl;
    nrrdNuke(nout);
    
    return 0;
}


















