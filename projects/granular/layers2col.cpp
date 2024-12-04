#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <graphics/colors.hpp>
#include <image/nrrd_wrapper.hpp>

char* in_name, *out_name;
int nlayers;
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
    hestOptAdd(&hopt, "n",      "# layers",         airTypeInt,     0, 1, &nlayers,     "5",    "number of layers");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute color associated with layers",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    init(argc, argv);
    
    Nrrd* nin = spurt::readNrrd(in_name);
    
    std::vector<float> coords;
    spurt::to_vector<float>(coords, nin);
    
    size_t npart = coords.size()/3;
    
    std::vector<float> y(npart);
    for (int i=0 ; i<npart ; ++i) {
        y[i] = coords[3*i+1];
    }
    
    float min = *std::min_element(y.begin(), y.end());
    float max = *std::max_element(y.begin(), y.end());
    float thickness = (max-min)/(float)(nlayers -1);
    
    
    
    std::vector<nvis::fvec3> color_scale;
    color_scale.push_back(nvis::fvec3(235, 59, 87));
    color_scale.push_back(nvis::fvec3(255, 103, 74));
    color_scale.push_back(nvis::fvec3(255, 151, 87));
    color_scale.push_back(nvis::fvec3(255, 285, 115));
    color_scale.push_back(nvis::fvec3(255, 234, 158));
    color_scale.push_back(nvis::fvec3(252, 254, 194));
    color_scale.push_back(nvis::fvec3(238, 252, 179));
    color_scale.push_back(nvis::fvec3(196, 240, 172));
    color_scale.push_back(nvis::fvec3(147, 224, 180));
    color_scale.push_back(nvis::fvec3(82, 207, 182));
    color_scale.push_back(nvis::fvec3(39, 168, 199));
    
    spurt::adaptive_color_map<float> cmap(y, color_scale);
    
    float* colors = (float*)calloc(3*npart, sizeof(float));
    for (int i=0 ; i<npart ; ++i) {
        nvis::fvec3 c = 1./255.*cmap(y[i]);
        colors[3*i  ] = c[0];
        colors[3*i+1] = c[1];
        colors[3*i+2] = c[2];
    }
    
    size_t size[2] = {3, npart};
    
    Nrrd* nout = nrrdNew();
    if (nrrdWrap_nva(nout, colors, nrrdTypeFloat, 2, size) ||
            nrrdSave(out_name, nout, NULL)) {
        std::cerr << biffGetDone(NRRD) << std::endl;
        exit(-1);
    }
    std::cerr << "exported:\t" << out_name << std::endl;
    nrrdNuke(nout);
    
    return 0;
}


















