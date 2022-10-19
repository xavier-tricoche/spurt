#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <image/nrrd_wrapper.hpp>
#include <graphics/colors.hpp>
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
    hestOptAdd(&hopt, "i",      "input file",       airTypeString,  1, 1, &in_name,     NULL,   "input file name (nrrd)");
    hestOptAdd(&hopt, "o",      "output base",      airTypeString,  1, 1, &out_name,    NULL,   "output file name (nrrd)");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute color associated with each particle based on its y-coordinate",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    init(argc, argv);
    
    Nrrd* nin = xavier::readNrrd(in_name);
    
    std::vector<float> coords;
    xavier::to_vector<float>(coords, nin);
    
    size_t npart = coords.size()/3;
    
    std::vector<float> y(npart);
    for (int i=0 ; i<npart ; ++i) {
        y[i] = coords[3*i+1];
    }
    
    std::vector<nvis::fvec3> color_scale;
    // color_scale.push_back(nvis::fvec3(235, 59, 87));
    // color_scale.push_back(nvis::fvec3(255, 103, 74));
    // color_scale.push_back(nvis::fvec3(255, 151, 87));
    // color_scale.push_back(nvis::fvec3(255, 185, 115));
    // color_scale.push_back(nvis::fvec3(255, 234, 158));
    // color_scale.push_back(nvis::fvec3(252, 254, 194));
    // color_scale.push_back(nvis::fvec3(238, 252, 179));
    // color_scale.push_back(nvis::fvec3(196, 240, 172));
    // color_scale.push_back(nvis::fvec3(147, 224, 180));
    // color_scale.push_back(nvis::fvec3(82, 207, 182));
    // color_scale.push_back(nvis::fvec3(39, 168, 199));
    color_scale.push_back(nvis::fvec3(128, 0, 0));
    color_scale.push_back(nvis::fvec3(192, 0, 0));
    color_scale.push_back(nvis::fvec3(255, 0, 0));
    color_scale.push_back(nvis::fvec3(255, 64, 0));
    color_scale.push_back(nvis::fvec3(255, 128, 0));
    color_scale.push_back(nvis::fvec3(255, 192, 0));
    color_scale.push_back(nvis::fvec3(255, 255, 0));
    // color_scale.push_back(nvis::fvec3(128, 255, 0));
    // color_scale.push_back(nvis::fvec3(0, 255, 0));
    // color_scale.push_back(nvis::fvec3(0, 128, 128));
    // color_scale.push_back(nvis::fvec3(0, 0, 255));
    // color_scale.push_back(nvis::fvec3(0, 128, 255));
    color_scale.push_back(nvis::fvec3(255, 255, 64));
    color_scale.push_back(nvis::fvec3(255, 255, 128));
    color_scale.push_back(nvis::fvec3(255, 255, 192));
    color_scale.push_back(nvis::fvec3(255, 255, 255));
    
    float min = *std::min_element(y.begin(), y.end());
    float max = *std::max_element(y.begin(), y.end());
    float step = (max - min)/(float)(color_scale.size()-1);
    
    float* colors = (float*)calloc(3*npart, sizeof(float));
    for (int i=0 ; i<npart ; ++i) {
    
        int j = floor((y[i]-min)/step);
        float u = (y[i]-min) - j*step;
        nvis::fvec3 c = (1.-u) * color_scale[j] + u*color_scale[j+1];
        
        std::cerr << "color = " << c << std::endl;
        
        colors[3*i  ] = c[0]/255.;
        colors[3*i+1] = c[1]/255.;
        colors[3*i+2] = c[2]/255.;
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


















