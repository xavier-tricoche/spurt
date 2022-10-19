#include <vector>
#include <map>
#include <string>
#include <math/fixed_vector.hpp>
#include <math/bounding_box.hpp>

#include <teem/nrrd.h>
#include <fstream>
#include <iostream>


char* ftle_in, *coord_in, *out;
float ratio;

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, (airMopper)hestParmFree, airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input coordinates",    airTypeString,  1, 1, &coord_in,            NULL,       "input coordinates (NRRD)");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &out,                 NULL,       "output file (VTK)");
    hestOptAdd(&hopt, "f",      "ftle values",          airTypeString,  1, 1, &ftle_in,             NULL,       "FTLE file");
    hestOptAdd(&hopt, "r",      "threshold ratio",      airTypeFloat,   1, 1, &ratio,               NULL,       "FTLE cutoff ratio");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Visualize particle systems filtered by FTLE value",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
    Nrrd* nin_coord = nrrdNew();
    if (nrrdLoad(nin_coord, coord_in, NULL)) {
        std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
        return -1;
    }
    float* coord = (float*)nin_coord->data;
    
    Nrrd* nin_ftle = nrrdNew();
    if (nrrdLoad(nin_ftle, ftle_in, NULL)) {
        std::cerr << "ERROR in " << argv[0] << ": " << biffGetDone(NRRD) << std::endl;
        return -1;
    }
    float* ftle = (float*)nin_ftle->data;
    
    int N = nin_coord->axis[1].size;
    
    std::vector<double> __tmp(ftle, &ftle[2*N]);
    std::sort(__tmp.begin(), __tmp.end());
    
    float threshold = ratio*__tmp[floor(ratio*2.*(float)N)];
    
    std::vector<nvis::fvec3>    pos;
    std::vector<float>          val;
    
    for (int i = 0 ; i < N ; ++i) {
        float ftle_p = ftle[2*i  ];
        float ftle_m = ftle[2*i+1];
        
        if (ftle_p > ftle_m && ftle_p > threshold) {
            pos.push_back(nvis::fvec3(coord[3*i], coord[3*i+1], coord[3*i+2]));
            val.push_back(ftle_p);
        }
        
        if (ftle_m > ftle_p && ftle_m > threshold) {
            pos.push_back(nvis::fvec3(coord[3*i], coord[3*i+1], coord[3*i+2]));
            val.push_back(-ftle_m);
        }
    }
    
    std::fstream tmp(out, std::ios::out);
    tmp << "# vtk DataFile Version 2.0\n"
        << "particles with FTLE values for " << coord_in << "\n"
        << "ASCII\n"
        << "DATASET POLYDATA\n"
        << "POINTS " << pos.size() << " float\n";
    for (int i = 0 ; i < pos.size() ; ++i) {
        tmp << pos[i][0] << " " << pos[i][1] << " " << pos[i][2] << '\n';
    }
    tmp << "POINT_DATA " << pos.size() << '\n'
        << "SCALARS ftle float 1\n"
        << "LOOKUP_TABLE default\n";
    for (int i = 0 ; i < val.size() ; ++i) {
        tmp << val[i] << '\n';
    }
    tmp.close();
    std::cerr << out << " was exported\n";
    return -1;
}




