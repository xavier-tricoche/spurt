#include <iostream>
#include <fstream>
#include <teem/hest.h>
#include <teem/unrrdu.h>

#include <math/fixed_vector.hpp>
#include <misc/time_helper.hpp>
#include <image/nrrd_field.hpp>
#include "ftle.hpp"
#include <tensor/double_point.hpp>
#include <tensor/eigenvector_field.hpp>
#include <math/inverse_transform.hpp>
#include <image/nrrd_wrapper.hpp>

char* name_in;
char* name_out;
double origin[3], e0[3], e1[3];

void initialize(int argc, const char* argv[], hestOpt* hopt)
{
    hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    const char* me = argv[0];
    
    mop = airMopNew();
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "i",      "input",            airTypeString,  1,  1,  &name_in,       NULL,   "input image file");
    hestOptAdd(&hopt, "o",      "output",           airTypeString,  1,  1,  &name_out,      NULL,   "output seed file");
    hestOptAdd(&hopt, "origin", "origin",           airTypeDouble,  3,  3,  origin,         NULL,   "image plane origin");
    hestOptAdd(&hopt, "e0",     "e0",               airTypeDouble,  3,  3,  e0,             NULL,   "x-edge vector");
    hestOptAdd(&hopt, "e1",     "e1",               airTypeDouble,  3,  3,  e1,             NULL,   "y-edge vector");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute seed points from binary image",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

int main(int argc, const char* argv[])
{
    hestOpt* hopt;
    initialize(argc, argv, hopt);
    
    Nrrd* nin = spurt::readNrrd(name_in);
    if (nin->dim != 2) {
        std::cerr << "Input NRRD has invalid dimension (" << nin->dim << " != 2)\n";
        exit(-1);
    }
    
    std::vector<float> values;
    spurt::to_vector(values, nin);
    
    std::cerr << "values contains " << values.size() << " elements\n";
    
    std::cerr << "after conversion to float array: min = " << *std::min_element(values.begin(), values.end())
              << " and max = " << *std::max_element(values.begin(), values.end()) << '\n';
              
    spurt::vec3 o(origin[0], origin[1], origin[2]);
    spurt::vec3 b0(e0[0], e0[1], e0[2]);
    spurt::vec3 b1(e1[0], e1[1], e1[2]);
    b0 /= (double)(nin->axis[0].size - 1);
    b1 /= (double)(nin->axis[1].size - 1);
    
    int N = values.size();
    std::vector<spurt::vec3> seeds;
    for (int n = 0 ; n < N ; ++n) {
        if (values[n] != 0) {
            int i = n % nin->axis[0].size;
            int j = n / nin->axis[0].size;
            spurt::vec3 x = o + (double)i * b0 + (double)j * b1;
            seeds.push_back(x);
        }
    }
    
    std::fstream out(name_out, std::ios::out);
    out << seeds.size() << '\n';
    for (int n = 0 ; n < seeds.size() ; ++n) {
        out << seeds[n][0] << " " << seeds[n][1] << " " << seeds[n][2] << '\n';
    }
    out.close();
    
    nrrdNuke(nin);
    
    return 0;
}




























































































































































































